# reference:
# https://www.protocols.io/view/label-free-quantification-lfq-proteomic-data-analy-5qpvobk7xl4o/v2

#------- Imports and data handling --------
diann_import <- function(file_path) {
  # Importing the report file
  raw_data <- readr::read_csv(file = file_path)
  db <- raw_data |>
    # removes ambiguous records (Proteotypic == 0, means peptide matched more than 1 ptn)
    dplyr::filter(Proteotypic != 0) |>
    # removes unneeded columns
    dplyr::select(-c(First.Protein.Description:Precursor.Charge)) |>
    # condenses the value cols into a single col
    # separates sample and replicate information into different cols
    tidyr::pivot_longer(cols = where(is.numeric),
                        names_to = c("Samples", "Replicate"),
                        names_sep = "-",
                        values_to = "peptide_intensities") |>
    # changes NAs to 0
    dplyr::mutate(peptide_intensities = tidyr::replace_na(peptide_intensities, 0))


  return(db)

  # note on the use of Protein.Names to group samples:
  # It seems to be the least ambiguous information we have available. Sometimes, Protein
  # I've noticed on a dataset that Protein.Ids could be ""Q88FA7", and  "Q88FA7;CP0001", even though names and group are uniquely identified as ""Q88FA7" "Q88FA7", and proteotypic != 0
}

diann_collapse <- function(database) {

  # this is needed as a separate function to allow the top3 method to integrate in our workflow

  # sum peptide intensities to make protein intensities
  db_collapse <- database |>
                  dplyr::group_by(Protein.Names, Samples, Replicate) |>
                  dplyr::summarise(value_sum = sum(peptide_intensities), .groups = "drop")

  return(db_collapse)

}

diann_summary <- function(database) {

  # checks if the data was already processed. If not, do it
  if (!"value_sum" %in% names(database)) {
    database <- diann_collapse(database)
  }

  # generating the summarised version of the db with some main statistics
  db_summary <- database |>
    dplyr::group_by(Protein.Names, Samples) |>
    dplyr::summarise(Counts_mean = mean(value_sum),
                     Counts_std = sd(value_sum),
                     `%CV` = Counts_std / Counts_mean * 100,
                     .groups = "drop")

  return(db_summary)

}

find_lod <- function(x) {
  input_sorted <- base::sort(base::unique(x))
  lod <- ifelse(input_sorted[1] == 0, input_sorted[2], input_sorted[1])

  return(lod)
}

# ------ T-tests ---------

diann_ttest <- function(data, group1, group2) {

  # change data so that the levels order match the user preference
  # (for determining fold-change)
  data <- data |>
    dplyr::mutate(Samples = forcats::fct_relevel(Samples, c(group1, group2)))

  # retain each row's information about the proteins being compared
  ptn_info <- data |>
    dplyr::select(-c(Samples:value_sum)) |>
    dplyr::distinct()

  # run Welch's t-test
  ttest <- try(
    t.test(data = data, value_sum ~ Samples), silent = TRUE)

  # retrieve t-test results, unless it's an error
  if (class(ttest) != "try-error") {

    ttest_results <- c(
      ttest$estimate[1],
      ttest$estimate[2],
      statistic = ttest$statistic,
      p.value = ttest$p.value)

  } else {

    ttest_results <- rep(NA, times = 4)
  }

  # define names for the ttest results vector
  names(ttest_results) <- c("log2_mean_A", "log2_mean_B", "statistic", "p.value")

  # assemble output
  output <- c(ptn_info, A = group1, B = group2, ttest_results)

  return(output)

  # TODO - sd info
  # since R's t.test() does not compute sd, we could do it ourselves
  # this is a bit of a hack, but it works for most cases
  # the issue is how to make sure the order is maintained for the user choice
  # sd <- data |>
  #   group_by(Samples) |>
  #   mutate(Samples = glue::glue("log2_std_{Samples}")) |>
  #   summarise(std = sd(value_sum)) |>
  #   deframe()
}
batch_ttest <- function(data, a, b) {
  message(glue::glue("running batch Welch's T-test for sample {a} over sample {b}..."))

  sample_selection <- c(a, b)

  # preparing the database for the t-test
  ttest_db <- data |>
    dplyr::filter(Samples %in% sample_selection) |>
    dplyr::group_by(Samples) |>
    dplyr::mutate(value_sum = ifelse(value_sum == 0, find_lod(value_sum), value_sum), #  inputing zeroes with the lowest of detected (LOD) value in each sample
           value_sum = log2(value_sum)) |> # log2 transforming Abundance values prior to the t-Test
    dplyr::ungroup() |>

    # separate the dataframe into a list with individual proteins for comparison
    dplyr::group_by(Protein.Names) |>
    dplyr::group_split() |>

    # running the t-test with our custom function
    purrr::map_dfr(.f = diann_ttest, group1 = sample_selection[1], group2 = sample_selection[2]) |>
    dplyr::mutate(`log2FC_A/B` = log2(2 ** log2_mean_A / 2 ** log2_mean_B))

    # including the BH correction considering all of the measured p-values
    ttest_db$`p.adjusted.BH` <- p.adjust(ttest_db$p.value, method = "BH")

  message(glue::glue("finished running analysis!"))

  return(ttest_db)

}

# ------ Gene ontology ---------

runEnrichment <- function(dataset, significance_threshold) {

  library(topGO)


  # 1) Loading GO annotations for P. putida KT2440
  GO_db_list <- rlist::list.load(file = "databases/GO_db_list.json") # returns the GO_db_list object from the JSON file created by createGOAnnot.R

  message("read GO_db_list")

  # 2) read the pre-built proteomic metadata and the user-supplied dataset table
  metadata <- readr::read_delim(file = "databases/PutidaMetadataDB.tsv")
  message("read metadata")

  universe <- dataset
  message("read universe")

  # 3) Attach the metadata column to the universe to select info
  db <- dplyr::left_join(universe, metadata, dplyr::join_by(Protein.Names == id))

  ## We need step 3) in order to convert the protein accessions to gene loci !
  universe_converted <- db |>
    dplyr::select(locusTag, p.adjusted.BH) |>
    tibble::deframe()

  # 4) Defining a selection criteria for the enrichment tests (i.e. our significant genes)
  target_genes <- function(universe_converted) universe_converted < significance_threshold

  ### Building the topGO data

  # ## TODO -> this section is extremely verbose. Look into ways of automating it....

  GOdata_BP <- new("topGOdata",
                   description = "Enrichment analysis - Biological process",
                   ontology = "BP",
                   allGenes = universe_converted,
                   geneSel = target_genes,
                   annot = annFUN.gene2GO,
                   gene2GO = GO_db_list)

  GOdata_CC <- new("topGOdata",
                   description = "Enrichment analysis - Cellular component",
                   ontology = "CC",
                   allGenes = universe_converted,
                   geneSel = target_genes,
                   annot = annFUN.gene2GO,
                   gene2GO = GO_db_list)

  GOdata_MF <- new("topGOdata",
                   description = "Enrichment analysis - Molecular Function",
                   ontology = "MF",
                   allGenes = universe_converted,
                   geneSel = target_genes,
                   annot = annFUN.gene2GO,
                   gene2GO = GO_db_list)

  test_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
  test_CC <- runTest(GOdata_CC, algorithm = "weight01", statistic = "fisher")
  test_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")


  results_BP <- GenTable(GOdata_BP , Fisher = test_BP) |>
    dplyr::mutate(GO.Domain = "Biological process")
  results_CC <- GenTable(GOdata_CC , Fisher = test_CC) |>
    dplyr::mutate(GO.Domain = "Cellular component")
  results_MF <- GenTable(GOdata_MF , Fisher = test_MF) |>
    dplyr::mutate(GO.Domain = "Molecular Function")



  result_table <- dplyr::bind_rows(results_BP, results_MF, results_CC) |>
    dplyr::mutate(Fisher = as.numeric(Fisher),
                  Annotated = as.numeric(Annotated),
                  Significant = as.numeric(Significant),
                  Expected = as.numeric(Expected)) |>
    dplyr::relocate(GO.Domain)

  message(print(result_table))

  return(result_table)

  detach("package:topGO", unload = TRUE)


}


#------ PCA -----------


diann_pca <- function(database) {

  # data needs to be wide for the clustering to work
  pca_data <- database |>
    dplyr::select(Protein.Names, Samples, Replicate, value_sum) |>
    dplyr::mutate(sample_fixed = glue::glue("{Samples}_{Replicate}")) |>
    dplyr::select(-Samples, -Replicate) |>
    tidyr::pivot_wider(names_from = Protein.Names,
                       values_from = value_sum) |>
    janitor::remove_constant() # important in case there's a column with all 0's

  pca_res <- pca_data |>
    tibble::column_to_rownames(var = "sample_fixed") |>
    prcomp(scale = TRUE)

  return(pca_res)
}

diann_pca_view <- function(pca_res) {

  # obtaining info on explained variance for axes labels
  explained_variance <- pca_res$sdev ** 2 / sum(pca_res$sdev ** 2)

  # plotting
  res1$x |>
    tibble::as_tibble(rownames = "Sample") |>
    ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::coord_fixed(1) +
      ggplot2::labs(x = glue::glue("PC1 ({round(explained_variance[1] * 100, 2)}%)"),
                    y = glue::glue("PC2 ({round(explained_variance[2] * 100, 2)}%)")) +
      ggplot2::theme_light()

  # todo: allow users to color using pre-defined groups as arguments
}

# ----- Top3 Calculations ---------------



# ----- TODO? ---------------
# File sanitation

