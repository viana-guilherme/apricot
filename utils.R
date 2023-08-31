# reference:
# https://www.protocols.io/view/label-free-quantification-lfq-proteomic-data-analy-5qpvobk7xl4o/v2

# ------ T-tests ---------
import_diann <- function(file_path) {
  # Importing the report file
  db <- readr::read_csv(file = file_path) |>
    # removes unneeded columns
    dplyr::select(-c(First.Protein.Description:Precursor.Id)) |>
    # condenses the value cols into a single col
    # separates sample and replicate information into different cols
    tidyr::pivot_longer(cols = where(is.numeric),
                        names_to = c("Samples", "Replicate"),
                        names_sep = "-",
                        values_to = "peptide_intensities") |>
    # changes NAs to 0
    dplyr::mutate(peptide_intensities = tidyr::replace_na(peptide_intensities, 0)) |>
    # sum peptide intensities to make protein intensities
    dplyr::group_by(Protein.Group, Protein.Ids, Protein.Names, Genes, Samples, Replicate) |>
    dplyr::summarise(value_sum = sum(peptide_intensities), .groups = "drop")

  # generating the summarised version of the db with some main statistics
  db_summary <- db |>
    dplyr::group_by(Protein.Group, Protein.Ids, Protein.Names, Genes, Samples) |>
    dplyr::summarise(Counts_mean = mean(value_sum),
                     Counts_std = sd(value_sum),
                     `%CV` = Counts_std / Counts_mean * 100,
                     .groups = "drop")

  return(db)
}
find_lod <- function(x) {
  input_sorted <- base::sort(base::unique(x))
  lod <- ifelse(input_sorted[1] == 0, input_sorted[2], input_sorted[1])

  return(lod)
}
format_ttest <- function(data, group1, group2) {

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
run_ttest <- function(data, a, b) {
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
    purrr::map_dfr(.f = format_ttest, group1 = sample_selection[1], group2 = sample_selection[2]) |>
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


# ----- TODO? ---------------
# PCA plot calculations
# File sanitation
#
