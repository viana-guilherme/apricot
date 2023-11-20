# reference:
# https://www.protocols.io/view/label-free-quantification-lfq-proteomic-data-analy-5qpvobk7xl4o/v2

#------- Imports and data handling --------
diann_import <- function(file_path, filter.type = "unique") {
  # Importing the report file
  raw_data <- readr::read_csv(file = file_path) |>
    # checks and corrects weirdly formatted column names
    sanitize_colnames()

  # TODO: Replace this with a sensibe case_when()
  if (filter.type == "proteotypic") { # removes ambiguous records (Proteotypic == 0, means peptide matched more than 1 ptn)
    proteo.filter <- parse(text = "Proteotypic != 0")
  } else {
    if (filter.type == "none") {# in practice does not filter anything
    proteo.filter <- parse(text = "Proteotypic >= 0")
    } else {
      if (filter.type == "unique") { # this seems to be the filtering used by the original analysis
        proteo.filter <- parse(text = "stringr::str_detect(Protein.Names, pattern = ';', negate = TRUE)")
      }
    }
  }

  db <- raw_data |>
      dplyr::filter(eval(proteo.filter)) |>
    # removes unneeded columns
    dplyr::select(-c(First.Protein.Description:Precursor.Charge, Protein.Ids)) |>
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
  # I've noticed on a dataset that Protein.Ids could be ""Q88FA7", and  "Q88FA7;CP0001",
  # even though names and group are uniquely identified as ""Q88FA7" "Q88FA7", and proteotypic != 0
}

diann_collapse <- function(diann_database, attach.locusTag = TRUE) {

  # this is needed as a separate function to allow the top3 method to integrate in our workflow

  # sum peptide intensities to make protein intensities
  db_collapse <- diann_database |>
                  dplyr::group_by(Protein.Names, Samples, Replicate) |>
                  dplyr::summarise(value_sum = sum(peptide_intensities), .groups = "drop") |>

    # Calculate percent abundance:
    # ((intensity of individual protein / sum of all protein intensities in a given sample) * 100)
    dplyr::mutate(percent_abundance = value_sum / sum(value_sum))

  if (attach.locusTag) {
    db_collapse <- dplyr::left_join(db_collapse, UniprotIds, by = c("Protein.Names" = "UNIPROT"), multiple = "first")

    # funnily enough, there might be a many-to-many relationship here (see Q877U6_PSEPK --> PP_1157 AND PP_3365)
    # a full list can be found with: UniprotIds[duplicated(UniprotIds$UNIPROT),] |> dplyr::arrange(UNIPROT)


  }

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

sanitize_colnames <- function(raw_data) {

  # checks if the column names in input file are a file path
  # and correct them
  # otherwise just return the input file original colname

  if (any(stringr::str_ends(colnames(raw_data), "\\\\.+raw"))) {
    a <- colnames(raw_data) |>  dplyr::if_else(
                                        stringr::str_ends(colnames(raw_data), "\\.+raw"),
                                        stringr::str_extract(colnames(raw_data), "(?<=\\-\\d{3}\\\\).+(?=.raw)"),
                                        colnames(raw_data)
                                        )
  }

  return(raw_data)
}

#------- Top3 Calculations ---------------

# top 3 method

# step 1) Filter the DIA-NN peptide report data (from import) to only proteins that have three or more peptides identified across all samples
# obs: the way we have the data better organized, we can it easily by using the "db" obj
# as returned by the diann_import function

find_eligible <- function(diann_database) {

  eligible_ptn <- diann_database |>
    dplyr::group_by(Protein.Names, Samples, Replicate) |> #here we use Protein.Names as an unique id
    dplyr::tally() |>
    dplyr::filter(n >= 3) |>
    dplyr::arrange(desc(n))

    return(eligible_ptn)

# (i'm still not sure if this is the best option, but it looks elegant for now)
}

# step 2) For each protein, rank the top 3 peptides by intensity (counts) in each of the samples
find_top3 <- function(diann_database, eligible) {

  eligible_list <- unique(eligible$Protein.Names)

  top3 <- diann_database |>
    dplyr::filter(Protein.Names %in% eligible_list) |>
    dplyr::group_by(Protein.Names, Samples, Replicate) |>
    dplyr::mutate(rank = rank(peptide_intensities, na.last = FALSE)) |>
    dplyr::slice_max(n = 3, order_by = peptide_intensities)

  return(top3)
}


# Calculate protein intensity  of the Top3 peptides by averaging the intensity (counts)
find_top3_intensities <- function(diann_database, eligible) {

  top3_mean_ranks <- find_top3(diann_database, eligible) |>
    dplyr::ungroup() |>
    dplyr::group_by(Protein.Names, Precursor.Id) |>
    dplyr::summarise(mean_rank = mean(rank)) |>

    # Filter the data to the three highest ranked peptides in each protein
    dplyr::slice_max(n = 3, order_by = mean_rank)

    # Calculate the percent of the total protein abundance
    diann_database |>
      dplyr::filter(Precursor.Id %in% top3_mean_ranks$Precursor.Id) |>
      dplyr::group_by(Protein.Names, Samples, Replicate) |>
      dplyr::summarise(protein_intensity = mean(peptide_intensities))
}

# make nice qc plots

diann_qc_counts <- function(diann_database, top3_database, output.nrow = 7) {

  # The idea of this function is to quickly generate a barplot with vital information
  # about the quality of the samples

  # 1) retrieve info on the total amount of proteins found across samples
  total_ptn <- diann_database |> dplyr::pull(Protein.Names) |> unique() |> length()

  # 2) retrieve information on the total amount of peptides in each sample
  peptides_per_sample <- diann_database |>
    dplyr::filter(peptide_intensities != 0) |>
    dplyr::group_by(Samples, Replicate) |>
    dplyr::summarise(ptn_count = length(unique(Protein.Names))) |>
    dplyr::mutate(unique_name = glue::glue("{Samples}_{Replicate}"), total = total_ptn)

  # 3) retrieve information on the top3 peptides in each sample
  peptides_per_sample_top3 <- top3_database |>
    dplyr::filter(protein_intensity != 0) |>
    dplyr::group_by(Samples, Replicate) |>
    dplyr::summarise(ptn_count = length(unique(Protein.Names))) |>
    dplyr::mutate(unique_name = glue::glue("{Samples}_{Replicate}"),
           total = total_ptn)

  # plot the information
  ggplot2::ggplot(peptides_per_sample,
                  ggplot2::aes(x = ptn_count, y = forcats::fct_rev(Replicate))) +

    # add bars
    ggplot2::geom_bar(ggplot2::aes(x = total, fill = "lightblue"),
                      stat = "identity", position = "dodge", width = 0.9) +

    ggplot2::geom_bar(ggplot2::aes(fill = "goldenrod2"),
                      stat = "identity", position = "dodge", width = 0.9) +

    ggplot2::geom_bar(data = peptides_per_sample_top3,
                      ggplot2::aes(x = ptn_count, y = forcats::fct_rev(Replicate), fill = "darkseagreen"),
                      stat = "identity", position = "dodge", width = 0.9) +

    # add the annotations
    ggplot2::geom_text(ggplot2::aes( x = ptn_count + 90, label = ptn_count), color = "darkgoldenrod4") +
    ggplot2::geom_text(ggplot2::aes( x = total - 100, label = total), color = "steelblue") +
    ggplot2::geom_text(data = peptides_per_sample_top3, ggplot2::aes( x = 100, label = ptn_count), color = "darkgreen") +

    # cosmetic adjustments
    ggplot2::facet_wrap(~Samples, nrow = output.nrow) +
    ggplot2::scale_fill_identity(guide = "legend",
                                 name = "Number of proteins",
                                 breaks = c("lightblue","goldenrod2", "darkseagreen"),
                                 labels = c("Total in dataset", "Total in sample", "Top3 method")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Counts", y = "") +
    ggplot2::theme_minimal(base_size = 14)
}

##### GO Enrichment and annotation

## t-test and enrichment

UniprotIds <- read.delim("./UP000000556_160488.idmapping.gz", sep = "\t", header = FALSE) |>
  dplyr::filter(V2 %in% c("UniProtKB-ID", "Gene_OrderedLocusName")) |>
  tidyr::pivot_wider(names_from = V2, values_from = V3, values_fn = list) |>
  tidyr::unnest(cols = c(Gene_OrderedLocusName, `UniProtKB-ID`)) |>
  dplyr::select(-1)

colnames(UniprotIds) <- c("GID", "UNIPROT")


#------- Simple bar plots ----

# TODO
### maybe highlight and add more info to the top protein in both functions below?

view_top_proteins <- function(dataset = data_wt_collapsed, sample = "IP5D_exp_ac", show.n = 15) {

  # this function works with a collapsed data set
  dataset |>

    # pre-processing the data frame
    dplyr::filter(Samples == sample) |>
    dplyr::group_by(Protein.Names, GID) |>
    dplyr::summarise(mean = mean(value_sum), sd = sd(value_sum), .groups = "drop") |>
    dplyr::slice_max(order_by = mean, n = show.n) |>

    # making a basic bar plot
    ggplot(aes(y = mean, x = forcats::fct_reorder(GID, mean, .desc = TRUE))) +
    geom_bar(stat = "identity", color = "black", fill = "steelblue") +
    geom_errorbar(aes(ymin = mean - sd , ymax = mean + sd), width = 0.25) +
    theme_light()

}

view_ttest <- function(ttest_table, direction = "up", n.show = 15,  p.threshold = 0.05, use.p.adj = FALSE) {

  # determines which p-value column to use for enrichment
  if (use.p.adj) {
    pvalue_choice <- "p.adjusted.BH"

  } else {
    pvalue_choice <- "p.value"
  }

  # determines whether to show up- or down- regulated genes
  if (direction == "up") {
    cutoff <- parse(text = "dplyr::slice_max(order_by = `log2FC_A/B`, n = n.show)")
    } else {
      if (direction == "down") {
        cutoff <- parse(text = "dplyr::slice_min(order_by = `log2FC_A/B`, n = n.show)")
        } else {
          message("invalid direction! Please use either 'up' or 'down'")
          break()
        }
      }

  # This function allows us to visualize the top 'n.show' up/downregulated genes
  ttest_table |>
    dplyr::filter(p.value <= p.threshold) |>
    eval(cutoff) |>
    ggplot(aes(y = `log2FC_A/B`, x = forcats::fct_reorder(GID, `log2FC_A/B`, .desc = TRUE))) +
    geom_bar(stat = "identity")

}


#------- PCA -----------

diann_pca <- function(diann_database) {

  # obs: tried writing a function to use the raw data as is, but the fact we
  # have several rows for the same protein in the same sample causes issues
  # Therefore, we start with the summed (collapsed) peptide values

  # checks if the data was already processed. If not, do it
  if (!"value_sum" %in% names(diann_database)) {
    diann_database <- diann_collapse(diann_database)
  }

  # data needs to be wide for
  pca_data <- diann_database |>
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

diann_pca_view <- function(diann_pca_results, color.filter = "samples") {

  # defines the logic for coloring
  if (color.filter == "samples") {
    mutate.call <- parse(text = "Sample")
  } else {
    mutate.call <- parse(text = "fix_names(Sample, color.filter)")
  }

  # obtaining info on explained variance for axes labels
  explained_variance <- diann_pca_results$sdev ** 2 / sum(diann_pca_results$sdev ** 2)

  # preparing the data frame

   lol <- diann_pca_results$x |>
    tibble::as_tibble(rownames = "Sample") |>
     dplyr::mutate(color = eval(mutate.call),
                   color = dplyr::if_else(color == "", "Others", color))


    # plotting
    ggplot2::ggplot(lol, ggplot2::aes(x = PC1, y = PC2, color = color)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::scale_color_manual(values = rep(
                                               RColorBrewer::brewer.pal(n = 9, "Set1"),
                                               times = length(unique(lol$Sample))),
                                  name = "Groups") +
      ggplot2::coord_fixed(1) +
      ggplot2::labs(x = glue::glue("PC1 ({round(explained_variance[1] * 100, 2)}%)"),
                    y = glue::glue("PC2 ({round(explained_variance[2] * 100, 2)}%)")) +
      ggplot2::theme_light()

  # todo: allow users to color using pre-defined groups as arguments
}

#------- T-tests ---------

format_ttest <- function(data, group1, group2) {

  # change data so that the levels order match the user preference
  # (for determining fold-change)
  data <- data |>
    dplyr::mutate(Samples = forcats::fct_relevel(Samples, c(group1, group2)))

  # retain each row's information about the proteins being compared
  ptn_info <- data |>
    dplyr::select(-c(Samples:percent_abundance)) |>
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
diann_ttest <- function(data, a, b) {
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

#------- Gene ontology ---------


prepare_go_enrichment <- function(ttest.table, direction = "up", direction.threshold = 4, p.threshold = 0.05, use.p.adj = TRUE)  {

# determines which pvalue column to use for enrichment
if (use.p.adj) {
  pvalue_choice <- "p.adjusted.BH"

} else {
  pvalue_choice <- "p.value"
}


if (direction == "up") {
  dir.filter <- parse(text = "logFC >= direction.threshold")
  } else {
  if (direction == "down") {
    dir.filter <- parse(text = "logFC <= -1 * direction.threshold")
    } else {
    message("invalid direction! Please use either 'up' or 'down'")
    break
    }
  }

# prepares the dataframe to match the clusterProfiler requisites
input_format <- data.frame(ttest.table[pvalue_choice],
                           ttest.table$`log2FC_A/B`,
                           ttest.table$Protein.Names) |>
  na.omit() |>
  dplyr::rename_with(~c("pvalue", "logFC","protein")) |>
  dplyr::left_join(UniprotIds, by = c("protein" = "UNIPROT"))

# extracts the significant proteins in the t-test table, performing the respective filters
significantList <- input_format |>
  dplyr::filter(pvalue <= p.threshold, eval(dir.filter)) |>
  dplyr::pull(GID) |>
  unique()


# extract all of the proteins present in the sample (aka universe)
universeList <- unique(input_format$GID)

return(list(significant = significantList, universe = universeList))

}

calculateFoldEnrichment <- function(go_results) {

  results <- go_results@result |>
                tibble::as_tibble() |>
                dplyr::mutate(sig.Match = stringr::str_extract(GeneRatio, pattern = "^.+(?=\\/)") |> as.numeric(),
                              sig.Total = stringr::str_extract(GeneRatio, pattern = "(?<=\\/).+") |> as.numeric(),
                              bg.Match = stringr::str_extract(BgRatio, pattern = "^.+(?=\\/)") |> as.numeric(),
                              bg.Total = stringr::str_extract(BgRatio, pattern = "(?<=\\/).+") |> as.numeric(),
                              FoldEnrichment = (sig.Match/sig.Total) / (bg.Match/bg.Total)) |>
                dplyr::select(ONTOLOGY, ID, Description, pvalue, p.adjust, qvalue, geneID, FoldEnrichment)

  return(results)
}


run_go_enrichment <- function(dataset, significance_threshold, direction, use.p.adj = FALSE, direction.threshold = 2) {

  # first we run prepare_go_enrichment to retrieve the list of significant / universe genes
  toEnrich <- prepare_go_enrichment(dataset,
                                    direction = direction,
                                    use.p.adj = use.p.adj,
                                    direction.threshold = direction.threshold)

  enriched <- clusterProfiler::enrichGO(
    gene = toEnrich$significant,
    universe = toEnrich$universe,
    org.Pputida.eg.db,
    keyType = "GID",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    pool = FALSE
  )


  output <- calculateFoldEnrichment(enriched)

  return(output)

}

view_go_enrichment <- function(dataset, fill.var = "p.adjust") {

  fill_var <- parse(text = fill.var)

  ggplot(dataset,
         aes(x = FoldEnrichment,
             y = forcats::fct_reorder(Description, FoldEnrichment),
             fill = eval(fill_var))) +
    geom_bar(stat = "identity") +
    facet_wrap(~ONTOLOGY, ncol = 1, scales = "free_y") +
    scale_fill_distiller(palette = "Blues", direction = -1) +
    theme_light()
}



#------- Volcano plots ---------


view_volcanoplot <- function(dataset = ttests_wt, index = 1,
                             up.threshold = 2, down.threshold = -2, p.threshold = 0.05, label.max = 5, label.by = "GID", use.p.adj = FALSE) {


  # determines which pvalue column to use for plotting and subsetting
  if (use.p.adj) {
    pvalue_choice <- parse(text = "p.adjusted.BH")

  } else {
    pvalue_choice <- parse(text = "p.value")
  }

  # retrieves information about the dataset
  dataset_name <- names(dataset[index])
  label.by <- parse(text = label.by)

  # now we apply some transformations to the data set for the data viz
  prepare_ttest <- dataset[[index]] |>
    na.omit() |>
    dplyr::mutate(
      # sets the groupings
      category = dplyr::case_when(`log2FC_A/B` >= up.threshold & eval(pvalue_choice) <= p.threshold ~ "up-regulated",
                                  `log2FC_A/B` <= down.threshold  & eval(pvalue_choice) <= p.threshold ~ "down-regulated",
                                  .default = "not significant"),

      # sets the manual color scale
         color = dplyr::case_when(`log2FC_A/B` >= up.threshold & eval(pvalue_choice) <= p.threshold~ "orange",
                                  `log2FC_A/B` <= down.threshold  & eval(pvalue_choice) <= p.threshold~ "steelblue",
                                  .default = "grey"),

      # prepares the column of choice for pvalues
                  pvalue_to_plot = -log10(eval(pvalue_choice))
      )

  # subseting the data further for cosmetic adjustments
  manual_colors <- prepare_ttest |> dplyr::select(category,color) |> tibble::deframe()
  to_annotate_up <-  prepare_ttest |> dplyr::filter(eval(pvalue_choice) <= p.threshold) |> dplyr::slice_max(n = label.max, order_by = `log2FC_A/B`)
  to_annotate_down <-  prepare_ttest |> dplyr::filter(eval(pvalue_choice) <= p.threshold) |> dplyr::slice_min(n = label.max, order_by = `log2FC_A/B`)

  # making the final plot
  ggplot(prepare_ttest, aes(x = `log2FC_A/B`, y = pvalue_to_plot, color = category)) +

    #geoms
    geom_point(alpha = 0.95, size = 4) +
    geom_hline(yintercept = -log10(p.threshold), linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = up.threshold, linetype = "dashed", linewidth = 0.2) +
    geom_vline(xintercept = down.threshold, linetype = "dashed", linewidth = 0.2) +
    ggrepel::geom_text_repel(data = to_annotate_up, aes(label = eval(label.by)), size = 4, force = 2, max.time = 1, force_pull = 0.1, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = to_annotate_down, aes(label = eval(label.by)), size = 4,max.time = 1, force = 2, force_pull = 0.1,  show.legend = FALSE) +

    # cosmetic options
    scale_color_manual(values = manual_colors, name = "Significance", breaks = c("up-regulated", "down-regulated", "not significant")) +
    xlim(-max(abs(prepare_ttest$`log2FC_A/B`)), max(abs(prepare_ttest$`log2FC_A/B`))) + # centers the x-axis around 0
    labs(x = "Log2 Fold Change", y = glue::glue("-log10({as.character(pvalue_choice)})"), title = dataset_name,
         subtitle = glue::glue("up.threshold = {up.threshold}, down.threshold = {down.threshold}, p.threshold = {p.threshold}, use.p.adj = {use.p.adj}")) +
    theme_light()

}
# ----- Little helpers ---------------

# identifies how many unique samples there are in a diann output
unique_samples <- function(dataset) {
  return(unique(dataset$Samples))
}

# this one is used in the diann_pca_view function
fix_names <- function(chr, color.filter) {

  matches <- stringr::str_extract_all(chr, color.filter, simplify = FALSE)

  fixed <- purrr::map_chr(.x = matches, ~ stringr::str_flatten(.x, collapse = '-'))

  return(fixed)
}

# ----- TODO ---------------

# 1. File sanitation

# 2. diann qc should also return a table (so that the result is not re-run every time)

# 3. why are results slightly different from the 'canonical' way? <= UNDER INVESTIGATION

# 4. check if in collapsed table the percent abundance sums up to 1
  # it does, but for the entire table --> need to do it PER SAMPLE

# 5. check if NA inputation is working well --> compare t-tests thomas


