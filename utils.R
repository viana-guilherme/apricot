# reference:
# https://www.protocols.io/view/label-free-quantification-lfq-proteomic-data-analy-5qpvobk7xl4o/v2

# ------
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
  input_sorted <- sort(unique(x))
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

# Significantly changing proteins are defined as:

#    a p-value (or adjusted p-value) < 0.05
#    a fold change of > 2 (UP) or < 0.5 (DOWN)

# ----- TODO?
# PCA plot calculations
# gene ontologies
# File sanitation
#
