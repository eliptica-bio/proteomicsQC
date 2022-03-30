#' Creates metadata from DIA-NN report
#'
#' @param \strong{report} DIA-NN report file
#' @param \strong{INTO} extracted columns in metadata output file, e.g.
#' \code{c("Batch", "sample_id", "run_order", "plate_pos",  "well")}
#' @param \strong{file_pattern} regular expression for extracting experimental metadata, e.g.
#' \code{"(.*?)_(.*?)_([0-9]+)_([0-9]+)_([A-Za-z0-9]+)\\\\.wiff\\\\.dia$"}
#'
#' @return returns metadata tibble
#' @export
#'
#' @examples
#' create_metadata(diann_report,
#'                 INTO = c("Batch", "sample_id", "run_order", "plate_pos",  "well"),
#'                 file_pattern = "(.*?)_(.*?)_([0-9]+)_([0-9]+)_([A-Za-z0-9]+)\\.wiff\\.dia$")
#'
create_metadata = function(report, INTO =  c("Batch", "sample_id", "run_order", "plate_pos",  "well"),
                           file_pattern = "(.*?)_(.*?)_([0-9]+)_([0-9]+)_([A-Za-z0-9]+)\\.wiff\\.dia$") {

  report %>%
    select(all_of("File.Name")) %>%
    distinct(File.Name) %>%
    mutate(base_File.Name = basename(as.character(fs::as_fs_path(File.Name)))) %>%
    extract(col = base_File.Name,
            regex = file_pattern,
            into = INTO, remove = T) %>%
    mutate(run_order = {if("run_order" %in% names(.)) as.numeric(run_order) else NULL}) -> ret

  if("well" %in% names(ret)) {
    ret %>% mutate(well = ifelse(well == "NA", NA, well)) -> ret
    ret %>%
      separate(well,
               into = c("row", "column"),
               sep = "(?<=[A-Za-z])(?=[0-9])",
               remove = F,
      ) -> ret
    ret %>% mutate(column = {if("column" %in% names(.)) as.numeric(column) else NULL}) -> ret
  }

  return(ret)
}

#' Counts prevalence of features (e.g. Precursor.Id) from DIA-NN report
#'
#' @param \strong{report} DIA-NN report file
#' @param \strong{feature_var} a single feature such as Precursor.Id
#' @param \strong{Q_THR} Q-value filtering threshold, default: 0.01
#' @return df with feature prevalence values
#' @export
#'
#' @examples
#' Precursor.Id_prevalence <- getPrevalence(diann_report, feature_var = "Precursor.Id")
#' @import dplyr
#' @importFrom magrittr %>%
#' @import tidyr

getPrevalence <- function(report, Q_THR = 0.01, feature_var = "Precursor.Id") {

  report %>%
    filter(Q.Value <= Q_THR) %>%
    select(all_of(c("File.Name", feature_var))) %>% distinct() -> data_tmp

  data_tmp %>% mutate(isPresent = 1) -> data_tmp

  data_tmp %>%
    ungroup() %>%
    complete(File.Name, !!as.name(feature_var), fill = list(isPresent = 0)) -> data_tmp_complete

  n_samples = data_tmp %>% distinct(File.Name) %>% nrow
  data_tmp_complete %>%
    ungroup() %>%
    group_by(!!as.name(feature_var)) %>%
    summarise(prevalence = sum(isPresent, na.rm = T)/n_samples) -> feature_prevalence

  return(feature_prevalence)
}

#' Removes trend from a batch of DIA-NN experiment
#'
#' @param \strong{report} DIA-NN report file
#' @param \strong{metadata} metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{QC_regex} regexp to extract QC values for fitting model
#' @param \strong{fit_model} model to fit, one of \code{lm}, \code{loess}, \code{rq} (robust regression using medians)
#' @param \strong{feature_value} feature value to remove trend from, default: "Precursor.Quantity"
#' @param \strong{feature_var} feature variable, default: "Precursor.Id"
#' @return returns tbl with corrected feature values from \code{feature_var}
#' @export
#' @examples
#' metadata <- create_metadata(diann_report)
#' diann_report magrittr::`%>%`()
#'    dplyr::filter(!grepl(pattern = "blank", x="File.Name", ignore.case = T)) magrittr::`%>%`()
#'    removeTrend(metadata = metadata, QC_regex = "MSQC", fit_model = "lm",
#'    feature_var = "Precursor.Id" , feature_value = "Precursor.Quantity" ) -> Precursor.Id_corrected
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import tidyr
#' @importFrom purrr map map2
#' @importFrom quantreg rq

removeTrend <- function(report, metadata, QC_regex = "MSQC", fit_model = "loess", feature_var = "Precursor.Id", feature_value = "Precursor.Quantity") {
  # report <- diann_report
  # feature_var = "Precursor.Id"
  # feature_value = "Precursor.Quantity"
  # QC_regex = "MSQC"
  #
  report %>%
    select(all_of(c("File.Name", feature_var, feature_value))) -> dataset

  if (!(fit_model %in% c("lm", "loess", "rq"))) {
    stop("Wrong fit model specified", call. = TRUE, domain = NULL)
  }

  dataset %>% filter(grepl(QC_regex, File.Name)) -> data_filtered_qc_pooled

  if (nrow(data_filtered_qc_pooled) < 4) {
    stop(paste("check regular expression match for QC samples, few samples returned", nrow(data_filtered_qc_pooled)), call. = TRUE, domain = NULL)
  }


  data_filtered_qc_pooled %>%
    left_join(metadata, by = "File.Name") %>%
    group_by(!!as.name(feature_var)) %>%
    nest() -> data_nested_qc

  model_formula = formula(paste(feature_value,  "~ run_order", sep = ""))
  qc_model = NULL
  if (fit_model == "lm") {
    try(data_nested_qc %>% group_by(!!as.name(feature_var)) %>% mutate(model = map(data, .f = function(x) { lm(model_formula, data = x, na.action=na.exclude)})) -> qc_model)
  }

  if (fit_model == "loess") {
    try(data_nested_qc %>% group_by(!!as.name(feature_var)) %>% mutate(model   = map(data, .f = function(x) { loess(model_formula, data = x, na.action=na.exclude, control = loess.control(surface = "direct"))})) -> qc_model)
  }

  if (fit_model == "rq") {
    try(data_nested_qc %>% group_by(!!as.name(feature_var)) %>% mutate(model   = map(data, .f = function(x) { quantreg::rq(model_formula, data = x, na.action=na.exclude)})) -> qc_model)
  }

  if (is.null(qc_model)) {
    stop("No model was fitted")
  }

  dataset %>% left_join(metadata, by = "File.Name") %>% ungroup() %>% group_by(!!as.name(feature_var)) %>% nest() -> dataset_nested

  qc_model %>%
    left_join(dataset_nested, by = c(feature_var), suffix = c(".train", ".test")) %>%

    mutate(predicted = map2(model, data.test,
                            function(x, y) {
                              pred = predict(x, y)
                              return(y %>% mutate(predicted = pred))
                            })) -> dataset_predicted


  correct_feature <- function(df, feature_value) {
    varname <- paste0(as_label(enquo(feature_value)), "_corrected", sep = "")
    df %>%
      mutate(!!varname := ({{feature_value}} - predicted) + median(predicted, na.rm = T)) %>%
      mutate(!!varname := !!as.name(varname) + 1.5*abs(min(!!as.name(varname), na.rm = T))) # to prevent corrected variables become <0
  }

  dataset_predicted %>%
    select(Precursor.Id, predicted) %>% unnest(predicted) %>%
    group_by(Precursor.Id) %>%
    correct_feature(!!as.name(feature_value)) -> dataset_corrected

  feature_value_corrected = paste0(feature_value, "_corrected", sep = "")
  dataset_corrected %>%
    select(all_of(c("File.Name", feature_var, feature_value, feature_value_corrected))) -> ret

  return(ret)
}

#
#' Counts features from DIA-NN report e.g. Precursor.Id in golden standard list
#'
#' @param \strong{report} DIA-NN report file
#' @param \strong{golden_standard} golden standard df with \code{feature_var}
#' @param \strong{feature_var} feature variable to count
#' @param QC_type comment, default: "linear peptides"
#'
#' @return
#' @export
#'
#' @examples
#' LIN_PROTOCOL = "5-min-sswath"
#' LIN_N = 10
#' LIN_COR = 0.8
#' golden_linear = linear_peptides_5minsswath %>% ungroup %>% filter(protocol == LIN_PROTOCOL, n >= LIN_N, cor >= LIN_COR)
#' countInGolden(report = diann_report,
#'              golden_standard = golden_linear,
#'              QC_type = paste(LIN_PROTOCOL, LIN_N, LIN_COR, sep = ":" )) -> golden_linear_counts
#'
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import tidyr

countInGolden = function(report, golden_standard,
                         feature_var = "Precursor.Id",
                         QC_type = "linear_peptides" ) {

  golden_standard %>% ungroup() %>% select(all_of(feature_var)) %>% distinct() -> golden_list
  golden_list %>% ungroup() %>% summarise(n()) %>% pull -> golden_n

  report %>%
    select(all_of(c("File.Name", feature_var))) %>%
    full_join(golden_list %>% mutate(isGolden = 1), by = feature_var) %>%
    filter(!is.na(File.Name)) %>% #filters out rows (files with goldem features that dont exist in entire batch)
    mutate(isGolden = ifelse(is.na(isGolden), 0, isGolden)) %>%
    ungroup() %>%
    group_by_at(.vars = c("File.Name", "isGolden")) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(golden_n = golden_n,
           golden_fraction = n/golden_n) %>% #fraction of standard
    mutate(QC_type = QC_type) -> dataset_summary

  return(dataset_summary)
}



#' Counts Z-score stats of precursors
#'
#' @param \strong{report} DIA-NN report file
#' @param \strong{Q_THR} Q-value filtering threshold, default: 0.01
#' @param \strong{feature_value} feature value to count
#'
#' @return returns tbl with Z-scores of feature values from \code{feature_value}
#' @export
#'
#' @examples
#' diann_report %>%
#'  filter(!grepl(pattern = ".*?BLANK.*?", perl = T, ignore.case = T, x = File.Name)) %>% #removes blanks
#'         countStats() -> report_stats
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @import tidyr
countStats <- function(report, Q_THR = 0.01, feature_value = "Precursor.Quantity") {
  report %>%
    filter(Q.Value <= Q_THR) %>%
    group_by(File.Name) %>%
    summarise(TIC = sum(!!as.name(feature_value)),
              n  = n()) %>%
    ungroup() %>%
    mutate(Z_n = (n - mean(n, na.rm = T))/sd(n, na.rm = T),
           Z_TIC = (TIC - mean(TIC, na.rm = T))/sd(TIC, na.rm = T),
           Zmod_n = (n - median(n, na.rm = T))/(1.486*mad(n, na.rm = T)), # modified Z-score using medians https://www.ibm.com/docs/en/cognos-analytics/11.1.0?topic=terms-modified-z-score
           Zmod_TIC = (TIC - median(TIC, na.rm = T))/(1.486*mad(TIC, na.rm = T))
           ) -> report_summary

  return(report_summary)
}

