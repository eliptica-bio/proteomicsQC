#' Plots feature counts
#'
#' @param \strong{report} DIA-NN report file
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{Q_THR} Q-value filtering threshold (Default)
#' @param \strong{features} list of features to plot, e.g. c("Protein.Group", "Protein.Ids", "Precursor.Id")
#' @param run_order_var variable for sample running order
#' @return
#' @export
#'
#' @examples
#' plotFeatureCounts(diann_report, Q_THR = 0.01, features = c("Protein.Ids", "Precursor.Id"))
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
plotFeatureCounts <- function(report, metadata, Q_THR = 0.01, run_order_var = "run_order", features = c("Protein.Group", "Protein.Ids", "Precursor.Id")) {
  report %>%
    filter(across(any_of("Q.Value"), ~.x < Q_THR)) %>%
    select(all_of(c("File.Name", features))) %>%
    pivot_longer(cols = all_of(features)) %>% distinct() -> features_selected

  if(!missing(metadata)) {
    features_selected %>%
      ungroup() %>%
      left_join(metadata %>% select(File.Name, !!as.name(run_order_var)), by = "File.Name") %>%
      mutate(File.Name = fct_reorder(File.Name, !!as.name(run_order_var))) -> features_selected
  }
  features_selected %>%
    ggplot(aes(x=File.Name)) +
    geom_bar() +
    stat_count() +
    facet_wrap(~name, ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") -> p
  return(p)
}




#' Plots missingness matrix
#'
#' @param \strong{report} DIA-NN report file
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{feature_var} feature for plotting missingness values
#' @param run_order_var variable for sample running order
#' @param subtitle string for a subtitle
#' @param prevalence_filter filter out non-prevalent proteins, default: 0.1
#' @param \strong{Q_THR} Q-value filtering threshold, default: 0.01
#'
#' @return
#' @export
#'
#' @examples
#' metadata <- create_metadata(diann_report)
#' plotMissingness(diann_report, metadata = metadata, prevalence_filter = 0.5)
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
#' @import forcats
plotMissingness <- function(report, metadata, Q_THR = 0.01,  feature_var = "Precursor.Id", run_order_var = "run_order", subtitle = "", prevalence_filter = 0.1) {

  feature_prevalence <- getPrevalence(report = report, Q_THR = Q_THR, feature_var = feature_var)

  report %>%
    filter(across(any_of("Q.Value"), ~.x < Q_THR)) %>%
    select(all_of(c("File.Name", feature_var))) %>%
    semi_join(feature_prevalence %>% filter(prevalence > prevalence_filter), by = feature_var) -> data_tmp

  data_tmp %>% mutate(isPresent = 1) -> data_tmp
  data_tmp %>%
    ungroup() %>%
    complete(File.Name, !!as.name(feature_var), fill = list(isPresent = 0)) -> data_tmp_complete

  completeness = round(data_tmp_complete %>% ungroup() %>% summarise(sum(isPresent)/length(isPresent)) %>% pull, 3)

  feature_per_sample <- data_tmp_complete %>% group_by(File.Name) %>% summarise(feature_count = sum(isPresent)) %>% pull %>% median

  if(!missing(metadata)) {
    data_tmp_complete %>%
      ungroup() %>%
      left_join(metadata, by = "File.Name") %>%
      mutate(File.Name = fct_reorder(File.Name, !!as.name(run_order_var))) -> data_tmp_complete
  }


  data_tmp_complete %>%
    group_by(!!as.name(feature_var)) %>%
    mutate(inN = sum(isPresent, na.rm = T)) %>% ungroup() %>%
    mutate(!!feature_var := fct_reorder(!!as.name(feature_var), inN),
           isPresent = as_factor(isPresent)) -> toPlot

  toPlot %>%
    ggplot() +
    geom_tile(aes(x = File.Name, y = !!as.name(feature_var), fill = isPresent)) +
    scale_fill_brewer() +
    labs(subtitle = subtitle,
         x = "Running order",
         caption = paste("Completeness:", completeness,
                         "Prevalence >:", prevalence_filter,
                         "Features (median):", feature_per_sample,
                         sep = " " )) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background = element_blank(), legend.position = "bottom", aspect.ratio = 5/8) -> p

  return(p)
}


#TODO
#' Plots combined DIA-NN experiment with Golden standarts metrics
#'
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param golden_thr threshold for QC line, default = 0.8
#' @param qc_pattern QC file pattern
#' @param \strong{data_counts} output from \code{\link{countInGolden}}
#' @return list of plot and samples bellow golden_thr
#' @export
#'
#' @examples
#' LIN_PROTOCOL <- "5-min-sswath"
#' LIN_N <- 10
#' LIN_COR <- 0.8
#' LIN_PROTOCOL = "5-min-sswath"
#' golden_linear <- linear_peptides_5minsswath %>% ungroup %>% filter(protocol == LIN_PROTOCOL, n >= LIN_N, cor >= LIN_COR)
#'
#' countInGolden(report = diann_report,
#'              golden_standard = golden_linear,
#'              QC_type = paste(LIN_PROTOCOL, LIN_N, LIN_COR, sep = ":" )) -> golden_linear_counts
#'
#' metadata <- create_metadata(diann_report)
#' ret = plotGoldenCounts(golden_linear_counts, metadata, golden_thr = 0.1)
#'
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
#' @import forcats
#' @importFrom fs as_fs_path
plotGoldenCounts <- function(data_counts, metadata, golden_thr = 0.9, qc_pattern = ".*?QC.*?") {

  if(!missing(metadata)) {
    data_counts %>%
      ungroup() %>%
      left_join(metadata, by = "File.Name") %>%
      mutate(File.Name = fct_reorder(File.Name, run_order)) -> data_counts
  }

  total_samples <- data_counts %>% distinct(File.Name) %>% nrow
  total_qc <- data_counts %>% filter(grepl(pattern = qc_pattern, perl = T, ignore.case = T,
                                           x = basename(as.character(fs::as_fs_path(as.character(File.Name)))))) %>% distinct(File.Name) %>% nrow

  lowQC <- data_counts %>%
    filter(isGolden == 1) %>%
    mutate(isLowQC = ifelse(round(golden_fraction, 2) < golden_thr, 1, 0)) %>%
    filter(isLowQC == 1) %>%
    select(File.Name, golden_n, golden_fraction, n) %>%
    mutate(golden_thr = golden_thr)

  total_QCpassed = data_counts %>% anti_join(lowQC, by = "File.Name") %>% distinct(File.Name) %>% nrow()

  data_counts %>%
    ungroup() %>%
    mutate(isGolden = as_factor(isGolden)) %>%
    ggplot() +
    geom_bar(stat = "identity", aes(x = File.Name, y = n, fill = isGolden)) +
    geom_point(data = lowQC, aes(x = File.Name, y = n), colour = "blue", size = 3, shape = 17) +
    labs(subtitle = paste("QC line: ", unique(data_counts$QC_type), "QC presence threshold: ", golden_thr),
         ylab = "Number of Precursor.ID",
         caption = paste("Total samples:", total_samples,
                         "QC passed:", round(total_QCpassed/total_samples, 2),
                         "QCs:", total_qc, sep = " ") ) +
    geom_hline(yintercept = unique(data_counts$golden_n)*golden_thr, linetype = 2) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") -> p

  return(list(plot = p, lowQC = lowQC, PASSED = 1 - nrow(lowQC)/total_samples))
}

#' Plots Experiment
#'
#' @param \strong{report} DIA-NN report file
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{Q_THR} Q-value filtering threshold, default: 0.01
#' @param \strong{feature_value} feature value to plot, default: "Precursor.Quantity"
#' @param \strong{feature_var} feature variable, default: "Precursor.Id"
#' @param run_order_var variable for sample running order
#' @param subtitle a string for a plot's subtitle
#' @param qc_pattern QC file pattern
#'
#' @return
#' @export
#'
#' @examples
#' metadata <- create_metadata(diann_report)
#' plotExperiment(diann_report, metadata = metadata, Q_THR - 0.01)
#'
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
#' @import forcats
plotExperiment <- function(report, metadata, Q_THR = 0.01, feature_var = "Precursor.Id", feature_value = "Precursor.Quantity",  run_order_var = "run_order", subtitle = "",  qc_pattern = ".*?QC.*?") {

  # report = diann_report
  #
  # metadata = metadata
  # Q_THR = 0.01
  # feature_var = "Precursor.Id"
  # run_order_var = "run_order"
  # subtitle = "TEST"
  # qc_pattern = ".*?QC.*?"
  # feature_value = "Precursor.Quantity"

  report %>%
    filter(across(any_of("Q.Value"), ~.x < Q_THR)) %>%
    select(all_of(c("File.Name", feature_var, feature_value))) -> dataset




  if(!missing(metadata)) {
    dataset %>%
      ungroup() %>%
      left_join(metadata, by = "File.Name") %>%
      mutate(File.Name = fct_reorder(File.Name, !!as.name(run_order_var))) -> dataset
  }

  total_samples <- dataset %>% distinct(File.Name) %>% nrow
  total_qc <- dataset %>% filter(grepl(pattern = qc_pattern, perl = T, ignore.case = T,
                                       x = basename(as.character(fs::as_fs_path(as.character(File.Name)))))) %>% distinct(File.Name) %>% nrow

  dataset %>%
    mutate(type = ifelse(str_detect(File.Name, qc_pattern), "QC", "sample")) -> toPlot

  toPlot %>% left_join(metadata) %>%  group_by(File.Name, type) %>% summarise(TIC = sum(!!as.name(feature_value), na.rm = T), n  = n()) -> toPlot.summary

  toPlot %>%
    left_join(metadata) %>%
    mutate(File.Name = fct_reorder(File.Name, !!as.name(run_order_var))) %>%

    ggplot(aes(x = File.Name, y = log(!!as.name(feature_value)), fill=type)) +
    stat_boxplot(geom ='errorbar', width = 0.6) +
    geom_boxplot() +
    geom_point(data = toPlot.summary, aes(x =File.Name, y = log(TIC), colour = type)) +
    labs(subtitle = subtitle,
         x = "Running order",
         caption = paste("Total Samples:", total_samples,
                         "Total QCs:", total_qc,
                         sep = " ")) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") -> p

  return(p)
}


#' Plots outliers for Total Ion Count (TIC) and number of features
#'
#' @param \strong{report_stats} TIC and n stats, output from \code{\link{countStats}}
#' @param \strong{Z_THR} Z-score filtering threshold, default: 3
#' @param \strong{stats} list of stats to plot, e.g. c("Z_n", "Z_TIC", "Zmod_n", "Zmod_TIC")
#'
#'
#' @return a list of plots of lenght 2
#' @export
#'
#' @examples
#' report_stats <- countStats(diann_report)
#' plots <- plotStats(report_stats, Z_THR = 3, stats = c("Zmod_n", "Z_TIC", "Zmod_TIC"))
#' plots[[1]]/plots[[2]]
#'
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
plotStats <- function(report_stats, Z_THR = 3, stats = c("Z_n", "Z_TIC", "Zmod_n", "Zmod_TIC") ) {
  #stats = c("Z_n", "Z_TIC", "Zmod_n", "Zmod_TIC")
  #Z_THR = 3

  TIC_flag = 0
  n_flag = 0

  if(any(!(stats %in% c("Z_n", "Z_TIC", "Zmod_n", "Zmod_TIC")))) {
    stop("Wrong stats specified")
  }

  if ( any(as_vector(lapply(stats, str_detect, "TIC"))) ) {
    TIC_flag = 1
  }

  if (any(as_vector(lapply(stats, str_detect, "n"))) ) {
    n_flag = 1
  }

  ret = list()

  if (TIC_flag) {
    report_stats %>%
      select(File.Name, TIC, n, all_of(stats)) %>%
      select(File.Name, matches("TIC")) %>%
      pivot_longer(cols = starts_with("Z")) -> toPlot

    toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR) %>%
      group_by(name) %>%
      summarise(x_val = max(abs(value), na.rm = T),
                n = n(),
                y_val = 5) -> outliers


    toPlot %>%
      ggplot(aes(x = TIC)) +
      geom_histogram(bins = 50, fill = "lightgrey", colour = "white") +
      geom_jitter(data = toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR), aes(x = TIC, y = 1 ), colour = "red") +
      facet_wrap(~name, ncol = 1, scales = "free") +
      ggrepel::geom_text_repel(data = toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR),
                               aes(x = TIC, y = 1, label = basename(fs::as_fs_path(as.character(File.Name)))), max.overlaps = 100, size = 3.5) +
      geom_label(data = outliers, aes(x = x_val, y = 5, label = paste("Total outliers:", n)), color = "red") +
      scale_x_log10() +
      labs(x = "Total signal per sample, log10(TIC)",
           subtitle = paste("Outlier abs Z-score threshold > " , Z_THR, sep ="")) +
      theme_bw() -> p

    ret = lappend(ret, p)

  }

  if (n_flag) {
    report_stats %>%
      select(File.Name, TIC, n, all_of(stats)) %>%
      select(File.Name, matches("n")) %>%
      pivot_longer(cols = starts_with("Z")) -> toPlot

    toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR) %>%
      group_by(name) %>%
      summarise(x_val = max(abs(value), na.rm = T),
                n = n(),
                y_val = 5) -> outliers

    toPlot %>%
      ggplot(aes(x = n)) +
      geom_histogram(bins = 50, fill = "lightgrey", colour = "white") +
      geom_jitter(data = toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR), aes(x = n, y = 1 ), colour = "red") +
      facet_wrap(~name, ncol = 1, scales = "free") +
      ggrepel::geom_text_repel(data = toPlot %>% filter(round(abs(value), digits = 1) >= Z_THR),
                               aes(x = n, y = 1, label = basename(fs::as_fs_path(as.character(File.Name)))), max.overlaps = 100, size = 3.5) +
      geom_label(data = outliers, aes(x = x_val, y = 5, label = paste("Total outliers:", n)), color = "red") +
      labs(x = "Number or identified precursors per sample",
           subtitle = paste("Outlier abs Z-score threshold > " , Z_THR, sep ="")) +
      theme_bw() -> p

    ret = lappend(ret, p)
  }
  return(ret)
}


