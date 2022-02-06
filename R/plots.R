#' Plots feature counts
#'
#' @param \strong{report} DIA-NN report file
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{Q_THR} Q-value filtering threshold (Default)
#' @param \strong{features} list of features to plot, e.g. c("Protein.Group", "Protein.Ids", "Precursor.Id")
#'
#' @return
#' @export
#'
#' @examples
#' plotFeatureCounts(diann_report, Q_THR = 0.01, features = c("Protein.Ids", "Precursor.Id"))
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import ggplot2
#' @import tidyr
plotFeatureCounts <- function(report, metadata, Q_THR = 0.01, features = c("Protein.Group", "Protein.Ids", "Precursor.Id")) {
  report %>%
    filter(Q.Value < Q_THR) %>%
    select(all_of(c("File.Name", features))) %>%
    pivot_longer(cols = all_of(features)) %>% distinct() -> features_selected

  if(!missing(metadata)) {
    features_selected %>%
      ungroup() %>%
      left_join(metadata %>% select(File.Name, run_order), by = "File.Name") %>%
      mutate(File.Name = fct_reorder(File.Name, run_order)) -> features_selected
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

#' Plots missingness values
#'
#' @param \strong{report} DIA-NN report file
#' @param metadata experiment description to arrange by run_order output from \code{\link{create_metadata}}
#' @param \strong{feature_var} feature for plotting missingness values
#' @param run_order_var sample running order
#' @param subtitle string for a subtitle
#' @param prevalence_filter filter out non-prevalent proteins, default: 0.1
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
plotMissingness <- function(report, metadata, feature_var = "Precursor.Id", run_order_var = "run_order", subtitle = "", prevalence_filter = 0.1) {

  feature_prevalence <- getPrevalence(report, feature_var = feature_var)

  report %>%
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
          panel.background = element_blank(), legend.position = "bottom", aspect.ratio = 5/8)
}

#TODO
#' Plots
#'
#' @param data_counts
#' @param metadata
#' @param golden_thr
#' @param qc_pattern
#'
#' @return list of plot and samples bellow golden_thr
#' @export
#'
#' @examples
plotGoldenCounts <- function(data_counts, metadata, golden_thr = 0.9, qc_pattern = ".*?QC.*?") {

  #data_counts = golden_linear_counts # from countInGolden()
  #golden_thr = 0.8
  #qc_pattern = ".*?QC.*?"

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

  return(list(plot = p, lowQC = lowQC))
}



