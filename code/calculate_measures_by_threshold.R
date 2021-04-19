calculate_measures_by_threshold <- function(prob, y_truth, positive_value) {
  library(dplyr)
  real_positives <- sum(y_truth == positive_value)
  
  as.data.frame(t(sapply(seq(0, 1, length.out = 1000), function(thresh) {
    true_positives <- sum((prob >= thresh) & (y_truth == positive_value))
    det_positives <- sum(prob >= thresh)
    c(thresh = thresh, 
      prec = true_positives / det_positives, 
      rec = true_positives / real_positives)
  }))) %>% 
    mutate(prec = ifelse(is.nan(prec), 1, prec))
}