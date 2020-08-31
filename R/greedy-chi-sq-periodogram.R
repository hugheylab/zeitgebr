#' @param time_resolution the resolution of periods to scan
#' @rdname periodogram_methods
#' @export
greedy_chi_sq_periodogram <- function(x,
                               period_range = c(hours(16), hours(32)),
                               sampling_rate = 1 / mins(1),
                               alpha = 0.05,
                               time_resolution = hours(0.1)){
  # lsp can handle time series, so lets use that feature!
  p_value = power = signif_threshold = period = .N = signif_threshold = NULL
  out <- data.table::data.table(
    period = seq(period_range[1],period_range[2],by=time_resolution)
  )
  corrected_alpha <- 1 - out[,(1 - alpha) ^ (1 / .N)]
  out[ , power := sapply(period, calc_greedy_Qp, x, sampling_rate)]
  out[ , signif_threshold := qchisq(corrected_alpha , round(period * sampling_rate), lower.tail = FALSE)]
  out[ , p_value := pchisq(power, round(period * sampling_rate), lower.tail = FALSE)]


}


#' Calculate Qp
#'
#' @param values activity values (each value represents the measured activity in a minute)
#' @param target_period a period at which the chi-squared statistics is to be calculated
#'
#' @return a numeric of the calculated chi-squared statistics at the given varPer
#' @noRd
calc_greedy_Qp <- function(target_period, values, sampling_rate){
  # . = NULL
  # col_num <- round(target_period * sampling_rate)
  # #row_num <- ceiling(length(values) / col_num)
  # row_num <- length(values) / col_num
  # dt <- data.table::data.table( col = (0 : (length(values) -1) %% col_num) + 1,
  #                               #row = ceiling(1:length(values) / col_num),
  #                               values = values,
  #                               key = "col")
  # avg_P <- dt[, .(avg_P = mean(values)),by=col]$avg_P
  # avg_all <- mean(values)
  # numerator <- sum((avg_P - avg_all) ^ 2) *  (nrow(dt) * row_num)
  #
  # denom <- sum((values - avg_all) ^ 2)
  # numerator / denom

  m = mean(values)
  n = length(values)

  k = n / (target_period * sampling_rate)
  xNow = c(x, rep.int(NA, ceiling(k) * (target_period * sampling_rate) - n))
  xMat = matrix(xNow, ncol = (target_period * sampling_rate), byrow = TRUE)
  xH = colMeans(xMat, na.rm = TRUE)
  qP = k * n * sum((xH - m)^2) / sum((x - m)^2)
}
