#' Obtain critical values and threshold
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#'
#' @return A list containing the critical values and
#' the threshold parameter \eqn{\omega}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
#' @examples
#' # generate data
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#' # construct no_nbd object
#' ts_no_nbd <- ts_hdchange(data_no_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' crit <- get_critical(ts_no_nbd)
#'
get_critical <- function(hdobj) {
  UseMethod("get_critical", hdobj)
}

#' Obtain the test statistics
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#'
#' @return A list containing the test statistics \eqn{\mathcal{\boldsymbol{Q}}_{n}}
#' and a sequence of standardised \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#' @examples
#' # generate data
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#' # construct no_nbd object
#' ts_no_nbd <- ts_hdchange(data_no_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' teststat <- get_teststats(ts_no_nbd)
#'
get_teststats <- function(hdobj) {
  UseMethod("get_teststats", hdobj)
}

#' Obtain the standardised gap vector
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#'
#' @return An array of \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#'
#' @examples
#' # generate data
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#' # construct no_nbd object
#' ts_no_nbd <- ts_hdchange(data_no_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' V_12_MAinf <- get_V_l2_MAinf(ts_no_nbd)
#'
get_V_l2_MAinf<-function(hdobj){
  UseMethod("get_V_l2_MAinf", hdobj)
}


#' Obtain the simulated standardised gap vector
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#'
#' @return An array of the simulated counterpart of \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#'
#' @examples
#' # generate data
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#' # construct no_nbd object
#' ts_no_nbd <- ts_hdchange(data_no_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' GS_MAinf <- get_GS_MAinf(ts_no_nbd)
#'
get_GS_MAinf <- function(hdobj) {
  UseMethod("get_GS_MAinf", hdobj)
}

#' Obtain the time-stamps and spatial locations with breaks
#'
#' @param estobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [est_hdchange()].
#'
#' @return A list containing the time-stamps and spatial locations with breaks.
#' For S3 class 'no_nbd', it returns the total number of breaks \eqn{\widehat{K}}
#' and the time-stamps \eqn{\hat{\tau}_{k}}. See Algorithm 1 of Li et al. (2023).
#' For S3 class 'nbd', it returns the total number of breaks \eqn{\widehat{R}} and the
#' spatial-temporal location of the break \eqn{(\hat{\tau}_{r},\hat{s}_{r})}. See Algorithm
#' 2 of Li et al. (2023).
#'
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#'
#' @examples
#' # generate data
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#' # construct no_nbd object
#' ts_no_nbd <- ts_hdchange(data_no_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' teststats <- get_teststats(ts_no_nbd)
#' V_12_MAinf <- get_V_l2_MAinf(ts_no_nbd)
#'
#' estobj <- est_hdchange(hdobj = ts_no_nbd, test_stats = teststats$stat_max,
#' threshold = 1e-5, stat_all = V_12_MAinf, critical_values = c(0.01, 0.05, 0.1))
#'
#' breaks <- get_breaks(estobj)
#'
get_breaks <- function(estobj) {
  UseMethod("get_breaks", estobj)
}

#' Plot the time series and change-points
#'
#' @param est_result An S3 object of class 'result_no_nbd' or 'result_nbd' created by [get_breaks()].
#' @param ... Additional arguments.
#'
#' @details See [hdchange()] for examples.
#' @return No return value. Presents the plot of the data and breaks.
#' @export
#'
#' @examples
#' \donttest{
#' # generate data
#' data_nbd <- sim_hdchange_nbd(n = 300,
#' p = 70,
#' nbd_info =
#'  list(
#'    (1:9), (2:31), (32:41), (42:70),
#'    (3:15), (16:35), (31:55)
#'  ),
#' sp_tp_break = rbind(c(2, 50), c(4, 150), c(2, 250)),
#' dist_info =
#'   list(dist = "t", dependence = "iid", param = 5),
#' jump_max = 1)
#'
#' # construct nbd object
#' ts_nbd <- ts_hdchange(data_nbd,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1),
#' nbd_info =
#'  list(
#'    (1:9), (2:31), (32:41), (42:70),
#'    (3:15), (16:35), (31:55)
#'  ))
#'
#' # Estimate the time-stamps of the breaks
#' est_result_nbd <- hdchange(ts_nbd)
#'
#' # Plot the results
#' plot_result(est_result_nbd, nbd_index = 2)
#' pairs <- est_result_nbd$nbd_and_stamps_pair
#' time_stamps <- pairs[pairs[, 1] == 2, 2]
#' axis(1,
#'   at = time_stamps,
#'   labels = c("break 1", "break 2")
#' )
#' title(main = "Change-points estimation for neibourhood 2")
#' }
plot_result <- function(est_result,...) {
  UseMethod("plot_result", est_result)
}
