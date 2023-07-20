#' 'no_nbd' or 'nbd' object construction
#'
#' @description
#' This function creates an S3 object of class 'no_nbd' or 'nbd' containing the
#' initialising information supplied to the main function [hdchange()].
#' 'no_nbd' or 'nbd' are constructed depending on whether the
#' neighbourhood information is provided. The resulting object will be used in the
#' test and estimation functions.
#'
#' @param data p by n data matrix, n = number of time series observations,
#' p = cross-sectional dimension.
#' @param window_size \eqn{window\_size = b \times n}, e.g. \eqn{n=100}, \eqn{b=30}.
#' @param m Number of blocks in long-run variance estimation, 8 by default.
#' @param h Parameter in long-run variance estimation, 1 by default.
#' @param N_rep Number of repetitions in MC simulation.
#' @param alpha A small positive number controlling for the threshold in break estimation.
#' @param quantiles An array of quantiles for critical values.
#' @param nbd_info A list containing the neighbourhood information, NULL by default
#' indicating no neighbourhoods.
#'
#' @details
#' 'nbd_info' indicates the location of individuals in the data matrix.
#' For example, 'nbd_info = list(c(1:10), c(25:35), c(7:18))' means that
#' there are three neighbourhoods. The first neighbourhood contains from the 1st
#' to 10th individuals and the same rule applies to the rest of neighbourhoods.
#' The neighbourhoods are allowed to be overlapped. See also the illustrating
#' example in [hdchange()].
#'
#' @return The return value is an S3 object of class 'no_nbd' or 'nbd'.
#' It contains a list of the following items:
#' \itemize{
#'
#' \item data, m, h, N_rep, alpha, quantiles, and nbd_info are the same as in
#' the arguments.
#'
#' \item n = number of time series observations.
#'
#' \item p = cross-sectional dimension.
#'
#' \item b = bandwith parameter \eqn{b = window\_size/n}.
#'
#' }
#'
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#'
#' @examples
#' data <- covid_data
#'
#' # No neighbourhood case
#' ts_no_nbd <- ts_hdchange(data,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1))
#'
#' # Neighbourhood case
#' ts_nbd <- ts_hdchange(data,
#' window_size = 30,
#' m = 8,
#' h = 1,
#' N_rep = 999,
#' alpha = 1e-5,
#' quantiles = c(0.01, 0.05, 0.1),
#' nbd_info = list(c(1:10), c(25:35), c(7:18)))
ts_hdchange <- function(data,
                        window_size = 30,
                        m = 8,
                        h = 1,
                        N_rep = 999,
                        alpha = 1e-5,
                        quantiles = c(0.01, 0.05, 0.1),
                        nbd_info = NULL) {
  n <- dim(data)[2] # number of observations
  p <- dim(data)[1] # cross-sectional dimension
  b <- window_size / n # bandwidth parameter

  if (length(alpha) > 1) {
    alpha <- alpha[1]
  }

  if (is.null(nbd_info)) {
    ts_hdchange_init <- structure(
      list(
        data,
        window_size,
        n,
        p,
        b,
        m,
        h,
        N_rep,
        alpha,
        quantiles,
        nbd_info
      ),
      class = c("no_nbd")
    )
  } else {
    check_nbd <- check_nbd(nbd_info)
    ts_hdchange_init <- structure(
      list(
        data,
        window_size,
        n,
        p,
        b,
        m,
        h,
        N_rep,
        alpha,
        quantiles,
        nbd_info
      ),
      class = c("nbd")
    )
  }


  names(ts_hdchange_init) <-
    c(
      "data",
      "window_size",
      "n",
      "p",
      "b",
      "m",
      "h",
      "N_rep",
      "alpha",
      "quantiles",
      "nbd_info"
    )

  return(ts_hdchange_init)
}

#' Test the existence of change-points in the data
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#' @param display A logical. If 'display = TRUE', the test statistics and
#' critical values will be printed.
#'
#' @details See [hdchange()] for examples.
#' @return A list containing the following elements:
#' \itemize{
#'
#' \item 'test_stats' The test statistics
#' \eqn{\mathcal{\boldsymbol{Q}}_{n}}.
#'
#' \item 'critical_values' The critical values.
#'
#' \item 'stat_all' An array of \eqn{|V_{i}|_{2}^{2}}.
#'
#' \item 'critical_value_alpha' The threshold value \eqn{\omega} depending on alpha.
#'
#' }
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
#' test <- test_existence(ts_no_nbd, display = TRUE)
#'
test_existence <- function(hdobj, display = TRUE) {
  #--------Test the existence of breaks---------#

  critical_vals <- get_critical(hdobj)

  critical_values <- critical_vals$critical_values

  critical_value_alpha <- critical_vals$critical_value_alpha

  get_stats <- get_teststats(hdobj)

  Tn_MAinf <- get_stats$stat_max

  V_l2_MAinf <- get_stats$stat_all

  if (display) {
    cat(
      "===================================================",
      "\n"
    )
    cat(
      "             l2 change-point detection             ",
      "\n"
    )
    cat(
      "===================================================",
      "\n"
    )
    cat(
      "H0: no breaks ",
      "\n"
    )
    cat(
      "Ha: at least one break",
      "\n"
    )
    cat(
      "Test statistics (Qn) = ", Tn_MAinf,
      "\n"
    )
    cat(
      "===================================================",
      "\n"
    )
    cat(
      "Critical values : ",
      "\n"
    )
    print(critical_values)
    cat(
      "===================================================",
      "\n"
    )
    return(invisible(NULL))
  }else{

    ret_list <-
      list(Tn_MAinf, critical_values, V_l2_MAinf, critical_value_alpha)
    names(ret_list) <- c(
      "test_stats",
      "critical_values",
      "stat_all",
      "critical_value_alpha"
    )
    return(ret_list)
  }


}



#' Estimate the time-stamps and spatial locations with breaks
#' @description
#' The main function of this package. It performs a test for existence of
#' breaks and estimates the time-stamps and locations of the breaks.
#'
#' @param hdobj An S3 object of class 'no_nbd' or 'nbd' generated
#' by [ts_hdchange()].
#'
#'
#' @return
#' The return value is an S3 object of class 'no_nbd' or 'nbd' containing a list
#' of the test results and change-point locations.
#'
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#' @export
#'
#' @examples
#' ############ No neighbourhood case ############
#'
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
#' \donttest{
#'
#' # Estimate the time-stamps of the breaks
#' est_result_no_nbd <- hdchange(ts_no_nbd)
#'
#' # Summarize the results
#' summary(est_result_no_nbd)
#'
#' # Plot the results
#' plot_result(est_result_no_nbd)
#' axis(1,
#'   at = est_result_no_nbd$time_stamps,
#'   labels = c("break 1", "break 2", "break 3")
#' )
#' title(main = "Change-points estimation")
#' }
#'
#' ############ Neighbourhood case ############
#'
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
#' \donttest{
#'
#' # Estimate the time-stamps of the breaks
#' est_result_nbd <- hdchange(ts_nbd)
#'
#' # Summarize the results
#' summary(est_result_nbd)
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
hdchange <- function(hdobj) {
  exist <- test_existence(hdobj, display = FALSE)

  estobj <- est_hdchange(
    hdobj,
    test_stats = exist$test_stats,
    threshold = exist$critical_value_alpha,
    stat_all = exist$stat_all,
    critical_values = exist$critical_values
  )

  est <- get_breaks(estobj)

  return(est)
}

