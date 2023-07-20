#' Obtain critical values and threshold
#' @inheritParams get_critical
#'
#' @return A list containing the critical values and
#' the threshold parameter \eqn{\omega}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
get_critical.no_nbd <- function(hdobj) {
  n <- hdobj$n
  p <- hdobj$p
  b <- hdobj$b
  quantiles <- hdobj$quantiles
  alpha <- hdobj$alpha
  N_rep <- hdobj$N_rep

  GS_MAinf <- get_GS_MAinf(hdobj)
  critical_values <- c(stats::quantile(GS_MAinf, 1 - quantiles))
  critical_value_alpha <- unlist(stats::quantile(GS_MAinf, 1 - alpha))

  ret_list <- list(critical_values, critical_value_alpha)
  names(ret_list) <- c("critical_values", "critical_value_alpha")

  return(ret_list)
}

#' Obtain the test statistics
#'
#' @inheritParams get_teststats
#'
#' @return A list containing the test statistics \eqn{\mathcal{\boldsymbol{Q}}_{n}}
#' and a sequence of standardised \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
get_teststats.no_nbd <- function(hdobj) {
  V_l2_MAinf <- get_V_l2_MAinf(hdobj)
  Tn_MAinf <- max(V_l2_MAinf)
  ret_list <- list(Tn_MAinf, V_l2_MAinf)
  names(ret_list) <- c("stat_max", "stat_all")

  return(ret_list)
}

#' Obtain the standardised gap vector
#'
#' @inheritParams get_V_l2_MAinf
#'
#' @return An array of \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
get_V_l2_MAinf.no_nbd <- function(hdobj) {
  data <- hdobj$data
  n <- hdobj$n
  p <- hdobj$p
  b <- hdobj$b
  weight <- 1 / (b * n)

  index <- expand.grid(j = c(1:p), i = (1:(n - 2 * b * n)))

  var_MAinf <- diag(get_lr_var(hdobj))

  V_MAinf_vec <- mapply(V_MAinf_i_j,
                        i = index$i,
                        j = index$j,
                        MoreArgs = list(n = n, p = p, b = b, data = data,
                                        var_MAinf = var_MAinf, weight)
  )

  V_MAinf_vec_mat <- matrix(V_MAinf_vec, nrow = p, ncol = n - 2 * b * n)

  V_l2_MAinf <- apply(V_MAinf_vec_mat, 2, function(v) sqrt(sum(v^2)))

  return(V_l2_MAinf)
}

#' Obtain the simulated standardised gap vector
#'
#' @inheritParams get_GS_MAinf
#'
#' @return An array of the simulated counterpart of \eqn{|V_{i}|_{2}^{2}}.
#'
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
get_GS_MAinf.no_nbd <- function(hdobj) {
  z_MAinf <- genZ(hdobj)
  n <- hdobj$n
  p <- hdobj$p
  b <- hdobj$b
  N_rep <- hdobj$N_rep

  index <- expand.grid(i = c(1:(n - 2 * b * n)), s = c(1:N_rep))

  GS_MAinf_vec <- mapply(GS_MAinf_s_i,
                         s = index$s, i = index$i,
                         MoreArgs = list(n = n, p = p, b = b, z_MAinf = z_MAinf)
  )

  GS_MAinf_mat <- matrix(GS_MAinf_vec, nrow = (n - 2 * b * n), ncol = N_rep)

  GS_MAinf <- apply(GS_MAinf_mat, 2, function(k) {
    sqrt(max(abs(k)))
  })

  return(GS_MAinf)
}

#' Obtain the time-stamps and spatial locations without break
#'
#' @inheritParams get_breaks
#'
#' @return A list containing the total number of breaks \eqn{\widehat{K}}
#' and the time-stamps \eqn{\hat{\tau}_{k}}. See Algorithm 1 of Li et al. (2023).
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
get_breaks.no_nbd <- function(estobj) {
  threshold <- estobj$threshold
  n <- estobj$hdobj$n
  b <- estobj$hdobj$b
  Tn_MAinf <- estobj$test_stats
  V_l2_MAinf <- estobj$stat_all

  if (Tn_MAinf > threshold) {
    K0_hat <- 0
    tau_hat <- rep(0, n - 2 * b * n)
    A1_MAinf <- cbind(
      subset(V_l2_MAinf, V_l2_MAinf > threshold),
      which(V_l2_MAinf > threshold)
    )
    while (length(A1_MAinf) > 0) {
      K0_hat <- K0_hat + 1
      if (length(A1_MAinf) > 2) {
        tau_hat[K0_hat] <- A1_MAinf[which.max(A1_MAinf[, 1]), 2]
        if (max(abs(A1_MAinf[, 2]) - tau_hat[K0_hat]) > 2 * b * n) {
          A1_MAinf <- A1_MAinf[abs(A1_MAinf[, 2] - tau_hat[K0_hat]) > b * n, ]
        } else {
          A1_MAinf <- c(0, 0)
          break
        }
      } else {
        tau_hat[K0_hat] <- A1_MAinf[2]
        break
      }
    }
    tau_hat <- tau_hat[tau_hat != 0] + b * n - 1

    ret_list <- list(
        K0_hat,
        tau_hat,
        Tn_MAinf,
        threshold,
        estobj$critical_values,
        estobj$hdobj
      )

  } else {
    ret_list <- list(
        0,
        NA,
        Tn_MAinf,
        threshold,
        estobj$critical_values,
        estobj$hdobj
      )

  }

  names(ret_list) <- c(
    "numbers",
    "time_stamps",
    "test_stats",
    "threshold",
    "critical_values",
    "hdobj"
  )
  class(ret_list) <- "result_no_nbd"

  return(ret_list)
}

#' Summarize the estimation results
#'
#' @param object An S3 object of class 'result_no_nbd' created by [get_breaks()].
#' @param ... Additional arguments.
#'
#' @details See [hdchange()] for examples.
#' @return No return value. Presents the summary of the test and estimation results.
#' @export
#'
summary.result_no_nbd <- function(object,...) {


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
    cat("H0: no breaks ", "\n")
    cat("Ha: at least one break", "\n")
    cat(
      "Test statistics (Qn) = ", object$test_stats,
      "\n"
    )
    cat(
      "===================================================",
      "\n"
    )
    cat("Critical values : ", "\n")
    print(object$critical_values)
    cat(
      "===================================================",
      "\n"
    )
    cat("\n")
    cat(
      "====================================================",
      "\n"
    )
    cat(
      "             l2 change-point estimation             ",
      "\n"
    )
    cat(
      "====================================================",
      "\n"
    )
    cat(
      "Number of breaks: ", object$numbers,
      "\n"
    )
    cat(
      "Time-stamps: ", object$time_stamps,
      "\n"
    )
    cat(
      "====================================================",
      "\n"
    )
    return(invisible(NULL))
}

#' Plot the time series and change-points
#'
#' @param est_result An S3 object of class 'result_no_nbd' created by [get_breaks()].
#' @param ... Additional arguments.
#'
#' @details See [hdchange()] for examples.
#' @return No return value. Presents the plot of the data and breaks.
#'
#' @export
#'
plot_result.result_no_nbd <- function(est_result,...) {
  data <- est_result$hdobj$data
  time_stamps <- est_result$time_stamps
  graphics::matplot(t(data), type = "l", pch = 1, xaxt = "n", xlab = "Time", ylab = " ")
  graphics::abline(v = time_stamps, col = "red")
  return(invisible(NULL))
}
