#' Obtain critical values and threshold
#'
#' @inheritParams get_critical
#'
#' @return A list containing the critical values and
#' the threshold parameter \eqn{\omega}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
get_critical.nbd <- function(hdobj) {

  quantiles <- hdobj$quantiles
  alpha <- hdobj$alpha

  GS_MAinf<-get_GS_MAinf(hdobj)

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
get_teststats.nbd <- function(hdobj) {

  V_l2_MAinf<-get_V_l2_MAinf(hdobj)
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
#'
get_V_l2_MAinf.nbd <- function(hdobj) {
  data <- hdobj$data
  n <- hdobj$n
  p <- hdobj$p
  b <- hdobj$b
  nbd_info <- hdobj$nbd_info
  weight <- 1 / (b * n) # weight in local linear estimation

  S_h <- length(nbd_info) # number of spatial neighborhoods
  nbd_size <- lengths(nbd_info) # size of each neighborhood

  var_MAinf <- diag(get_lr_var(hdobj))

  j <- rep(
    unlist(
      sapply(
        c(1:S_h),
        function(r) {
          nbd_info[[r]]
        }
      )
    ), (n - 2 * b * n)
  )
  r <- rep(1:S_h, nbd_size)
  index <- expand.grid(r = r, i = c(1:(n - 2 * b * n)))
  index$j <- j

  V_MAinf_vec <- mapply(V_MAinf_i_j, index$i, index$j,
                        MoreArgs = list(
                          n = n, b = b, data = data,
                          var_MAinf = var_MAinf, weight = weight
                        )
  )

  nbd_size_rep <- rep(nbd_size, (n - 2 * b * n))
  row_take_start <- c(1, cumsum(nbd_size_rep) + 1)
  row_take_start <- row_take_start[-length(row_take_start)]
  row_take_end <- cumsum(nbd_size_rep)
  row_take <- cbind(row_take_start, row_take_end, rep(1:S_h, (n - 2 * b * n)))

  V_l2_MAinf <- apply(row_take, 1, function(index_) {
    rstart <- index_[1]
    rend <- index_[2]
    nbd_index <- index_[3]
    nbd_factorK <- (p / nbd_size[nbd_index])^(1 / 4)
    nbd_factorK * sqrt(sum(V_MAinf_vec[rstart:rend]^2))
  })

  V_l2_MAinf <- matrix(V_l2_MAinf, nrow = S_h, ncol = n - 2 * b * n)

  return(V_l2_MAinf)
}

#' Obtain the simulated standardised gap vector
#'
#' @inheritParams get_GS_MAinf
#'
#' @return An array of the simulated counterpart of \eqn{|V_{i}|_{2}^{2}}.
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
#'
get_GS_MAinf.nbd <- function(hdobj) {
  z_MAinf <- genZ(hdobj)
  n <- hdobj$n
  p <- hdobj$p
  b <- hdobj$b
  N_rep <- hdobj$N_rep
  nbd_info <- hdobj$nbd_info
  S_h <- length(nbd_info)
  nbd_size <- lengths(nbd_info)

  index <- expand.grid(r = c(1:S_h), i = c(1:(n - 2 * b * n)), s = c(1:N_rep))

  GS_MAinf_vec <- mapply(GS_MAinf_s_i_r,
                         index$s, index$i, index$r,
                         MoreArgs = list(
                           n = n, p = p, b = b, nbd_info = nbd_info,
                           nbd_size = nbd_size, z_MAinf = z_MAinf
                         )
  )

  GS_MAinf_mat <- matrix(GS_MAinf_vec, ncol = N_rep)

  GS_MAinf <- apply(GS_MAinf_mat, 2, function(k) {
    sqrt(max(abs(k)))
  })

  return(GS_MAinf)
}

#' Obtain the time-stamps and spatial locations with breaks
#'
#' @inheritParams get_breaks
#'
#' @return A list containing the total number of breaks \eqn{\widehat{R}} and the
#' spatial-temporal location of the break \eqn{(\hat{\tau}_{r},\hat{s}_{r})}. See Algorithm
#' 2 of Li et al. (2023).
#'
#' @references Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
#' \emph{arXiv preprint arXiv:2208.13074}.
#'
#' @export
get_breaks.nbd <- function(estobj) {
  threshold <- estobj$threshold
  n <- estobj$hdobj$n
  b <- estobj$hdobj$b
  p <- estobj$hdobj$p
  Tn_MAinf <- estobj$test_stats
  V_l2_MAinf <- estobj$stat_all
  nbd_info <- estobj$hdobj$nbd_info
  S_h <- length(nbd_info)
  nbd_size <- lengths(nbd_info)

  if (Tn_MAinf > threshold) {
    R0_hat <- 0
    K0_hat <- 0
    S0_hat <- 0
    sp_tp_break_hat <-
      matrix(0, nrow = p * (n - 2 * b * n), ncol = 2)
    size_A1_MAinf <- length(V_l2_MAinf[V_l2_MAinf > threshold])
    A1_MAinf <- cbind(
      V_l2_MAinf[V_l2_MAinf > threshold],
      which(V_l2_MAinf > threshold, arr.ind = TRUE),
      rep(1, size_A1_MAinf)
    )
    while (dim(A1_MAinf)[1] > 0) {
      R0_hat <- R0_hat + 1
      sp_tp_break_hat[R0_hat, 1] <-
        A1_MAinf[which.max(A1_MAinf[, 1]), 2]
      sp_tp_break_hat[R0_hat, 2] <-
        A1_MAinf[which.max(A1_MAinf[, 1]), 3]
      if (dim(A1_MAinf)[1] > 2 * b * n) {
        nbd_overlap <- rep(0, S_h)
        for (l in 1:S_h) {
          nbd_overlap[l] <-
            length(intersect(nbd_info[[sp_tp_break_hat[R0_hat, 1]]], nbd_info[[l]]))
          if (nbd_overlap[l] != 0) {
            A1_MAinf[which(A1_MAinf[, 2] == l &
                             abs(A1_MAinf[, 3] - sp_tp_break_hat[R0_hat, 2]) <= 2 * b * n), 4] <- 0
          }
          nbd_overlap_band <- rep(0, S_h)
          for (ll in 1:S_h) {
            nbd_overlap_band[l] <-
              length(intersect(nbd_info[[l]], nbd_info[[ll]]))
            if (nbd_overlap_band[ll] != 0) {
              A1_MAinf[which(A1_MAinf[, 2] == ll &
                               abs(A1_MAinf[, 3] - sp_tp_break_hat[R0_hat, 2]) <= 2 * b * n), 4] <- 0
            }
          }
        }
        A1_MAinf <- A1_MAinf[which(A1_MAinf[, 4] != 0), ]
      } else {
        A1_MAinf <- c(0, 0, 0)
        break
      }
    }
    sp_tp_break_hat <-
      matrix(sp_tp_break_hat[apply(
        sp_tp_break_hat, 1,
        function(x) {
          !all(x == 0)
        }
      ), ], ncol = 2)
    if (nrow(sp_tp_break_hat) > 0) {
      sp_tp_break_hat[, 2] <- sp_tp_break_hat[, 2] + b * n - 1
      tau_hat <- unique(sp_tp_break_hat[, 2])
      l_hat <- unique(sp_tp_break_hat[, 1])
      K0_hat <- length(unique(sp_tp_break_hat[, 2]))
      S0_hat <- length(unique(sp_tp_break_hat[, 1]))
    } else {
      tau_hat <- sp_tp_break_hat[2]
      l_hat <- sp_tp_break_hat[1]
      K0_hat <- 1
      S0_hat <- 1
    }

    ret_list <- list(
        R0_hat,
        K0_hat,
        S0_hat,
        sp_tp_break_hat,
        tau_hat,
        l_hat,
        Tn_MAinf,
        threshold,
        estobj$critical_values,
        estobj$hdobj
      )


  } else {
    ret_list <- list(
        0,
        0,
        0,
        NA,
        NA,
        NA,
        Tn_MAinf,
        threshold,
        estobj$critical_values,
        estobj$hdobj
      )

  }

  names(ret_list) <- c(
    "total_num_breaks",
    "num_time_stamps",
    "num_of_nbd_break",
    "nbd_and_stamps_pair",
    "time_stamps",
    "nbd_break",
    "test_stats",
    "threshold",
    "critical_values",
    "hdobj"
  )

  class(ret_list) <- "result_nbd"

  return(ret_list)
}

#' Summarize the estimation results
#'
#' @param object An S3 object of class 'result_nbd' created by [get_breaks()].
#' @param ... Additional arguments.
#'
#' @details See [hdchange()] for examples.
#' @return No return value. Presents the summary of the test and estimation results.
#'
#' @export
summary.result_nbd <- function(object, ...) {

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
      "Test statistics (Qn) = ", object$test_stats,
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
    print(object$critical_values)
    cat(
      "===================================================",
      "\n"
    )
    cat("\n")
    cat(
      "================================================================",
      "\n"
    )
    cat(
      "                   l2 change-point estimation                   ",
      "\n"
    )
    cat(
      "================================================================",
      "\n"
    )
    cat(
      "Total number of breaks: ",
      object$total_num_breaks,
      "\n"
    )
    cat(
      "Number of time-stamps: ", object$num_time_stamps,
      "\n"
    )
    cat(
      "Number of neighbourhoods with break: ",
      object$num_of_nbd_break,
      "\n"
    )
    cat("Neighbourhood-time pair: ")
    cat(apply(object$nbd_and_stamps_pair, 1, function(x) {
      paste0("(", x[1], ",", x[2], ")")
    }), "\n")
    cat(
      "Time-stamps: ", object$time_stamps,
      "\n"
    )
    cat(
      "Neighbourhoods with break: ", object$nbd_break,
      "\n"
    )
    cat(
      "================================================================",
      "\n"
    )
    return(invisible(NULL))
}

#' Plot the time series and change-points
#'
#' @param est_result An S3 object of class 'result_nbd' created by [get_breaks()].
#' @param nbd_index An integer indicating which neighbourhood to be plotted.
#' @param ... Additional arguments.
#'
#' @details See [hdchange()] for examples.
#' @return No return value. Presents the plot of the data and breaks.
#'
#' @export
#'
plot_result.result_nbd <- function(est_result, ..., nbd_index) {

  nbd_break <- est_result$nbd_break

  if (!(nbd_index %in% nbd_break)) {
    stop("Breaks not detected in this neighbourhood.")
  } else {
    data <- est_result$hdobj$data
    nbd_info <- est_result$hdobj$nbd_info
    id <- nbd_info[[nbd_index]]
    pairs <- est_result$nbd_and_stamps_pair
    time_stamps <- pairs[pairs[, 1] == nbd_index, 2]
    graphics::matplot(t(data[id, ]),
      xaxt = "n", type = "l", xlab = "Time",
      ylab = " "
    )
    graphics::abline(v = time_stamps, col = "red")
  }
  return(invisible(NULL))
}
