#' Simulate data with neighbourhood
#'
#' @param n Number of time series observations.
#' @param p Number of individual.
#' @param nbd_info A list containing the neighbourhood information. See [ts_hdchange()].
#' @param sp_tp_break A \eqn{K \times 2} matrix indicating the spatial-temporal break location.
#' @param dist_info A list specifying the distribution of the innovation.
#' @param jump_max Maximum jump size of the breaks.
#'
#' @details
#' 'sp_tp_break' should be a \eqn{K \times 2} matrix with first column
#' indicating the neighbourhoods and the second column indicating the time stamps.
#' For example, 'sp_tp_break = rbind(c(2, 50), c(4, 150), c(2, 250))' means that
#' the second neighbourhood has two breaks taking place at \eqn{i = 50, 250} and
#' the fourth neighbourhood has one break taking place at \eqn{i = 150}.
#'
#' 'dist_info' should be a list containing the following items:
#'  \itemize{
#'
#' \item dist: distribution of the innovations, either "normal" or "t".
#'
#' \item dependence: iid or \eqn{MA(\infty)}, either "iid" or "MA_inf".
#'
#' \item param = parameter of the distribution, standard deviation for normal distribution
#' and degree of freedom for t distribution
#'
#' }
#'
#' 'jump_max' is set equal in nbd case for convenience.
#'
#' See [ts_hdchange()] for example.
#'
#' @return A \eqn{p \times n} simulated data matrix.
#' @export
#'
#' @examples
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
#'
sim_hdchange_nbd <- function(n = 300,
                             p = 70,
                             nbd_info =
                               list(
                                 (1:9), (2:31), (32:41), (42:70),
                                 (3:15), (16:35), (31:55)
                               ),
                             sp_tp_break = rbind(c(2, 50), c(4, 150), c(2, 250)),
                             dist_info =
                               list(dist = "normal", dependence = "iid", param = 1),
                             jump_max = 1) {
  S_h <- length(nbd_info) # number of spatial neighborhoods
  nbd_size <- lengths(nbd_info) # size of each neighborhood
  ave_nbd_size <- mean(nbd_size)

  R0 <- nrow(sp_tp_break) # number of spatial-temporal windows with break
  K0 <- length(unique(sp_tp_break[, 2])) # number of true time points with breaks (any neighborhood)
  S0 <- length(unique(sp_tp_break[, 1])) # number of true spatial neighborhoods with breaks
  tau <- sort(unique(sp_tp_break[, 2])) # break time stamps, length(tau)=K0, min(tau)>bn & max(tau)<n-bn
  l <- sort(unique(sp_tp_break[, 1])) # break spatial locations, length(l)=S0, distance>1
  S <- nbd_size[l] # number of sectional dimensions with jumps, smaller than the corresponding nbd size, length(S)=S0,


  #------------------------------------------------------------------------------------------------

  # 1. Generate observations Y = trend + epsilon
  beta <- 2
  q <- n * 2
  scale_a <- seq(0.5, 0.9, (0.9 - 0.5) / max(1, (p - 1)))
  a <- scale_a %*% t((1:q)^(-beta))
  # var_MAinf <- (rowSums(a))^2 * (nu / (nu - 2))

  # 1.1 jumps in the trends

  jump_size <- matrix(0, nrow = p, ncol = K0) # jump sizes of all p coordinates
  break_size_l2 <- rep(0, K0) # break sizes of all p coordinates, i.e. l2-norm of jump size
  trend <- matrix(0, nrow = p, ncol = n) # different trends due to breaks

  for (k in 1:K0) {
    for (h in 1:S0) {
      if (any(apply(sp_tp_break, 1, function(x, want) isTRUE(all.equal(x, want)), c(l[h], tau[k])))) {
        jump_size[nbd_info[[l[h]]], k] <- c(rep(jump_max, S[h]), rep(0, nbd_size[l[h]] - S[h]))
      }
    }
    break_size_l2[k] <- norm(jump_size[, k], type = "2")
    if (K0 == 1) {
      trend[, (tau + 1):n] <- rep(jump_size[, k], n - tau)
    } else {
      if (k == 1) {
        trend[, (tau[k] + 1):tau[k + 1]] <- rep(jump_size[, k], tau[k + 1] - tau[k])
      } else if (k == K0) {
        trend[, (tau[k] + 1):n] <- matrix(
          rep(trend[, tau[k]], n - tau[K0]),
          p, n - tau[K0]
        ) + jump_size[, k]
      } else {
        trend[, (tau[k] + 1):tau[k + 1]] <- rep(trend[, tau[k]], tau[k + 1] - tau[k]) + jump_size[, k]
      }
    }
  }

  # 1.2 innovations

  epsilon <- matrix(0, nrow = p, ncol = n)
  distr <- dist_info$dist
  depend <- dist_info$dependence
  if (distr == "t" & depend == "MA_inf") {
    df <- dist_info[[3]]
    for (j in 1:p) {
      eta_j <- stats::rt(q, df = df)
      epsilon[j, ] <- (Re(stats::fft(stats::fft(a[j, ]) * stats::fft(eta_j), inverse = TRUE)) / q)[(q - n + 1):q]
    }
  } else if (distr == "t" & depend == "iid") {
    df <- dist_info[[3]]
    for (j in 1:p) {
      epsilon[j, ] <- stats::rt(n, df = df)
    }
  } else if (distr == "normal" & depend == "MA_inf") {
    sd <- dist_info[[3]]

    for (j in 1:p) {
      eta_j <- stats::rnorm(q, mean = 0, sd = sd)
      epsilon[j, ] <- (Re(stats::fft(stats::fft(a[j, ]) * stats::fft(eta_j), inverse = TRUE)) / q)[(q - n + 1):q]
    }
  } else if (distr == "normal" & depend == "iid") {
    sd <- dist_info[[3]]
    for (j in 1:p) {
      epsilon[j, ] <- stats::rnorm(n, mean = 0, sd = sd)
    }
  } else {
    stop("distribution info not recognised.")
  }

  data <- trend + epsilon
}

#' Simulate data without neighbourhood
#'
#' @param n Number of time series observations.
#' @param p Number of individuals.
#' @param S Number of individuals with jumps.
#' @param tau An array of length \eqn{K} for time stamps for breaks.
#' @param dist_info A list specifying the distribution of the innovation.
#' @param jump_max An array of length \eqn{K} for jump sizes of the breaks.
#'
#' @details
#' 'dist_info' should be a list containing the following items:
#'  \itemize{
#'
#' \item dist: distribution of the innovations, either "normal" or "t".
#'
#' \item dependence: iid or MA(\eqn{\infty}), either "iid" or "MA_inf".
#'
#' \item param = parameter of the distribution, standard deviation for normal distribution
#' and degree of freedom for t distribution
#'
#' }
#' See [ts_hdchange()] for example.
#'
#' @return A \eqn{p \times n} simulated data matrix.
#' @export
#'
#' @examples
#' data_no_nbd <- sim_hdchange_no_nbd(n = 200,
#' p = 30,
#' S = 30,
#' tau = c(40, 100, 160),
#' dist_info =
#'   list(dist = "normal", dependence = "MA_inf", param = 1),
#' jump_max = c(2, 2, 1.5))
#'
#'
sim_hdchange_no_nbd <- function(n = 200,
                                p = 30,
                                S = 30,
                                tau = c(40, 100, 160),
                                dist_info =
                                  list(dist = "normal", dependence = "iid", param = 1),
                                jump_max = c(2, 2, 1.5)) {
  # 1. Generate observations Y = trend + epsilon

  # 1.1 jumps in the trends

  K0 <- length(tau)
  jump_size <- matrix(0, nrow = p, ncol = K0)
  break_size_l2 <- rep(0, K0)
  trend <- matrix(0, nrow = p, ncol = n)
  q <- n * 2

  beta <- 2
  scale_a <- seq(0.5, 0.9, (0.9 - 0.5) / max(1, (p - 1)))
  a <- scale_a %*% t((1:q)^(-beta))

  for (k in 1:K0) {
    jump_size[, k] <- c(rep(jump_max[k], S), rep(0, p - S))
    break_size_l2[k] <- norm(jump_size[, k], type = "2")
    if (K0 == 1) {
      trend[, (tau + 1):n] <- rep(jump_size[, k], n - tau)
    } else {
      if (k == 1) {
        trend[, (tau[k] + 1):tau[k + 1]] <- rep(jump_size[, k], tau[k + 1] - tau[k])
      } else if (k == K0) {
        trend[, (tau[k] + 1):n] <- matrix(
          rep(trend[, tau[k]], n - tau[K0]),
          p, n - tau[K0]
        ) + jump_size[, k]
      } else {
        trend[, (tau[k] + 1):tau[k + 1]] <- rep(trend[, tau[k]], tau[k + 1] - tau[k]) + jump_size[, k]
      }
    }
  }

  # 1.2 innovations

  epsilon <- matrix(0, nrow = p, ncol = n)
  distr <- dist_info$dist
  depend <- dist_info$dependence
  if (distr == "t" & depend == "MA_inf") {
    df <- dist_info[[3]]
    for (j in 1:p) {
      eta_j <- stats::rt(q, df = df)
      epsilon[j, ] <- (Re(stats::fft(stats::fft(a[j, ]) * stats::fft(eta_j), inverse = TRUE)) / q)[(q - n + 1):q]
    }
  } else if (distr == "t" & depend == "iid") {
    df <- dist_info[[3]]
    for (j in 1:p) {
      epsilon[j, ] <- stats::rt(n, df = df)
    }
  } else if (distr == "normal" & depend == "MA_inf") {
    sd <- dist_info[[3]]

    for (j in 1:p) {
      eta_j <- stats::rnorm(q, mean = 0, sd = sd)
      epsilon[j, ] <- (Re(stats::fft(stats::fft(a[j, ]) * stats::fft(eta_j), inverse = TRUE)) / q)[(q - n + 1):q]
    }
  } else if (distr == "normal" & depend == "iid") {
    sd <- dist_info[[3]]
    for (j in 1:p) {
      epsilon[j, ] <- stats::rnorm(n, mean = 0, sd = sd)
    }
  } else {
    stop("distribution info not recognised.")
  }

  data <- trend + epsilon
}
