# vec_move.R - vectorized hourly move step. A single-step twin of simulate_burn_in's inner step
#   that additionally returns the realized step length and turn angle, so the season loop can write
#   the same five movement columns per agent that simulate_move returns. Operates on the living
#   subset handed to it; held agents (cornered or zero-weight) keep position and heading with a
#   zero-length step, matching simulate_move's guards.

#' Advance a set of agents one movement step at once
#'
#' @param x Numeric vector of current x coordinates, one element per agent.
#' @param y Numeric vector of current y coordinates.
#' @param heading Numeric vector of current headings (radians).
#' @param prm List of per-agent movement parameters already subset to these agents and with the
#'   night coefficient offsets already folded in when the hour is night: \code{shape},
#'   \code{scale}, \code{kappa} (numeric vectors), \code{mu} (list), and the six working
#'   coefficient vectors \code{b_sl}, \code{b_lsl}, \code{b_cta}, \code{b_fo}, \code{b_es},
#'   \code{b_ca}.
#' @param move The movement data list described in \code{simulate_move}.
#'
#' @return A numeric matrix with one row per agent and columns \code{x}, \code{y},
#'   \code{step_length}, \code{turn_angle}, \code{heading}.
#'
#' @importFrom stats rgamma runif
#' @keywords internal
simulate_move_step_vec <- function(x, y, heading, prm, move) {

  n <- length(x)
  nc <- move$n_candidates
  ext <- move$ext

  # ----------------------------------------------------------------------------------------------------------------------
  # candidate draws
  # ----------------------------------------------------------------------------------------------------------------------

  #1) step lengths for every agent in one call; matrix() fills by column, so repeating the
  #   per-agent vector places agent i's draws in row i

  sl <- matrix(
    stats::rgamma(
      n = n * nc,
      shape = rep(prm$shape, times = nc),
      scale = rep(prm$scale, times = nc)
    ),
    nrow = n
  )

  #2) turn angles one agent at a time (circular::rvonmises accepts only a scalar kappa), wrapped to
  #   (-pi, pi] as amt wraps them

  ta <- matrix(0, nrow = n, ncol = nc)

  for (i in seq_len(n)) {
    angles <- suppressWarnings(
      as.numeric(
        circular::rvonmises(
          n = nc,
          mu = prm$mu[[i]],
          kappa = prm$kappa[i]
        )
      )
    )
    angles <- angles %% (2 * pi)
    ta[i, ] <- ifelse(angles > pi, angles - (2 * pi), angles)
  }

  #3) candidate endpoints from each agent's current position and heading

  bearing <- heading + ta
  x2 <- x + sl * cos(bearing)
  y2 <- y + sl * sin(bearing)

  # ----------------------------------------------------------------------------------------------------------------------
  # guards and weights
  # ----------------------------------------------------------------------------------------------------------------------

  #1) outside-extent guard per agent

  frac_outside <- rowMeans(
    x2 < ext["xmin"] | x2 > ext["xmax"] | y2 < ext["ymin"] | y2 > ext["ymax"]
  )

  #2) end-of-step covariates by cell; off-raster candidates return NA and fall out as zero weight

  cells <- terra::cellFromXY(move$geom, cbind(as.vector(x2), as.vector(y2)))
  fo <- matrix(move$forage[cells], nrow = n)
  es <- matrix(move$escape[cells], nrow = n)
  ca <- matrix(move$canopy[cells], nrow = n)

  #3) linear predictor with each agent's coefficients recycled down its own row, then per-row
  #   centering and exponentiation exactly as simulate_move does per call

  w <- sl * prm$b_sl +
    log(sl) * prm$b_lsl +
    cos(ta) * prm$b_cta +
    fo * prm$b_fo +
    es * prm$b_es +
    ca * prm$b_ca

  finite <- is.finite(w)
  w_bar <- rowSums(ifelse(finite, w, 0)) / rowSums(finite)

  w <- exp(w - w_bar)
  w[!is.finite(w)] <- 0

  #4) the two hold-position conditions

  row_w <- rowSums(w)
  hold <- frac_outside > 0.5 | row_w <= 0

  # ----------------------------------------------------------------------------------------------------------------------
  # select one candidate per agent and assemble the step
  # ----------------------------------------------------------------------------------------------------------------------

  #1) inverse-CDF sample across each agent's own row of weights

  u <- stats::runif(n) * row_w
  acc <- numeric(n)
  sel <- rep.int(nc, n)
  found <- logical(n)

  for (j in seq_len(nc)) {
    acc <- acc + w[, j]
    newly <- !found & acc >= u
    sel[newly] <- j
    found <- found | newly
  }

  pick <- cbind(seq_len(n), sel)

  #2) moved agents take the selected candidate (step length and turn angle are the drawn values
  #   that generated the endpoint); held agents keep position and heading with a zero-length step

  out_x <- ifelse(hold, x, x2[pick])
  out_y <- ifelse(hold, y, y2[pick])
  out_sl <- ifelse(hold, 0, sl[pick])
  out_ta <- ifelse(hold, 0, ta[pick])
  out_heading <- ifelse(hold, heading, bearing[pick])

  cbind(
    x = out_x,
    y = out_y,
    step_length = out_sl,
    turn_angle = out_ta,
    heading = out_heading
  )
}
