#' Advance a set of agents one movement step at once
#'
#' @description Draws one movement step for every agent supplied, from each agent's own step-length
#'   and turn-angle distributions and selection coefficients. For each agent, \code{n_candidates}
#'   candidate steps are drawn (gamma step lengths and von Mises turn angles), the candidate
#'   endpoints are scored by the habitat and movement covariates of the cell they land in, and one
#'   endpoint is sampled with probability proportional to its selection weight. An agent holds
#'   position, keeping its prior location and heading with a zero-length step, when more than half
#'   of its candidate endpoints fall outside the raster extent or when no candidate carries positive
#'   weight (every endpoint on a no-data cell).
#'
#' @param x Numeric vector of current x coordinates, one element per agent.
#' @param y Numeric vector of current y coordinates.
#' @param heading Numeric vector of current headings (radians).
#' @param prm List of per-agent movement parameters, subset to these agents, with the night
#'   coefficient offsets already added when the hour is night: \code{shape}, \code{scale},
#'   \code{kappa} (numeric vectors), \code{mu} (list), and the six selection coefficient vectors
#'   \code{b_sl}, \code{b_lsl}, \code{b_cta}, \code{b_fo}, \code{b_es}, \code{b_ca}.
#' @param move A list of movement data: \code{geom} (a single-layer \code{terra::SpatRaster}
#'   supplying the grid for \code{cellFromXY}), \code{forage}, \code{escape}, \code{canopy}
#'   (numeric covariate vectors indexed by cell number), \code{ext} (the named grid extent, xmin /
#'   xmax / ymin / ymax), and \code{n_candidates}.
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

  #1) candidate step lengths for every agent in one gamma draw. matrix() fills by column, so repeating each agent's
  #   shape and scale nc times places agent i's draws in row i

  sl <- matrix(
    stats::rgamma(
      n = n * nc,
      shape = rep(prm$shape, times = nc),
      scale = rep(prm$scale, times = nc)
    ),
    nrow = n
  )

  #2) candidate turn angles, drawn one agent at a time because circular::rvonmises() accepts only a single kappa,
  #   then wrapped to (-pi, pi]

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

  #3) candidate endpoints: each turn angle is added to the agent's current heading to give a bearing, and the step
  #   length is laid out along it

  bearing <- heading + ta
  x2 <- x + sl * cos(bearing)
  y2 <- y + sl * sin(bearing)

  # ----------------------------------------------------------------------------------------------------------------------
  # guards and weights
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the fraction of each agent's candidate endpoints that fall outside the grid extent

  frac_outside <- rowMeans(
    x2 < ext["xmin"] | x2 > ext["xmax"] | y2 < ext["ymin"] | y2 > ext["ymax"]
  )

  #2) the forage, escape-terrain, and canopy values of the cell each endpoint lands in. endpoints off the raster
  #   return NA, which becomes a zero weight below

  cells <- terra::cellFromXY(move$geom, cbind(as.vector(x2), as.vector(y2)))
  fo <- matrix(move$forage[cells], nrow = n)
  es <- matrix(move$escape[cells], nrow = n)
  ca <- matrix(move$canopy[cells], nrow = n)

  #3) the selection linear predictor for every candidate, with each agent's coefficients recycled along its own
  #   row. each row is centered on its own finite mean before exponentiation, and non-finite weights (from NA
  #   covariates or a zero step length) are set to zero so those candidates are never sampled

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

  #4) the two hold-position conditions: more than half of the candidates outside the extent, or no candidate with
  #   positive weight

  row_w <- rowSums(w)
  hold <- frac_outside > 0.5 | row_w <= 0

  # ----------------------------------------------------------------------------------------------------------------------
  # select one candidate per agent and assemble the step
  # ----------------------------------------------------------------------------------------------------------------------

  #1) sample one candidate per agent by inverse cumulative distribution across its row of weights: a uniform draw
  #   scaled to the row total is compared with the running sum of weights, and the first candidate whose running
  #   sum reaches it is selected

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

  #2) agents that move take the selected candidate's endpoint, step length, turn angle, and bearing as their new
  #   heading; held agents keep their position and heading with a zero-length step and zero turn angle

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
