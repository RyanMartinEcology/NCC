#' Validate the iSSF coefficient structure
#'
#' @description The movement functions (\code{simulate_move_step_vec} and
#'   \code{simulate_burn_in}) compute the selection linear predictor by name
#'   rather than through \code{model.matrix}, so they depend on a fixed term
#'   structure: six movement and habitat covariates, each with a day main effect
#'   and a \code{:tod_end_night_end} night offset. This check confirms that an iSSF
#'   carries exactly those twelve named coefficients and stops the run otherwise,
#'   so a change to the term structure fails at the start of a run rather than
#'   being miscomputed. Only the names are checked, not the values. The names are
#'   shared across agents, so the caller checks one agent once per run.
#'
#' @param coefs The named coefficient vector \code{issf$coefficients}.
#'
#' @return Invisibly \code{TRUE}; called for its side effect (error on mismatch).
#'
#' @keywords internal
validate_move_coefs <- function(coefs) {

  # ----------------------------------------------------------------------------------------------------------------------
  # check the coefficient names
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the twelve expected names: step length, log step length, cosine of turn angle, forage biomass, distance to
  #   escape terrain, and canopy cover, each as a day main effect and as a night offset

  expected <- c(
    "sl_", "log(sl_)", "cos(ta_)",
    "forage_biomass_end", "escape_terrain_end", "canopy_cover_end",
    "sl_:tod_end_night_end", "log(sl_):tod_end_night_end", "cos(ta_):tod_end_night_end",
    "forage_biomass_end:tod_end_night_end", "escape_terrain_end:tod_end_night_end",
    "canopy_cover_end:tod_end_night_end"
  )

  #2) stop the run if the supplied names are not exactly this set, in any order

  if (!setequal(names(coefs), expected)) {
    stop(
      "the movement functions expect an iSSF with exactly these coefficients:\n  ",
      paste(expected, collapse = ", "),
      "\nThe selection predictor must be updated if the term ",
      "structure changes."
    )
  }
  invisible(TRUE)
}

#' Simulate the movement burn-in for all agents at once
#'
#' @description Moves every agent forward \code{n_steps} movement-only steps under its own iSSF and
#'   returns only the endpoint. Each step draws \code{n_candidates} candidate steps per agent
#'   (gamma step lengths and von Mises turn angles), scores the candidate endpoints by the habitat
#'   and movement covariates of the cell they land in, and samples one endpoint per agent with
#'   probability proportional to its selection weight. An agent holds position for a step when more
#'   than half of its candidate endpoints fall outside the raster extent or when no candidate
#'   carries positive weight. All agents advance together one step at a time.
#'
#'   Step lengths and turn angles do not depend on position, so they are drawn ahead of the step
#'   loop in blocks of \code{block} steps: one \code{rgamma} call covering every agent at once, and
#'   one \code{circular::rvonmises} call per agent (which accepts only a single kappa). The size of
#'   each pre-drawn array is \code{n * block * n_candidates}.
#'
#'   No state other than position and heading is changed: there is no foraging, no forage
#'   depletion, and no energetics, and the intermediate path is discarded.
#'
#' @param x0 Numeric vector of starting x coordinates, one element per agent.
#' @param y0 Numeric vector of starting y coordinates, one element per agent.
#' @param heading0 Numeric vector of starting headings (radians), one element per agent.
#' @param issf A list of per-agent iSSF model objects from \code{amt::make_issf_model()}.
#' @param move A list of movement data: \code{geom} (a single-layer \code{terra::SpatRaster}
#'   supplying the grid for \code{cellFromXY}), \code{forage}, \code{escape}, \code{canopy}
#'   (numeric covariate vectors indexed by cell number), \code{tod} (0 for day, 1 for night;
#'   constant for the whole burn-in), \code{ext} (the named grid extent, xmin / xmax / ymin /
#'   ymax), and \code{n_candidates}.
#' @param n_steps Integer. Number of burn-in steps taken by each agent.
#' @param block Integer. Number of steps to pre-draw candidate steps for at a time.
#'
#' @return A numeric matrix with one row per agent and columns \code{x}, \code{y}, \code{heading}.
#'
#' @importFrom stats rgamma runif
#' @keywords internal
simulate_burn_in <- function(x0, y0, heading0, issf, move, n_steps, block = 50L) {

  # ----------------------------------------------------------------------------------------------------------------------
  # per-agent parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the agent count, the candidate count, and each agent's movement-distribution parameters as vectors across
  #   agents: the gamma step-length shape and scale and the von Mises turn-angle concentration (kappa) and mean
  #   (mu)

  n <- length(x0)
  nc <- move$n_candidates

  shape <- vapply(issf, function(m) m$sl_$params$shape, numeric(1))
  scale <- vapply(issf, function(m) m$sl_$params$scale, numeric(1))
  kappa <- vapply(issf, function(m) m$ta_$params$kappa, numeric(1))
  mu <- lapply(issf, function(m) m$ta_$params$mu)

  #2) the twelve selection coefficients of every agent as one matrix, agents in rows and coefficients in columns

  coef_names <- c(
    "sl_",
    "log(sl_)",
    "cos(ta_)",
    "forage_biomass_end",
    "escape_terrain_end",
    "canopy_cover_end",
    "sl_:tod_end_night_end",
    "log(sl_):tod_end_night_end",
    "cos(ta_):tod_end_night_end",
    "forage_biomass_end:tod_end_night_end",
    "escape_terrain_end:tod_end_night_end",
    "canopy_cover_end:tod_end_night_end"
  )

  coef_mat <- t(vapply(issf, function(m) m$coefficients[coef_names], numeric(length(coef_names))))
  colnames(coef_mat) <- coef_names

  #3) the six working coefficients: the day main effects, with each night offset added when the burn-in runs at
  #   night (tod = 1). time of day is constant for the whole burn-in, so this is done once

  b_sl <- coef_mat[, "sl_"]
  b_lsl <- coef_mat[, "log(sl_)"]
  b_cta <- coef_mat[, "cos(ta_)"]
  b_fo <- coef_mat[, "forage_biomass_end"]
  b_es <- coef_mat[, "escape_terrain_end"]
  b_ca <- coef_mat[, "canopy_cover_end"]

  if (move$tod == 1) {
    b_sl <- b_sl + coef_mat[, "sl_:tod_end_night_end"]
    b_lsl <- b_lsl + coef_mat[, "log(sl_):tod_end_night_end"]
    b_cta <- b_cta + coef_mat[, "cos(ta_):tod_end_night_end"]
    b_fo <- b_fo + coef_mat[, "forage_biomass_end:tod_end_night_end"]
    b_es <- b_es + coef_mat[, "escape_terrain_end:tod_end_night_end"]
    b_ca <- b_ca + coef_mat[, "canopy_cover_end:tod_end_night_end"]
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # walk every agent forward
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the working position and heading of every agent, and the grid extent used by the outside-extent guard

  x <- x0
  y <- y0
  heading <- heading0
  ext <- move$ext

  steps_done <- 0L

  while (steps_done < n_steps) {

    #2) the number of steps in this block: block, or the steps remaining if fewer

    nb <- min(block, n_steps - steps_done)

    #3) candidate step lengths for every agent and every step in the block, in one gamma draw. matrix() fills by
    #   column, so repeating each agent's shape and scale places agent i's draws in row i

    sl_blk <- matrix(
      stats::rgamma(
        n = n * nb * nc,
        shape = rep(shape, times = nb * nc),
        scale = rep(scale, times = nb * nc)
      ),
      nrow = n
    )

    #4) candidate turn angles, drawn one agent at a time because circular::rvonmises() accepts only a single kappa,
    #   then wrapped to (-pi, pi]

    ta_blk <- matrix(0, nrow = n, ncol = nb * nc)

    for (i in seq_len(n)) {

      angles <- suppressWarnings(
        as.numeric(
          circular::rvonmises(
            n = nb * nc,
            mu = mu[[i]],
            kappa = kappa[i]
          )
        )
      )

      angles <- angles %% (2 * pi)
      ta_blk[i, ] <- ifelse(angles > pi, angles - (2 * pi), angles)
    }

    #5) take the block one step at a time

    for (s in seq_len(nb)) {

      #6) this step's candidate step lengths and turn angles (the block columns for step s), and the endpoints
      #   they imply: each turn angle is added to the agent's current heading to give a bearing, and the step
      #   length is laid out along it

      cols <- ((s - 1) * nc + 1):(s * nc)
      sl <- sl_blk[, cols, drop = FALSE]
      ta <- ta_blk[, cols, drop = FALSE]

      bearing <- heading + ta
      x2 <- x + sl * cos(bearing)
      y2 <- y + sl * sin(bearing)

      #7) the fraction of each agent's candidate endpoints that fall outside the grid extent

      frac_outside <- rowMeans(
        x2 < ext["xmin"] | x2 > ext["xmax"] | y2 < ext["ymin"] | y2 > ext["ymax"]
      )

      #8) the forage, escape-terrain, and canopy values of the cell each endpoint lands in. endpoints off the
      #   raster return NA, which becomes a zero weight below

      cells <- terra::cellFromXY(move$geom, cbind(as.vector(x2), as.vector(y2)))
      fo <- matrix(move$forage[cells], nrow = n)
      es <- matrix(move$escape[cells], nrow = n)
      ca <- matrix(move$canopy[cells], nrow = n)

      #9) the selection linear predictor for every candidate, with each agent's coefficients recycled along its
      #   own row

      w <- sl * b_sl +
        log(sl) * b_lsl +
        cos(ta) * b_cta +
        fo * b_fo +
        es * b_es +
        ca * b_ca

      #10) each row is centered on its own finite mean before exponentiation, and non-finite weights (from NA
      #    covariates or a zero step length) are set to zero so those candidates are never sampled

      finite <- is.finite(w)
      w_bar <- rowSums(ifelse(finite, w, 0)) / rowSums(finite)

      w <- exp(w - w_bar)
      w[!is.finite(w)] <- 0

      #11) the two hold-position conditions: more than half of the candidates outside the extent, or no candidate
      #    with positive weight

      row_w <- rowSums(w)
      hold <- frac_outside > 0.5 | row_w <= 0

      #12) sample one candidate per agent by inverse cumulative distribution across its row of weights: a uniform
      #    draw scaled to the row total is compared with the running sum of weights, and the first candidate whose
      #    running sum reaches it is selected

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

      #13) agents that move take the selected endpoint and the bearing to it as their new heading; held agents keep
      #    their position and heading

      pick <- cbind(seq_len(n), sel)
      x_new <- x2[pick]
      y_new <- y2[pick]

      moved <- !hold
      heading[moved] <- atan2(y_new[moved] - y[moved], x_new[moved] - x[moved])
      x[moved] <- x_new[moved]
      y[moved] <- y_new[moved]
    }

    steps_done <- steps_done + nb
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # return the endpoints
  # ----------------------------------------------------------------------------------------------------------------------

  #1) one row per agent, in the order the agents were supplied

  cbind(
    x = x,
    y = y,
    heading = heading
  )
}
