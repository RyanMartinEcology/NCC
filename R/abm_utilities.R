#' Validate the iSSF coefficient structure for the hand-rolled predictor
#'
#' @description \code{simulate_move} computes the selection linear predictor by
#'   hand (rather than via \code{model.matrix}) and therefore hard-codes the
#'   model's term structure: six habitat / movement covariates, each with a
#'   day main effect and a \code{:tod_end_night_end} night offset. This guard
#'   checks that a fitted/constructed iSSF carries exactly those twelve named
#'   coefficients and \code{stop()}s otherwise, so a change to the term
#'   structure fails loudly instead of being silently miscomputed. Coefficient
#'   \emph{values} are irrelevant here; only the names are checked. The names are
#'   shared across agents, so the caller validates one agent once per run.
#'
#' @param coefs The named coefficient vector \code{issf$coefficients}.
#'
#' @return Invisibly \code{TRUE}; called for its side effect (error on mismatch).
#'
#' @keywords internal
validate_move_coefs <- function(coefs) {
  expected <- c(
    "sl_", "log(sl_)", "cos(ta_)",
    "forage_biomass_end", "escape_terrain_end", "canopy_cover_end",
    "sl_:tod_end_night_end", "log(sl_):tod_end_night_end", "cos(ta_):tod_end_night_end",
    "forage_biomass_end:tod_end_night_end", "escape_terrain_end:tod_end_night_end",
    "canopy_cover_end:tod_end_night_end"
  )
  if (!setequal(names(coefs), expected)) {
    stop(
      "simulate_move expects an iSSF with exactly these coefficients:\n  ",
      paste(expected, collapse = ", "),
      "\nThe hand-rolled selection predictor must be updated if the term ",
      "structure changes."
    )
  }
  invisible(TRUE)
}

#' Simulate one movement step for a single agent
#'
#' @description Draws one movement step for a living agent from its individual
#'   iSSF. This is a stripped, in-memory reimplementation of the continuous-case
#'   path of \code{amt::redistribution_kernel()}, specialized to this model:
#'   it draws candidate steps with amt's \code{random_steps_simple}, looks up the
#'   end-of-step covariates by cell index against precomputed in-memory value
#'   vectors (no \code{terra::extract} per call), scores them with amt's
#'   continuous-case weight formula inlined (no movement compensation, since
#'   \code{landscape = "continuous"}), and samples one endpoint. Step length,
#'   turn angle, and the updated heading are derived from the drawn endpoint.
#'
#'   The candidate draw and the final sample consume the RNG in the same order
#'   and amount as \code{amt::redistribution_kernel()}, and \code{cellFromXY}
#'   returns the same nearest cell \code{terra::extract} would, so for a given
#'   seed the result matches the amt path.
#'
#' @param x0 Numeric. Previous-row x coordinate (step start).
#' @param y0 Numeric. Previous-row y coordinate (step start).
#' @param heading0 Numeric. Previous-row heading (radians); the turn-angle reference.
#' @param issf The agent's iSSF model object from \code{amt::make_issf_model()}.
#' @param move A list of per-hour movement data: \code{geom} (a single-layer
#'   SpatRaster for \code{cellFromXY} geometry), \code{forage}, \code{escape},
#'   \code{canopy} (numeric covariate vectors indexed by cell), \code{tod}
#'   (0 day / 1 night), \code{ext} (named bbox xmin/xmax/ymin/ymax),
#'   \code{t_delta}, and \code{n_candidates}. The selection predictor is
#'   hand-rolled and assumes the twelve-term day/night structure enforced by
#'   \code{validate_move_coefs}.
#' @param time The current time step (\code{POSIXct}).
#'
#' @return A named numeric vector with elements \code{x}, \code{y},
#'   \code{step_length}, \code{turn_angle}, \code{heading}; or \code{NULL} if the
#'   step cannot be drawn, either because more than half the candidate endpoints
#'   fall outside the raster extent or because no candidate carries positive
#'   selection weight (all endpoints on no-data cells).
#'
#' @importFrom amt make_start
#' @keywords internal
simulate_move <- function(x0, y0, heading0, issf, move, time) {

  # ----------------------------------------------------------------------------------------------------------------------
  # draw candidate steps
  # ----------------------------------------------------------------------------------------------------------------------

  #1) tentative start, then draw n_candidates candidates from the agent's pooled
  #   Gamma / von Mises (amt internal; same RNG draw as redistribution_kernel)

  start <- amt::make_start(
    c(x0, y0),
    ta_ = heading0,
    time = time,
    dt = move$t_delta
  )

  xy <- amt:::random_steps_simple(
    start,
    sl_model = issf$sl_,
    ta_model = issf$ta_,
    n.control = move$n_candidates
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # outside-extent guard
  # ----------------------------------------------------------------------------------------------------------------------

  #1) abort if more than half the candidate endpoints spill past the rectangular
  #   extent (tolerance.outside = 0.5), matching the prior kernel behavior

  ext <- move$ext
  fraction_outside <- mean(
    xy$x2_ < ext["xmin"] | xy$x2_ > ext["xmax"] |
      xy$y2_ < ext["ymin"] | xy$y2_ > ext["ymax"]
  )
  if (fraction_outside > 0.5) {
    warning(
      round(fraction_outside * 100, 3),
      "% of candidate steps fall outside the extent (> 50% allowed); step not drawn"
    )
    return(NULL)
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # end-of-step covariates by cell index
  # ----------------------------------------------------------------------------------------------------------------------

  #1) cell each endpoint falls in (NA off-raster); look up the frozen day-start
  #   forage and the static escape / canopy as local vectors. tod is constant for
  #   the hour and enters the predictor below as a whole-call scalar, so it is not
  #   stored per candidate. off-raster NA covariates become zero-weight below,
  #   the same as terra::extract + ssf_weights

  cells <- terra::cellFromXY(move$geom, cbind(xy$x2_, xy$y2_))
  fo <- move$forage[cells]
  es <- move$escape[cells]
  ca <- move$canopy[cells]

  # ----------------------------------------------------------------------------------------------------------------------
  # selection weights (hand-rolled linear predictor, continuous case)
  # ----------------------------------------------------------------------------------------------------------------------

  #1) linear predictor X %*% beta, formed term by term by name rather than via
  #   model.matrix. day main effects are always present; the six
  #   :tod_end_night_end offsets are added only at night (tod == 1), which the
  #   scalar tod lets us branch on instead of carrying a candidate column.
  #   structure is validated once per run by validate_move_coefs(). continuous
  #   landscape => no movement-kernel compensation term

  coefs <- issf$coefficients
  sl <- xy$sl_
  lsl <- log(sl)
  cta <- cos(xy$ta_)

  w <- sl  * coefs[["sl_"]] +
    lsl * coefs[["log(sl_)"]] +
    cta * coefs[["cos(ta_)"]] +
    fo  * coefs[["forage_biomass_end"]] +
    es  * coefs[["escape_terrain_end"]] +
    ca  * coefs[["canopy_cover_end"]]

  if (move$tod == 1) {
    w <- w +
      sl  * coefs[["sl_:tod_end_night_end"]] +
      lsl * coefs[["log(sl_):tod_end_night_end"]] +
      cta * coefs[["cos(ta_):tod_end_night_end"]] +
      fo  * coefs[["forage_biomass_end:tod_end_night_end"]] +
      es  * coefs[["escape_terrain_end:tod_end_night_end"]] +
      ca  * coefs[["canopy_cover_end:tod_end_night_end"]]
  }

  #2) exponentiate and center; non-finite weights (off-raster NA covariates) are
  #   set to zero so they are never sampled

  w <- exp(w - mean(w[is.finite(w)], na.rm = TRUE))
  w[!is.finite(w)] <- 0

  # ----------------------------------------------------------------------------------------------------------------------
  # sample one endpoint and derive movement quantities
  # ----------------------------------------------------------------------------------------------------------------------

  #1) if no candidate carries positive weight, every endpoint fell on a no-data cell inside the
  #   extent; return NULL so the caller holds the agent in place, exactly as for a cornered agent.
  #   this fires only where sample.int would otherwise error on an all-zero probability vector, so
  #   any step that can be drawn is drawn as before

  if (!any(w > 0)) return(NULL)

  #2) sample a single endpoint with probability proportional to the weights

  idx <- sample.int(nrow(xy), size = 1, prob = w)
  x1 <- xy$x2_[idx]
  y1 <- xy$y2_[idx]

  #3) step length, absolute bearing, and turn angle relative to the prior heading

  step_length <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  bearing <- atan2(y1 - y0, x1 - x0)
  turn_angle <- atan2(sin(bearing - heading0), cos(bearing - heading0))

  #3) return the step, stripping any names amt attached to the candidate columns so
  #   the five elements carry exactly the names the caller writes by (x, y,
  #   step_length, turn_angle, heading)

  c(
    x = unname(x1),
    y = unname(y1),
    step_length = unname(step_length),
    turn_angle = unname(turn_angle),
    heading = unname(bearing)
  )
}


#' Simulate the movement burn-in for all agents at once
#'
#' @description Walks every agent forward \code{n_steps} movement-only steps under its own iSSF and
#'   returns only the endpoint. This is a vectorized twin of \code{simulate_move}: it draws the same
#'   gamma step lengths and von Mises turn angles, scores candidates with the same weight formula,
#'   and applies the same two hold-position guards, but advances all agents together one step at a
#'   time rather than running one agent to completion at a time.
#'
#'   Step lengths and turn angles do not depend on position, so they are drawn ahead of the step
#'   loop in blocks of \code{block} steps: one \code{rgamma} call covering every agent at once
#'   (\code{rgamma} recycles shape and scale), and one \code{circular::rvonmises} call per agent
#'   (which accepts only a scalar kappa). Blocking bounds the pre-draw at
#'   \code{n * block * n_candidates} doubles per array.
#'
#'   No state other than position and heading is touched: no foraging, no depletion, no energetics,
#'   and the intermediate path is discarded. Because the random numbers are consumed in a different
#'   order than a per-agent loop would consume them, the result is distributionally identical to
#'   repeated \code{simulate_move} calls but does not reproduce them at a given seed.
#'
#' @param x0 Numeric vector of starting x coordinates, one element per agent.
#' @param y0 Numeric vector of starting y coordinates, one element per agent.
#' @param heading0 Numeric vector of starting headings (radians), one element per agent.
#' @param issf A list of per-agent iSSF model objects from \code{amt::make_issf_model()}.
#' @param move The movement data list described in \code{simulate_move}. \code{tod} is constant for
#'   the whole burn-in.
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

  #1) agent count, candidate count, and the movement-distribution parameters lifted out of each
  #   agent's iSSF into vectors, so the draws below can be made for every agent at once

  n <- length(x0)
  nc <- move$n_candidates

  shape <- vapply(issf, function(m) m$sl_$params$shape, numeric(1))
  scale <- vapply(issf, function(m) m$sl_$params$scale, numeric(1))
  kappa <- vapply(issf, function(m) m$ta_$params$kappa, numeric(1))
  mu <- lapply(issf, function(m) m$ta_$params$mu)

  #2) the twelve selection coefficients as one matrix with agents in rows, so a whole step's weights
  #   form with vector arithmetic instead of a per-agent lookup

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

  #3) the six working coefficients, with the night offsets folded in once if the burn-in runs at
  #   night. tod is constant for the whole burn-in, so this happens here rather than per step

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

  #1) working position and heading, and the extent used by the outside guard

  x <- x0
  y <- y0
  heading <- heading0
  ext <- move$ext

  steps_done <- 0L

  while (steps_done < n_steps) {

    #2) the current block of steps

    nb <- min(block, n_steps - steps_done)

    #3) step lengths for every agent and every step in the block in a single call. matrix() fills by
    #   column, so repeating the per-agent vector places agent i's draws in row i

    sl_blk <- matrix(
      stats::rgamma(
        n = n * nb * nc,
        shape = rep(shape, times = nb * nc),
        scale = rep(scale, times = nb * nc)
      ),
      nrow = n
    )

    #4) turn angles one agent at a time, because circular::rvonmises accepts only a scalar kappa,
    #   then wrapped to (-pi, pi] exactly as amt's random_numbers.vonmises_distr wraps them

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

    #5) consume the block one step at a time

    for (s in seq_len(nb)) {

      #6) this step's candidate step lengths and turn angles, and the endpoints they imply. the
      #   turn angle is taken relative to each agent's current heading

      cols <- ((s - 1) * nc + 1):(s * nc)
      sl <- sl_blk[, cols, drop = FALSE]
      ta <- ta_blk[, cols, drop = FALSE]

      bearing <- heading + ta
      x2 <- x + sl * cos(bearing)
      y2 <- y + sl * sin(bearing)

      #7) outside-extent guard, evaluated per agent

      frac_outside <- rowMeans(
        x2 < ext["xmin"] | x2 > ext["xmax"] | y2 < ext["ymin"] | y2 > ext["ymax"]
      )

      #8) end-of-step covariates by cell index; off-raster candidates return NA and fall out as zero
      #   weight below

      cells <- terra::cellFromXY(move$geom, cbind(as.vector(x2), as.vector(y2)))
      fo <- matrix(move$forage[cells], nrow = n)
      es <- matrix(move$escape[cells], nrow = n)
      ca <- matrix(move$canopy[cells], nrow = n)

      #9) linear predictor, with each agent's coefficients recycled down its own row

      w <- sl * b_sl +
        log(sl) * b_lsl +
        cos(ta) * b_cta +
        fo * b_fo +
        es * b_es +
        ca * b_ca

      #10) exponentiate and center WITHIN each agent's own candidate set, matching simulate_move's
      #    per-call centering; non-finite weights become zero so they are never sampled

      finite <- is.finite(w)
      w_bar <- rowSums(ifelse(finite, w, 0)) / rowSums(finite)

      w <- exp(w - w_bar)
      w[!is.finite(w)] <- 0

      #11) the two hold-position conditions: more than half the candidates outside the extent, or no
      #    candidate carrying positive weight

      row_w <- rowSums(w)
      hold <- frac_outside > 0.5 | row_w <= 0

      #12) sample one candidate per agent by inverse CDF across its own row of weights

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

      #13) advance the agents that drew a step; held agents keep their position and heading

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


#' Simulate one foraging step for a single agent
#'
#' @description Foraging for a living agent at its current cell. During daylight,
#'   reads the forage density at the agent's position from the working values
#'   vector \code{vals}, computes consumed dry matter intake via \code{calc_dmi}
#'   (floored at available biomass), and converts the consumed intake to energy via
#'   \code{calc_energy_i}. At night no foraging occurs and intake is zero. Reads
#'   \code{vals} but does not modify it (to avoid copying the whole vector); the
#'   caller performs the in-place cell depletion using the returned \code{cell}.
#'
#' @param x Numeric. Agent's current x coordinate.
#' @param y Numeric. Agent's current y coordinate.
#' @param vals Numeric vector of current forage biomass (g/cell), indexed by cell
#'   number. Read only.
#' @param geom A single-layer \code{terra::SpatRaster} supplying the grid geometry
#'   for \code{cellFromXY} (the day-start forage layer). Not modified.
#' @param is_day Logical. Whether the current hour is daylight.
#' @param rep_status Numeric. The agent's reproductive status (1 = lactating,
#'   0 = not). Passed through to \code{calc_dmi} to select the agent's
#'   \code{intake_multiplier} slot.
#' @param bm Numeric. The agent's current body mass (kg), passed to \code{calc_dmi}.
#' @param forage_hours Numeric. Number of foraging (daylight) hours in the current
#'   day, passed to \code{calc_dmi}.
#' @param cell_area Numeric. Area of one raster cell (m^2), passed to
#'   \code{calc_dmi} for the g/cell to kg/ha conversion.
#'
#' @return A list with elements \code{forage_consumed} (g), \code{energy_i} (kJ),
#'   and \code{cell} (the depleted cell index, or NA at night).
#'
#' @importFrom terra cellFromXY
#' @keywords internal
simulate_forage <- function(x, y, vals, geom, is_day, rep_status, bm, forage_hours, cell_area) {

  # ----------------------------------------------------------------------------------------------------------------------
  # night gate
  # ----------------------------------------------------------------------------------------------------------------------

  #1) at night no foraging occurs

  if (!is_day) {
    return(list(forage_consumed = 0, energy_i = 0, cell = NA_integer_))
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # read density at the current cell
  # ----------------------------------------------------------------------------------------------------------------------

  #1) locate the cell and read its current biomass from the vector

  cell <- terra::cellFromXY(geom, cbind(x, y))
  density <- vals[cell]

  # ----------------------------------------------------------------------------------------------------------------------
  # compute intake and energy
  # ----------------------------------------------------------------------------------------------------------------------

  #1) consumed dry matter intake, floored at available biomass inside calc_dmi. the agent's
  #   reproductive status selects its intake_multiplier slot

  consumed <- calc_dmi(density, rep_status, bm, forage_hours, cell_area)

  # ----------------------------------------------------------------------------------------------------------------------
  # return intake, energy, and the cell for the caller to deplete
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the caller writes forage_consumed and energy_i and depletes vals[cell]

  list(
    forage_consumed = consumed,
    energy_i = calc_energy_i(consumed),
    cell = cell
  )
}
