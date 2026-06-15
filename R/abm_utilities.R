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
#'   \code{step_length}, \code{turn_angle}, \code{heading}; or \code{NULL} if more
#'   than half the candidate endpoints fall outside the raster extent.
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

  #1) sample a single endpoint with probability proportional to the weights

  idx <- sample.int(nrow(xy), size = 1, prob = w)
  x1 <- xy$x2_[idx]
  y1 <- xy$y2_[idx]

  #2) step length, absolute bearing, and turn angle relative to the prior heading

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
#' @param consumed_today Numeric. Dry matter (g) the agent has already consumed
#'   earlier in the current day; used to hold the day's intake under
#'   \code{max_daily_intake}, and ignored at night.
#'
#' @return A list with elements \code{forage_consumed} (g), \code{energy_i} (kJ),
#'   and \code{cell} (the depleted cell index, or NA at night).
#'
#' @importFrom terra cellFromXY
#' @keywords internal
simulate_forage <- function(x, y, vals, geom, is_day, consumed_today) {

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

  #1) consumed dry matter intake, floored at available biomass inside calc_dmi

  consumed <- calc_dmi(density)

  #2) cap the day's cumulative intake at max_daily_intake; NA leaves intake uncapped

  cap <- get_param("max_daily_intake")
  if (!is.na(cap)) {
    consumed <- min(consumed, max(0, cap - consumed_today))
  }

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
