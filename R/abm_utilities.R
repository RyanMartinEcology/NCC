#' Simulate one movement step for a single agent
#'
#' @description Draws one movement step for a living agent from its individual
#'   iSSF using a redistribution kernel (Signer et al. 2024). The kernel is built
#'   from the agent's stored \code{make_issf_model} object and the supplied
#'   per-hour covariate map. One endpoint is sampled; step length, turn angle, and
#'   the updated heading are derived from the start and the drawn endpoint.
#'   Endpoints falling outside the map or on NA cells are rejected and the kernel
#'   is redrawn. Operates on scalars only; the caller writes the results to the
#'   agent's row.
#'
#' @param x0 Numeric. Previous-row x coordinate (step start).
#' @param y0 Numeric. Previous-row y coordinate (step start).
#' @param heading0 Numeric. Previous-row heading (radians); the turn-angle reference.
#' @param issf The agent's iSSF model object from \code{amt::make_issf_model()}.
#' @param map A \code{terra::SpatRaster} for the current hour with layers named to
#'   match the model coefficients: \code{forage_biomass}, \code{escape_terrain},
#'   \code{canopy_cover}, and a constant \code{tod_end_night} layer (0 day, 1 night).
#' @param time The current time step (\code{POSIXct}).
#'
#' @return A named numeric vector with elements \code{x}, \code{y},
#'   \code{step_length}, \code{turn_angle}, and \code{heading}.
#'
#' @importFrom amt make_start redistribution_kernel
#' @keywords internal
simulate_move <- function(x0, y0, heading0, issf, map, time) {

  # ----------------------------------------------------------------------------------------------------------------------
  # build the start
  # ----------------------------------------------------------------------------------------------------------------------

  #1) construct the sim_start from current position, heading, and time

  start <- amt::make_start(
    c(x0, y0),
    ta_ = heading0,
    time = time,
    dt = get_param("t_delta")
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # draw a step
  # ----------------------------------------------------------------------------------------------------------------------

  #1) draw one step from the kernel
  # candidates that land on NA cells (masked, or past the raster extent) receive zero weight
  #   inside ssf_weights, so a valid cell is sampled without a rejection step.
  # tolerance.outside = 0.5 lets up to half the candidate endpoints spill past the rectangular
  #   extent before the kernel aborts; the default of 0 aborts on a single out-of-extent
  #   candidate, which corrupts the agent position and triggers a downstream "missing value"
  #   error on the following step. note that if more than 50% spill the kernel still returns
  #   NULL, so a deeply cornered agent can still terminate the run (no stay-put guard yet)

  n_candidates <- get_param("n_candidates")

  rdk <- amt::redistribution_kernel(
    issf,
    start = start,
    map = map,
    n.control = n_candidates,
    n.sample = 1,
    landscape = "continuous",
    as.rast = FALSE,
    tolerance.outside = 0.5
  )

  x1 <- rdk$redistribution.kernel$x_
  y1 <- rdk$redistribution.kernel$y_

  # ----------------------------------------------------------------------------------------------------------------------
  # derive movement quantities
  # ----------------------------------------------------------------------------------------------------------------------

  #1) step length is the euclidean distance from start to endpoint

  step_length <- sqrt((x1 - x0)^2 + (y1 - y0)^2)

  #2) absolute bearing of the realized step

  bearing <- atan2(y1 - y0, x1 - x0)

  #3) turn angle is the bearing relative to the previous heading, wrapped to (-pi, pi]

  turn_angle <- atan2(sin(bearing - heading0), cos(bearing - heading0))

  # ----------------------------------------------------------------------------------------------------------------------
  # return the movement quantities
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the caller writes these to the agent's row

  c(
    x = x1,
    y = y1,
    step_length = step_length,
    turn_angle = turn_angle,
    heading = bearing
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
#'
#' @return A list with elements \code{forage_consumed} (g), \code{energy_i} (kJ),
#'   and \code{cell} (the depleted cell index, or NA at night).
#'
#' @importFrom terra cellFromXY
#' @keywords internal
simulate_forage <- function(x, y, vals, geom, is_day) {

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
