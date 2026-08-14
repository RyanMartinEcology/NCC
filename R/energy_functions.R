#' Calculate lean body mass
#'
#' @description Returns the lean (fat-free) component of body mass as the
#'   complement of the ingesta-free body fat fraction.
#'
#' @param bm Numeric. Body mass (kg).
#' @param ifbfat Numeric. Ingesta-free body fat as a fraction of body mass.
#'
#' @return Numeric. Lean body mass (kg).
#'
#' @keywords internal
calc_lean_mass <- function(bm, ifbfat) {
  bm * (1 - ifbfat)
}

#' Calculate fat mass
#'
#' @description Returns the fat component of body mass from the ingesta-free
#'   body fat fraction.
#'
#' @param bm Numeric. Body mass (kg).
#' @param ifbfat Numeric. Ingesta-free body fat as a fraction of body mass.
#'
#' @return Numeric. Fat mass (kg).
#'
#' @keywords internal
calc_fat_mass <- function(bm, ifbfat) {
  bm * ifbfat
}

#' Calculate basal metabolic energy expenditure
#'
#' @description Returns the daily basal metabolic energy expenditure from the
#'   Chappel and Hudson (1980) basal metabolic rate (460 kJ * kg^-0.75 * day^-1),
#'   scaled by metabolic body mass.
#'
#' @param bm Numeric. Body mass (kg).
#'
#' @return Numeric. Basal metabolic energy expenditure (kJ/day).
#'
#' @keywords internal
calc_energy_bmr <- function(bm) {
  out <- 460 * bm^0.75
  attr(out, 'unit') <- 'kJ * day^-1'
  attr(out, 'source') <- 'Chappel and Hudson 1980'
  attr(out, 'full_name') <- 'Basal metabolic energy expenditure'
  out
}

#' Calculate heat increment of feeding energy expenditure
#'
#' @description Returns the daily heat increment of feeding from the Chappel and
#'   Hudson (1978b) per-body-mass rate, multiplied by the proportion of the diurnal
#'   period spent feeding (Courtemanch et al. 2014) and day length. Day length is
#'   read from the precomputed per-time-step \code{day_length} vector (set in
#'   \code{.set_defaults()} from study latitude and day-of-year).
#'
#' @param bm Numeric. Body mass (kg).
#' @param t Integer. Time-step index into the precomputed \code{day_length} vector
#'   (aligned with the hourly time sequence).
#'
#' @return Numeric. Heat increment of feeding energy expenditure (kJ/day).
#'
#' @keywords internal
calc_energy_hif <- function(bm, t) {
  day_length <- get_param("day_length")[t]
  get_param("HIF") * bm * get_param("prop_day_forage") * day_length
}

#' Calculate locomotion energy expenditure
#'
#' @description Returns the daily locomotion energy expenditure, summed over the
#'   day's 24 steps. Each step's cost is the Dailey and Hobbs (1989) distance cost
#'   factor (selected by the signed slope of the step) times body mass times
#'   horizontal step length. Slope is computed from DEM elevations at the step's
#'   start and end coordinates. The caller supplies the day's coordinate and step
#'   length windows and the body mass.
#'
#' @param x_window Numeric vector of x coordinates for the 25 rows spanning the
#'   day's steps (the 24 step rows plus the row immediately preceding the first).
#' @param y_window Numeric vector of y coordinates, aligned with \code{x_window}.
#' @param sl_window Numeric vector of the 24 step lengths for the day.
#' @param bm Numeric. Body mass (kg) for the day.
#' @param dem_vals Numeric vector of DEM elevations indexed by cell number,
#'   precomputed once via \code{terra::values(dem)} by the caller.
#'
#' @return Numeric. Locomotion energy expenditure (kJ/day).
#'
#' @keywords internal
calc_energy_loc <- function(x_window, y_window, sl_window, bm, dem_vals) {

  # ----------------------------------------------------------------------------------------------------------------------
  # locomotion energy summed over the day's steps
  # ----------------------------------------------------------------------------------------------------------------------

  #1) step start and end coordinates: start = rows 1:24, end = rows 2:25 of the window

  n <- length(sl_window)
  start <- cbind(x_window[1:n], y_window[1:n])
  end <- cbind(x_window[2:(n + 1)], y_window[2:(n + 1)])

  #2) start and end elevations: all rasters share the model CRS, so index the in-memory dem
  #   values by cell (cheaper than terra::extract) rather than reprojecting

  elev_start <- dem_vals[terra::cellFromXY(dem, start)]
  elev_end <- dem_vals[terra::cellFromXY(dem, end)]

  #3) signed slope of each step in degrees: descent negative, incline positive

  slope <- atan2(elev_end - elev_start, sl_window) * 180 / pi

  #4) distance cost factor per step, selected by slope bin: read the five constants once and
  #   fill by bin via which() (NA-safe: a step with NA slope keeps an NA factor, as the prior
  #   nested ifelse produced, so the day's energy_loc propagates NA)

  f_d10  <- get_param("distance_cost_factor_d_10")
  f_d110 <- get_param("distance_cost_factor_d_1_10")
  f_f    <- get_param("distance_cost_factor_f")
  f_i110 <- get_param("distance_cost_factor_i_1_10")
  f_i10  <- get_param("distance_cost_factor_i_10")

  factor <- rep(NA_real_, n)
  factor[which(slope < -10)]                <- f_d10
  factor[which(slope >= -10 & slope < -1)]  <- f_d110
  factor[which(slope >= -1  & slope <= 1)]  <- f_f
  factor[which(slope >  1   & slope <= 10)] <- f_i110
  factor[which(slope > 10)]                 <- f_i10

  #5) daily cost: factor times body mass times horizontal step length, summed and converted
  #   from joules to kilojoules

  sum(factor * bm * sl_window) / 1000
}

#' Calculate lactation energy expenditure
#'
#' @description Returns the daily lactation energy cost as basal metabolic energy
#'   scaled by the lactation modifier indexed at days post partum. Returns 0 for
#'   non-reproductive individuals.
#'
#' @param bmr Numeric. Basal metabolic energy expenditure (kJ/day).
#' @param j_post_partum Numeric. Days post partum; indexes the lactation modifier.
#' @param rep_status Numeric. Reproductive status (1 = lactating, 0 = not).
#'
#' @return Numeric. Lactation energy expenditure (kJ/day).
#'
#' @keywords internal
calc_energy_rep <- function(bmr, j_post_partum, rep_status) {
  if (rep_status == 0) return(0)
  bmr * get_param("lactation_modifier")[as.integer(j_post_partum)]
}

#' Calculate hourly dry matter intake
#'
#' @description Returns the hourly dry matter intake at the agent's current cell
#'   from a negative-exponential functional response of standing vegetation
#'   biomass, mass-specific in metabolic body mass:
#'   \deqn{I_{day} = intake\_max \times bm^{0.75} \times (1 - e^{-V / intake\_decay})}
#'   where \eqn{V} is biomass in kg/ha. The published response is a daily intake,
#'   so it is divided by the number of hours the agent can forage that day to give
#'   an hourly rate; summed across the day's foraging hours this returns the daily
#'   curve evaluated at the mean density the agent actually experienced. Intake is
#'   floored at the available biomass so consumption cannot exceed what the cell
#'   holds; the returned value is the actual consumed amount and is used both to
#'   deplete the cell and to compute energy intake. The whole response is scaled by
#'   the agent's \code{intake_multiplier} slot, which is 1 for both reproductive
#'   states by default.
#'
#' @param density Numeric. Forage density at the agent's current cell (g/cell).
#' @param rep_status Numeric. Reproductive status (1 = lactating, 0 = not); selects
#'   the agent's slot of the \code{intake_multiplier} vector
#'   (\code{rep_status + 1}: index 1 = nonrepro, index 2 = repro).
#' @param bm Numeric. The agent's current body mass (kg).
#' @param forage_hours Numeric. Number of foraging (daylight) hours in the current
#'   day, used to spread the daily intake across the hours the agent can eat.
#' @param cell_area Numeric. Area of one raster cell (m^2), used to convert the
#'   cell's g/cell biomass to the kg/ha units the response is defined in.
#'
#' @return Numeric. Consumed dry matter intake (g/hour), floored at \code{density}.
#'
#' @keywords internal
calc_dmi <- function(density, rep_status, bm, forage_hours, cell_area) {

  #1) convert the cell's biomass from g/cell to kg/ha. one kg/ha spread over a cell of cell_area
  #   square metres is cell_area / 10 grams, so dividing by that factor inverts the conversion

  biomass_kg_ha <- density / (cell_area / 10)

  #2) the daily intake for this body mass at this biomass, scaled by the agent's reproductive-state
  #   multiplier, then the even hourly share of it

  daily <- get_param("intake_max") * bm^0.75 *
    (1 - exp(-biomass_kg_ha / get_param("intake_decay"))) *
    get_param("intake_multiplier")[[rep_status + 1L]]

  dmi <- daily / forage_hours

  #3) an agent cannot eat more than the cell holds

  min(dmi, density)
}

#' Calculate hourly energy intake
#'
#' @description Converts consumed dry matter intake to metabolizable energy via
#'   digestible energy and the digestible-to-metabolizable conversion factor.
#'
#' @param dmi Numeric. Consumed dry matter intake (g/hour), from \code{calc_dmi}.
#'
#' @return Numeric. Metabolizable energy intake (kJ/hour).
#'
#' @keywords internal
calc_energy_i <- function(dmi) {
  dmi * get_param("DE") * get_param("DE_to_ME_conversion_factor")
}

#' Calculate daily net energy balance
#'
#' @description Sums the day's hourly energy intake and subtracts the daily
#'   expenditure terms (basal metabolism, heat increment of feeding, locomotion,
#'   lactation) to give net energy balance. Operates on scalars and the day's
#'   intake vector; the caller writes the returned values to the agent's row.
#'
#' @param energy_i_window Numeric vector of the day's 24 hourly energy intakes (kJ).
#' @param energy_bmr Numeric. Daily basal metabolic energy expenditure (kJ/day).
#' @param energy_hif Numeric. Daily heat increment of feeding (kJ/day).
#' @param energy_loc Numeric. Daily locomotion energy expenditure (kJ/day).
#' @param energy_rep Numeric. Daily lactation energy expenditure (kJ/day).
#'
#' @return A named numeric vector with elements \code{daily_intake} and
#'   \code{energy_net} (kJ/day).
#'
#' @keywords internal
calc_energy_net <- function(energy_i_window, energy_bmr, energy_hif, energy_loc, energy_rep) {

  # ----------------------------------------------------------------------------------------------------------------------
  # net daily energy balance
  # ----------------------------------------------------------------------------------------------------------------------

  #1) total metabolizable energy taken in over the day

  daily_intake <- sum(energy_i_window, na.rm = TRUE)

  #2) net balance: intake minus basal metabolism, heat increment, locomotion, and lactation

  net <- daily_intake - energy_bmr - energy_hif - energy_loc - energy_rep

  #3) return both, unnamed so the two elements carry exactly daily_intake and energy_net; a
  #   named input would otherwise compound the names

  c(daily_intake = unname(daily_intake), energy_net = unname(net))
}

#' Calculate daily fat-mass change
#'
#' @description Converts the day's net energy balance into a change in fat mass at
#'   the energy density of fat (\code{E_fat}), with no conversion loss in either
#'   direction: a surplus deposits and a deficit mobilizes at the same rate.
#'   Operates on a scalar; the caller writes the returned value to the agent's row.
#'
#' @param net Numeric. Daily net energy balance (kJ/day).
#'
#' @return Numeric. Daily fat-mass change (kg).
#'
#' @keywords internal
calc_fat_change <- function(net) {

  # ----------------------------------------------------------------------------------------------------------------------
  # convert net energy balance to a change in fat mass
  # ----------------------------------------------------------------------------------------------------------------------

  #1) surplus and deficit both convert at the energy density of fat, so the sign of net carries
  #   straight through and no branch is needed

  fat_change_g <- net / get_param("E_fat")

  #2) convert grams to kilograms

  fat_change_g / 1000
}
