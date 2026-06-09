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
#' @description Returns the daily heat increment of feeding from the Dailey and
#'   Hobbs (1989) per-body-mass rate, multiplied by the proportion of the diurnal
#'   period spent feeding (Courtemanch et al. 2014) and day length. Day length is
#'   computed from solar geometry using the study latitude and the day of year.
#'
#' @param bm Numeric. Body mass (kg).
#' @param time POSIXct. Current simulation time; used to derive day of year.
#'
#' @return Numeric. Heat increment of feeding energy expenditure (kJ/day).
#'
#' @importFrom lubridate yday
#' @keywords internal
calc_energy_hif <- function(bm, time) {
  doy <- lubridate::yday(time)
  decl <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  lat_rad <- get_param("study_lat") * pi / 180
  day_length <- (24 / pi) * acos(-tan(lat_rad) * tan(decl))

  get_param("HIF") * bm * get_param("prop_day_forage") * day_length
}

#' Calculate locomotion energy expenditure
#'
#' @description Returns the daily locomotion energy expenditure, summed over the
#'   24 hourly steps ending at and including \code{time}. Each step's cost is the
#'   Dailey and Hobbs (1989) distance cost factor (selected by the signed slope
#'   of the step) times body mass times horizontal step length. Slope is computed
#'   from DEM elevations at the step's start and end coordinates.
#'
#' @param agent A tibble for a single agent (one row per time step), containing
#'   \code{datetime}, \code{x}, \code{y}, \code{step_length}, and \code{bm}.
#' @param time POSIXct. The time step at which the daily call is made; the window
#'   is the 24 rows ending at and including this time.
#'
#' @return Numeric. Locomotion energy expenditure (kJ/day).
#'
#' @keywords internal
calc_energy_loc <- function(agent, time) {

  # window: the 24 rows ending at and including the row at `time`,
  #   plus the immediately preceding row to start the earliest step
  end_idx <- which(agent$datetime == time)
  rows <- (end_idx - 23):end_idx
  bm <- agent$bm[end_idx]

  # step start/end cells: start = row t-1, end = row t
  start <- agent[rows - 1, c("x", "y")]
  end <- agent[rows, c("x", "y")]
  sl <- agent$step_length[rows]

  # project points from the agent CRS to the DEM CRS and extract elevations
  pts_start <- terra::project(terra::vect(as.matrix(start), crs = get_param("epsg")), terra::crs(dem))
  pts_end <- terra::project(terra::vect(as.matrix(end), crs = get_param("epsg")), terra::crs(dem))
  elev_start <- terra::extract(dem, pts_start)[, 2]
  elev_end <- terra::extract(dem, pts_end)[, 2]

  # signed slope of each step (degrees): descent negative, incline positive
  slope <- atan2(elev_end - elev_start, sl) * 180 / pi

  # select the distance cost factor by slope bin
  factor <- ifelse(
    slope < -10,
    get_param("distance_cost_factor_d_10"),
    ifelse(
      slope < -1,
      get_param("distance_cost_factor_d_1_10"),
      ifelse(
        slope <= 1,
        get_param("distance_cost_factor_f"),
        ifelse(
          slope <= 10,
          get_param("distance_cost_factor_i_1_10"),
          get_param("distance_cost_factor_i_10")
        )
      )
    )
  )

  sum(factor * bm * sl) / 1000
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

#' Calculate hourly energy intake
#'
#' @description Returns the hourly metabolizable energy intake at the agent's
#'   current cell. Dry matter intake is a Michaelis-Menten functional response of
#'   forage density (g/cell) with a 372 g/hour asymptote, then converted to energy
#'   via digestible energy and the digestible-to-metabolizable conversion factor.
#'   The forage layer is selected by matching the date of \code{time}.
#'
#' @param x Numeric. Agent x coordinate (in the forage raster CRS).
#' @param y Numeric. Agent y coordinate (in the forage raster CRS).
#' @param time POSIXct. Current simulation time; selects the forage layer by date.
#'
#' @return Numeric. Metabolizable energy intake (kJ/hour).
#'
#' @keywords internal
calc_energy_i <- function(x, y, time) {
  target_date <- as.Date(time, tz = "America/Denver")
  layer_dates <- as.Date(terra::time(forage), tz = "America/Denver")
  layer_idx <- which(layer_dates == target_date)

  density <- terra::extract(forage[[layer_idx]], cbind(x, y))[, 1]
  dmi <- (372 * density) / (get_param("half_saturation") + density)

  dmi * get_param("DE") * get_param("DE_to_ME_conversion_factor")
}

#' Calculate daily net energy balance
#'
#' @description Sums the day's hourly energy intake and subtracts the daily
#'   expenditure terms (basal metabolism, heat increment of feeding, locomotion,
#'   lactation) to give net energy balance. The day is the 24 hourly rows ending
#'   at and including \code{time}. Writes \code{daily_intake} and \code{energy_net}
#'   to that day's 23:00 row and returns the updated tibble.
#'
#' @param agent A tibble for a single agent (one row per time step), containing
#'   \code{datetime}, \code{energy_i}, and the daily expenditure columns
#'   \code{energy_bmr}, \code{energy_hif}, \code{energy_loc}, \code{energy_rep}.
#' @param time POSIXct. The 23:00 time step at which the daily balance is computed.
#'
#' @return The agent tibble with \code{daily_intake} and \code{energy_net} set on
#'   the row at \code{time}.
#'
#' @keywords internal
calc_energy_net <- function(agent, time) {
  end_idx <- which(agent$datetime == time)
  day_rows <- (end_idx - 23):end_idx

  daily_intake <- sum(agent$energy_i[day_rows], na.rm = TRUE)

  net <- daily_intake -
    agent$energy_bmr[end_idx] -
    agent$energy_hif[end_idx] -
    agent$energy_loc[end_idx] -
    agent$energy_rep[end_idx]

  agent$daily_intake[end_idx] <- daily_intake
  agent$energy_net[end_idx] <- net

  agent
}

#' Calculate daily fat-mass change
#'
#' @description Converts the day's net energy balance into a change in fat mass.
#'   A surplus is deposited at the fat deposition efficiency; a deficit is covered
#'   by mobilizing fat at the fat catabolism efficiency. Reads \code{energy_net}
#'   from the row at \code{time}, writes \code{fat_change} (kg) to that row, and
#'   returns the updated tibble.
#'
#' @param agent A tibble for a single agent (one row per time step), containing
#'   \code{datetime} and \code{energy_net}.
#' @param time POSIXct. The 23:00 time step at which the daily change is computed.
#'
#' @return The agent tibble with \code{fat_change} (kg) set on the row at
#'   \code{time}.
#'
#' @keywords internal
calc_fat_change <- function(agent, time) {
  end_idx <- which(agent$datetime == time)
  net <- agent$energy_net[end_idx]

  if (net >= 0) {
    fat_change_g <- net * get_param("fat_dep") / get_param("E_fat")
  } else {
    fat_change_g <- net / (get_param("fat_eff") * get_param("E_fat"))
  }

  agent$fat_change[end_idx] <- fat_change_g / 1000

  agent
}
