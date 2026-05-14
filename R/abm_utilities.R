#' Update agent states at the end of a simulation day
#' @description Updates all agents in the list at the end of each simulation
#'   day using forage consumed across the 24 hourly steps of that day.
#'   Called once per day when hour == 23. Looks back 24 time steps to retrieve
#'   the previous day's end-of-day body condition values.
#' @param agents A named list of agent tibbles.
#' @param t Integer. Current time step index (hour 23 of the current day).
#' @keywords internal
update_agents_test <- function(agents, t) {

  for (i in seq_along(agents)) {
    prev <- t - 24L

    # check last known status by walking back to last non-NA status
    last_known <- tail(which(!is.na(agents[[i]]$status[1:prev])), 1)

    if (length(last_known) > 0 && agents[[i]]$status[last_known] == "DEAD") {
      agents[[i]]$status[t] <- "DEAD"
      next
    }

    # pull previous end-of-day values
    bm_prev         <- agents[[i]]$bm[prev]
    fat_mass_prev   <- agents[[i]]$fat_mass[prev]
    fat_change_prev <- agents[[i]]$fat_change[prev]


    # check if agent dies this time step
    if (length(fat_mass_prev) == 0 || length(fat_change_prev) == 0 ||
        is.na(fat_mass_prev) || is.na(fat_change_prev)) {
      agents[[i]]$status[t] <- "ALIVE"
    } else if ((fat_mass_prev + fat_change_prev) < 0) {
      agents[[i]]$status[t] <- "DEAD"
    } else {
      agents[[i]]$status[t] <- "ALIVE"
    }

    if (agents[[i]]$status[t] == "DEAD") next

    # bm
    agents[[i]]$bm[t] <- bm_prev + fat_change_prev

    # ifbfat
    agents[[i]]$ifbfat[t] <- (fat_mass_prev + fat_change_prev) / agents[[i]]$bm[t]

    # lean_mass
    agents[[i]]$lean_mass[t] <- calc_lean_mass(agents[[i]]$bm[t], agents[[i]]$ifbfat[t])

    # fat_mass
    agents[[i]]$fat_mass[t] <- fat_mass_prev + fat_change_prev

    # j_post_partum
    agents[[i]]$j_post_partum[t] <- agents[[i]]$j_post_partum[prev] + 1L

    # dmi: sum of forage consumed across all 24 hourly steps of the day
    agents[[i]]$dmi[t] <- agents[[i]]$forage_consumed[t]

    # energy_i
    agents[[i]]$energy_i[t] <- calc_energy_i(agents[[i]]$forage_consumed[t])

    # energy_rmr
    agents[[i]]$energy_rmr[t] <- calc_energy_rmr(agents[[i]]$bm[t])

    # energy_hif
    agents[[i]]$energy_hif[t] <- calc_energy_hif(agents[[i]]$bm[t])

    # energy_loc
    agents[[i]]$energy_loc[t] <- calc_energy_loc(agents[[i]]$bm[t])

    # energy_rep
    agents[[i]]$energy_rep[t] <- calc_energy_rep(
      agents[[i]]$energy_rmr[t],
      agents[[i]]$j_post_partum[t],
      agents[[i]]$rep_status[t]
    )

    # energy_net
    agents[[i]]$energy_net[t] <- calc_energy_net(
      agents[[i]]$energy_i[t],
      agents[[i]]$energy_rmr[t],
      agents[[i]]$energy_hif[t],
      agents[[i]]$energy_loc[t],
      agents[[i]]$energy_rep[t]
    )

    # fat_change
    agents[[i]]$fat_change[t] <- calc_fat_change(agents[[i]]$energy_net[t])
  }

  agents
}

#' Determine whether a given datetime is during daylight hours
#'
#' @description Computes sunrise and sunset times at the study area location
#'   using a standard solar position algorithm (NOAA), then returns TRUE if
#'   the supplied datetime falls within that window. Latitude and longitude are
#'   read from globals (\code{study_lat}, \code{study_lon}).
#'
#' @param datetime POSIXct. The datetime to evaluate. Must be in UTC.
#'
#' @return Logical. TRUE if datetime is between sunrise and sunset, FALSE
#'   otherwise.
#'
#' @keywords internal
is_daylight <- function(datetime) {

  lat <- get_param("study_lat")
  lon <- get_param("study_lon")

  # julian day of year
  jd <- as.integer(format(datetime, "%j"))

  # fractional hour in local solar time
  # longitude correction: 4 minutes per degree, converted to hours
  hour_utc   <- as.numeric(format(datetime, "%H")) +
    as.numeric(format(datetime, "%M")) / 60
  hour_local <- hour_utc + lon / 15

  # solar declination (radians)
  declination <- 23.45 * (pi / 180) *
    sin(2 * pi * (284 + jd) / 365)

  # latitude in radians
  lat_rad <- lat * (pi / 180)

  # hour angle at sunrise/sunset (radians)
  cos_hour_angle <- -tan(lat_rad) * tan(declination)

  # clamp to [-1, 1] to avoid NaN at extreme latitudes
  cos_hour_angle <- max(-1, min(1, cos_hour_angle))
  hour_angle     <- acos(cos_hour_angle) * (180 / pi) / 15  # hours from solar noon

  sunrise <- 12 - hour_angle
  sunset  <- 12 + hour_angle

  hour_local >= sunrise & hour_local <= sunset
}

#' Compute hourly dry matter intake via a Type II functional response
#'
#' @description Computes hourly DMI in grams per hour as a Michaelis-Menten
#'   functional response to forage availability. The maximum intake rate of
#'   372 g/hour is taken from Spalinger and Hobbs (1992). The half-saturation
#'   constant is set globally via \code{half_saturation}.
#'
#' @param extracted Numeric vector. Forage availability in grams per cell at
#'   each agent's current location.
#'
#' @return Numeric vector of hourly DMI values in grams per hour, one per agent.
#'
#' @references Spalinger, D.E. and Hobbs, N.T. (1992). Mechanisms of foraging
#'   in mammalian herbivores: new models of functional response. \emph{The
#'   American Naturalist}, 140(2), 325--348.
#'
#' @keywords internal
dmi <- function(extracted) {
  half_saturation <- get_param("half_saturation")
  (372 * extracted) / (half_saturation + extracted)
}

#' Simulate random movement for a single agent
#'
#' @description Moves a single agent one step by drawing a step length from a
#'   Gamma distribution and a direction from Uniform(0, 2*pi). The proposed
#'   position is clipped to the raster extent if it falls outside.
#'
#' @param x Numeric. Current x coordinate in UTM Zone 12N metres.
#' @param y Numeric. Current y coordinate in UTM Zone 12N metres.
#' @param ext_vec Numeric vector of length 4. Raster extent as
#'   c(xmin, xmax, ymin, ymax).
#' @param res_vec Numeric vector of length 2. Raster resolution as c(x, y).
#' @param shape Numeric. Shape parameter for the Gamma step length distribution.
#'   Default 2.
#' @param scale Numeric. Scale parameter for the Gamma step length distribution.
#'   Default 25.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{x}{New x coordinate.}
#'     \item{y}{New y coordinate.}
#'   }
#'
#' @importFrom stats rgamma runif
#' @keywords internal
simulate_move <- function(x, y, ext_vec, res_vec, shape = 2, scale = 25) {

  # -------------------------------------------------------------------------
  # draw step length and direction
  # -------------------------------------------------------------------------
  step_length <- stats::rgamma(1, shape = shape, scale = scale)
  direction   <- stats::runif(1, 0, 2 * pi)

  # -------------------------------------------------------------------------
  # propose new position and clip to raster extent
  # -------------------------------------------------------------------------
  x_new <- x + step_length * cos(direction)
  y_new <- y + step_length * sin(direction)

  x_new <- max(ext_vec[1], min(ext_vec[2] - res_vec[1], x_new))
  y_new <- max(ext_vec[3], min(ext_vec[4] - res_vec[2], y_new))

  list(x = x_new, y = y_new)
}

#' Simulate forage consumption for a single agent
#'
#' @description Extracts forage availability at the agent's current position,
#'   computes hourly DMI via a Type II functional response using \code{dmi()},
#'   consumes the lesser of hourly DMI and available forage, and decrements
#'   the cell values vector in place.
#'
#' @param x Numeric. Current x coordinate in UTM Zone 12N metres.
#' @param y Numeric. Current y coordinate in UTM Zone 12N metres.
#' @param cell_vals Numeric vector. Forage availability in grams per cell for
#'   the current day, extracted from the raster as a plain vector.
#' @param ext_vec Numeric vector of length 4. Raster extent as
#'   c(xmin, xmax, ymin, ymax).
#' @param res_vec Numeric vector of length 2. Raster resolution as c(x, y).
#' @param ncols Integer. Number of columns in the raster.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{forage_consumed}{Grams of forage consumed at this step.}
#'     \item{cell_vals}{Updated forage values vector after consumption.}
#'   }
#'
#' @keywords internal
simulate_forage <- function(x, y, cell_vals, ext_vec, res_vec, ncols) {

  # -------------------------------------------------------------------------
  # compute cell id from coordinates without calling terra
  # -------------------------------------------------------------------------
  col_id  <- floor((x - ext_vec[1]) / res_vec[1]) + 1L
  row_id  <- floor((ext_vec[4] - y) / res_vec[2]) + 1L
  cell_id <- (row_id - 1L) * ncols + col_id

  # guard against NA cells — treat as empty
  forage_avail <- ifelse(is.na(cell_vals[cell_id]), 0, cell_vals[cell_id])

  # -------------------------------------------------------------------------
  # compute hourly dmi and consume lesser of dmi and available forage
  # -------------------------------------------------------------------------
  hourly_dmi      <- dmi(forage_avail)
  forage_consumed <- min(hourly_dmi, forage_avail)

  # -------------------------------------------------------------------------
  # decrement cell values vector in place
  # -------------------------------------------------------------------------
  cell_vals[cell_id] <- forage_avail - forage_consumed

  list(
    forage_consumed = forage_consumed,
    cell_vals       = cell_vals
  )
}

#' Plot IFB Fat Through Time and End-of-Season Distribution
#' @description Two-panel plot: left panel shows ingesta-free body fat through
#' time for all agents colored by reproductive status; right panel shows a
#' density plot of final ifbfat values with rep_status-specific means.
#' @param agents A named list of agent tibbles from a completed simulation.
#' @return A patchwork ggplot object.
#' @export
plot_result_test <- function(agents) {

  plot_data <- dplyr::bind_rows(agents) |>
    dplyr::filter(status == "ALIVE", !is.na(ifbfat)) |>
    dplyr::mutate(rep_status = factor(rep_status,
                                      levels = c(0, 1),
                                      labels = c("Non-Lactating", "Lactating")))

  final_data <- plot_data |>
    dplyr::filter(!is.na(ifbfat)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(datetime, n = 1) |>
    dplyr::ungroup()

  means_data <- final_data |>
    dplyr::group_by(rep_status) |>
    dplyr::summarise(mean_ifbfat = mean(ifbfat, na.rm = TRUE), .groups = "drop")

  colors <- c("Non-Lactating" = "#2166AC", "Lactating" = "#D6604D")

  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = datetime, y = ifbfat,
                                                group = id, color = rep_status)) +
    ggplot2::geom_line(alpha = 0.6) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Date", y = "Ingesta-Free Body Fat (Percent)",
                  color = "Reproductive Status") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")

  p2 <- ggplot2::ggplot(final_data, ggplot2::aes(x = ifbfat, fill = rep_status,
                                                 color = rep_status)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(data = means_data,
                        ggplot2::aes(xintercept = mean_ifbfat, color = rep_status),
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "Early Winter Ingesta-Free Body Fat (Percent)", y = "Density",
                  color = "Reproductive Status", fill = "Reproductive Status") +
    ggplot2::theme_classic()

  p1 + p2 + patchwork::plot_layout(guides = "collect")
}

##' Plot histogram of end-of-season forage depletion
#' @description Plots a histogram of cell-wise differences between the last
#'   layer of the before and after forage rasters, representing end-of-season
#'   forage depletion per cell due to agent consumption.
#' @param forage_before A \code{terra::SpatRaster} with one layer per simulation
#'   day, with values in grams per cell, before simulation.
#' @param forage_after A \code{terra::SpatRaster} with one layer per simulation
#'   day, with values in grams per cell, after simulation.
#' @return A ggplot object.
#' @export
plot_forage <- function(forage_before, forage_after) {

  n_layers <- terra::nlyr(forage_before)

  diff_data <- dplyr::tibble(
    depletion = as.vector(terra::values(forage_before[[n_layers]])) -
      as.vector(terra::values(forage_after[[n_layers]]))
  )

  ggplot2::ggplot(diff_data, ggplot2::aes(x = depletion)) +
    ggplot2::geom_histogram(fill = "#D6604D", color = "white", bins = 40) +
    ggplot2::labs(x     = "Forage depletion (g/cell)",
                  y     = "Count",
                  title = "End-of-season forage depletion per cell") +
    ggplot2::theme_classic()
}

#' Plot Daily Dry Matter Intake Through Time
#' @description Plots daily dry matter intake through time for all agents
#'   colored by reproductive status.
#' @param agents A named list of agent tibbles from a completed simulation.
#' @return A ggplot object.
#' @export
plot_dmi <- function(agents) {

  plot_data <- dplyr::bind_rows(agents) |>
    dplyr::filter(status == "ALIVE", !is.na(dmi)) |>
    dplyr::mutate(rep_status = factor(rep_status,
                                      levels = c(0, 1),
                                      labels = c("Non-Lactating", "Lactating")))

  means_data <- plot_data |>
    dplyr::group_by(rep_status, datetime) |>
    dplyr::summarise(mean_dmi = mean(dmi, na.rm = TRUE), .groups = "drop")

  colors <- c("Non-Lactating" = "#2166AC", "Lactating" = "#D6604D")

  ggplot2::ggplot() +
    ggplot2::geom_line(data = plot_data,
                       ggplot2::aes(x = datetime, y = dmi,
                                    group = id, color = rep_status),
                       alpha = 0.3) +
    ggplot2::geom_line(data = means_data,
                       ggplot2::aes(x = datetime, y = mean_dmi,
                                    color = rep_status),
                       linewidth = 1.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Date", y = "Dry Matter Intake (g/day)",
                  color = "Reproductive Status") +
    ggplot2::theme_classic()
}

#' Plot global mean ifbfat as a function of number of agents
#' @description Plots the global mean end-of-season ifbfat across replicates
#'   as a function of number of agents, with a horizontal line indicating the
#'   carrying capacity threshold. If the number of replicates is >= 3, a 95%
#'   confidence interval is plotted around the mean.
#' @param output A named list returned by \code{ncc_replicator()}.
#' @return A ggplot object.
#' @export
plot_replicator <- function(output) {

  carrying_capacity <- as.numeric(get_param("carrying_capacity"))
  replicates        <- as.integer(get_param("replicates"))

  # compute summary statistics per agent level
  summary_data <- output$results |>
    dplyr::group_by(n_agents) |>
    dplyr::summarise(
      sd_ifbfat   = sd(mean_ifbfat, na.rm = TRUE),
      mean_ifbfat = mean(mean_ifbfat, na.rm = TRUE),
      n           = dplyr::n(),
      .groups     = "drop"
    ) |>
    dplyr::mutate(
      se       = sd_ifbfat / sqrt(n),
      ci_lower = ifelse(n >= 3, mean_ifbfat - qt(0.975, df = n - 1) * se, NA_real_),
      ci_upper = ifelse(n >= 3, mean_ifbfat + qt(0.975, df = n - 1) * se, NA_real_)
    )

  p <- ggplot2::ggplot(summary_data,
                       ggplot2::aes(x = n_agents, y = mean_ifbfat)) +
    ggplot2::geom_hline(yintercept = carrying_capacity,
                        linetype   = "dashed",
                        color      = "#D6604D",
                        linewidth  = 0.8) +
    ggplot2::annotate("text",
                      x      = min(summary_data$n_agents),
                      y      = carrying_capacity,
                      label  = paste("Carrying capacity threshold:", carrying_capacity),
                      hjust  = 0,
                      vjust  = -0.5,
                      color  = "#D6604D",
                      size   = 3.5)

  # add confidence interval ribbon if replicates >= 3
  if (replicates >= 3 && any(!is.na(summary_data$ci_lower))) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           fill  = "#2166AC",
                           alpha = 0.2,
                           na.rm = TRUE)
  }

  p <- p +
    ggplot2::geom_line(color     = "#2166AC", linewidth = 1) +
    ggplot2::geom_point(color    = "#2166AC", size      = 3) +
    ggplot2::labs(x     = "Number of agents",
                  y     = "Global mean ifbfat",
                  title = "Nutritional carrying capacity") +
    ggplot2::theme_classic()

  p
}
