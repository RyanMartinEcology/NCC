#' Simulate a spatiotemporally autocorrelated forage raster stack
#'
#' @description Generates a daily \code{SpatRaster} of forage biomass (kg/ha)
#'   over the simulation window. Spatial autocorrelation is imposed by Gaussian
#'   smoothing of random noise. Temporal decline follows a logistic decay curve
#'   with a mean of \code{peak_mean} on the first day of the simulation (July 1)
#'   and a mean of \code{end_mean} by the end of the season (Oct 31).
#'   Values are clipped to \[0, max_val\] kg/ha.
#'
#' @param nrow Integer. Number of raster rows. Default 100.
#' @param ncol Integer. Number of raster columns. Default 100.
#' @param res Numeric. Cell resolution in metres. Default 30.
#' @param xmin Numeric. Left edge of raster extent in UTM Zone 12N metres.
#'   Default 500000.
#' @param ymin Numeric. Bottom edge of raster extent in UTM Zone 12N metres.
#'   Default 4800000.
#' @param smooth_sigma Numeric. Standard deviation of the Gaussian smoothing
#'   kernel in cells, controlling spatial autocorrelation. Default 5.
#' @param peak_mean Numeric. Target mean forage at peak on day 1 (kg/ha).
#'   Default 700.
#' @param end_mean Numeric. Target mean forage on the final simulation day
#'   (kg/ha). Default 200.
#' @param max_val Numeric. Hard upper ceiling on forage values (kg/ha).
#'   Default 1200.
#' @param seed Integer or NULL. Random seed for reproducibility. Default NULL.
#'
#' @return A \code{terra::SpatRaster} with one layer per simulation day, named
#'   by date (e.g. \code{"2025-07-01"}). Values are in kg/ha. CRS is set to
#'   the value of \code{get_param("epsg")}.
#'
#' @importFrom terra rast ext crs clamp focal focalMat values
#' @importFrom stats rnorm
#' @export
simulate_forage_raster <- function(nrow         = 100L,
                                   ncol         = 100L,
                                   res          = 30,
                                   xmin         = 500000,
                                   ymin         = 4800000,
                                   smooth_sigma = 5,
                                   peak_mean    = 150,
                                   end_mean     = 5,
                                   max_val      = 200,
                                   seed         = NULL) {

  if (!is.null(seed)) set.seed(seed)

  t_start <- get_param("t_start")
  t_end   <- get_param("t_end")
  epsg    <- get_param("epsg")

  dates   <- seq(as.Date(t_start), as.Date(t_end), by = "day")
  n_days  <- length(dates)
  n_cells <- nrow * ncol

  # -------------------------------------------------------------------------
  # Gaussian smoothing kernel for spatial autocorrelation
  # -------------------------------------------------------------------------
  template <- terra::rast(nrows = nrow, ncols = ncol)
  terra::ext(template) <- terra::ext(xmin, xmin + ncol * res, ymin, ymin + nrow * res)
  terra::crs(template) <- epsg

  kernel <- terra::focalMat(
    template,
    d    = smooth_sigma * res,
    type = "Gauss"
  )

  # -------------------------------------------------------------------------
  # temporal decline: logistic decay anchored to peak_mean on day 1
  # and end_mean on day n_days
  # -------------------------------------------------------------------------
  day_index   <- seq_len(n_days)
  midpoint    <- n_days * 0.6
  growth_rate <- 0.05
  temporal_wt <- 1 / (1 + exp(growth_rate * (day_index - midpoint)))

  # rescale so temporal_wt[1] = peak_mean and temporal_wt[n_days] = end_mean
  wt_min      <- temporal_wt[n_days]
  wt_max      <- temporal_wt[1]
  temporal_wt <- (temporal_wt - wt_min) / (wt_max - wt_min)      # [0, 1]
  temporal_wt <- temporal_wt * (peak_mean - end_mean) + end_mean  # [end_mean, peak_mean]

  # -------------------------------------------------------------------------
  # base spatial layer: smoothed random noise rescaled so mean = 1
  # -------------------------------------------------------------------------
  base_noise <- matrix(stats::rnorm(n_cells), nrow = nrow, ncol = ncol)
  base_r     <- terra::rast(base_noise)
  terra::ext(base_r) <- terra::ext(xmin, xmin + ncol * res, ymin, ymin + nrow * res)
  terra::crs(base_r) <- epsg

  base_r <- terra::focal(base_r, w = kernel, fun = "sum", na.policy = "omit")

  # fill edge NAs introduced by focal using iterative mean filter
  while (any(is.na(as.vector(terra::values(base_r))))) {
    base_r <- terra::focal(base_r, w = 3, fun = "mean", na.rm = TRUE)
  }

  # shift to >= 0 and rescale so mean = 1 (per-day scaling applied below)
  v      <- as.vector(terra::values(base_r))
  v      <- v - min(v, na.rm = TRUE)
  v      <- v / mean(v, na.rm = TRUE)
  terra::values(base_r) <- v

  # -------------------------------------------------------------------------
  # build daily layers by scaling spatial base by the temporal weight for
  # that day, then clamp to [0, max_val]
  # -------------------------------------------------------------------------
  layers <- vector("list", n_days)

  for (d in seq_len(n_days)) {
    layer_d     <- base_r * temporal_wt[d]
    layer_d     <- terra::clamp(layer_d, lower = 0, upper = max_val)
    layers[[d]] <- layer_d
  }

  forage_stack        <- terra::rast(layers)
  names(forage_stack) <- as.character(dates)

  forage_stack
}

#' Convert forage raster from kg/ha to grams per cell
#'
#' @description Converts a \code{SpatRaster} of forage biomass in kg/ha to
#'   absolute forage mass in grams per cell, based on cell resolution.
#'
#' @param forage_stack A \code{terra::SpatRaster} with values in kg/ha.
#' @param res Numeric. Cell resolution in metres. Default 30.
#'
#' @return A \code{terra::SpatRaster} with values in grams per cell.
#'
#' @keywords internal
convert_forage_units <- function(forage_stack, res = 30) {
  cell_area_ha <- (res * res) / 10000
  forage_stack * cell_area_ha * 1000
}
