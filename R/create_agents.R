#' Create and initialize agent population
#' @description Creates a named list where each element represents an individual
#' sheep agent as a tibble with one row per simulation time step.
#' @param forage_raster A \code{terra::SpatRaster}. Used to define the spatial
#'   extent for initializing agent starting positions.
#' @return A named list of tibbles, one per agent.
#' @importFrom stats rnorm runif
#' @keywords internal
create_agents <- function(forage_raster) {
  n          <- get_param("n_agents")
  t_start    <- get_param("t_start")
  t_end      <- get_param("t_end")
  t_delta    <- get_param("t_delta")
  times      <- seq(t_start, t_end, by = as.numeric(t_delta, units = "secs"))
  bm         <- get_param("bm")
  ifbfat     <- sapply(seq_len(n), function(i) draw_param("ifbf"))
  rep_status <- sapply(seq_len(n), function(i) draw_param("rep_status"))
  ids        <- paste0("BHS_", sprintf("%03d", seq_len(n)))
  agents     <- vector("list", n)
  names(agents) <- ids

  # draw random starting positions from raster extent
  e      <- terra::ext(forage_raster)
  x_init <- stats::runif(n, e$xmin, e$xmax)
  y_init <- stats::runif(n, e$ymin, e$ymax)

  for (i in seq_len(n)) {
    empty_rows <- dplyr::tibble(
      datetime        = times,
      id              = ids[i],
      status          = NA_character_,
      x               = NA_real_,
      y               = NA_real_,
      bm              = NA_real_,
      ifbfat          = NA_real_,
      lean_mass       = NA_real_,
      fat_mass        = NA_real_,
      rep_status      = rep_status[i],
      j_post_partum   = NA_real_,
      forage_consumed = NA_real_,
      dmi             = NA_real_,
      energy_i        = NA_real_,
      energy_rmr      = NA_real_,
      energy_hif      = NA_real_,
      energy_loc      = NA_real_,
      energy_rep      = NA_real_,
      energy_net      = NA_real_,
      fat_change      = NA_real_
    )

    empty_rows$status[1]          <- "ALIVE"
    empty_rows$x[1]               <- x_init[i]
    empty_rows$y[1]               <- y_init[i]
    empty_rows$bm[1]              <- bm
    empty_rows$ifbfat[1]          <- ifbfat[i]
    empty_rows$lean_mass[1]       <- calc_lean_mass(bm, ifbfat[i])
    empty_rows$fat_mass[1]        <- calc_fat_mass(bm, ifbfat[i])
    empty_rows$j_post_partum[1]   <- get_param("j_post_partum")
    empty_rows$energy_rmr[1]      <- calc_energy_rmr(bm)
    empty_rows$energy_hif[1]      <- calc_energy_hif(bm)
    empty_rows$energy_loc[1]      <- calc_energy_loc(bm)
    empty_rows$forage_consumed[1] <- 0
    empty_rows$fat_change[1]      <- 0

    agents[[i]] <- empty_rows
  }

  agents
}
