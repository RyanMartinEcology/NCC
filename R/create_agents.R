#' Create and initialize agent population
#'
#' @description Creates a named list where each element represents an individual
#'   sheep agent as a tibble with one row per simulation time step containing
#'   only time-varying state variables. Fixed individual parameters (reproductive
#'   status, initial body condition, movement parameters) are stored in a
#'   separate \code{agent_params} dataframe. Both are returned as a named list.
#'
#' @param forage_raster A \code{terra::SpatRaster}. Used to define the spatial
#'   extent for initializing agent starting positions.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{agents}{A named list of tibbles, one per agent, containing
#'       time-varying state variables.}
#'     \item{agent_params}{A dataframe with one row per agent containing
#'       fixed individual parameters.}
#'   }
#'
#' @importFrom stats rnorm runif rlnorm
#' @keywords internal
create_agents <- function(forage_raster) {

  n          <- get_param("n_agents")
  t_start    <- get_param("t_start")
  t_end      <- get_param("t_end")
  t_delta    <- get_param("t_delta")
  times      <- seq(t_start, t_end, by = as.numeric(t_delta, units = "secs"))
  bm         <- get_param("bm")
  ids        <- paste0("BHS_", sprintf("%03d", seq_len(n)))

  # draw fixed individual parameters
  ifbfat     <- sapply(seq_len(n), function(i) draw_param("ifbf"))
  rep_status <- sapply(seq_len(n), function(i) draw_param("rep_status"))

  # draw individual movement parameters from population-level log-normal distributions
  day_shape  <- stats::rlnorm(n,
                              meanlog = get_param("day_shape_meanlog"),
                              sdlog   = get_param("day_shape_sdlog"))
  day_scale  <- stats::rlnorm(n,
                              meanlog = get_param("day_scale_meanlog"),
                              sdlog   = get_param("day_scale_sdlog"))
  night_shape <- stats::rlnorm(n,
                               meanlog = get_param("night_shape_meanlog"),
                               sdlog   = get_param("night_shape_sdlog"))
  night_scale <- stats::rlnorm(n,
                               meanlog = get_param("night_scale_meanlog"),
                               sdlog   = get_param("night_scale_sdlog"))
  day_kappa   <- stats::rlnorm(n,
                               meanlog = get_param("day_kappa_meanlog"),
                               sdlog   = get_param("day_kappa_sdlog"))
  night_kappa <- stats::rlnorm(n,
                               meanlog = get_param("night_kappa_meanlog"),
                               sdlog   = get_param("night_kappa_sdlog"))

  # build agent_params dataframe — one row per agent, all fixed parameters
  agent_params <- data.frame(
    id            = ids,
    rep_status    = rep_status,
    j_post_partum = get_param("j_post_partum"),
    bm_init       = bm,
    ifbfat_init   = ifbfat,
    day_shape     = day_shape,
    day_scale     = day_scale,
    night_shape   = night_shape,
    night_scale   = night_scale,
    day_kappa     = day_kappa,
    night_kappa   = night_kappa,
    stringsAsFactors = FALSE
  )

  # draw random starting positions from raster extent
  e      <- terra::ext(forage_raster)
  x_init <- stats::runif(n, e$xmin, e$xmax)
  y_init <- stats::runif(n, e$ymin, e$ymax)

  # build agent tibbles — time-varying state only
  agents        <- vector("list", n)
  names(agents) <- ids

  for (i in seq_len(n)) {

    bm_i      <- bm
    ifbfat_i  <- ifbfat[i]

    empty_rows <- dplyr::tibble(
      datetime        = times,
      id              = ids[i],
      status          = NA_character_,
      x               = NA_real_,
      y               = NA_real_,
      heading         = NA_real_,
      bm              = NA_real_,
      ifbfat          = NA_real_,
      lean_mass       = NA_real_,
      fat_mass        = NA_real_,
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
    empty_rows$heading[1]         <- stats::runif(1, 0, 2 * pi)
    empty_rows$bm[1]              <- bm_i
    empty_rows$ifbfat[1]          <- ifbfat_i
    empty_rows$lean_mass[1]       <- calc_lean_mass(bm_i, ifbfat_i)
    empty_rows$fat_mass[1]        <- calc_fat_mass(bm_i, ifbfat_i)
    empty_rows$j_post_partum[1]   <- get_param("j_post_partum")
    empty_rows$energy_rmr[1]      <- calc_energy_rmr(bm_i)
    empty_rows$energy_hif[1]      <- calc_energy_hif(bm_i)
    empty_rows$energy_loc[1]      <- calc_energy_loc(bm_i)
    empty_rows$forage_consumed[1] <- 0
    empty_rows$fat_change[1]      <- 0

    agents[[i]] <- empty_rows
  }

  list(agents = agents, agent_params = agent_params)
}
