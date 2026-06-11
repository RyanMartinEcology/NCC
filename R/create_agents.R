#' Draw an individual movement model from the population MVN
#'
#' @description Draws one agent's 15-dimensional movement vector from the
#'   population-level multivariate normal distribution (\code{mvn_mu},
#'   \code{mvn_sigma}), exponentiates the three log-scale pooled tentative
#'   distribution dimensions, and builds a single iSSF model via
#'   \code{amt::make_issf_model()}. The model carries \code{tod_end_}
#'   interactions, so day and night are encoded in one model rather than two.
#'
#' @return A single iSSF model object from \code{amt::make_issf_model()}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom amt make_issf_model make_gamma_distr make_vonmises_distr
#' @keywords internal
draw_movement_params <- function() {

  # ----------------------------------------------------------------------------------------------------------------------
  # draw one movement model from the population MVN
  # ----------------------------------------------------------------------------------------------------------------------

  #1) population mean vector and covariance

  mu <- get_param("mvn_mu")
  sigma <- get_param("mvn_sigma")

  #2) draw one 15-dimensional movement vector

  draw <- MASS::mvrnorm(
    1,
    mu = mu,
    Sigma = sigma
  )

  #3) exponentiate the three log-scale pooled tentative distribution parameters

  shape <- exp(draw[["log_shape"]])
  scale <- exp(draw[["log_scale"]])
  kappa <- exp(draw[["log_kappa"]])

  #4) the remaining twelve dimensions are the selection coefficients

  coefs <- draw[setdiff(names(draw), c("log_shape", "log_scale", "log_kappa"))]

  #5) build the single iSSF model, carrying day and night via tod_end_ interactions

  amt::make_issf_model(
    coefs = coefs,
    sl = amt::make_gamma_distr(shape = shape, scale = scale),
    ta = amt::make_vonmises_distr(kappa = kappa)
  )
}

#' Create and initialize agent population
#'
#' @description Creates a named list where each element represents an individual
#'   sheep agent as a tibble with one row per simulation time step containing
#'   only time-varying state variables. Fixed individual parameters (reproductive
#'   status, initial body condition, movement parameters) are stored in a
#'   separate \code{agent_params} dataframe. Both are returned as a named list.
#'
#' @param forage_raster A \code{terra::SpatRaster}. Used to define the set of
#'   valid (non-NA) cells available for initializing agent starting positions.
#' @param dem A \code{terra::SpatRaster} of elevation in meters. Used to restrict
#'   starting positions to a target elevation band.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{agents}{A named list of tibbles, one per agent, containing
#'       time-varying state variables.}
#'     \item{agent_params}{A tibble with one row per agent containing
#'       fixed individual parameters, including a single iSSF model object
#'       (\code{issf}) with tod_end_ interactions.}
#'   }
#'
#' @importFrom stats rnorm runif
#' @keywords internal
create_agents <- function(forage_raster, dem) {

  # ----------------------------------------------------------------------------------------------------------------------
  # set up dimensions, time axis, and identifiers
  # ----------------------------------------------------------------------------------------------------------------------

  #1) agent count, the hourly time sequence, starting body mass, and agent ids

  n <- get_param("n_agents")
  t_start <- get_param("t_start")
  t_end <- get_param("t_end")
  t_delta <- get_param("t_delta")
  times <- seq(
    t_start,
    t_end,
    by = as.numeric(t_delta, units = "secs")
  )
  bm <- get_param("bm")
  ids <- paste0("BHS_", sprintf("%03d", seq_len(n)))

  # ----------------------------------------------------------------------------------------------------------------------
  # draw fixed individual parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) per-agent body condition and reproductive status

  ifbfat <- sapply(seq_len(n), function(i) draw_param("ifbf"))
  rep_status <- sapply(seq_len(n), function(i) draw_param("rep_status"))

  #2) per-agent movement model drawn from the population-level MVN

  issf <- lapply(seq_len(n), function(i) draw_movement_params())

  #3) the fixed-parameter tibble, one row per agent

  agent_params <- dplyr::tibble(
    id = ids,
    rep_status = rep_status,
    j_post_partum = ifelse(rep_status == 1, get_param("j_post_partum"), NA_real_),
    bm_init = bm,
    ifbfat_init = ifbfat,
    issf = issf
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # place agents in the starting elevation band
  # ----------------------------------------------------------------------------------------------------------------------

  #1) target band of 10000-11500 ft converted to meters (1 ft = 0.3048 m); agents are
  #   restricted to this elevation band at initialization

  elev_min <- 10000 * 0.3048
  elev_max <- 11500 * 0.3048

  #2) candidate starting cells: valid (non-NA) forage cells whose dem elevation is in band;
  #   cells outside the study polygon are NA in forage and so are excluded

  ref <- forage_raster[[1]]
  ref_cells <- which(!is.na(terra::values(ref, mat = FALSE)))
  ref_xy <- terra::xyFromCell(ref, ref_cells)

  elev <- terra::extract(dem, ref_xy)[, 1]
  in_band <- !is.na(elev) & elev >= elev_min & elev <= elev_max
  start_xy <- ref_xy[in_band, , drop = FALSE]

  #3) sample one starting cell per agent, with replacement only if too few cells are in band

  start_idx <- sample(
    nrow(start_xy),
    n,
    replace = nrow(start_xy) < n
  )
  x_init <- start_xy[start_idx, 1]
  y_init <- start_xy[start_idx, 2]

  # ----------------------------------------------------------------------------------------------------------------------
  # build per-agent state tibbles
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the named list that will hold one state tibble per agent

  agents <- vector("list", n)
  names(agents) <- ids

  for (i in seq_len(n)) {

    #2) the agent's first-row body condition

    bm_i <- bm
    ifbfat_i <- ifbfat[i]

    #3) an all-NA hourly state frame, one row per time step and one column per state variable

    empty_rows <- dplyr::tibble(
      datetime = times,
      id = ids[i],
      status = NA_character_,
      x = NA_real_,
      y = NA_real_,
      heading = NA_real_,
      step_length = NA_real_,
      turn_angle = NA_real_,
      bm = NA_real_,
      ifbfat = NA_real_,
      lean_mass = NA_real_,
      fat_mass = NA_real_,
      forage_consumed = NA_real_,
      energy_i = NA_real_,
      daily_intake = NA_real_,
      energy_bmr = NA_real_,
      energy_hif = NA_real_,
      energy_loc = NA_real_,
      energy_rep = NA_real_,
      energy_net = NA_real_,
      fat_change = NA_real_
    )

    #4) initialize the first row: status, start position, heading, body state, and the
    #   first-step energetics (locomotion and intake begin at zero)

    empty_rows$status[1] <- "ALIVE"
    empty_rows$x[1] <- x_init[i]
    empty_rows$y[1] <- y_init[i]
    empty_rows$heading[1] <- stats::runif(
      1,
      0,
      2 * pi
    )
    empty_rows$step_length[1] <- 0
    empty_rows$turn_angle[1] <- pi / 2
    empty_rows$bm[1] <- bm_i
    empty_rows$ifbfat[1] <- ifbfat_i
    empty_rows$lean_mass[1] <- calc_lean_mass(bm_i, ifbfat_i)
    empty_rows$fat_mass[1] <- calc_fat_mass(bm_i, ifbfat_i)
    empty_rows$energy_bmr[1] <- calc_energy_bmr(bm_i)
    empty_rows$energy_hif[1] <- calc_energy_hif(bm_i, 1L)
    empty_rows$energy_loc[1] <- 0
    empty_rows$energy_i[1] <- 0
    empty_rows$energy_rep[1] <- calc_energy_rep(
      empty_rows$energy_bmr[1],
      agent_params$j_post_partum[i],
      rep_status[i]
    )
    empty_rows$forage_consumed[1] <- 0
    empty_rows$fat_change[1] <- 0

    agents[[i]] <- empty_rows
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # return
  # ----------------------------------------------------------------------------------------------------------------------

  #1) return the agents and their fixed individual parameters

  list(agents = agents, agent_params = agent_params)
}
