#' Draw an individual movement model from the population MVN
#'
#' @description Draws one agent's 15-dimensional movement vector from the
#'   population-level multivariate normal distribution (\code{mvn_mu},
#'   \code{mvn_sigma}), exponentiates the three log-scale step-length and
#'   turn-angle distribution parameters, and builds a single iSSF model via
#'   \code{amt::make_issf_model()}. The model carries the night offsets as
#'   \code{:tod_end_night_end} interaction coefficients, so day and night are
#'   encoded in one model rather than two.
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

  #1) the population mean vector and covariance matrix

  mu <- get_param("mvn_mu")
  sigma <- get_param("mvn_sigma")

  #2) one 15-dimensional draw

  draw <- MASS::mvrnorm(
    1,
    mu = mu,
    Sigma = sigma
  )

  #3) the gamma step-length shape and scale and the von Mises turn-angle concentration are stored on the log scale,
  #   so they are exponentiated back to their natural scale

  shape <- exp(draw[["log_shape"]])
  scale <- exp(draw[["log_scale"]])
  kappa <- exp(draw[["log_kappa"]])

  #4) the remaining twelve dimensions are the selection coefficients: six day main effects and six night offsets

  coefs <- draw[setdiff(names(draw), c("log_shape", "log_scale", "log_kappa"))]

  #5) build the single iSSF model from the coefficients and the two movement distributions

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
#'   separate \code{agent_params} tibble. Both are returned as a named list.
#'
#' @param forage_reference A \code{terra::SpatRaster} of daily potential forage
#'   biomass. Agent starting locations are drawn uniformly at random from the
#'   non-NA cells of its first layer, so agents may begin anywhere the forage
#'   surface carries data.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{agents}{A named list of tibbles, one per agent, containing
#'       time-varying state variables.}
#'     \item{agent_params}{A tibble with one row per agent containing
#'       fixed individual parameters, including a single iSSF model object
#'       (\code{issf}) carrying the day coefficients and night offsets.}
#'   }
#'
#' @importFrom stats rnorm runif
#' @keywords internal
create_agents <- function(forage_reference) {

  # ----------------------------------------------------------------------------------------------------------------------
  # set up dimensions, time axis, and identifiers
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the agent count, the hourly time sequence from season start to season end, the starting body mass shared by
  #   all agents, and the agent ids (BHS_001, BHS_002, ...)

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

  #1) each agent's starting body fat fraction and reproductive status (1 = lactating, 0 = not), drawn from the
  #   parameter environment's draw functions

  ifbfat <- sapply(seq_len(n), function(i) draw_param("ifbf"))
  rep_status <- sapply(seq_len(n), function(i) draw_param("rep_status"))

  #2) each agent's movement model, drawn from the population-level multivariate normal distribution

  issf <- lapply(seq_len(n), function(i) draw_movement_params())

  #3) the fixed-parameter tibble, one row per agent. days post partum starts at j_post_partum for lactating agents
  #   and is NA for the rest

  agent_params <- dplyr::tibble(
    id = ids,
    rep_status = rep_status,
    j_post_partum = ifelse(rep_status == 1, get_param("j_post_partum"), NA_real_),
    bm_init = bm,
    ifbfat_init = ifbfat,
    issf = issf
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # place agents on the forage surface
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the cells eligible to hold an agent are the non-NA cells of the first forage layer (terra::cells() returns
  #   non-NA cell numbers only). the NA mask is taken to be the same in every daily layer

  valid_cells <- terra::cells(forage_reference[[1]])

  stopifnot("forage_reference has no non-NA cells to place agents in" = length(valid_cells) > 0)

  #2) one start cell per agent, drawn uniformly and with replacement, so every cell carrying forage data is
  #   equally likely regardless of its biomass or distance to escape terrain

  start_cells <- sample(
    valid_cells,
    size = n,
    replace = TRUE
  )

  #3) the start coordinates: the centers of the drawn cells, in the raster's coordinate reference system

  start_xy <- terra::xyFromCell(forage_reference, start_cells)
  x_init <- start_xy[, 1]
  y_init <- start_xy[, 2]

  # ----------------------------------------------------------------------------------------------------------------------
  # build per-agent state tibbles
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the named list that holds one state tibble per agent

  agents <- vector("list", n)
  names(agents) <- ids

  for (i in seq_len(n)) {

    #2) the agent's starting body mass and body fat fraction

    bm_i <- bm
    ifbfat_i <- ifbfat[i]

    #3) an all-NA hourly state tibble, one row per time step and one column per state variable

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

    #4) fill the first row: the agent starts alive at its start cell with a uniformly random heading, a
    #   zero-length step, and a turn angle of pi / 2; body mass, body fat, lean mass, and fat mass from its
    #   starting values; basal metabolism, heat increment of feeding, and lactation cost evaluated at those
    #   values; and locomotion, intake, forage consumed, and fat change at zero

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

  #1) the agents and their fixed individual parameters

  list(agents = agents, agent_params = agent_params)
}
