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
#' @param start_range An sf polygon or \code{terra::SpatVector} delimiting the
#'   population range (the 95\% KDE home range) within which agent starting
#'   locations are drawn uniformly at random.
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
create_agents <- function(start_range) {

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
  # place agents within the starting range
  # ----------------------------------------------------------------------------------------------------------------------

  #1) accept the starting range as an sf polygon or a SpatVector (the 95% KDE population home range)

  if (inherits(start_range, "sf")) start_range <- terra::vect(start_range)

  #2) draw one uniform-random location per agent from anywhere inside the range polygon(s); random
  #   points are area-weighted across the disjoint sub-ranges, with no cell or elevation constraint

  start_pts <- terra::spatSample(
    start_range,
    size = n,
    method = "random"
  )

  stopifnot("spatSample returned fewer points than agents" = nrow(start_pts) >= n)

  #3) pull the coordinates (exactly n, in the range's CRS)

  start_xy <- terra::crds(start_pts)[seq_len(n), , drop = FALSE]
  x_init <- start_xy[, 1]
  y_init <- start_xy[, 2]

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
