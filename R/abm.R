#' Run the NCC Agent-Based Model
#'
#' @export
ncc_abm <- function(forage, forage_reference, dem, canopy, escape) {

  # ------------------------------------------------------------------------------------------------------------------------
  # resolve time parameters
  # ------------------------------------------------------------------------------------------------------------------------

  #1) pull time parameters from globals

  t_start <- get_param("t_start")
  t_end <- get_param("t_end")
  t_delta <- get_param("t_delta")

  #2) build hourly time sequence

  times <- seq(
    t_start,
    t_end,
    by = as.numeric(t_delta, units = "secs")
  )

  #3) build daily date sequence in America/Denver

  dates <- seq(
    as.Date(t_start, tz = "America/Denver"),
    as.Date(t_end, tz = "America/Denver"),
    by = "day"
  )

  # ------------------------------------------------------------------------------------------------------------------------
  # assign input rasters and check rasters
  # ------------------------------------------------------------------------------------------------------------------------

  message("1) Importing and checking rasters")

  #1) check forage and forage_reference have one layer per simulation day

  stopifnot(
    "forage must have one layer per simulation day" = terra::nlyr(forage) == length(dates),
    "forage_reference must have one layer per simulation day" = terra::nlyr(forage_reference) == length(dates)
  )

  #2) check canopy, escape, and dem each have a single layer

  stopifnot(
    "canopy must have a single layer" = terra::nlyr(canopy) == 1,
    "escape must have a single layer" = terra::nlyr(escape) == 1,
    "dem must have a single layer" = terra::nlyr(dem) == 1
  )

  #3) assign input rasters to the global environment

  list2env(
    list(
      forage = forage,
      forage_reference = forage_reference,
      dem = dem,
      canopy = canopy,
      escape = escape
    ),
    envir = .GlobalEnv
  )

  # ------------------------------------------------------------------------------------------------------------------------
  # create agents
  # ------------------------------------------------------------------------------------------------------------------------

  message("2) Creating agents")

  #1) initialize agents and fixed individual parameters

  init <- create_agents(forage)
  agents <- init$agents
  agent_params <- init$agent_params

}
