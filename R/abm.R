#' Run the NCC Agent-Based Model
#'
#' @param verbose Logical. If \code{TRUE} (the default), progress messages are
#'   printed as the simulation runs; if \code{FALSE}, all messages are suppressed.
#'
#' @export
ncc_abm <- function(forage_reference, dem, canopy, escape, verbose = TRUE) {

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

  if (verbose) message("1) Importing and checking rasters")

  #1) check forage_reference has one layer per simulation day

  stopifnot(
    "forage_reference must have one layer per simulation day" = terra::nlyr(forage_reference) == length(dates)
  )

  #2) check canopy, escape, and dem each have a single layer

  stopifnot(
    "canopy must have a single layer" = terra::nlyr(canopy) == 1,
    "escape must have a single layer" = terra::nlyr(escape) == 1,
    "dem must have a single layer" = terra::nlyr(dem) == 1
  )

  #3) force all rasters into memory so per-day and per-step access is a RAM read
  # the season pulls one reference layer per day and the kernel extracts escape/canopy
  #   at every candidate step; reading these from disk is slow, so load all rasters
  #   into memory once up front. arithmetic forcing strips names and time, so both
  #   are restored immediately after

  reference_dates <- terra::time(forage_reference)

  forage_reference <- forage_reference * 1
  dem <- dem * 1
  canopy <- canopy * 1
  escape <- escape * 1

  terra::time(forage_reference) <- reference_dates
  names(escape) <- "escape_terrain"
  names(canopy) <- "canopy_cover"

  #4) assign input rasters to the global environment

  list2env(
    list(
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

  if (verbose) message("2) Creating agents")

  #1) initialize agents and fixed individual parameters

  init <- create_agents(forage_reference, dem)
  agents <- init$agents
  agent_params <- init$agent_params

  # ------------------------------------------------------------------------------------------------------------------------
  # simulate daily and hourly dynamics
  # ------------------------------------------------------------------------------------------------------------------------

  if (verbose) message("3) Simulating movement and foraging")

  #1) initialize per-run state: the persistent grazing deficit (g/cell, zero
  #   everywhere), the per-time-step daylight schedule, and the dem values vector
  #   (for cell-index elevation lookups in calc_energy_loc), all hoisted out of the loop
  # realized biomass on any day is potential minus this deficit; grazing adds to it
  #   and it decays each day toward zero. starts at zero (ungrazed at season start)

  deficit <- rep(0, terra::ncell(forage_reference))
  is_daylight <- get_param("is_daylight")
  dem_vals <- terra::values(dem, mat = FALSE)

  #2) loop over each simulation day

  for (d in seq_along(dates)) {

    #3) announce the current day

    if (verbose) message("Day ", d, " of ", length(dates), ": ", dates[d])

    #4) build this day's realized biomass: potential minus the carried deficit
    # forage_reference holds potential biomass; realized is potential_d - deficit.
    #   vals is the day-start realized biomass and depletes within the day as agents
    #   forage

    potential_d <- terra::values(forage_reference[[d]], mat = FALSE)
    vals <- potential_d - deficit

    forage_layer <- forage_reference[[d]]
    names(forage_layer) <- "forage_biomass"
    terra::values(forage_layer) <- vals

    #5) prebuild the day and night covariate maps once for the day
    # forage_layer is frozen for the day (movement sees day-start realized forage),
    #   so both maps stay valid all day; each hour selects between them by daylight

    tod_day <- forage_layer
    names(tod_day) <- "tod_end_night"
    values(tod_day) <- 0

    tod_night <- forage_layer
    names(tod_night) <- "tod_end_night"
    values(tod_night) <- 1

    map_day <- c(forage_layer, escape, canopy, tod_day)
    map_night <- c(forage_layer, escape, canopy, tod_night)

    #6) loop over each hour within the day

    for (h in 1:24) {

      #7) absolute hour index across the season; skip the initialized first row

      t <- (d - 1) * 24 + h
      if (t < 2 || t > length(times)) next

      time_t <- times[t]

      #8) select the hour's covariate map by daylight

      is_day <- is_daylight[t]
      map_t <- if (is_day) map_day else map_night

      #9) filter to living agents and shuffle their processing order

      living <- which(
        vapply(
          agents,
          function(a) a$status[t - 1] == "ALIVE",
          logical(1)
        )
      )
      living <- sample(living)

      #10) carry status forward for living agents
      # at the first hour of a day, status comes from the previous day's last row
      #   (where the survival check was recorded); within a day it is constant

      for (i in living) {
        agents[[i]]$status[t] <- agents[[i]]$status[t - 1]
      }

      #11) carry body state forward within the day
      # the first hour of each day already holds the body state: set by
      #   create_agents on day 1, and by the previous day's energetics
      #   propagation thereafter, so this runs only from the second hour onward

      if (h > 1) {
        for (i in living) {
          agents[[i]]$bm[t] <- agents[[i]]$bm[t - 1]
          agents[[i]]$ifbfat[t] <- agents[[i]]$ifbfat[t - 1]
          agents[[i]]$lean_mass[t] <- agents[[i]]$lean_mass[t - 1]
          agents[[i]]$fat_mass[t] <- agents[[i]]$fat_mass[t - 1]
        }
      }

      #12) move phase — every living agent moves

      for (i in living) {
        mv <- simulate_move(
          agents[[i]]$x[t - 1],
          agents[[i]]$y[t - 1],
          agents[[i]]$heading[t - 1],
          agent_params$issf[[i]],
          map_t,
          time_t
        )
        agents[[i]]$x[t] <- mv[["x"]]
        agents[[i]]$y[t] <- mv[["y"]]
        agents[[i]]$step_length[t] <- mv[["step_length"]]
        agents[[i]]$turn_angle[t] <- mv[["turn_angle"]]
        agents[[i]]$heading[t] <- mv[["heading"]]
      }

      #13) forage phase — every living agent forages in shuffled order, depleting vals

      for (i in living) {
        fg <- simulate_forage(
          agents[[i]]$x[t],
          agents[[i]]$y[t],
          vals,
          forage_layer,
          is_day
        )
        agents[[i]]$forage_consumed[t] <- fg$forage_consumed
        agents[[i]]$energy_i[t] <- fg$energy_i
        if (!is.na(fg$cell)) vals[fg$cell] <- vals[fg$cell] - fg$forage_consumed
      }
    }

    #14) end of day: add the day's consumption to the deficit, then cap and decay it
    # day-start realized was potential_d - deficit; vals is what remains after
    #   grazing, so (potential_d - deficit) - vals is the day's consumption per
    #   cell. add it to the deficit, then cap at the next day's potential and decay.
    #   skipped on the last day (no next day uses the result)

    if (d < length(dates)) {
      deficit <- deficit + ((potential_d - deficit) - vals)
      potential_next <- terra::values(forage_reference[[d + 1]], mat = FALSE)
      deficit <- update_forage(deficit, potential_next)
    }

    #15) end of day: update each living agent's energy balance, mass, and survival
    # daily values are written to the day's last hourly row (t_day); new body mass
    #   and body fat are propagated to the first row of the next day so day d+1
    #   reads the updated state

    t_day <- d * 24
    if (t_day > length(times)) t_day <- length(times)

    day_rows <- ((d - 1) * 24 + 1):t_day

    day_living <- which(
      vapply(
        agents,
        function(a) a$status[t_day] == "ALIVE",
        logical(1)
      )
    )

    for (i in day_living) {

      bm_i <- agents[[i]]$bm[t_day]

      #1) daily energy expenditure terms
      # the locomotion window spans the day's position rows and the steps between
      #   them; on day 1 it starts at the initial row (row 1), giving one fewer step

      loc_start <- max(1, t_day - 24)

      energy_bmr <- calc_energy_bmr(bm_i)
      energy_hif <- calc_energy_hif(bm_i, t_day)
      energy_loc <- calc_energy_loc(
        agents[[i]]$x[loc_start:t_day],
        agents[[i]]$y[loc_start:t_day],
        agents[[i]]$step_length[(loc_start + 1):t_day],
        bm_i,
        dem_vals
      )
      energy_rep <- calc_energy_rep(
        energy_bmr,
        agent_params$j_post_partum[i],
        agent_params$rep_status[i]
      )

      #2) daily net energy balance from the day's hourly intake

      net <- calc_energy_net(
        agents[[i]]$energy_i[day_rows],
        energy_bmr,
        energy_hif,
        energy_loc,
        energy_rep
      )

      #3) fat-mass change and propagated body mass and body fat

      fat_change <- calc_fat_change(net[["energy_net"]])
      mass <- update_mass(bm_i, agents[[i]]$ifbfat[t_day], fat_change)

      #4) write the day's energetics to the day's last row

      agents[[i]]$energy_bmr[t_day] <- energy_bmr
      agents[[i]]$energy_hif[t_day] <- energy_hif
      agents[[i]]$energy_loc[t_day] <- energy_loc
      agents[[i]]$energy_rep[t_day] <- energy_rep
      agents[[i]]$daily_intake[t_day] <- net[["daily_intake"]]
      agents[[i]]$energy_net[t_day] <- net[["energy_net"]]
      agents[[i]]$fat_change[t_day] <- fat_change

      #5) propagate the updated body state to the next day's first row

      if (d < length(dates)) {
        t_next <- t_day + 1
        agents[[i]]$bm[t_next] <- mass[["bm"]]
        agents[[i]]$ifbfat[t_next] <- mass[["ifbfat"]]
        agents[[i]]$lean_mass[t_next] <- calc_lean_mass(mass[["bm"]], mass[["ifbfat"]])
        agents[[i]]$fat_mass[t_next] <- calc_fat_mass(mass[["bm"]], mass[["ifbfat"]])
      }

      #6) survival and days-post-partum: mark dead if fat reserves are exhausted,
      #   and advance the days-post-partum counter by one day
      # death is recorded at t_day so the next day's living filter (which reads the
      #   previous row) excludes the agent from all further processing

      upd <- update_status(
        calc_fat_mass(mass[["bm"]], mass[["ifbfat"]]),
        agent_params$j_post_partum[i]
      )
      agents[[i]]$status[t_day] <- upd$status
      agent_params$j_post_partum[i] <- upd$j_post_partum
    }
  }

  #1) return the agent population and final parameters

  list(agents = agents, agent_params = agent_params)
}
