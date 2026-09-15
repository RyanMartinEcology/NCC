# ----------------------------------------------------------------------------------------------------------------------
# internal agent-state representation
# ----------------------------------------------------------------------------------------------------------------------

# during a run, each agent's time-varying state is stored in a numeric matrix with one row per hourly time step and
#   one column per state variable. writing single values into a matrix is much faster than writing them into a
#   tibble, and the hourly movement, foraging, and carry-forward writes account for most of the run time. the
#   matrices are converted back to one tibble per agent (agents_to_tibble) before ncc_abm() returns, which is the
#   format summary_abm() and plot_abm() read.
#
# status is stored as a numeric code in the matrix and as character ('ALIVE' / 'DEAD') in the returned tibbles.
#   datetime and id are the same on every row of an agent, so they are left out of the matrix and added back when
#   the tibbles are rebuilt.

#1) numeric status codes used inside the agent matrices

.status_alive <- 1
.status_dead <- 0

#2) matrix column names: every per-agent tibble column except datetime and id. this order defines the matrix
#   layout only; the column order of the returned tibbles is set in agents_to_tibble()

.agent_cols <- c(
  "status", "x", "y", "heading", "step_length", "turn_angle",
  "bm", "ifbfat", "lean_mass", "fat_mass", "forage_consumed", "energy_i",
  "daily_intake", "energy_bmr", "energy_hif", "energy_loc", "energy_rep",
  "energy_net", "fat_change"
)

#' Convert agent tibbles to numeric state matrices
#'
#' @description Converts each per-agent tibble from \code{create_agents} into a numeric matrix of
#'   time-varying state for fast in-loop writes. The character \code{status} column is encoded as a
#'   numeric code (\code{.status_alive} / \code{.status_dead}, \code{NA} for unset rows); the
#'   constant \code{datetime} and \code{id} columns are dropped and restored on conversion back.
#'
#' @param agents A named list of per-agent tibbles.
#'
#' @return A named list of numeric matrices, one per agent, with columns \code{.agent_cols}.
#'
#' @keywords internal

agents_to_matrix <- function(agents) {

  # ----------------------------------------------------------------------------------------------------------------------
  # build one state matrix per agent
  # ----------------------------------------------------------------------------------------------------------------------

  lapply(agents, function(a) {

    #1) an empty numeric matrix with one row per time step and one named column per state variable

    m <- matrix(
      NA_real_,
      nrow = nrow(a),
      ncol = length(.agent_cols),
      dimnames = list(NULL, .agent_cols)
    )

    #2) encode status numerically: 'ALIVE' becomes .status_alive, 'DEAD' becomes .status_dead, and rows with no

    status_code <- rep(NA_real_, nrow(a))
    status_code[which(a$status == "ALIVE")] <- .status_alive
    status_code[which(a$status == "DEAD")] <- .status_dead
    m[, "status"] <- status_code

    #3) copy every remaining state column from the tibble into the matrix

    for (col in .agent_cols[-1]) {
      m[, col] <- a[[col]]
    }

    m
  })
}

#' Convert agent state matrices back to tibbles
#'
#' @description Rebuilds the per-agent output tibbles from the internal numeric matrices, restoring
#'   the constant \code{datetime} and \code{id} columns and decoding the numeric \code{status} code
#'   back to \code{"ALIVE"} / \code{"DEAD"} / \code{NA}. The columns and their order match the
#'   tibbles built by \code{create_agents}.
#'
#' @param agents_m A named list of per-agent state matrices.
#' @param times The hourly \code{POSIXct} sequence shared by all agents (the \code{datetime} column).
#'
#' @return A named list of per-agent tibbles.
#'
#' @keywords internal

agents_to_tibble <- function(agents_m, times) {

  # ----------------------------------------------------------------------------------------------------------------------
  # rebuild one tibble per agent
  # ----------------------------------------------------------------------------------------------------------------------

  #1) agent ids, taken from the names of the matrix list

  ids <- names(agents_m)

  #2) rebuild each agent's tibble from its matrix

  out <- lapply(seq_along(agents_m), function(i) {

    m <- agents_m[[i]]

    #3) decode the numeric status code back to 'ALIVE' / 'DEAD', leaving NA where no status was recorded

    status <- rep(NA_character_, nrow(m))
    status[which(m[, "status"] == .status_alive)] <- "ALIVE"
    status[which(m[, "status"] == .status_dead)] <- "DEAD"

    #4) assemble the tibble in create_agents() column order, adding back the shared hourly datetime sequence and the agent's id

    dplyr::tibble(
      datetime = times,
      id = ids[i],
      status = status,
      x = m[, "x"],
      y = m[, "y"],
      heading = m[, "heading"],
      step_length = m[, "step_length"],
      turn_angle = m[, "turn_angle"],
      bm = m[, "bm"],
      ifbfat = m[, "ifbfat"],
      lean_mass = m[, "lean_mass"],
      fat_mass = m[, "fat_mass"],
      forage_consumed = m[, "forage_consumed"],
      energy_i = m[, "energy_i"],
      daily_intake = m[, "daily_intake"],
      energy_bmr = m[, "energy_bmr"],
      energy_hif = m[, "energy_hif"],
      energy_loc = m[, "energy_loc"],
      energy_rep = m[, "energy_rep"],
      energy_net = m[, "energy_net"],
      fat_change = m[, "fat_change"]
    )
  })

  #5) name each tibble by its agent id

  names(out) <- ids
  out
}

#' Run the NCC Agent-Based Model
#'
#' @description Simulates one season of movement, foraging, energetics, and survival for a
#'   population of female bighorn sheep agents on a daily forage landscape. Agents are created by
#'   \code{create_agents} and moved through a burn-in under their own integrated step-selection
#'   functions (iSSFs). Each hour, every living agent then takes one movement step, with separate
#'   day and night selection coefficients, and during daylight eats from the cell it occupies,
#'   following a negative-exponential (type II) intake response to forage biomass. Grazing lowers
#'   realized biomass below the potential biomass in \code{forage_reference}; the shortfall is
#'   carried between days as a deficit and, when \code{forage_regrowth} is \code{TRUE}, recovers
#'   through geometric regrowth. At the end of each day, each living agent's intake energy minus
#'   its basal metabolism, heat increment of feeding, locomotion, and lactation costs sets its
#'   change in fat mass, which updates body mass and body fat. An agent dies when its fat mass
#'   reaches zero.
#'
#' @param forage_reference A multi-layer \code{terra::SpatRaster} of daily potential forage
#'   biomass (g/cell), one layer per day, with layer dates set through \code{terra::time()}.
#'   Layers are matched to the simulated days by date, so the raster may cover a longer window
#'   than the season. The first layer also defines the model grid: agent start cells, the cell
#'   lookups for movement and foraging, the grid extent, and the cell area.
#' @param dem A single-layer \code{terra::SpatRaster} of elevation (m), used to compute the slope
#'   of each step for locomotion costs.
#' @param canopy A single-layer \code{terra::SpatRaster} of canopy cover, used as a habitat
#'   covariate in movement. Values are looked up by the cell numbers of \code{forage_reference},
#'   so this raster must share its grid.
#' @param escape A single-layer \code{terra::SpatRaster} of distance to escape terrain (m), used
#'   as a habitat covariate in movement. Values are looked up by the cell numbers of
#'   \code{forage_reference}, so this raster must share its grid.
#' @param year Numeric. Calendar year the simulated season runs in, defaulting to 2024.
#'   Sets \code{t_start} and \code{t_end} from the stored season month-day bounds and
#'   refreshes the daylight schedule, then selects the \code{forage_reference} layers whose
#'   \code{terra::time()} dates match the simulated days. Point this at the year of the
#'   forage rasters being supplied.
#' @param forage_regrowth Logical. If \code{TRUE} (the default), depleted forage
#'   recovers each day through the geometric regrowth function
#'   (\code{update_forage}); if \code{FALSE}, geometric regrowth is disabled so
#'   grazed biomass is not recovered. The end-of-day deficit is still capped so
#'   realized biomass stays non-negative, and the seasonal phenology of potential
#'   biomass (the daily \code{forage_reference} layers) is unaffected either way.
#' @param verbose Logical. If \code{TRUE} (the default), progress messages are
#'   printed as the simulation runs; if \code{FALSE}, all messages are suppressed.
#' @param report_time Logical. If \code{TRUE}, the wall-clock run time and the
#'   \code{n_candidates} setting are printed when the run finishes; defaults to
#'   \code{FALSE}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{agents}{A named list of tibbles, one per agent, with one row per hourly time step
#'       holding the agent's status, position, movement, body condition, hourly intake, and
#'       daily energetics.}
#'     \item{agent_params}{A tibble of fixed individual parameters, one row per agent, as built
#'       by \code{create_agents}, with \code{j_post_partum} advanced through each agent's last
#'       day alive.}
#'   }
#'
#' @export

ncc_abm <- function(forage_reference, dem, canopy, escape, year = 2024, forage_regrowth = TRUE, verbose = TRUE, report_time = FALSE) {

  # ----------------------------------------------------------------------------------------------------------------------
  # point the clock at the requested year
  # ----------------------------------------------------------------------------------------------------------------------

  #1) rebuild t_start, t_end, and the daylight vectors for the requested year before anything reads them

  .set_season(year)

  # ----------------------------------------------------------------------------------------------------------------------
  # resolve time parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) start the run clock, read at the end for the optional run-time report

  start_time <- Sys.time()

  #2) pull the season start, season end, and time step from the global parameters

  t_start <- get_param("t_start")
  t_end <- get_param("t_end")
  t_delta <- get_param("t_delta")

  #3) build the hourly time sequence; each element is one row of every agent's state

  times <- seq(
    t_start,
    t_end,
    by = as.numeric(t_delta, units = "secs")
  )

  #4) build the daily date sequence in the America/Denver time zone; each element is one simulated day

  dates <- seq(
    as.Date(t_start, tz = "America/Denver"),
    as.Date(t_end, tz = "America/Denver"),
    by = "day"
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # assign input rasters and check rasters
  # ----------------------------------------------------------------------------------------------------------------------

  if (verbose) message("1) Importing and checking rasters")

  #1) select the forage_reference layers whose dates match the simulated days. matching on date rather than on
  #   layer count lets a raster covering a longer window serve a shorter season, and stops the run if a layer is
  #   missing for any simulated day

  reference_time <- terra::time(forage_reference)

  stopifnot(
    "forage_reference must carry layer dates (terra::time)" = !anyNA(reference_time)
  )

  layer_idx <- match(dates, as.Date(reference_time))

  stopifnot(
    "forage_reference has no layer for one or more simulated days" = !anyNA(layer_idx)
  )

  forage_reference <- forage_reference[[layer_idx]]

  #2) check that canopy, escape, and dem each have a single layer

  stopifnot(
    "canopy must have a single layer" = terra::nlyr(canopy) == 1,
    "escape must have a single layer" = terra::nlyr(escape) == 1,
    "dem must have a single layer" = terra::nlyr(dem) == 1
  )

  #3) load every raster into memory so each daily layer read and each per-step covariate lookup reads from RAM
  #   rather than disk. multiplying by 1 forces the values into memory but drops layer names and dates, so the
  #   forage dates and the escape and canopy layer names are restored immediately afterward

  reference_dates <- terra::time(forage_reference)

  forage_reference <- forage_reference * 1
  dem <- dem * 1
  canopy <- canopy * 1
  escape <- escape * 1

  terra::time(forage_reference) <- reference_dates
  names(escape) <- "escape_terrain"
  names(canopy) <- "canopy_cover"

  #4) copy the four rasters into the global environment. calc_energy_loc() reads the dem raster from there for its
  #   cell geometry, so this assignment is required. objects with these names in the global environment are
  #   overwritten

  list2env(
    list(
      forage_reference = forage_reference,
      dem = dem,
      canopy = canopy,
      escape = escape
    ),
    envir = .GlobalEnv
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # create agents
  # ----------------------------------------------------------------------------------------------------------------------

  if (verbose) message("2) Creating agents")

  #1) create the agents (one state tibble each) and their fixed individual parameters

  init <- create_agents(forage_reference)
  agents <- init$agents
  agent_params <- init$agent_params

  #2) convert each agent's tibble to a numeric state matrix for fast writes inside the simulation loop; the
  #   matrices are converted back to tibbles before the function returns

  agents <- agents_to_matrix(agents)

  # ----------------------------------------------------------------------------------------------------------------------
  # set up per-run state
  # ----------------------------------------------------------------------------------------------------------------------

  if (verbose) message("3) Simulating movement and foraging")

  #1) grazing deficit, daylight schedule, and elevation values. realized biomass in a cell on any day is that
  #   day's potential biomass minus the deficit. the deficit starts at zero in every cell (no grazing before the
  #   season), grows as agents eat, and is capped and decayed at the end of each day. the hourly daylight flags
  #   come from the season clock, and the dem values are extracted once so calc_energy_loc() can index
  #   elevations by cell number

  deficit <- rep(0, terra::ncell(forage_reference))
  is_daylight <- get_param("is_daylight")
  dem_vals <- terra::values(dem, mat = FALSE)

  #2) static lookups for movement and foraging: the escape and canopy values indexed by cell number, the first
  #   forage layer as the reference grid for cellFromXY(), the area of one cell (m^2, used to convert forage
  #   from g/cell to kg/ha), the grid extent (an agent whose candidate steps mostly leave it holds position),
  #   and the number of candidate steps per move. validate_move_coefs() then checks that the agents' iSSF
  #   carries exactly the twelve coefficients extracted by name below; the first agent is checked because all
  #   agents share the same term structure

  escape_vals <- terra::values(escape, mat = FALSE)
  canopy_vals <- terra::values(canopy, mat = FALSE)
  geom_ref <- forage_reference[[1]]
  cell_area <- prod(terra::res(forage_reference))
  move_ext <- as.vector(terra::ext(forage_reference))
  move_n_candidates <- get_param("n_candidates")
  validate_move_coefs(agent_params$issf[[1]]$coefficients)

  #3) per-agent movement parameters, stored as vectors across agents: the gamma step-length shape and scale,
  #   the von Mises turn-angle concentration (kappa) and mean (mu), and the six selection coefficients (step
  #   length, log step length, cosine of turn angle, forage biomass, distance to escape terrain, and canopy
  #   cover). prm_day holds the day coefficients. prm_night adds each coefficient's night offset (its
  #   :tod_end_night_end term) to the day value; the step-length and turn-angle distributions are the same by
  #   day and by night

  vec_coef_names <- c(
    "sl_",
    "log(sl_)",
    "cos(ta_)",
    "forage_biomass_end",
    "escape_terrain_end",
    "canopy_cover_end",
    "sl_:tod_end_night_end",
    "log(sl_):tod_end_night_end",
    "cos(ta_):tod_end_night_end",
    "forage_biomass_end:tod_end_night_end",
    "escape_terrain_end:tod_end_night_end",
    "canopy_cover_end:tod_end_night_end"
  )
  vec_coef <- t(vapply(agent_params$issf, function(m) m$coefficients[vec_coef_names], numeric(12L)))
  colnames(vec_coef) <- vec_coef_names

  prm_day <- list(
    shape = vapply(agent_params$issf, function(m) m$sl_$params$shape, numeric(1)),
    scale = vapply(agent_params$issf, function(m) m$sl_$params$scale, numeric(1)),
    kappa = vapply(agent_params$issf, function(m) m$ta_$params$kappa, numeric(1)),
    mu = lapply(agent_params$issf, function(m) m$ta_$params$mu),
    b_sl = vec_coef[, "sl_"],
    b_lsl = vec_coef[, "log(sl_)"],
    b_cta = vec_coef[, "cos(ta_)"],
    b_fo = vec_coef[, "forage_biomass_end"],
    b_es = vec_coef[, "escape_terrain_end"],
    b_ca = vec_coef[, "canopy_cover_end"]
  )
  prm_night <- prm_day
  prm_night$b_sl <- prm_day$b_sl + vec_coef[, "sl_:tod_end_night_end"]
  prm_night$b_lsl <- prm_day$b_lsl + vec_coef[, "log(sl_):tod_end_night_end"]
  prm_night$b_cta <- prm_day$b_cta + vec_coef[, "cos(ta_):tod_end_night_end"]
  prm_night$b_fo <- prm_day$b_fo + vec_coef[, "forage_biomass_end:tod_end_night_end"]
  prm_night$b_es <- prm_day$b_es + vec_coef[, "escape_terrain_end:tod_end_night_end"]
  prm_night$b_ca <- prm_day$b_ca + vec_coef[, "canopy_cover_end:tod_end_night_end"]

  #4) intake parameters for the forage phase: the intake asymptote per unit metabolic mass (intake_max), the
  #   biomass scale of the intake response (intake_decay), the reproductive-status multipliers (a length-2
  #   vector indexed by rep_status + 1, so single brackets pick one multiplier per agent), and the factor that
  #   converts consumed forage to metabolizable energy (digestible energy content times the DE-to-ME
  #   conversion factor)

  fg_intake_max <- get_param("intake_max")
  fg_intake_decay <- get_param("intake_decay")
  fg_intake_multiplier <- get_param("intake_multiplier")
  fg_energy_factor <- get_param("DE") * get_param("DE_to_ME_conversion_factor")

  #5) integer positions of the state columns in the agent matrices, resolved once so the hour loop writes by
  #   position rather than by name. body_idx groups the four body-condition columns, move_idx the five movement
  #   columns, forage_idx the two hourly intake columns, and energy_idx the seven daily energetics columns

  status_idx <- match("status", .agent_cols)
  x_idx <- match("x", .agent_cols)
  y_idx <- match("y", .agent_cols)
  heading_idx <- match("heading", .agent_cols)
  step_length_idx <- match("step_length", .agent_cols)
  turn_angle_idx <- match("turn_angle", .agent_cols)
  bm_idx <- match("bm", .agent_cols)
  ifbfat_idx <- match("ifbfat", .agent_cols)
  energy_i_idx <- match("energy_i", .agent_cols)
  body_idx <- match(c("bm", "ifbfat", "lean_mass", "fat_mass"), .agent_cols)
  move_names <- c("x", "y", "step_length", "turn_angle", "heading")
  move_idx <- match(move_names, .agent_cols)
  forage_idx <- match(c("forage_consumed", "energy_i"), .agent_cols)
  energy_idx <- match(
    c("energy_bmr", "energy_hif", "energy_loc", "energy_rep", "daily_intake", "energy_net", "fat_change"),
    .agent_cols
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # burn-in
  # ----------------------------------------------------------------------------------------------------------------------

  #1) burn-in movement data: day 1 potential biomass (no grazing has happened yet, so it is also the realized
  #   biomass), the static escape and canopy values, and the time of day fixed at day (tod = 0), all held
  #   constant for every burn-in step

  n_burn_in <- get_param("burn_in")

  if (verbose) message("   Burn-in: ", n_burn_in, " movement steps per agent")

  move_burn <- list(
    geom = geom_ref,
    forage = terra::values(forage_reference[[1]], mat = FALSE),
    escape = escape_vals,
    canopy = canopy_vals,
    tod = 0,
    ext = move_ext,
    n_candidates = move_n_candidates
  )

  #2) when burn_in is above zero, move every agent n_burn_in steps under its own iSSF and overwrite its first-row
  #   x, y, and heading with where it ends up. simulate_burn_in() advances all agents together one step at a
  #   time, with the same hold-position rules as the hourly movement. only the endpoint is kept; there is no
  #   foraging, forage depletion, or energetics during the burn-in

  if (n_burn_in > 0) {

    burn_time <- system.time(
      burn_end <- simulate_burn_in(
        x0 = vapply(agents, function(a) a[1, x_idx], numeric(1)),
        y0 = vapply(agents, function(a) a[1, y_idx], numeric(1)),
        heading0 = vapply(agents, function(a) a[1, heading_idx], numeric(1)),
        issf = agent_params$issf,
        move = move_burn,
        n_steps = n_burn_in
      )
    )

    for (i in seq_along(agents)) {
      agents[[i]][1, c(x_idx, y_idx, heading_idx)] <- burn_end[i, ]
    }

    #3) when report_time is TRUE, report the burn-in's elapsed time with its number of steps and agents

    if (report_time) {
      message(
        "Burn-in time: ", round(burn_time[["elapsed"]], 2), " sec  (",
        n_burn_in, " steps x ", length(agents), " agents)"
      )
    }
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # simulate daily and hourly dynamics
  # ----------------------------------------------------------------------------------------------------------------------

  #1) each agent's current x, y, and heading, held in vectors so the hour loop reads and updates positions
  #   without indexing into every agent's matrix. they start from row 1, which holds the burn-in endpoint when
  #   burn-in is on. each agent's matrix still receives every hourly position as the permanent record. entries
  #   for dead agents are no longer updated or read

  pos_x <- vapply(agents, function(a) a[1, x_idx], numeric(1))
  pos_y <- vapply(agents, function(a) a[1, y_idx], numeric(1))
  pos_h <- vapply(agents, function(a) a[1, heading_idx], numeric(1))

  #2) day 1 potential biomass (g/cell) for every cell. each later day's vector is carried over from the end-of-day
  #   forage update (step 17), so each forage_reference layer is extracted once per run

  potential_d <- terra::values(forage_reference[[1]], mat = FALSE)

  #3) loop over each simulated day

  for (d in seq_along(dates)) {

    #4) when verbose is TRUE, announce the current day every seventh day, starting with day 1

    if (verbose && (d - 1) %% 7 == 0) message("Day ", d, " of ", length(dates), ": ", dates[d])

    #5) this day's realized biomass (vals, g/cell): potential biomass minus the deficit carried from earlier days.
    #   vals is depleted in place as agents eat during the day

    vals <- potential_d - deficit

    #6) the number of daylight hourly steps in this day (forage_hours_d). each hourly intake is the daily intake
    #   response divided by this count. it is used only in daylight hours, so a day with no daylight steps never
    #   divides by it

    day_steps <- ((d - 1) * 24 + 1):min(d * 24, length(times))
    forage_hours_d <- sum(is_daylight[day_steps])

    #7) the movement data for this day. the list holds vals as it stands at the start of the day; because R copies
    #   vals when it is later modified, movement sees day-start forage all day while foraging depletes vals. the
    #   static escape, canopy, grid, and extent values come from the per-run setup. day and night differ in the
    #   selection coefficients (prm_day and prm_night), not in the movement data

    move_day <- list(
      geom = geom_ref,
      forage = vals,
      escape = escape_vals,
      canopy = canopy_vals,
      ext = move_ext,
      n_candidates = move_n_candidates
    )

    #8) the agents alive at the start of the day and their body masses. status and body mass change only in the
    #   end-of-day block (survival check and mass update), so both are computed once per day. status is read
    #   from row 1 on day 1 and from the previous day's last row after that, where the survival check is
    #   recorded; body mass is read from the day's first row, where the previous day's update was written

    day_status_row <- if (d == 1) 1L else (d - 1L) * 24L
    day_living <- which(
      vapply(
        agents,
        function(a) a[day_status_row, status_idx] == .status_alive,
        logical(1)
      )
    )
    day_bm_row <- (d - 1L) * 24L + 1L
    bm_day <- vapply(agents, function(a) a[day_bm_row, bm_idx], numeric(1))

    #9) loop over each hour within the day

    for (h in 1:24) {

      #10) the absolute hour index across the season (t). row 1 holds the initial state and is skipped, as are
      #    indices past the last time step (the season's final day can have fewer than 24 hourly rows)

      t <- (d - 1) * 24 + h
      if (t < 2 || t > length(times)) next

      #11) this hour's daylight flag

      is_day <- is_daylight[t]

      #12) shuffle the order in which living agents move and eat this hour. sample.int() permutes positions, so
      #    a single surviving agent is returned as itself (sample() on a length-1 numeric vector would instead
      #    draw from 1 to that value)

      living <- day_living[sample.int(length(day_living))]

      #13) copy each living agent's status from the previous row. on the first hour of a day this brings forward
      #    the survival check recorded on the previous day's last row; within a day, status does not change

      for (i in living) {
        agents[[i]][t, status_idx] <- agents[[i]][t - 1, status_idx]
      }

      #14) from the second hour of the day on, copy each living agent's body mass, body fat, lean mass, and fat
      #    mass from the previous row. the first hour of each day already holds these values: create_agents()
      #    sets them on day 1, and the previous day's energetics update (step 23) writes them afterward

      if (h > 1) {
        for (i in living) {
          agents[[i]][t, body_idx] <- agents[[i]][t - 1, body_idx]
        }
      }

      #15) move phase: every living agent takes one step. prm_sub holds the day or night movement parameters of
      #    the living agents, and simulate_move_step_vec() draws a step for all of them at once. an agent holds
      #    position (prior position and heading, zero step length) when more than half of its candidate steps
      #    leave the grid or none of its candidates has positive selection weight. the new positions update
      #    pos_x, pos_y, and pos_h and are written to each agent's matrix row for this hour

      if (length(living) > 0) {

        prm_t <- if (is_day) prm_day else prm_night
        prm_sub <- list(
          shape = prm_t$shape[living],
          scale = prm_t$scale[living],
          kappa = prm_t$kappa[living],
          mu = prm_t$mu[living],
          b_sl = prm_t$b_sl[living],
          b_lsl = prm_t$b_lsl[living],
          b_cta = prm_t$b_cta[living],
          b_fo = prm_t$b_fo[living],
          b_es = prm_t$b_es[living],
          b_ca = prm_t$b_ca[living]
        )

        mv_mat <- simulate_move_step_vec(
          pos_x[living],
          pos_y[living],
          pos_h[living],
          prm_sub,
          move_day
        )

        pos_x[living] <- mv_mat[, "x"]
        pos_y[living] <- mv_mat[, "y"]
        pos_h[living] <- mv_mat[, "heading"]

        for (k in seq_along(living)) {
          agents[[living[k]]][t, move_idx] <- mv_mat[k, move_names]
        }
      }

      #16) forage phase. at night, hourly intake and intake energy are recorded as zero. in daylight, each living
      #    agent eats from the cell it now occupies. forage density is converted from g/cell to kg/ha (divided by
      #    cell_area / 10); daily intake is intake_max * bm^0.75 * (1 - exp(-density / intake_decay)) times the
      #    agent's reproductive-status multiplier, using its day-start body mass; and the hour's consumption is
      #    that daily intake divided by the day's foraging hours, capped at the biomass left in the cell. agents
      #    alone in their cell are computed together. agents sharing a cell are computed one at a time in the
      #    shuffled order, each seeing the biomass left by the agents before it. consumption is removed from
      #    vals (positions with no cell are skipped), and intake energy is consumption times the metabolizable
      #    energy factor

      if (!is_day) {

        for (i in living) {
          agents[[i]][t, forage_idx] <- c(0, 0)
        }

      } else if (length(living) > 0) {

        fmult <- fg_intake_multiplier[agent_params$rep_status[living] + 1L]
        fbm <- bm_day[living]

        fcells <- terra::cellFromXY(geom_ref, cbind(pos_x[living], pos_y[living]))
        fdens <- vals[fcells]

        fdaily <- fg_intake_max * fbm^0.75 *
          (1 - exp(-(fdens / (cell_area / 10)) / fg_intake_decay)) *
          fmult
        fconsumed <- pmin(fdaily / forage_hours_d, fdens)

        contested <- duplicated(fcells) | duplicated(fcells, fromLast = TRUE)

        if (any(contested)) {

          # agents sharing a cell: recompute each one's consumption in the shuffled order, removing it from the
          #   cell before the next agent eats

          for (k in which(contested)) {
            if (is.na(fcells[k])) next
            d_k <- vals[fcells[k]]
            daily_k <- fg_intake_max * fbm[k]^0.75 *
              (1 - exp(-(d_k / (cell_area / 10)) / fg_intake_decay)) *
              fmult[k]
            c_k <- min(daily_k / forage_hours_d, d_k)
            fconsumed[k] <- c_k
            vals[fcells[k]] <- d_k - c_k
          }

          # agents alone in their cell: remove each one's consumption

          unc <- which(!contested & !is.na(fcells))
          vals[fcells[unc]] <- vals[fcells[unc]] - fconsumed[unc]

        } else {

          # no shared cells this hour: remove every agent's consumption at once

          ok <- which(!is.na(fcells))
          vals[fcells[ok]] <- vals[fcells[ok]] - fconsumed[ok]
        }

        fenergy <- fconsumed * fg_energy_factor

        for (k in seq_along(living)) {
          agents[[living[k]]][t, forage_idx] <- c(fconsumed[k], fenergy[k])
        }
      }
    }

    #17) end of day, forage: add the day's consumption to the deficit, then cap and decay it. day-start realized
    #    biomass was potential_d - deficit and vals is what remained after grazing, so their difference is each
    #    cell's consumption for the day. update_forage() caps the deficit at the next day's potential biomass,
    #    so realized biomass cannot go negative, and, when forage_regrowth is TRUE, shrinks it by one day of
    #    geometric regrowth. skipped on the final day because no later day uses the result

    if (d < length(dates)) {
      deficit <- deficit + ((potential_d - deficit) - vals)
      potential_next <- terra::values(forage_reference[[d + 1]], mat = FALSE)
      deficit <- update_forage(deficit, potential_next, regrowth = forage_regrowth)

      # the next day's potential biomass becomes potential_d for the next pass through the day loop

      potential_d <- potential_next
    }

    #18) end of day, energetics: daily energy terms are written to the day's last hourly row (t_day, which is the
    #    season's last row on the final day), and the updated body state goes to the next day's first row.
    #    day_living is refreshed from t_day, so the loop below covers every agent alive at the end of the day

    t_day <- d * 24
    if (t_day > length(times)) t_day <- length(times)

    day_rows <- ((d - 1) * 24 + 1):t_day

    day_living <- which(
      vapply(
        agents,
        function(a) a[t_day, status_idx] == .status_alive,
        logical(1)
      )
    )

    for (i in day_living) {

      bm_i <- agents[[i]][t_day, bm_idx]

      #19) the day's expenditure terms at the agent's end-of-day body mass: basal metabolism, heat increment of
      #    feeding, locomotion over the day's steps, and lactation (zero for non-reproductive agents). the
      #    locomotion window runs from the row before the day's first step to t_day; on day 1 it starts at
      #    row 1, so day 1 has one fewer step

      loc_start <- max(1, t_day - 24)

      energy_bmr <- calc_energy_bmr(bm_i)
      energy_hif <- calc_energy_hif(bm_i, t_day)
      energy_loc <- calc_energy_loc(
        agents[[i]][loc_start:t_day, x_idx],
        agents[[i]][loc_start:t_day, y_idx],
        agents[[i]][(loc_start + 1):t_day, step_length_idx],
        bm_i,
        dem_vals
      )
      energy_rep <- calc_energy_rep(
        energy_bmr,
        agent_params$j_post_partum[i],
        agent_params$rep_status[i]
      )

      #20) the day's net energy balance: the summed hourly intake energy minus the four expenditure terms.
      #    calc_energy_net() returns both the daily intake and the net balance

      net <- calc_energy_net(
        agents[[i]][day_rows, energy_i_idx],
        energy_bmr,
        energy_hif,
        energy_loc,
        energy_rep
      )

      #21) convert the net balance to a change in fat mass, then update body mass and body fat with lean mass held
      #    constant

      fat_change <- calc_fat_change(net[["energy_net"]])
      mass <- update_mass(bm_i, agents[[i]][t_day, ifbfat_idx], fat_change)

      #22) write the day's seven energetics values to the day's last row

      agents[[i]][t_day, energy_idx] <-
        c(energy_bmr, energy_hif, energy_loc, energy_rep, net[["daily_intake"]], net[["energy_net"]], fat_change)

      #23) write the updated body mass, body fat, lean mass, and fat mass to the next day's first row, where the
      #    next day's hours read them. skipped on the final day

      if (d < length(dates)) {
        t_next <- t_day + 1
        agents[[i]][t_next, body_idx] <-
          c(
            mass[["bm"]],
            mass[["ifbfat"]],
            calc_lean_mass(mass[["bm"]], mass[["ifbfat"]]),
            calc_fat_mass(mass[["bm"]], mass[["ifbfat"]])
          )
      }

      #24) survival and days post partum: update_status() marks the agent dead when its fat mass has reached zero
      #    and advances its days-post-partum counter by one day. the status is written to t_day, the row the next
      #    day's living set is read from, so a dead agent is excluded from all later processing

      upd <- update_status(
        calc_fat_mass(mass[["bm"]], mass[["ifbfat"]]),
        agent_params$j_post_partum[i]
      )
      agents[[i]][t_day, status_idx] <- if (upd$status == "ALIVE") .status_alive else .status_dead
      agent_params$j_post_partum[i] <- upd$j_post_partum
    }
  }

  # ----------------------------------------------------------------------------------------------------------------------
  # return
  # ----------------------------------------------------------------------------------------------------------------------

  #1) convert the agent matrices back to one tibble per agent

  agents <- agents_to_tibble(agents, times)

  #2) when report_time is TRUE, report the total run time (minutes) and the n_candidates setting

  if (report_time) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(
      "Run time: ", round(elapsed, 2), " min  (n_candidates = ",
      get_param("n_candidates"), ")"
    )
  }

  #3) return the agent tibbles and the fixed agent parameters

  list(agents = agents, agent_params = agent_params)
}
