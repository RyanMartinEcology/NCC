# ------------------------------------------------------------------------------------------------------------------------
# internal agent-state representation
# ------------------------------------------------------------------------------------------------------------------------

# during simulation each agent's time-varying state is held in a numeric matrix (one row per
#   hourly step, one column per state variable) rather than a tibble. single-cell writes into a
#   matrix are far cheaper than into a tbl_df, and the hourly move / forage / carry-forward writes
#   dominate the run. agents are converted back to tibbles before ncc_abm returns, so every
#   downstream consumer (summary_abm, plot_abm) is unaffected.
#
# status is character in the output but is encoded as a numeric code in the matrix; datetime and
#   id are constant down the rows and are restored on conversion back to tibbles.

.status_alive <- 1
.status_dead <- 0

# the matrix columns: every tibble column except datetime and id, with status encoded numerically.
#   ordering here is the matrix layout only; the output tibble column order is set in agents_to_tibble
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
  lapply(agents, function(a) {

    m <- matrix(
      NA_real_,
      nrow = nrow(a),
      ncol = length(.agent_cols),
      dimnames = list(NULL, .agent_cols)
    )

    status_code <- rep(NA_real_, nrow(a))
    status_code[which(a$status == "ALIVE")] <- .status_alive
    status_code[which(a$status == "DEAD")] <- .status_dead
    m[, "status"] <- status_code

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
#'   back to \code{"ALIVE"} / \code{"DEAD"} / \code{NA}. The column set and order match
#'   \code{create_agents} exactly, so the result is value-identical to running the model on tibbles.
#'
#' @param agents_m A named list of per-agent state matrices.
#' @param times The hourly \code{POSIXct} sequence shared by all agents (the \code{datetime} column).
#'
#' @return A named list of per-agent tibbles.
#'
#' @keywords internal
agents_to_tibble <- function(agents_m, times) {

  ids <- names(agents_m)

  out <- lapply(seq_along(agents_m), function(i) {

    m <- agents_m[[i]]

    status <- rep(NA_character_, nrow(m))
    status[which(m[, "status"] == .status_alive)] <- "ALIVE"
    status[which(m[, "status"] == .status_dead)] <- "DEAD"

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

  names(out) <- ids
  out
}

#' Run the NCC Agent-Based Model
#'
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
#' @export
ncc_abm <- function(forage_reference, dem, canopy, escape, year = 2024, forage_regrowth = TRUE, verbose = TRUE, report_time = FALSE) {

  # ------------------------------------------------------------------------------------------------------------------------
  # point the clock at the requested year
  # ------------------------------------------------------------------------------------------------------------------------

  #1) rebuild t_start / t_end and the daylight vectors together, before anything reads them

  .set_season(year)

  # ------------------------------------------------------------------------------------------------------------------------
  # resolve time parameters
  # ------------------------------------------------------------------------------------------------------------------------

  #1) start the run clock (read at the end for the optional run-time report)

  start_time <- Sys.time()

  #2) pull time parameters from globals

  t_start <- get_param("t_start")
  t_end <- get_param("t_end")
  t_delta <- get_param("t_delta")

  #3) build hourly time sequence

  times <- seq(
    t_start,
    t_end,
    by = as.numeric(t_delta, units = "secs")
  )

  #4) build daily date sequence in America/Denver

  dates <- seq(
    as.Date(t_start, tz = "America/Denver"),
    as.Date(t_end, tz = "America/Denver"),
    by = "day"
  )

  # ------------------------------------------------------------------------------------------------------------------------
  # assign input rasters and check rasters
  # ------------------------------------------------------------------------------------------------------------------------

  if (verbose) message("1) Importing and checking rasters")

  #1) select the forage_reference layers whose dates match the simulated days. matching on date
  #   rather than counting layers lets a raster covering a longer window serve a shorter season,
  #   and catches a raster that is misaligned rather than merely the wrong length

  reference_time <- terra::time(forage_reference)

  stopifnot(
    "forage_reference must carry layer dates (terra::time)" = !anyNA(reference_time)
  )

  layer_idx <- match(dates, as.Date(reference_time))

  stopifnot(
    "forage_reference has no layer for one or more simulated days" = !anyNA(layer_idx)
  )

  forage_reference <- forage_reference[[layer_idx]]

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

  init <- create_agents(forage_reference)
  agents <- init$agents
  agent_params <- init$agent_params

  #2) convert each agent to a numeric state matrix for fast in-loop writes; the matrices are
  #   converted back to tibbles just before the function returns

  agents <- agents_to_matrix(agents)

  # ------------------------------------------------------------------------------------------------------------------------
  # simulate daily and hourly dynamics
  # ------------------------------------------------------------------------------------------------------------------------

  if (verbose) message("3) Simulating movement and foraging")

  #1) initialize per-run state hoisted out of the loop: the grazing deficit, the
  #   daylight schedule, the dem values (for calc_energy_loc), and the movement
  #   lookups — the static escape/canopy covariate vectors, a reference layer for
  #   cellFromXY geometry, the extent for the outside guard, and the t_delta /
  #   n_candidates constants. the agents' shared iSSF term structure is validated
  #   once here, since simulate_move's hand-rolled predictor assumes it. the matrix column
  #   indices used in the hour loop are resolved just below
  # realized biomass on any day is potential minus this deficit; grazing adds to it
  #   and it decays each day toward zero. starts at zero (ungrazed at season start)

  deficit <- rep(0, terra::ncell(forage_reference))
  is_daylight <- get_param("is_daylight")
  dem_vals <- terra::values(dem, mat = FALSE)

  escape_vals <- terra::values(escape, mat = FALSE)
  canopy_vals <- terra::values(canopy, mat = FALSE)
  geom_ref <- forage_reference[[1]]
  cell_area <- prod(terra::res(forage_reference))
  move_ext <- as.vector(terra::ext(forage_reference))
  move_t_delta <- get_param("t_delta")
  move_n_candidates <- get_param("n_candidates")
  validate_move_coefs(agent_params$issf[[1]]$coefficients)

  # per-agent movement parameters lifted into vectors once, mirroring simulate_burn_in, with the
  #   night coefficient offsets folded into a second parameter set so the hour loop just picks day
  #   or night. consumed by simulate_move_step_vec in the vectorized move phase below
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

  # intake parameters hoisted once for the vectorized forage phase. the multiplier is kept as the
  #   full length-2 vector and indexed per agent with [ rather than [[
  fg_intake_max <- get_param("intake_max")
  fg_intake_decay <- get_param("intake_decay")
  fg_intake_multiplier <- get_param("intake_multiplier")
  fg_energy_factor <- get_param("DE") * get_param("DE_to_ME_conversion_factor")

  #2) column indices into the agent matrices, resolved once. integer indexing drops the per-write
  #   column-name lookup in the hour loop and removes any chance of a name-subscript error
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

  # ------------------------------------------------------------------------------------------------------------------------
  # burn-in
  # ------------------------------------------------------------------------------------------------------------------------

  #1) the burn-in movement data: day 1 potential biomass, held frozen for every burn-in step, with
  #   tod fixed at day. no grazing has occurred yet, so potential and realized biomass are identical

  n_burn_in <- get_param("burn_in")

  if (verbose) message("   Burn-in: ", n_burn_in, " movement steps per agent")

  move_burn <- list(
    geom = geom_ref,
    forage = terra::values(forage_reference[[1]], mat = FALSE),
    escape = escape_vals,
    canopy = canopy_vals,
    tod = 0,
    ext = move_ext,
    t_delta = move_t_delta,
    n_candidates = move_n_candidates
  )

  #2) walk every agent forward n_burn_in steps under its own iSSF and overwrite its first-row
  #   position and heading with the endpoint. nothing else is updated — no foraging, no depletion,
  #   no energetics — and the intermediate path is discarded. simulate_burn_in advances all agents
  #   together one step at a time, applying the same hold-position guards as the hour loop

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

    #3) report the burn-in's own elapsed time when report_time is on, so its cost is measured
    #   directly rather than inferred by differencing two whole-run times

    if (report_time) {
      message(
        "Burn-in time: ", round(burn_time[["elapsed"]], 2), " sec  (",
        n_burn_in, " steps x ", length(agents), " agents)"
      )
    }
  }

  # ------------------------------------------------------------------------------------------------------------------------
  # simulate daily and hourly dynamics, continued
  # ------------------------------------------------------------------------------------------------------------------------

  #3) persistent position state: every agent's current x, y, and heading held as plain vectors so
  #   the hour loop reads and advances positions without per-agent matrix lookups. initialized from
  #   the first row (which the burn-in has already overwritten); the agent matrices still receive
  #   every hourly write as the permanent record. dead agents' entries go stale and are never read

  pos_x <- vapply(agents, function(a) a[1, x_idx], numeric(1))
  pos_y <- vapply(agents, function(a) a[1, y_idx], numeric(1))
  pos_h <- vapply(agents, function(a) a[1, heading_idx], numeric(1))

  #4) day 1's potential biomass, extracted once here; each subsequent day's vector is carried over
  #   from the end-of-day deficit update below, so every layer is extracted exactly once per run

  potential_d <- terra::values(forage_reference[[1]], mat = FALSE)

  #5) loop over each simulation day

  for (d in seq_along(dates)) {

    #6) announce the current day

    if (verbose && (d - 1) %% 7 == 0) message("Day ", d, " of ", length(dates), ": ", dates[d])

    #7) build this day's realized biomass: potential minus the carried deficit
    # forage_reference holds potential biomass; realized is potential_d - deficit.
    #   vals is the day-start realized biomass and depletes within the day as agents
    #   forage. potential_d was extracted at loop entry (day 1) or carried from the
    #   previous day's deficit update

    vals <- potential_d - deficit

    # the day's foraging-hour count, used by calc_dmi to spread the daily intake response across
    #   the hours an agent can eat. counting the daylight steps rather than using the continuous
    #   day_length keeps the realized daily total equal to the curve. a day with no daylight steps
    #   never reaches calc_dmi, since simulate_forage returns zero at night before it is called

    day_steps <- ((d - 1) * 24 + 1):min(d * 24, length(times))
    forage_hours_d <- sum(is_daylight[day_steps])

    #6) prebuild the day and night movement-data lists once for the day
    # movement sees day-start realized forage: the lists capture vals here, and the
    #   hour loop's in-place depletion of vals copy-on-writes, so the captured forage
    #   stays frozen all day. day and night differ only in tod (0 vs 1); the static
    #   escape/canopy/geometry/extent come from the per-run setup above

    move_day <- list(
      geom = geom_ref,
      forage = vals,
      escape = escape_vals,
      canopy = canopy_vals,
      tod = 0,
      ext = move_ext,
      t_delta = move_t_delta,
      n_candidates = move_n_candidates
    )
    move_night <- move_day
    move_night$tod <- 1

    # the day's living set and body masses, computed once: status and body mass change only at day
    #   boundaries (the survival check and mass propagation both run in the end-of-day block), so
    #   scanning every agent every hour rediscovers the same answer 24 times. day 1 reads row 1;
    #   later days read the previous day's last row, where the survival check was recorded

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

    #7) loop over each hour within the day

    for (h in 1:24) {

      #8) absolute hour index across the season; skip the initialized first row

      t <- (d - 1) * 24 + h
      if (t < 2 || t > length(times)) next

      time_t <- times[t]

      #9) select the hour's movement data by daylight

      is_day <- is_daylight[t]
      move_t <- if (is_day) move_day else move_night

      #10) shuffle the day's living set into this hour's processing order
      # permute the vector's own elements: sample(x) on a length-1 x would draw from seq_len(x)
      #   instead of returning x, resurrecting dead agents whenever one animal is the sole survivor

      living <- day_living[sample.int(length(day_living))]

      #11) carry status forward for living agents
      # at the first hour of a day, status comes from the previous day's last row
      #   (where the survival check was recorded); within a day it is constant

      for (i in living) {
        agents[[i]][t, status_idx] <- agents[[i]][t - 1, status_idx]
      }

      #12) carry body state forward within the day
      # the first hour of each day already holds the body state: set by
      #   create_agents on day 1, and by the previous day's energetics
      #   propagation thereafter, so this runs only from the second hour onward

      if (h > 1) {
        for (i in living) {
          agents[[i]][t, body_idx] <- agents[[i]][t - 1, body_idx]
        }
      }

      #13) move phase — all living agents advance in one vectorized step. hold-position guards are
      #   applied inside simulate_move_step_vec, which returns held agents with their prior
      #   position and heading and a zero-length step. current positions come from the persistent
      #   pos_ vectors; the agent matrices receive the same writes as before as the record

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
          move_t
        )

        pos_x[living] <- mv_mat[, "x"]
        pos_y[living] <- mv_mat[, "y"]
        pos_h[living] <- mv_mat[, "heading"]

        for (k in seq_along(living)) {
          agents[[living[k]]][t, move_idx] <- mv_mat[k, move_names]
        }
      }

      #14) forage phase — every living agent forages in shuffled order, depleting vals. batched:
      #   night writes zeros; day batches the cell lookup and the intake curve across all living
      #   agents, reading positions straight from the move step's output. foraging consumes no RNG
      #   and the arithmetic is elementwise identical to calc_dmi (vector [ for the multiplier,
      #   pmin for the biomass floor), so uncontested cells reproduce the per-agent loop
      #   bit-for-bit. agents sharing a cell this hour see sequential depletion, so those are
      #   recomputed one at a time in the same shuffled order

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

          # sequential recompute for co-grazing agents, in the shuffled living order the loop used
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

          unc <- which(!contested & !is.na(fcells))
          vals[fcells[unc]] <- vals[fcells[unc]] - fconsumed[unc]

        } else {
          ok <- which(!is.na(fcells))
          vals[fcells[ok]] <- vals[fcells[ok]] - fconsumed[ok]
        }

        fenergy <- fconsumed * fg_energy_factor

        for (k in seq_along(living)) {
          agents[[living[k]]][t, forage_idx] <- c(fconsumed[k], fenergy[k])
        }
      }
    }

    #15) end of day: add the day's consumption to the deficit, then cap and decay it
    # day-start realized was potential_d - deficit; vals is what remains after
    #   grazing, so (potential_d - deficit) - vals is the day's consumption per
    #   cell. add it to the deficit, then cap at the next day's potential and decay.
    #   skipped on the last day (no next day uses the result)

    if (d < length(dates)) {
      deficit <- deficit + ((potential_d - deficit) - vals)
      potential_next <- terra::values(forage_reference[[d + 1]], mat = FALSE)
      deficit <- update_forage(deficit, potential_next, regrowth = forage_regrowth)
      # carry tomorrow's potential forward as the next iteration's potential_d, so each layer is
      #   extracted exactly once per run instead of twice
      potential_d <- potential_next
    }

    #16) end of day: update each living agent's energy balance, mass, and survival
    # daily values are written to the day's last hourly row (t_day); new body mass
    #   and body fat are propagated to the first row of the next day so day d+1
    #   reads the updated state

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

      #1) daily energy expenditure terms
      # the locomotion window spans the day's position rows and the steps between
      #   them; on day 1 it starts at the initial row (row 1), giving one fewer step

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

      #2) daily net energy balance from the day's hourly intake

      net <- calc_energy_net(
        agents[[i]][day_rows, energy_i_idx],
        energy_bmr,
        energy_hif,
        energy_loc,
        energy_rep
      )

      #3) fat-mass change and propagated body mass and body fat

      fat_change <- calc_fat_change(net[["energy_net"]])
      mass <- update_mass(bm_i, agents[[i]][t_day, ifbfat_idx], fat_change)

      #4) write the day's energetics to the day's last row

      agents[[i]][t_day, energy_idx] <-
        c(energy_bmr, energy_hif, energy_loc, energy_rep, net[["daily_intake"]], net[["energy_net"]], fat_change)

      #5) propagate the updated body state to the next day's first row

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

      #6) survival and days-post-partum: mark dead if fat reserves are exhausted,
      #   and advance the days-post-partum counter by one day
      # death is recorded at t_day so the next day's living filter (which reads the
      #   previous row) excludes the agent from all further processing

      upd <- update_status(
        calc_fat_mass(mass[["bm"]], mass[["ifbfat"]]),
        agent_params$j_post_partum[i]
      )
      agents[[i]][t_day, status_idx] <- if (upd$status == "ALIVE") .status_alive else .status_dead
      agent_params$j_post_partum[i] <- upd$j_post_partum
    }
  }

  # ------------------------------------------------------------------------------------------------------------------------
  # return
  # ------------------------------------------------------------------------------------------------------------------------

  #1) convert agents back to per-agent tibbles, then return with final parameters

  agents <- agents_to_tibble(agents, times)

  #2) optional run-time report: elapsed wall clock and the n_candidates setting

  if (report_time) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(
      "Run time: ", round(elapsed, 2), " min  (n_candidates = ",
      get_param("n_candidates"), ")"
    )
  }

  list(agents = agents, agent_params = agent_params)
}
