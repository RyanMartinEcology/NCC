#' Run the NCC Agent-Based Model
#' @description Runs the NCC agent-based model from t_start to t_end,
#'   iterating by t_delta. Agents move hourly and forage consumption is
#'   accumulated across 24 hourly steps before daily state updates are applied.
#'   Agent order is randomized at each hourly step. Forage depletion carries
#'   forward across days with seasonal decline applied multiplicatively.
#' @param forage_raster A \code{terra::SpatRaster} of forage availability in
#'   grams per cell, with one layer per simulation day. Generated outside the
#'   loop via \code{simulate_forage_raster()} and \code{convert_forage_units()}.
#' @param agents Optional. A named list of agent tibbles from a previous run.
#'   If NULL, agents are initialized from scratch.
#' @param start_time Optional. POSIXct. Time to resume from. If NULL, starts
#'   from t_start.
#' @param n_days Optional. Integer. Number of days to run the simulation. If
#'   NULL, runs the full simulation from t_start to t_end.
#' @param keep_time Logical. If TRUE, prints elapsed time alongside the day
#'   completion message. Default FALSE.
#' @param return_forage Logical. If TRUE, returns a named list with elements
#'   \code{agents} and \code{forage_raster} (the depleted forage raster after
#'   simulation). If FALSE, returns only the agents list. Default FALSE.
#' @return If \code{return_forage = FALSE}, a named list of agent tibbles. If
#'   \code{return_forage = TRUE}, a named list with elements \code{agents} and
#'   \code{forage_raster}.
#' @importFrom lubridate hour
#' @importFrom terra ext res ncol values setValues
#' @export
ncc_abm <- function(forage_raster, agents = NULL, start_time = NULL,
                    n_days = NULL, keep_time = FALSE, return_forage = FALSE) {

  t_start <- get_param("t_start")
  t_end   <- get_param("t_end")
  t_delta <- get_param("t_delta")
  times   <- seq(t_start, t_end, by = as.numeric(t_delta, units = "secs"))
  dates   <- seq(as.Date(t_start), as.Date(t_end), by = "day")

  # initialize agents if not resuming
  if (is.null(agents)) {
    agents <- create_agents(forage_raster)
  }

  # determine starting time step
  if (is.null(start_time)) {
    t_index <- 2L
  } else {
    t_index <- which(times == start_time)
    if (length(t_index) == 0) stop("start_time does not match a simulation time step.")
  }

  # determine ending time step
  if (!is.null(n_days)) {
    t_end_index <- min(t_index + (n_days * 24L) - 1L, length(times))
  } else {
    t_end_index <- length(times)
  }

  # pre-extract raster geometry — used in simulate_step() every hour
  e       <- terra::ext(forage_raster)
  ext_vec <- c(e$xmin, e$xmax, e$ymin, e$ymax)
  res_vec <- terra::res(forage_raster)
  ncols   <- terra::ncol(forage_raster)

  # pre-compute daily means from original raster for temporal decline ratios
  daily_means <- sapply(seq_len(terra::nlyr(forage_raster)), function(d) {
    mean(as.vector(terra::values(forage_raster[[d]])), na.rm = TRUE)
  })

  # initialize daily forage accumulator for each agent
  daily_forage <- rep(0, length(agents))

  # record start time for elapsed time tracking
  if (keep_time) run_start <- proc.time()

  # track current day to avoid redundant layer extraction
  current_day_layer <- NA_integer_
  cell_vals         <- NULL

  # main simulation loop
  for (t in t_index:t_end_index) {

    current_time <- times[t]
    current_date <- as.Date(current_time)
    day_layer    <- which(dates == current_date)

    # extract or update cell values at the start of each new day
    if (!identical(day_layer, current_day_layer)) {
      if (is.na(current_day_layer)) {
        # first day — extract fresh from raster
        cell_vals <- as.vector(terra::values(forage_raster[[day_layer]]))
      } else {
        # subsequent days — apply temporal decline ratio to depleted cell_vals
        temporal_ratio <- daily_means[day_layer] / daily_means[current_day_layer]
        cell_vals      <- pmax(cell_vals * temporal_ratio, 0)
      }
      current_day_layer <- day_layer
    }

    # reset daily forage accumulator at start of each new day
    if ((t - 1L) %% 24L == 1L) {
      daily_forage <- rep(0, length(agents))
    }

    # randomize agent order for this step
    agent_order <- sample(seq_along(agents))

    # movement phase — all agents move in randomized order
    for (i in agent_order) {

      # skip dead agents — look back to last valid status
      prev_status <- t - 24L
      if (prev_status >= 1L &&
          !is.na(agents[[i]]$status[prev_status]) &&
          agents[[i]]$status[prev_status] == "DEAD") next

      move_result <- simulate_move(
        x       = agents[[i]]$x[t - 1L],
        y       = agents[[i]]$y[t - 1L],
        ext_vec = ext_vec,
        res_vec = res_vec
      )

      agents[[i]]$x[t] <- move_result$x
      agents[[i]]$y[t] <- move_result$y
    }

    # foraging phase — agents forage in the same randomized order
    for (i in agent_order) {

      # skip dead agents
      prev_status <- t - 24L
      if (prev_status >= 1L &&
          !is.na(agents[[i]]$status[prev_status]) &&
          agents[[i]]$status[prev_status] == "DEAD") next

      forage_result <- simulate_forage(
        x         = agents[[i]]$x[t],
        y         = agents[[i]]$y[t],
        cell_vals = cell_vals,
        ext_vec   = ext_vec,
        res_vec   = res_vec,
        ncols     = ncols
      )

      cell_vals       <- forage_result$cell_vals
      daily_forage[i] <- daily_forage[i] + forage_result$forage_consumed
    }

    # run daily state update at the end of each day
    if ((t - 1L) %% 24L == 0L) {

      # write depleted cell values back to raster once per day
      forage_raster[[day_layer]] <- terra::setValues(forage_raster[[day_layer]], cell_vals)

      # write accumulated daily forage consumed to agent tibble
      for (i in seq_along(agents)) {
        agents[[i]]$forage_consumed[t] <- daily_forage[i]
      }

      # print progress message
      if (keep_time) {
        elapsed     <- proc.time() - run_start
        elapsed_sec <- elapsed["elapsed"]
        elapsed_str <- sprintf("%02d:%02d:%02d",
                               floor(elapsed_sec / 3600),
                               floor((elapsed_sec %% 3600) / 60),
                               floor(elapsed_sec %% 60))
        cat("Completed:", format(current_date, "%B %d"),
            "| Elapsed:", elapsed_str, "\n")
        cat("Completed:", format(current_date, "%B %d"), "\n")
      }

      agents <- update_agents_test(agents, t)
    }
  }

  if (return_forage) {
    list(agents = agents, forage_raster = forage_raster)
  } else {
    agents
  }
}

#' Run the NCC Agent-Based Model with population-level replication and carrying capacity detection
#' @description Wrapper around \code{ncc_abm()} that runs \code{replicates}
#'   iterations of the model at a given agent population size, calculates the
#'   global mean end-of-season ifbfat across replicates, and incrementally adds
#'   \code{delta_n} agents until the global mean falls below
#'   \code{carrying_capacity} or \code{n_max} agents are reached. Dead agents
#'   contribute 0 to the mean. Replicates can optionally be parallelized using
#'   the future backend.
#' @param parallel Logical. If TRUE, replicates are run in parallel using the
#'   current \code{future} plan. Default FALSE.
#' @param find_carrying_capacity Logical. If TRUE, stops when the global mean
#'   ifbfat falls below \code{carrying_capacity}. If FALSE, runs all agent
#'   levels from \code{n_agents} to \code{n_max} regardless of threshold.
#'   Default TRUE.
#' @param ... Additional arguments passed to \code{ncc_abm()}.
#' @return A named list with elements:
#'   \describe{
#'     \item{results}{A tibble with columns \code{n_agents}, \code{replicate},
#'       and \code{mean_ifbfat} for each run.}
#'     \item{global_means}{A tibble with columns \code{n_agents} and
#'       \code{global_mean_ifbfat} for each agent level.}
#'     \item{carrying_capacity_n}{Integer. The number of agents at which the
#'       global mean first fell below \code{carrying_capacity}, or \code{NA}
#'       if not reached or not searched for.}
#'   }
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom future plan
#' @export
ncc_replicator <- function(parallel = FALSE, find_carrying_capacity = TRUE, ...) {

  replicates        <- as.integer(get_param("replicates"))
  carrying_capacity <- as.numeric(get_param("carrying_capacity"))
  delta_n           <- as.integer(get_param("delta_n"))
  n_max             <- as.integer(get_param("n_max"))
  n_start           <- as.integer(get_param("n_agents"))

  results      <- dplyr::tibble(n_agents    = integer(),
                                replicate   = integer(),
                                mean_ifbfat = numeric())
  global_means <- dplyr::tibble(n_agents           = integer(),
                                global_mean_ifbfat = numeric())

  carrying_capacity_n <- NA_integer_
  n_current           <- n_start

  # capture additional arguments to pass to ncc_abm
  abm_args <- list(...)

  # define single replicate function
  run_replicate <- function(i) {
    forage_i <- simulate_forage_raster()
    forage_i <- convert_forage_units(forage_i)

    agents_i <- do.call(ncc_abm, c(list(forage_raster = forage_i), abm_args))

    end_ifbfat <- sapply(agents_i, function(a) {
      last_val <- tail(a$ifbfat[!is.na(a$ifbfat)], 1)
      if (length(last_val) == 0 || tail(a$status[!is.na(a$status)], 1) == "DEAD") {
        0
      } else {
        last_val
      }
    })

    dplyr::tibble(
      replicate   = i,
      mean_ifbfat = mean(end_ifbfat)
    )
  }

  # start global timer
  global_start <- proc.time()

  repeat {

    cat(crayon::cyan(paste0("Running ", replicates, " replicates with ",
                            n_current, " agents\n")))

    set_param("n_agents", n_current)

    if (parallel) {
      replicate_results <- furrr::future_map_dfr(
        seq_len(replicates),
        run_replicate,
        .options = furrr::furrr_options(seed = TRUE)
      )
    } else {
      replicate_results <- purrr::map_dfr(seq_len(replicates), run_replicate)
    }

    replicate_results <- dplyr::mutate(replicate_results, n_agents = n_current)

    results <- dplyr::bind_rows(results, replicate_results)

    global_mean <- mean(replicate_results$mean_ifbfat)

    global_means <- dplyr::bind_rows(global_means, dplyr::tibble(
      n_agents           = n_current,
      global_mean_ifbfat = global_mean
    ))

    # compute cumulative elapsed time
    global_elapsed     <- proc.time() - global_start
    global_elapsed_sec <- global_elapsed["elapsed"]
    global_elapsed_str <- sprintf("%02d:%02d:%02d",
                                  floor(global_elapsed_sec / 3600),
                                  floor((global_elapsed_sec %% 3600) / 60),
                                  floor(global_elapsed_sec %% 60))

    cat(crayon::green(paste0("Global mean ifbfat at n = ", n_current, ": ",
                             round(global_mean, 4),
                             " | Cumulative elapsed: ", global_elapsed_str, "\n")))

    # check carrying capacity threshold only if find_carrying_capacity is TRUE
    if (find_carrying_capacity && global_mean < carrying_capacity) {
      carrying_capacity_n <- n_current
      cat(crayon::red(paste0("Carrying capacity reached at n = ",
                             carrying_capacity_n, "\n")))
      break
    }

    n_current <- n_current + delta_n

    if (n_current > n_max) {
      if (find_carrying_capacity) {
        cat(crayon::yellow("n_max reached without meeting carrying capacity threshold\n"))
      } else {
        cat(crayon::yellow("n_max reached\n"))
      }
      break
    }
  }

  list(
    results             = results,
    global_means        = global_means,
    carrying_capacity_n = carrying_capacity_n
  )
}
