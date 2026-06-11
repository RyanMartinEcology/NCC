#' Summarize an NCC agent-based model run
#'
#' @description Produces agent-level summary statistics from an \code{ncc_abm}
#'   result: survival, ending body condition, and monthly means of net energy,
#'   dry matter intake, and distance moved, each broken out for the whole
#'   population and by reproductive status. The returned object has a
#'   \code{print} method that formats these as a readable summary.
#'
#' @param result A list returned by \code{ncc_abm}, with elements \code{agents}
#'   (a named list of per-agent hourly tibbles) and \code{agent_params} (fixed
#'   per-agent parameters, including \code{rep_status}).
#'
#' @return An object of class \code{summary_abm}: a list with elements
#'   \code{n_agents}, \code{n_survived}, \code{prop_survived}, \code{ifbfat}
#'   (named vector: all, rep, non_rep), and \code{net_energy}, \code{dmi},
#'   \code{distance} (each a months-by-group matrix).
#'
#' @export
summary_abm <- function(result) {

  agents <- result$agents
  agent_params <- result$agent_params

  # ----------------------------------------------------------------------------------------------------------------------
  # assemble a long agent-day table
  # ----------------------------------------------------------------------------------------------------------------------

  #1) reproductive status per agent, aligned to the agents list order

  rep_status <- agent_params$rep_status[match(names(agents), agent_params$id)]

  #2) collapse each agent's hourly tibble to one row per day: daily values (energy_net,
  #   daily_intake) sit on the day's last row, distance is the day's summed hourly step
  #   length, and a day counts only while the agent is alive

  day_table <- do.call(rbind, lapply(seq_along(agents), function(i) {

    a <- agents[[i]]
    day <- as.Date(a$datetime)
    alive <- a$status == "ALIVE" & !is.na(a$status)

    dist_day <- tapply(a$step_length[alive], day[alive], sum)
    net_day <- tapply(a$energy_net[alive], day[alive], function(x) x[!is.na(x)][1])
    dmi_day <- tapply(a$forage_consumed[alive], day[alive], sum)

    days <- as.Date(names(dist_day))

    data.frame(
      month = factor(
        format(days, "%b"),
        levels = c("Jul", "Aug", "Sep", "Oct")
      ),
      rep_status = rep_status[i],
      net_energy = as.numeric(net_day),
      dmi = as.numeric(dmi_day),
      distance = as.numeric(dist_day)
    )
  }))

  # ----------------------------------------------------------------------------------------------------------------------
  # survival and ending body condition
  # ----------------------------------------------------------------------------------------------------------------------

  #1) survival counts

  n_agents <- length(agents)
  n_survived <- sum(vapply(
    agents,
    function(a) {
      alive <- which(a$status == "ALIVE")
      length(alive) > 0 && a$status[max(alive)] == "ALIVE"
    },
    logical(1)
  ))

  #2) ending ifbfat per agent (last recorded value), dead animals included

  end_ifbfat <- vapply(
    agents,
    function(a) {
      recorded <- which(!is.na(a$ifbfat))
      a$ifbfat[max(recorded)]
    },
    numeric(1)
  )

  ifbfat <- c(
    all = mean(end_ifbfat),
    rep = mean(end_ifbfat[rep_status == 1]),
    non_rep = mean(end_ifbfat[rep_status == 0])
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # monthly means by group
  # ----------------------------------------------------------------------------------------------------------------------

  #1) mean of a metric over living agent-days, by month and group

  monthly <- function(metric) {
    all_grp <- tapply(
      day_table[[metric]],
      day_table$month,
      mean,
      na.rm = TRUE
    )
    rep_grp <- tapply(
      day_table[[metric]][day_table$rep_status == 1],
      day_table$month[day_table$rep_status == 1],
      mean,
      na.rm = TRUE
    )
    non_grp <- tapply(
      day_table[[metric]][day_table$rep_status == 0],
      day_table$month[day_table$rep_status == 0],
      mean,
      na.rm = TRUE
    )
    out <- cbind(all = all_grp, rep = rep_grp, non_rep = non_grp)
    out[c("Jul", "Aug", "Sep", "Oct"), , drop = FALSE]
  }

  #2) build the three monthly tables

  out <- list(
    n_agents = n_agents,
    n_survived = n_survived,
    prop_survived = n_survived / n_agents,
    ifbfat = ifbfat,
    net_energy = monthly("net_energy"),
    dmi = monthly("dmi"),
    distance = monthly("distance")
  )

  class(out) <- "summary_abm"
  out
}

#' Print an NCC agent-based model summary
#'
#' @description Formats a \code{summary_abm} object as a readable console summary:
#'   agent count and survival, ending body condition broken out by reproductive
#'   status, and monthly tables of net energy, dry matter intake, and distance moved.
#'
#' @param x An object of class \code{summary_abm} from \code{summary_abm}.
#' @param ... Ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @export
print.summary_abm <- function(x, ...) {

  #1) header and survival

  cat("NCC agent-based model summary\n")
  cat(sprintf(
    "Agents: %d    Survived: %d (%.0f%%)\n\n",
    x$n_agents,
    x$n_survived,
    100 * x$prop_survived
  ))

  #2) ending body condition

  cat("Ending body condition (ifbfat)\n")
  cat(sprintf(
    "  All: %.3f    Reproductive: %.3f    Non-reproductive: %.3f\n\n",
    x$ifbfat[["all"]],
    x$ifbfat[["rep"]],
    x$ifbfat[["non_rep"]]
  ))

  #3) monthly tables

  print_month_table <- function(mat, title) {
    cat(title, "\n")
    cat(sprintf("  %-5s %9s %9s %9s\n", "", "All", "Rep", "Non-rep"))
    fmt <- function(v) if (is.nan(v) || is.na(v)) "-" else sprintf("%.1f", v)
    for (m in rownames(mat)) {
      cat(sprintf(
        "  %-5s %9s %9s %9s\n",
        m,
        fmt(mat[m, "all"]),
        fmt(mat[m, "rep"]),
        fmt(mat[m, "non_rep"])
      ))
    }
    cat("\n")
  }

  print_month_table(x$net_energy, "Net energy (kJ/day)")
  print_month_table(x$dmi, "Dry matter intake (g/day)")
  print_month_table(x$distance, "Distance moved (m/day)")

  invisible(x)
}

#' Plot season-long trajectories and tracks from an NCC run
#'
#' @description Visual companion to \code{summary_abm}: individual daily
#'   trajectories of body condition and net energy, a survival curve, and
#'   hourly movement tracks over the forage landscape. Trajectory lines are
#'   colored by reproductive status; dead individuals are retained and their
#'   lines simply end at death.
#'
#' @param result A list returned by \code{ncc_abm}.
#' @param forage_reference The forage SpatRaster passed to \code{ncc_abm}; its
#'   first layer is used as the backdrop for the movement-tracks panel.
#'
#' @return A named list of ggplot objects: \code{condition}, \code{net_energy},
#'   \code{expenditure}, \code{survival}, and \code{tracks}.
#'
#' @export
plot_abm <- function(result, forage_reference) {

  # ----------------------------------------------------------------------------------------------------------------------
  # setup
  # ----------------------------------------------------------------------------------------------------------------------

  #1) unpack

  agents <- result$agents
  agent_params <- result$agent_params

  #2) daily row index and date axis

  n_times <- nrow(agents[[1]])
  daily_rows <- seq(24, n_times, by = 24)
  dates_axis <- agents[[1]]$datetime[daily_rows]

  #3) reproductive status aligned to the agents list order

  rep_status <- agent_params$rep_status[match(names(agents), agent_params$id)]

  #4) shared color scale for reproductive status

  status_scale <- scale_color_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("0" = "non-reproductive", "1" = "reproductive"),
    name = "reproductive status"
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # assemble plotting data
  # ----------------------------------------------------------------------------------------------------------------------

  #1) daily trajectories (one row per agent per day)

  traj <- do.call(rbind, lapply(seq_along(agents), function(i) {
    a <- agents[[i]]
    data.frame(
      id = names(agents)[i],
      date = a$datetime[daily_rows],
      rep_status = factor(rep_status[i], levels = c(0, 1)),
      ifbfat = a$ifbfat[daily_rows],
      energy_net = a$energy_net[daily_rows],
      total_expenditure = a$energy_bmr[daily_rows] +
        a$energy_hif[daily_rows] +
        a$energy_loc[daily_rows] +
        a$energy_rep[daily_rows]
    )
  }))

  #2) survival curve (agents alive at each daily row)

  n_alive <- vapply(
    daily_rows,
    function(t) sum(vapply(agents, function(a) isTRUE(a$status[t] == "ALIVE"), logical(1))),
    numeric(1)
  )
  surv <- data.frame(date = dates_axis, n_alive = n_alive)

  #3) hourly movement tracks (dead agents truncate at death)

  trk <- do.call(rbind, lapply(seq_along(agents), function(i) {
    a <- agents[[i]]
    data.frame(id = names(agents)[i], x = a$x, y = a$y)
  }))

  #4) forage backdrop (first layer)

  backdrop <- as.data.frame(forage_reference[[1]], xy = TRUE, na.rm = TRUE)
  names(backdrop)[3] <- "forage"

  # ----------------------------------------------------------------------------------------------------------------------
  # build panels
  # ----------------------------------------------------------------------------------------------------------------------

  #1) body condition with the nutritional carrying-capacity threshold

  cc <- get_param("carrying_capacity")

  p_condition <- ggplot(traj, aes(date, ifbfat, group = id, color = rep_status)) +
    geom_hline(yintercept = cc, linetype = "dotdash", linewidth = 1.1, color = "grey25") +
    geom_line(linewidth = 1.1, alpha = 0.5) +
    status_scale +
    labs(title = "Body condition", x = "date", y = "ifbfat") +
    theme_martin(base_size = 14)

  #2) daily net energy with a zero reference

  p_net <- ggplot(traj, aes(date, energy_net, group = id, color = rep_status)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.1, color = "grey50") +
    geom_line(linewidth = 1.1, alpha = 0.5) +
    status_scale +
    labs(title = "Net energy", x = "date", y = "net energy (kJ/day)") +
    theme_martin(base_size = 14)

  #3) survival curve

  p_survival <- ggplot(surv, aes(date, n_alive)) +
    geom_step(linewidth = 1.1, color = "grey25") +
    labs(title = "Survival", x = "date", y = "individuals alive") +
    theme_martin(base_size = 14)

  #4) movement tracks over the forage backdrop

  p_tracks <- ggplot() +
    geom_raster(data = backdrop, aes(x, y, fill = forage)) +
    geom_path(data = trk, aes(x, y, group = id), color = "white", linewidth = 0.7) +
    scale_fill_viridis_c(name = "forage") +
    coord_fixed() +
    labs(title = "Movement", x = NULL, y = NULL) +
    theme_martin(base_size = 14)

  #5) total metabolic expenditure (basal metabolism + heat increment + locomotion + lactation)

  p_expenditure <- ggplot(traj, aes(date, total_expenditure, group = id, color = rep_status)) +
    geom_line(linewidth = 1.1, alpha = 0.5) +
    status_scale +
    labs(title = "Metabolic expenditure", x = "date", y = "expenditure (kJ/day)") +
    theme_martin(base_size = 14)

  #6) return the panels

  list(
    condition = p_condition,
    net_energy = p_net,
    expenditure = p_expenditure,
    survival = p_survival,
    tracks = p_tracks
  )
}

