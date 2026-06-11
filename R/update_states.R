#' Cap and decay the grazing deficit at the end of a day
#'
#' @description Advances the grazing deficit from the end of one day to the start
#'   of the next. The deficit is the amount (g/cell) by which grazing has pushed
#'   realized biomass below potential; realized biomass on any day is
#'   \code{potential - deficit}. The deficit is first capped at the next day's
#'   potential so that, after decay, realized biomass on the next day cannot fall
#'   below zero (the deficit can never represent more missing biomass than the
#'   next day's potential supports). It is then decayed by a fixed fraction
#'   (\code{plant_regrowth_rate}, a 21-day half-life; Oesterheld 1992), shrinking
#'   geometrically toward zero. An ungrazed cell has zero deficit and therefore
#'   sits exactly at potential every day, regardless of the seasonal phenology of
#'   potential.
#'
#' @param deficit Numeric vector of the day's grazing deficit (g/cell), indexed by
#'   cell number, after the day's consumption has been added.
#' @param potential_next Numeric vector of the next day's potential biomass
#'   (g/cell), indexed by cell number (the next layer of forage_reference).
#'
#' @return Numeric vector of the next day's start-of-day deficit (g/cell).
#'
#' @keywords internal
update_forage <- function(deficit, potential_next) {

  # ----------------------------------------------------------------------------------------------------------------------
  # cap then decay the deficit
  # ----------------------------------------------------------------------------------------------------------------------

  #1) cap the deficit at the next day's potential so realized biomass stays >= 0

  deficit <- pmin(deficit, potential_next)

  #2) recover a fixed fraction of the capped deficit

  rate <- get_param("plant_regrowth_rate")

  (1 - rate) * deficit
}

#' Update body mass and body fat from the day's fat-mass change
#'
#' @description Propagates an agent's body mass and ingesta-free body fat fraction
#'   from one day to the next given the day's change in fat mass. Lean mass is held
#'   constant: only the fat compartment changes with energy balance. The new fat
#'   mass is the current fat mass plus \code{fat_change}, the new body mass is lean
#'   mass plus the new fat mass, and the new body fat fraction is the new fat mass
#'   over the new body mass. Fat mass is floored at zero so the fat compartment
#'   cannot go negative.
#'
#' @param bm Numeric. Current body mass (kg).
#' @param ifbfat Numeric. Current ingesta-free body fat as a fraction of body mass.
#' @param fat_change Numeric. The day's change in fat mass (kg), from
#'   \code{calc_fat_change}.
#'
#' @return A named numeric vector with elements \code{bm} (new body mass, kg) and
#'   \code{ifbfat} (new body fat fraction).
#'
#' @keywords internal
update_mass <- function(bm, ifbfat, fat_change) {

  # ----------------------------------------------------------------------------------------------------------------------
  # propagate the fat compartment, holding lean mass constant
  # ----------------------------------------------------------------------------------------------------------------------

  #1) current lean and fat compartments

  lean_mass <- calc_lean_mass(bm, ifbfat)
  fat_mass <- calc_fat_mass(bm, ifbfat)

  #2) apply the day's fat-mass change, floored at zero

  fat_mass <- max(fat_mass + fat_change, 0)

  #3) recombine into body mass and body fat fraction

  bm_new <- lean_mass + fat_mass

  c(bm = bm_new, ifbfat = fat_mass / bm_new)
}

#' Update survival status from fat reserves and advance days post partum
#'
#' @description Determines whether an agent survives the day and advances its
#'   days-post-partum counter by one day. An agent dies when its fat reserves are
#'   exhausted: if fat mass is at or below zero it is marked \code{"DEAD"},
#'   otherwise \code{"ALIVE"}. The days-post-partum counter is incremented by one;
#'   for non-reproductive agents it is \code{NA} and stays \code{NA}.
#'
#' @param fat_mass Numeric. The agent's current fat mass (kg).
#' @param j_post_partum Numeric. The agent's current days post partum (\code{NA}
#'   for non-reproductive agents).
#'
#' @return A named list with elements \code{status} (\code{"ALIVE"} or
#'   \code{"DEAD"}) and \code{j_post_partum} (the counter advanced by one day).
#'
#' @keywords internal
update_status <- function(fat_mass, j_post_partum) {
  status <- if (fat_mass <= 0) "DEAD" else "ALIVE"
  list(status = status, j_post_partum = j_post_partum + 1)
}
