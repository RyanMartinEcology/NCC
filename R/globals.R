#' Package parameter environment
#'
#' @description Internal environment that holds every global model parameter. Storing
#'   parameters in a dedicated environment (rather than as package objects) lets them be
#'   read, overridden, and redrawn at run time through \code{get_param()},
#'   \code{set_param()}, and \code{draw_param()} without reassigning into the package
#'   namespace. Populated by \code{.set_defaults()} at load.
#'
#' @keywords internal
.ncc_env <- new.env(parent = emptyenv())

#' Populate the parameter environment with model defaults
#'
#' @description Writes every default model parameter into \code{.ncc_env}: simulation and
#'   time settings, the study-area latitude, the precomputed daylight schedule, the
#'   energetic and body-condition constants, the lactation schedule, the plant regrowth
#'   rate, the intake response, and the population-level movement distribution. Run once
#'   at load by \code{.onLoad()}. Most scalars carry \code{unit}, \code{source}, and
#'   \code{full_name} attributes that document the value.
#'
#' @return Called for its side effect of populating \code{.ncc_env}; the return value is
#'   not used.
#'
#' @keywords internal
.set_defaults <- function() {

  # ----------------------------------------------------------------------------------------------------------------------
  # simulation parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the ingesta-free body fat fraction drawn as the reference line on the body-condition panel of plot_abm()

  .ncc_env$carrying_capacity <- 0.1466
  attr(.ncc_env$carrying_capacity, 'unit') <- 'Percent'
  attr(.ncc_env$carrying_capacity, 'source') <- NA
  attr(.ncc_env$carrying_capacity, 'full_name') <- 'Global mean ifbfat threshold for carrying capacity'

  #2) the number of agents created by create_agents(); set per run with set_param()

  .ncc_env$n_agents <- 125 # denotes current number of agents, can be updated within the model

  #3) the number of movement-only steps each agent takes before the season starts, so that starting positions
  #   reflect the agents' own movement models rather than the cells they were dropped in

  .ncc_env$burn_in <- 1000L
  attr(.ncc_env$burn_in, 'unit') <- 'steps'
  attr(.ncc_env$burn_in, 'source') <- NA
  attr(.ncc_env$burn_in, 'full_name') <- 'Number of burn-in movement steps drawn before the simulation starts'

  # ----------------------------------------------------------------------------------------------------------------------
  # time parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the hourly time step (t_delta) and the season window in the study time zone. the season is stored as
  #   month-day strings, and t_start and t_end are built from them here for the default year and rebuilt by
  #   .set_season() for the year ncc_abm() is given. every hourly and daily sequence in the model runs from
  #   t_start to t_end

  .ncc_env$t_delta <- lubridate::hours(1)
  .ncc_env$season_start_md <- "07-01"
  .ncc_env$season_end_md <- "10-15"
  default_year <- 2024

  .ncc_env$t_start <- as.POSIXct(
    paste0(default_year, "-", .ncc_env$season_start_md),
    tz = "America/Denver"
  )
  .ncc_env$t_end <- as.POSIXct(
    paste0(default_year, "-", .ncc_env$season_end_md),
    tz = "America/Denver"
  )

  # ----------------------------------------------------------------------------------------------------------------------
  # spatial parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the study-area latitude, used only to compute the daylight schedule. the rasters carry their own projection,
  #   so no coordinate reference system is stored here

  .ncc_env$study_lat <- 43.74075
  attr(.ncc_env$study_lat, 'unit') <- 'decimal degrees'
  attr(.ncc_env$study_lat, 'source') <- 'Grand Teton summit — placeholder'
  attr(.ncc_env$study_lat, 'full_name') <- 'Study area latitude for solar position calculation'

  # ----------------------------------------------------------------------------------------------------------------------
  # daylight schedule
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the hourly daylight flags and day lengths for the season, built by .refresh_daylight() from the time and
  #   latitude parameters above

  .refresh_daylight()

  # ----------------------------------------------------------------------------------------------------------------------
  # energy parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) forage energy: the digestible energy content of suitable forage (DE) and the fraction of digestible energy
  #   that is metabolizable. their product converts grams of forage eaten to kJ of metabolizable energy

  .ncc_env$DE <- 12.98544 # this is the weighted mean of digestible energy of suitable forage biomass in vegetation transects.
  attr(.ncc_env$DE, 'unit') <- 'kJ/g'
  attr(.ncc_env$DE, 'source') <- NA
  attr(.ncc_env$DE, 'full_name') <- 'Digestible Energy'

  .ncc_env$DE_to_ME_conversion_factor <- 0.82
  attr(.ncc_env$DE_to_ME_conversion_factor, 'unit') <- 'Percent'
  attr(.ncc_env$DE_to_ME_conversion_factor, 'source') <- 'NRC 2007'
  attr(.ncc_env$DE_to_ME_conversion_factor, 'full_name') <- NA

  #2) the locomotion cost per kg of body mass per m travelled, by slope class: descent steeper than 10 degrees
  #   (d_10), descent of 1 to 10 degrees (d_1_10), flat (f), incline of 1 to 10 degrees (i_1_10), and incline
  #   steeper than 10 degrees (i_10). calc_energy_loc() selects one per step by the signed slope of the step

  .ncc_env$distance_cost_factor_d_10 <- 5.34
  attr(.ncc_env$distance_cost_factor_d_10, 'unit') <- 'J * kg^-1 * m^-1'
  attr(.ncc_env$distance_cost_factor_d_10, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$distance_cost_factor_d_10, 'full_name') <- NA

  .ncc_env$distance_cost_factor_d_1_10 <- 2.00
  attr(.ncc_env$distance_cost_factor_d_1_10, 'unit') <- 'J * kg^-1 * m^-1'
  attr(.ncc_env$distance_cost_factor_d_1_10, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$distance_cost_factor_d_1_10, 'full_name') <- NA

  .ncc_env$distance_cost_factor_f <- 4.95
  attr(.ncc_env$distance_cost_factor_f, 'unit') <- 'J * kg^-1 * m^-1'
  attr(.ncc_env$distance_cost_factor_f, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$distance_cost_factor_f, 'full_name') <- NA

  .ncc_env$distance_cost_factor_i_1_10 <- 7.44
  attr(.ncc_env$distance_cost_factor_i_1_10, 'unit') <- 'J * kg^-1 * m^-1'
  attr(.ncc_env$distance_cost_factor_i_1_10, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$distance_cost_factor_i_1_10, 'full_name') <- NA

  .ncc_env$distance_cost_factor_i_10 <- 21.08
  attr(.ncc_env$distance_cost_factor_i_10, 'unit') <- 'J * kg^-1 * m^-1'
  attr(.ncc_env$distance_cost_factor_i_10, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$distance_cost_factor_i_10, 'full_name') <- NA

  #3) the heat increment of feeding (HIF), an hourly cost per kg of body mass, and the proportion of the diurnal
  #   period spent feeding. calc_energy_hif() multiplies them by body mass and day length for the daily cost

  .ncc_env$HIF <- 1.799 # 0.43 kcal * kg^-1 * h^-1 (females) x 4.184 kJ/kcal; Chappel and Hudson 1978b
  attr(.ncc_env$HIF, 'unit') <- 'kJ * kg^-1 * h^-1'
  attr(.ncc_env$HIF, 'source') <- 'Chappel and Hudson 1978b'
  attr(.ncc_env$HIF, 'full_name') <- 'Heat Increment of Feeding'

  .ncc_env$prop_day_forage <- 0.72
  attr(.ncc_env$prop_day_forage, 'unit') <- 'proportion'
  attr(.ncc_env$prop_day_forage, 'source') <- 'Courtemanch et al. 2014'
  attr(.ncc_env$prop_day_forage, 'full_name') <- 'Proportion of the diurnal period spent feeding'

  #4) the energy density of fat, which converts a daily net energy balance into a change in fat mass. the
  #   conversion is lossless in both directions: a surplus deposits and a deficit mobilizes at the same rate

  .ncc_env$E_fat <- 39.5
  attr(.ncc_env$E_fat, 'unit') <- 'kJ/g'
  attr(.ncc_env$E_fat, 'source') <- 'Robbins 1993'
  attr(.ncc_env$E_fat, 'full_name') <- 'Energy Value of Fat Reserves'

  # ----------------------------------------------------------------------------------------------------------------------
  # body condition and body mass parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the starting body mass shared by all agents (bm) and the draw function for each agent's starting
  #   ingesta-free body fat fraction (ifbf): a normal draw in percent, converted to a fraction and floored at 0.01

  .ncc_env$bm <- 56.91 # calculated from Teton capture data

  .ncc_env$ifbf <- function() max(rnorm(1, mean = 8.39, sd = 2.91) / 100, 0.01) # taken from Smiley et al. 2022, floored at 0.01 (1% body fat) to bar impossible negative/near-zero draws

  # ----------------------------------------------------------------------------------------------------------------------
  # pregnancy and lactation parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the draw function for each agent's reproductive status (1 = lactating, 0 = not), the days post partum at
  #   which lactating agents start the season, and the lactation modifier: the multiple of basal metabolism that
  #   lactation costs on each day post partum, which calc_energy_rep() indexes by the agent's days post partum.
  #   the vector rises from 0.65 at day 1 to 1 at day 21, holds at 1 through day 42, then declines to 0.379 at
  #   day 168 and is 0 at day 169

  .ncc_env$rep_status <- function() rbinom(1, 1, 0.678571429) #this is the proportion of captured females that showed some evidence of lactation
  .ncc_env$j_post_partum <- 25
  .ncc_env$lactation_modifier <- c(0.65, 0.664152312, 0.678612759, 0.693388051, 0.708485042,
                                   0.723910736, 0.739672291, 0.755777018, 0.772232391, 0.789046043,
                                   0.806225775, 0.823779557, 0.841715535, 0.860042028, 0.878767541,
                                   0.89790076, 0.917450563, 0.937426019, 0.957836397, 0.978691165,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   0.992320762, 0.984700494, 0.977138744, 0.969635063, 0.962189004,
                                   0.954800126, 0.947467988, 0.940192156, 0.932972196, 0.92580768,
                                   0.918698182, 0.91164328, 0.904642554, 0.897695588, 0.89080197,
                                   0.883961289, 0.87717314, 0.870437118, 0.863752824, 0.857119861,
                                   0.850537833, 0.84400635, 0.837525024, 0.83109347, 0.824711305,
                                   0.818378151, 0.81209363, 0.805857369, 0.799668998, 0.79352815,
                                   0.787434458, 0.781387561, 0.7753871, 0.769432717, 0.76352406,
                                   0.757660777, 0.751842519, 0.746068941, 0.7403397, 0.734654455,
                                   0.729012869, 0.723414605, 0.717859332, 0.712346719, 0.706876439,
                                   0.701448166, 0.696061579, 0.690716356, 0.68541218, 0.680148737,
                                   0.674925713, 0.669742797, 0.664599683, 0.659496063, 0.654431636,
                                   0.6494061, 0.644419155, 0.639470507, 0.634559861, 0.629686924,
                                   0.624851408, 0.620053025, 0.61529149, 0.610566521, 0.605877835,
                                   0.601225154, 0.596608203, 0.592026707, 0.587480392, 0.582968991,
                                   0.578492233, 0.574049853, 0.569641587, 0.565267174, 0.560926353,
                                   0.556618865, 0.552344457, 0.548102872, 0.543893859, 0.539717169,
                                   0.535572552, 0.531459763, 0.527378557, 0.523328691, 0.519309925,
                                   0.515322021, 0.51136474, 0.507437848, 0.503541112, 0.4996743,
                                   0.495837182, 0.49202953, 0.488251118, 0.484501721, 0.480781117,
                                   0.477089084, 0.473425404, 0.469789857, 0.466182229, 0.462602304,
                                   0.459049871, 0.455524718, 0.452026635, 0.448555415, 0.445110851,
                                   0.441692738, 0.438300875, 0.434935058, 0.431595088, 0.428280766,
                                   0.424991896, 0.421728282, 0.41848973, 0.415276048, 0.412087044,
                                   0.408922529, 0.405782316, 0.402666217, 0.399574047, 0.396505623,
                                   0.393460761, 0.390439282, 0.387441006, 0.384465754, 0.38151335,
                                   0.378583618, 0)

  # ----------------------------------------------------------------------------------------------------------------------
  # plant trait parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the fraction of the grazing deficit recovered each day, set from a 42-day half-life. update_forage() applies
  #   it to the deficit at the end of each day

  .ncc_env$plant_regrowth_rate <- 1 - 0.5^(1/42) # Oesterheld 1992

  # ----------------------------------------------------------------------------------------------------------------------
  # intake parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the dry-matter-intake functional response, a negative exponential in standing vegetation biomass and
  #   mass-specific in metabolic body mass:
  #     daily intake (g) = intake_max * bm^0.75 * (1 - exp(-V / intake_decay))
  #   with V the biomass in kg/ha. intake_max is the asymptote in g per kg^0.75 per day and intake_decay the
  #   biomass scale over which intake approaches it (half-maximal intake occurs at intake_decay * log(2) kg/ha).
  #   the same response applies to both reproductive states, and ncc_abm() spreads the daily total evenly across
  #   the day's foraging hours

  .ncc_env$intake_max <- 83.4
  attr(.ncc_env$intake_max, 'unit') <- 'g * kg^-0.75 * day^-1'
  attr(.ncc_env$intake_max, 'source') <- NA
  attr(.ncc_env$intake_max, 'full_name') <- 'Mass-specific daily dry matter intake asymptote'

  .ncc_env$intake_decay <- 111
  attr(.ncc_env$intake_decay, 'unit') <- 'kg/ha'
  attr(.ncc_env$intake_decay, 'source') <- NA
  attr(.ncc_env$intake_decay, 'full_name') <- 'Biomass scale constant of the intake functional response'

  #2) the multiplier applied to the intake response, a named length-2 vector indexed by rep_status + 1 (nonrepro =
  #   non-lactating, repro = lactating). non-lactating agents take the response as is; lactating agents take 1.5
  #   times it, raising their asymptote at bm = 56.91 kg from 1728 to 2592 g/day

  .ncc_env$intake_multiplier <- c(
    nonrepro = 1,
    repro = 1.5
  )
  attr(.ncc_env$intake_multiplier, 'unit') <- 'multiplier'
  attr(.ncc_env$intake_multiplier, 'source') <- NA
  attr(.ncc_env$intake_multiplier, 'full_name') <- 'Intake functional response multiplier, by reproductive status'

  # ----------------------------------------------------------------------------------------------------------------------
  # movement parameters
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the number of candidate steps drawn for each agent at each movement step

  .ncc_env$n_candidates <- 25L
  attr(.ncc_env$n_candidates, 'unit') <- 'candidate steps'
  attr(.ncc_env$n_candidates, 'source') <- NA
  attr(.ncc_env$n_candidates, 'full_name') <- 'Number of candidate steps drawn per redistribution kernel'

  #2) the population-level multivariate normal distribution of movement parameters, from the 22 per-individual
  #   iSSF fits. each agent's movement model is one draw from MVN(mvn_mu, mvn_sigma) in draw_movement_params().
  #   the 15 dimensions, in order:
  #     log_shape, log_scale: the gamma step-length distribution parameters, on the log scale
  #     log_kappa: the von Mises turn-angle concentration, on the log scale
  #     sl_, log(sl_), cos(ta_): the day coefficients that adjust the step-length and turn-angle distributions
  #     forage_biomass_end, escape_terrain_end, canopy_cover_end: the day habitat selection coefficients, on
  #       forage in g/cell, distance to escape terrain in m, and canopy cover
  #     the six :tod_end_night_end terms: the night offsets, added to the matching day coefficient at night
  #   the geometry terms are written in call form (log(sl_), cos(ta_)) and the habitat terms with the _end suffix,
  #   matching the fitted models; the movement functions extract the twelve coefficients by these names, and
  #   validate_move_coefs() stops a run whose model does not carry exactly this set

  mvn_names <- c("log_shape", "log_scale", "log_kappa",
                 "sl_", "log(sl_)", "cos(ta_)",
                 "forage_biomass_end", "escape_terrain_end", "canopy_cover_end",
                 "sl_:tod_end_night_end", "log(sl_):tod_end_night_end", "cos(ta_):tod_end_night_end",
                 "forage_biomass_end:tod_end_night_end", "escape_terrain_end:tod_end_night_end",
                 "canopy_cover_end:tod_end_night_end")

  #3) the mean vector

  .ncc_env$mvn_mu <- c(
    "log_shape" = -6.905901e-01,
    "log_scale" = 5.752240e+00,
    "log_kappa" = -2.401935e+00,
    "sl_" = 2.708199e-04,
    "log(sl_)" = 2.547561e-01,
    "cos(ta_)" = 2.411655e-01,
    "forage_biomass_end" = -1.314966e-05,
    "escape_terrain_end" = -2.759348e-03,
    "canopy_cover_end" = -1.499584e-02,
    "sl_:tod_end_night_end" = -6.798309e-04,
    "log(sl_):tod_end_night_end" = -3.099616e-01,
    "cos(ta_):tod_end_night_end" = -4.453259e-01,
    "forage_biomass_end:tod_end_night_end" = 2.202214e-06,
    "escape_terrain_end:tod_end_night_end" = -9.140983e-03,
    "canopy_cover_end:tod_end_night_end" = -1.817424e-02
  )

  #4) the among-individual covariance matrix, filled by row in the order of mvn_names; the matrix is symmetric

  .ncc_env$mvn_sigma <- matrix(
    c(
      2.882366e-03, -4.464718e-03, 8.722793e-03, -8.740639e-06, 7.661434e-04, -8.658777e-04, 2.402033e-07, 4.015119e-05, 2.828031e-04, 3.965397e-06, -1.174691e-03, 9.540974e-04, -2.387599e-09, 3.378749e-05, -2.795464e-04,
      -4.464718e-03, 2.188443e-02, 3.953503e-02, -1.503428e-05, 2.730466e-03, 4.803810e-03, -6.424334e-07, 2.570127e-05, -6.163560e-04, 3.125365e-05, -3.545337e-03, -2.995137e-03, 1.022565e-07, 6.021152e-05, 1.024416e-03,
      8.722793e-03, 3.953503e-02, 3.878611e-01, -1.520050e-04, 1.840269e-02, 1.506794e-02, -2.863984e-07, 4.646876e-04, 4.274598e-04, 1.799358e-04, -2.586744e-02, -1.799202e-02, -5.238914e-08, 6.929005e-04, 6.535663e-04,
      -8.740639e-06, -1.503428e-05, -1.520050e-04, 1.471616e-07, -1.038839e-05, -1.108390e-05, -1.489394e-09, -5.894355e-07, -1.005024e-06, -1.147480e-07, 1.522050e-05, 1.189478e-06, 1.139270e-09, -7.536502e-07, 1.162117e-07,
      7.661434e-04, 2.730466e-03, 1.840269e-02, -1.038839e-05, 3.600731e-03, 2.473636e-03, 2.442349e-08, -2.955850e-05, -2.777794e-04, 1.634906e-05, -4.301191e-03, -3.552449e-03, -5.159011e-08, 3.595848e-05, 3.651190e-05,
      -8.658777e-04, 4.803810e-03, 1.506794e-02, -1.108390e-05, 2.473636e-03, 5.733813e-03, -7.615615e-08, -2.505558e-06, -1.948975e-04, 1.062269e-05, -2.993055e-03, -6.036426e-03, -1.567597e-08, 9.747783e-05, -8.080815e-06,
      2.402033e-07, -6.424334e-07, -2.863984e-07, -1.489394e-09, 2.442349e-08, -7.615615e-08, 9.984019e-11, 6.135135e-09, 3.359245e-08, -1.013181e-10, -2.852680e-08, 3.068313e-07, -7.074911e-11, 1.130349e-08, -5.066978e-08,
      4.015119e-05, 2.570127e-05, 4.646876e-04, -5.894355e-07, -2.955850e-05, -2.505558e-06, 6.135135e-09, 5.284780e-06, 1.544404e-05, -5.966437e-08, 1.969635e-05, 9.263255e-05, -6.425070e-10, 4.218420e-06, -4.320912e-06,
      2.828031e-04, -6.163560e-04, 4.274598e-04, -1.005024e-06, -2.777794e-04, -1.948975e-04, 3.359245e-08, 1.544404e-05, 1.052512e-04, -1.370212e-06, 2.690999e-04, 3.768379e-04, 1.749792e-08, 7.852534e-06, -2.410731e-05,
      3.965397e-06, 3.125365e-05, 1.799358e-04, -1.147480e-07, 1.634906e-05, 1.062269e-05, -1.013181e-10, -5.966437e-08, -1.370212e-06, 2.818054e-07, -2.338597e-05, -1.142270e-05, -9.489072e-10, 1.616836e-07, 2.302001e-06,
      -1.174691e-03, -3.545337e-03, -2.586744e-02, 1.522050e-05, -4.301191e-03, -2.993055e-03, -2.852680e-08, 1.969635e-05, 2.690999e-04, -2.338597e-05, 5.298479e-03, 4.124893e-03, 3.964541e-08, -5.424057e-05, -6.580214e-05,
      9.540974e-04, -2.995137e-03, -1.799202e-02, 1.189478e-06, -3.552449e-03, -6.036426e-03, 3.068313e-07, 9.263255e-05, 3.768379e-04, -1.142270e-05, 4.124893e-03, 1.026174e-02, -1.205752e-07, -2.201733e-05, -4.733059e-05,
      -2.387599e-09, 1.022565e-07, -5.238914e-08, 1.139270e-09, -5.159011e-08, -1.567597e-08, -7.074911e-11, -6.425070e-10, 1.749792e-08, -9.489072e-10, 3.964541e-08, -1.205752e-07, 1.094024e-10, -8.291312e-09, 2.031977e-08,
      3.378749e-05, 6.021152e-05, 6.929005e-04, -7.536502e-07, 3.595848e-05, 9.747783e-05, 1.130349e-08, 4.218420e-06, 7.852534e-06, 1.616836e-07, -5.424057e-05, -2.201733e-05, -8.291312e-09, 6.472494e-06, -7.952111e-06,
      -2.795464e-04, 1.024416e-03, 6.535663e-04, 1.162117e-07, 3.651190e-05, -8.080815e-06, -5.066978e-08, -4.320912e-06, -2.410731e-05, 2.302001e-06, -6.580214e-05, -4.733059e-05, 2.031977e-08, -7.952111e-06, 1.149536e-04
    ),
    nrow = 15,
    ncol = 15,
    byrow = TRUE,
    dimnames = list(mvn_names, mvn_names)
  )
}

#' Recompute the daylight schedule from the current clock and location
#'
#' @description Rebuilds \code{is_daylight} and \code{day_length} in \code{.ncc_env} from the
#'   current \code{t_start}, \code{t_end}, \code{t_delta}, and \code{study_lat}. Both vectors have
#'   one element per hourly time step of \code{seq(t_start, t_end, by = t_delta)} and are indexed by
#'   time-step position, so they are valid only for the season they were built from. Called at load
#'   by \code{.set_defaults()} and again by \code{.set_season()} whenever the season window moves.
#'
#'   Sunrise and sunset are geometric, from study latitude and day of year, with clock noon taken
#'   as solar noon.
#'
#' @return Called for its side effect; the return value is not used.
#'
#' @keywords internal
.refresh_daylight <- function() {

  # ----------------------------------------------------------------------------------------------------------------------
  # daylight schedule, one entry per hourly time step
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the hourly time sequence of the model clock, from t_start to t_end by t_delta

  daylight_times <- seq(
    .ncc_env$t_start,
    .ncc_env$t_end,
    by = as.numeric(.ncc_env$t_delta, units = "secs")
  )

  #2) solar geometry from latitude and day of year: the solar declination and the resulting half-day length
  #   (hours from solar noon to sunset), plus the local clock hour of each time step

  daylight_lat_rad <- .ncc_env$study_lat * pi / 180
  daylight_doy <- lubridate::yday(daylight_times)
  daylight_decl <- 0.409 * sin(2 * pi / 365 * daylight_doy - 1.39)
  daylight_half_day <- (12 / pi) * acos(-tan(daylight_lat_rad) * tan(daylight_decl))
  daylight_hour <- lubridate::hour(daylight_times) + lubridate::minute(daylight_times) / 60

  #3) the daylight flag per time step: TRUE when the clock hour falls within the half-day on either side of noon

  .ncc_env$is_daylight <- daylight_hour >= (12 - daylight_half_day) &
    daylight_hour <= (12 + daylight_half_day)

  #4) the day length (hours) per time step, twice the half-day. it depends only on day of year, so it is constant
  #   within a day

  .ncc_env$day_length <- 2 * daylight_half_day

  invisible(NULL)
}

#' Point the simulation clock at a given year
#'
#' @description Rebuilds \code{t_start} and \code{t_end} for \code{year} from the stored
#'   \code{season_start_md} and \code{season_end_md}, then refreshes the daylight schedule so the
#'   clock and the daylight vectors always describe the same season. \code{ncc_abm()} calls it
#'   before reading the clock, so the model can be pointed at a different year's forage rasters
#'   through its \code{year} argument.
#'
#'   Because \code{t_start} and \code{t_end} are rebuilt here, assigning them directly through
#'   \code{set_param()} does not survive the next \code{ncc_abm()} call. Change the season through
#'   \code{season_start_md} / \code{season_end_md} and the \code{year} argument instead.
#'
#' @param year Numeric or integer. Calendar year the season runs in.
#'
#' @return Called for its side effect; the return value is not used.
#'
#' @keywords internal
.set_season <- function(year) {

  # ----------------------------------------------------------------------------------------------------------------------
  # rebuild the season window and the daylight schedule
  # ----------------------------------------------------------------------------------------------------------------------

  #1) the season window for this year, from the stored month-day bounds

  .ncc_env$t_start <- as.POSIXct(
    paste0(year, "-", .ncc_env$season_start_md),
    tz = "America/Denver"
  )

  .ncc_env$t_end <- as.POSIXct(
    paste0(year, "-", .ncc_env$season_end_md),
    tz = "America/Denver"
  )

  #2) rebuild the daylight vectors for the new window

  .refresh_daylight()

  invisible(NULL)
}

#' Package load hook
#'
#' @description Run automatically by R when the package is loaded. Populates the
#'   parameter environment with model defaults by calling \code{.set_defaults()}.
#'
#' @param libname Library path; supplied by R.
#' @param pkgname Package name; supplied by R.
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .set_defaults()
}

#' Get a global model parameter
#'
#' @description Returns a parameter from the package parameter environment exactly as
#'   stored, including any \code{unit}, \code{source}, and \code{full_name} attributes.
#'   Parameters stored as draw functions (for example \code{ifbf} and \code{rep_status})
#'   are returned as the function itself; use \code{draw_param()} to obtain a drawn value.
#'
#' @param param Character string. Name of the parameter.
#'
#' @return The stored parameter value, of whatever type it was set to, or \code{NULL} if
#'   the parameter is not set.
#'
#' @importFrom lubridate hours yday
#' @export
get_param <- function(param) {
  .ncc_env[[param]]
}

#' Set or override a global model parameter
#'
#' @description Writes a value into the package parameter environment, creating the
#'   parameter if it does not exist and replacing it if it does. Used to override a
#'   default at run time, for example raising \code{n_agents} or \code{n_candidates} for
#'   a particular run.
#'
#' @param param Character string. Name of the parameter.
#' @param value New value to store; any type.
#'
#' @return The assigned \code{value}, invisibly.
#'
#' @export
set_param <- function(param, value) {
  .ncc_env[[param]] <- value
}

#' Draw a parameter value
#'
#' @description Returns a parameter value, first evaluating it if it is stored as a
#'   function. Several individual-level parameters (\code{ifbf}, \code{rep_status}) are
#'   stored as random-draw functions so that each agent receives its own value; this
#'   evaluates the function once per call. Parameters that are not functions are returned
#'   unchanged, so the same call works for fixed and stochastic parameters alike.
#'
#' @param param Character string. Name of the parameter.
#'
#' @return For a parameter stored as a function, one draw from it; otherwise the stored
#'   value unchanged.
#'
#' @export
draw_param <- function(param) {
  val <- .ncc_env[[param]]
  if (is.function(val)) val() else val
}
