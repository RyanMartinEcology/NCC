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
#'   time settings, the study-area location, the precomputed daylight schedule, the
#'   energetic and body-condition constants, the pregnancy and lactation schedule, plant
#'   and calibration parameters, and the population-level movement distribution. Run once
#'   at load by \code{.onLoad()}; call again by hand to refresh the derived daylight
#'   vectors after changing any time or location parameter. Most scalars carry
#'   \code{unit}, \code{source}, and \code{full_name} attributes that document the value.
#'
#' @return Called for its side effect of populating \code{.ncc_env}; the return value is
#'   not used.
#'
#' @keywords internal
.set_defaults <- function() {

  # -----------------------------------------------------------------------------------------------------
  # simulation parameters
  # -----------------------------------------------------------------------------------------------------

  #1) agent-scaling loop constants: the model is run in replicate (replicates) at agent counts
  #   rising by delta_n up to n_max, and the level at which mean body condition reaches the
  #   carrying_capacity threshold is taken as nutritional carrying capacity

  .ncc_env$replicates <- 100L
  attr(.ncc_env$replicates, 'unit') <- 'iterations'
  attr(.ncc_env$replicates, 'source') <- NA
  attr(.ncc_env$replicates, 'full_name') <- 'Number of model replicates per agent level'

  .ncc_env$carrying_capacity <- 0.1466
  attr(.ncc_env$carrying_capacity, 'unit') <- 'Percent'
  attr(.ncc_env$carrying_capacity, 'source') <- NA
  attr(.ncc_env$carrying_capacity, 'full_name') <- 'Global mean ifbfat threshold for carrying capacity'

  .ncc_env$delta_n <- 25L
  attr(.ncc_env$delta_n, 'unit') <- 'Agents'
  attr(.ncc_env$delta_n, 'source') <- NA
  attr(.ncc_env$delta_n, 'full_name') <- 'Number of agents to add per iteration'

  .ncc_env$n_max <- 400L #CHANGE to match manuscript
  attr(.ncc_env$n_max, 'unit') <- 'Agents'
  attr(.ncc_env$n_max, 'source') <- NA
  attr(.ncc_env$n_max, 'full_name') <- 'Maximum number of agents'

  #2) the mutable agent count for the run currently in progress

  .ncc_env$n_agents <- 50 # denotes current number of agents, can be updated within the model

  # -----------------------------------------------------------------------------------------------------
  # time parameters
  # -----------------------------------------------------------------------------------------------------

  #1) hourly time step (t_delta) across a fixed summer-to-autumn window in the study time zone;
  #   t_start and t_end bound every hourly and daily sequence built elsewhere in the model

  .ncc_env$t_delta <- lubridate::hours(1)
  .ncc_env$t_start <- as.POSIXct("2025-07-01", tz = "America/Denver") # start of simulation
  .ncc_env$t_end <- as.POSIXct("2025-10-31", tz = "America/Denver") # end of simulation

  # -----------------------------------------------------------------------------------------------------
  # spatial parameters
  # -----------------------------------------------------------------------------------------------------

  #1) the projected coordinate system shared by all model rasters

  .ncc_env$epsg <- "EPSG:32612"

  #2) study-area coordinates, used only to derive the daylight schedule below

  .ncc_env$study_lat <- 43.74075
  attr(.ncc_env$study_lat, 'unit') <- 'decimal degrees'
  attr(.ncc_env$study_lat, 'source') <- 'Grand Teton summit — placeholder'
  attr(.ncc_env$study_lat, 'full_name') <- 'Study area latitude for solar position calculation'

  .ncc_env$study_lon <- -110.80252
  attr(.ncc_env$study_lon, 'unit') <- 'decimal degrees'
  attr(.ncc_env$study_lon, 'source') <- 'Grand Teton summit — placeholder'
  attr(.ncc_env$study_lon, 'full_name') <- 'Study area longitude for solar position calculation'

  # -----------------------------------------------------------------------------------------------------
  # daylight schedule — precomputed logical, one entry per hourly time step
  # geometric sunrise/sunset from study latitude and day-of-year (clock noon ~ solar noon)
  # positionally aligned with seq(t_start, t_end, by = t_delta); index by loop position t
  # NOTE: computed at load from current t_start/t_end/t_delta/study_lat; if those change,
  #   rerun .set_defaults() to refresh this vector
  # -----------------------------------------------------------------------------------------------------

  #1) hourly time sequence matching the model clock (t_start to t_end by t_delta)

  daylight_times <- seq(
    .ncc_env$t_start,
    .ncc_env$t_end,
    by = as.numeric(.ncc_env$t_delta, units = "secs")
  )

  #2) solar geometry from latitude and day of year: declination and the resulting half-day
  #   length (hours from solar noon to sunset), plus the local clock hour of each step

  daylight_lat_rad <- .ncc_env$study_lat * pi / 180
  daylight_doy <- lubridate::yday(daylight_times)
  daylight_decl <- 0.409 * sin(2 * pi / 365 * daylight_doy - 1.39)
  daylight_half_day <- (12 / pi) * acos(-tan(daylight_lat_rad) * tan(daylight_decl))
  daylight_hour <- lubridate::hour(daylight_times) + lubridate::minute(daylight_times) / 60

  #3) daylight flag per step: TRUE when the clock hour falls within the half-day of solar noon

  .ncc_env$is_daylight <- daylight_hour >= (12 - daylight_half_day) &
    daylight_hour <= (12 + daylight_half_day)

  #4) day length (hours) per step, equal to twice the half-day; depends only on day of year so
  #   it is constant within a day. this is exactly the quantity calc_energy_hif() formerly
  #   recomputed per call: (24/pi)*acos(...) = 2 * half-day

  .ncc_env$day_length <- 2 * daylight_half_day

  # -----------------------------------------------------------------------------------------------------
  # energy parameters
  # -----------------------------------------------------------------------------------------------------

  #1) forage energy: digestible energy of suitable forage (DE), the digestible-to-metabolizable
  #   conversion, and the resulting metabolizable energy (ME)

  .ncc_env$DE <- 12.98544 # this is the weighted mean of digestible energy of suitable forage biomass in vegetation transects.
  attr(.ncc_env$DE, 'unit') <- 'kJ/g'
  attr(.ncc_env$DE, 'source') <- NA
  attr(.ncc_env$DE, 'full_name') <- 'Digestible Energy'

  .ncc_env$DE_to_ME_conversion_factor <- 0.82
  attr(.ncc_env$DE_to_ME_conversion_factor, 'unit') <- 'Percent'
  attr(.ncc_env$DE_to_ME_conversion_factor, 'source') <- 'NRC 2007'
  attr(.ncc_env$DE_to_ME_conversion_factor, 'full_name') <- NA

  .ncc_env$ME <- .ncc_env$DE * .ncc_env$DE_to_ME_conversion_factor
  attr(.ncc_env$ME, 'unit') <- 'kJ/g'
  attr(.ncc_env$ME, 'source') <- NA
  attr(.ncc_env$ME, 'full_name') <- 'Metabolizable Energy'

  #2) slope-binned locomotion cost factors, from steep descent (d_10) through flat (f) to steep
  #   incline (i_10); calc_energy_loc() selects one per step by the signed slope of the step

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

  #3) heat increment of feeding (HIF) and the proportion of the diurnal period spent foraging,
  #   which together scale the daily feeding-related heat cost

  .ncc_env$HIF <- 1.799 # 0.43 kcal * kg^-1 * h^-1 (ewes) x 4.184 kJ/kcal; Chappel and Hudson 1978b
  attr(.ncc_env$HIF, 'unit') <- 'kJ * kg^-1 * h^-1'
  attr(.ncc_env$HIF, 'source') <- 'Chappel and Hudson 1978b'
  attr(.ncc_env$HIF, 'full_name') <- 'Heat Increment of Feeding'

  .ncc_env$prop_day_forage <- 0.72
  attr(.ncc_env$prop_day_forage, 'unit') <- 'proportion'
  attr(.ncc_env$prop_day_forage, 'source') <- 'Courtemanch et al. 2014'
  attr(.ncc_env$prop_day_forage, 'full_name') <- 'Proportion of the diurnal period spent feeding'

  #4) fat-reserve energetics: the energy density of fat (E_fat) and the catabolism and
  #   deposition efficiencies that convert a net daily energy balance into a change in fat mass

  .ncc_env$E_fat <- 39.5
  attr(.ncc_env$E_fat, 'unit') <- 'kJ/g'
  attr(.ncc_env$E_fat, 'source') <- 'Robbins 1993'
  attr(.ncc_env$E_fat, 'full_name') <- 'Energy Value of Fat Reserves'

  .ncc_env$fat_eff <- 0.65
  attr(.ncc_env$fat_eff, 'unit') <- 'Percent'
  attr(.ncc_env$fat_eff, 'source') <- 'Robbins 1993'
  attr(.ncc_env$fat_eff, 'full_name') <- 'Fat Catabolism Efficiency'

  .ncc_env$fat_dep <- 0.65
  attr(.ncc_env$fat_dep, 'unit') <- 'Percent'
  attr(.ncc_env$fat_dep, 'source') <- 'Robbins 1993'
  attr(.ncc_env$fat_dep, 'full_name') <- 'Fat Deposition Efficiency'

  # -----------------------------------------------------------------------------------------------------
  # body condition and body mass parameters
  # -----------------------------------------------------------------------------------------------------

  #1) energy threshold below which an agent is treated as dead

  .ncc_env$death <- 0
  attr(.ncc_env$death, 'unit') <- 'kJ'
  attr(.ncc_env$death, 'source') <- NA
  attr(.ncc_env$death, 'full_name') <- 'Energy Threshold for Death'

  #2) starting body mass (bm) and the among-individual ingesta-free body fat draw (ifbf),
  #   sampled once per agent at initialization

  .ncc_env$bm <- 56.91 # calculated from Teton capture data

  .ncc_env$ifbf <- function() max(rnorm(1, mean = 8.39, sd = 2.91) / 100, 0.01) # taken from Smiley et al. 2022, floored at 0.01 (1% body fat) to bar impossible negative/near-zero draws

  # -----------------------------------------------------------------------------------------------------
  # pregnancy and lactation parameters
  # -----------------------------------------------------------------------------------------------------

  #1) reproductive state and its lactation energy cost: rep_status is drawn per agent,
  #   j_post_partum sets the starting day post partum, lactation_modifier gives the daily cost
  #   multiplier indexed by days post partum, and calc_lactation_modifier performs the lookup

  .ncc_env$rep_status <- function() rbinom(1, 1, 0.678571429) #this is the proportion of captured ewes that showed some evidence of lactation
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
  .ncc_env$calc_lactation_modifier <- function(j_post_partum) {
    get_param("lactation_modifier")[as.integer(j_post_partum)]
  }

  # -----------------------------------------------------------------------------------------------------
  # plant trait parameters
  # -----------------------------------------------------------------------------------------------------

  #1) daily fractional regrowth of grazed forage, set from a 21-day half-life so a depleted
  #   cell recovers toward its reference biomass over the season

  .ncc_env$plant_regrowth_rate <- 1 - 0.5^(1/21) #21 day half-life on 42 day recovery period Osterheild 1992

  # -----------------------------------------------------------------------------------------------------
  # model calibration parameters
  # -----------------------------------------------------------------------------------------------------

  #1) the dry-matter-intake functional response: max_dmi is the intake asymptote and
  #   half_saturation the forage density at half-maximal intake, both tuned during calibration

  .ncc_env$max_dmi <- 372
  attr(.ncc_env$max_dmi, 'unit') <- 'g/hour'
  attr(.ncc_env$max_dmi, 'source') <- NA
  attr(.ncc_env$max_dmi, 'full_name') <- 'Maximum dry matter intake (DMI functional response asymptote)'

  .ncc_env$half_saturation <- 5000
  attr(.ncc_env$half_saturation, 'unit') <- 'g/cell'
  attr(.ncc_env$half_saturation, 'source') <- 'Spalinger and Hobbs 1992'
  attr(.ncc_env$half_saturation, 'full_name') <- 'Half-saturation constant for DMI functional response'

  #2) daily dry matter intake cap applied within simulate_forage: once an agent's cumulative
  #   intake for the day reaches this value, further foraging that day yields zero. set to NA
  #   to remove the cap entirely

  .ncc_env$max_daily_intake <- 1914.6 #1914.6 # this is 1512 + 2 sd (201.3)
  attr(.ncc_env$max_daily_intake, 'unit') <- 'g/day'
  attr(.ncc_env$max_daily_intake, 'source') <- 'Kraussman et al. 1988'
  attr(.ncc_env$max_daily_intake, 'full_name') <- 'Daily dry matter intake cap (NA disables the cap)'

  # -----------------------------------------------------------------------------------------------------
  # movement parameters — population-level multivariate normal distribution
  # Empirical values from the 22-animal iSSF fit (issf_mvn_params.rds, fit_issf_v1.R, 2026-07-23)
  # One movement model per agent is drawn from MVN(mvn_mu, mvn_sigma) in create_agents()
  # The model is a single make_issf_model() with day/night carried by tod_end_night interactions;
  #   the tentative Gamma and von Mises are pooled (one set), and the day/night difference is
  #   carried by the :tod_end_night_end interaction coefficients, not by separate distributions
  #
  # Naming convention (required for redistribution_kernel() to evaluate the model):
  #   geometry transforms use call form — log(sl_), cos(ta_) — so the kernel computes them
  #     from the sl_ and ta_ columns it generates for each candidate step; a bare log_sl_ /
  #     cos_ta_ name fails because the kernel looks for a column of that literal name
  #   habitat terms use the _end suffix (forage_biomass_end, etc.) to bind to the end-of-step
  #     covariate extraction (fun = extract_covariates(where = "both"))
  #   day/night is supplied as a constant raster layer named tod_end_night in the per-hour map
  #     (NOT via the covars argument, which does not reach the formula evaluation in this amt
  #     version); interactions therefore reference its end-of-step extraction, tod_end_night_end
  #   interactions are geometry-first (e.g. cos(ta_):tod_end_night_end) to match the kernel-tested
  #     form; the fitting code must produce these exact names or the MVN means map to wrong terms
  #
  # Dimensions (15), in order:
  #   log_shape, log_scale  : pooled tentative Gamma step-length params, LOG scale (exp after draw)
  #   log_kappa             : pooled tentative von Mises concentration, LOG scale (exp after draw)
  #   main effects (day, tod_end_night = 0):
  #     sl_, log(sl_)       : Gamma step-length correction coefficients
  #     cos(ta_)            : von Mises turn-angle correction coefficient
  #     forage_biomass_end, escape_terrain_end, canopy_cover_end : habitat selection coefficients
  #   night interaction offsets (added to the main effect when tod_end_night_end == 1):
  #     sl_:tod_end_night_end, log(sl_):tod_end_night_end, cos(ta_):tod_end_night_end,
  #     forage_biomass_end:tod_end_night_end, escape_terrain_end:tod_end_night_end,
  #     canopy_cover_end:tod_end_night_end
  # create_agents() exponentiates the three log dims to build the pooled tentative distributions;
  #   the remaining twelve coefficients are passed to make_issf_model() coefs as one vector
  # -----------------------------------------------------------------------------------------------------

  #1) number of candidate steps drawn per redistribution kernel at each move

  .ncc_env$n_candidates <- 25L
  attr(.ncc_env$n_candidates, 'unit') <- 'candidate steps'
  attr(.ncc_env$n_candidates, 'source') <- NA
  attr(.ncc_env$n_candidates, 'full_name') <- 'Number of candidate steps drawn per redistribution kernel'

  #2) empirical population MVN from the 22 converged per-individual iSSF fits
  #   (back-transformed coefficients, raw g/cell forage and metre distance-to-escape units).
  #   the three log dims are on the log scale; day dims are the day slopes and the
  #   :tod_end_night_end dims are night minus day offsets, matching the kernel's coding.
  #   names and dimnames must stay exactly as below or validate_move_coefs() stops the run

  mvn_names <- c("log_shape", "log_scale", "log_kappa",
                 "sl_", "log(sl_)", "cos(ta_)",
                 "forage_biomass_end", "escape_terrain_end", "canopy_cover_end",
                 "sl_:tod_end_night_end", "log(sl_):tod_end_night_end", "cos(ta_):tod_end_night_end",
                 "forage_biomass_end:tod_end_night_end", "escape_terrain_end:tod_end_night_end",
                 "canopy_cover_end:tod_end_night_end")

  #3) empirical mean vector

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

  #4) empirical among-individual covariance (raw, no ridge), filled by row; the matrix is
  #   symmetric, so row order matches mvn_names and the lower triangle mirrors the upper

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
