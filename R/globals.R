#' @keywords internal
.ncc_env <- new.env(parent = emptyenv())

#' Set default model parameters
#' @keywords internal
.set_defaults <- function() {

  # -----------------------------------------------------------------------------------------------------
  # simulation parameters
  # -----------------------------------------------------------------------------------------------------

  .ncc_env$replicates <- 100L
  attr(.ncc_env$replicates, 'unit') <- 'iterations'
  attr(.ncc_env$replicates, 'source') <- NA
  attr(.ncc_env$replicates, 'full_name') <- 'Number of model replicates per agent level'

  .ncc_env$carrying_capacity <- 0.1466
  attr(.ncc_env$carrying_capacity, 'unit') <- 'Percent'
  attr(.ncc_env$carrying_capacity, 'source') <- NA
  attr(.ncc_env$carrying_capacity, 'full_name') <- 'Global mean ifbfat threshold for carrying capacity'

  .ncc_env$delta_n <- 5L
  attr(.ncc_env$delta_n, 'unit') <- 'Agents'
  attr(.ncc_env$delta_n, 'source') <- NA
  attr(.ncc_env$delta_n, 'full_name') <- 'Number of agents to add per iteration'

  .ncc_env$n_max <- 400L #CHANGE to match manuscript
  attr(.ncc_env$n_max, 'unit') <- 'Agents'
  attr(.ncc_env$n_max, 'source') <- NA
  attr(.ncc_env$n_max, 'full_name') <- 'Maximum number of agents'

  .ncc_env$n_agents <- 50 # denotes current number of agents, can be updated within the model

  # -----------------------------------------------------------------------------------------------------
  # time parameters
  # -----------------------------------------------------------------------------------------------------

  .ncc_env$t_delta <- lubridate::hours(1)
  .ncc_env$t_start <- as.POSIXct("2025-07-01", tz = "America/Denver") # start of simulation
  .ncc_env$t_end <- as.POSIXct("2025-10-31", tz = "America/Denver") # end of simulation

  # -----------------------------------------------------------------------------------------------------
  # spatial parameters
  # -----------------------------------------------------------------------------------------------------

  .ncc_env$epsg <- "EPSG:32612"

  .ncc_env$study_lat <- 43.74075
  attr(.ncc_env$study_lat, 'unit') <- 'decimal degrees'
  attr(.ncc_env$study_lat, 'source') <- 'Grand Teton summit — placeholder'
  attr(.ncc_env$study_lat, 'full_name') <- 'Study area latitude for solar position calculation'

  .ncc_env$study_lon <- -110.80252
  attr(.ncc_env$study_lon, 'unit') <- 'decimal degrees'
  attr(.ncc_env$study_lon, 'source') <- 'Grand Teton summit — placeholder'
  attr(.ncc_env$study_lon, 'full_name') <- 'Study area longitude for solar position calculation'

  # -----------------------------------------------------------------------------------------------------
  # energy parameters
  # -----------------------------------------------------------------------------------------------------

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

  .ncc_env$HIF <- 0.0018
  attr(.ncc_env$HIF, 'unit') <- 'kJ * kg^-1 * h^-1'
  attr(.ncc_env$HIF, 'source') <- 'Dailey and Hobbs 1989'
  attr(.ncc_env$HIF, 'full_name') <- 'Heat Increment of Feeding'

  .ncc_env$prop_day_forage <- 0.72
  attr(.ncc_env$prop_day_forage, 'unit') <- 'proportion'
  attr(.ncc_env$prop_day_forage, 'source') <- 'Courtemanch et al. 2014'
  attr(.ncc_env$prop_day_forage, 'full_name') <- 'Proportion of the diurnal period spent feeding'

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

  .ncc_env$death <- 0
  attr(.ncc_env$death, 'unit') <- 'kJ'
  attr(.ncc_env$death, 'source') <- NA
  attr(.ncc_env$death, 'full_name') <- 'Energy Threshold for Death'

  .ncc_env$bm <- 56.91 # calculated from Teton capture data

  .ncc_env$ifbf <- function() rnorm(1, mean = 8.39, sd = 2.91) / 100 # taken from Smiley et al. 2022

  # -----------------------------------------------------------------------------------------------------
  # pregnancy and lactation parameters
  # -----------------------------------------------------------------------------------------------------

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

  .ncc_env$plant_regrowth_rate <- 1 - 0.5^(1/21) #21 day half-life on 42 day recovery period Osterheild 1992

  # -----------------------------------------------------------------------------------------------------
  # model calibration parameters
  # -----------------------------------------------------------------------------------------------------

  .ncc_env$half_saturation <- 60
  attr(.ncc_env$half_saturation, 'unit') <- 'g/cell'
  attr(.ncc_env$half_saturation, 'source') <- 'Spalinger and Hobbs 1992'
  attr(.ncc_env$half_saturation, 'full_name') <- 'Half-saturation constant for DMI functional response'

  # -----------------------------------------------------------------------------------------------------
  # movement parameters — population-level multivariate normal distribution
  # PLACEHOLDER values pending empirical iSSF estimates
  # One movement model per agent is drawn from MVN(mvn_mu, mvn_sigma) in create_agents()
  # The model is a single make_issf_model() with tod_end_ interactions (day = reference level);
  #   the tentative Gamma and von Mises are pooled (one set), and the day/night difference is
  #   carried by the :tod_end_night interaction coefficients, not by separate distributions
  #
  # Dimensions (15), in order:
  #   log_shape, log_scale  : pooled tentative Gamma step-length params, LOG scale (exp after draw)
  #   log_kappa             : pooled tentative von Mises concentration, LOG scale (exp after draw)
  #   main effects (day, reference level):
  #     sl_, log_sl_        : Gamma step-length correction coefficients
  #     cos_ta_             : von Mises turn-angle correction coefficient
  #     forage_biomass, escape_terrain, canopy_cover : habitat selection coefficients
  #   night interaction offsets (added to the main effect when tod_end_ == night):
  #     tod_end_night:sl_, tod_end_night:log_sl_, tod_end_night:cos_ta_,
  #     tod_end_night:forage_biomass, tod_end_night:escape_terrain, tod_end_night:canopy_cover
  # Names are factor-first to match a fit_issf() formula of tod_end_ * (movement + habitat terms);
  #   they must match names(coefficients) of the fitted model exactly or the kernel update lookups fail
  # create_agents() exponentiates the three log dims to build the pooled tentative distributions;
  #   the remaining twelve coefficients are passed to make_issf_model() coefs as one vector
  # -----------------------------------------------------------------------------------------------------

  mvn_names <- c("log_shape", "log_scale", "log_kappa",
                 "sl_", "log_sl_", "cos_ta_",
                 "forage_biomass", "escape_terrain", "canopy_cover",
                 "tod_end_night:sl_", "tod_end_night:log_sl_", "tod_end_night:cos_ta_",
                 "tod_end_night:forage_biomass", "tod_end_night:escape_terrain",
                 "tod_end_night:canopy_cover")

  # mean vector (placeholder); the three log dims are on the log scale
  # main effects reproduce the prior daytime values; offsets reproduce the prior night values
  #   (night = main + offset), e.g. cos_ta_ day 0.50, night 0.20, so offset -0.30
  .ncc_env$mvn_mu <- c(
    "log_shape"                    =  log(2.0),
    "log_scale"                    =  log(22),
    "log_kappa"                    =  log(0.35),
    "sl_"                          =  0.00,
    "log_sl_"                      =  0.00,
    "cos_ta_"                      =  0.50,
    "forage_biomass"               =  0.40,
    "escape_terrain"               =  0.30,
    "canopy_cover"                 = -0.20,
    "tod_end_night:sl_"            =  0.00,
    "tod_end_night:log_sl_"        =  0.00,
    "tod_end_night:cos_ta_"        = -0.30,
    "tod_end_night:forage_biomass" = -0.30,
    "tod_end_night:escape_terrain" =  0.20,
    "tod_end_night:canopy_cover"   =  0.30
  )

  # covariance matrix (placeholder) — diagonal, off-diagonals zero pending empirical structure
  # log dims: variance 0.09 (= prior sdlog 0.3, squared); coefficient dims: variance 0.10
  mvn_var <- c(0.09, 0.09, 0.09, rep(0.10, 12))
  .ncc_env$mvn_sigma <- diag(mvn_var)
  dimnames(.ncc_env$mvn_sigma) <- list(mvn_names, mvn_names)


}

.onLoad <- function(libname, pkgname) {
  .set_defaults()
}

#' Get a global model parameter
#' @param param Character string. Name of the parameter.
#' @importFrom lubridate hours yday
#' @export
get_param <- function(param) {
  .ncc_env[[param]]
}

#' Set or override a global model parameter
#' @param param Character string. Name of the parameter.
#' @param value New value.
#' @export
set_param <- function(param, value) {
  .ncc_env[[param]] <- value
}

#' Draw a parameter value
#' @param param Character string. Name of the parameter.
#' @export
draw_param <- function(param) {
  val <- .ncc_env[[param]]
  if (is.function(val)) val() else val
}
