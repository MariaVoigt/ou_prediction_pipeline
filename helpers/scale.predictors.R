suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

scale.predictors.observation <- function(predictor_names_for_scaling,
                             predictor_names_add,
                             predictors,
                             geography){

  # SCALE PREDICTORS
  for (predictor_name in predictor_names_for_scaling){
    predictors[predictors$predictor == predictor_name,
               "value" ] <-
      as.numeric(as.vector(scale(predictors[predictors$predictor == predictor_name,
                                            "unscaled_value"])))
  }

  # Rename here
  predictors_obs <- predictors %>%
    dplyr::filter(predictor %in% predictor_names_for_scaling) %>%
    dcast(id + unscaled_year ~ predictor,  value.var = "value") %>%
    inner_join(geography, by = "id")



  predictors_obs_unscaled <- predictors %>%
    dplyr::filter(predictor %in% predictor_names_for_scaling)%>%
    dcast(id + unscaled_year ~ predictor,  value.var = "unscaled_value") %>%
    dplyr::select(-unscaled_year)


  # all predictors minus id
  names(predictors_obs_unscaled)[-1] <- paste0("unscaled_",
                                               names(predictors_obs_unscaled)[-1])

  predictors_obs <- left_join(predictors_obs, predictors_obs_unscaled, by = "id")


  # scale the additional predictors
  # manually adding the additional columns
  # NOT SO NICE HARD-CODING
  predictors_obs$year <- NA
  predictors_obs$x_center <- NA
  predictors_obs$y_center <- NA

  for (predictor_name in predictor_names_add){
    predictors_obs[ , predictor_name] <- as.numeric(as.vector(scale(predictors_obs[ ,
                                                                                    paste0("unscaled_", predictor_name) ])))
  }


  predictors_obs <- predictors_obs %>%
    dplyr::select(dplyr::all_of(c("id",
                                  predictor_names_for_scaling,
                                  predictor_names_add,
                                  paste0("unscaled_", predictor_names_for_scaling),
                                  paste0("unscaled_", predictor_names_add)))) %>%
    as.data.frame()
  return(predictors_obs)
}





scale.predictors.grid <- function(predictor_names_for_scaling,
                                         predictor_names_add,
                                         predictors,
                                         predictors_obs,
                                         geography){


  for (predictor_name in predictor_names_for_scaling){
    mean_predictor_obs <- mean(predictors_obs[ , paste0("unscaled_", predictor_name)])
    sd_predictor_obs <- sd(predictors_obs[ , paste0("unscaled_", predictor_name)])
    predictors[predictors$predictor == predictor_name,
               "value" ] <- (predictors[predictors$predictor == predictor_name,
                                        "unscaled_value" ] - mean_predictor_obs) /
      sd_predictor_obs

  }


  # cast it to wide
  predictors_grid <- dplyr::filter(predictors, predictor %in% predictor_names_for_scaling)  %>%
    dcast(id + year ~ predictor,  value.var = "value") %>%
    dplyr::select(-year)


  predictors_grid_unscaled <- dplyr::filter(predictors,
                                            predictor %in% predictor_names_for_scaling ) %>%
    dcast(id + year ~ predictor,  value.var = "unscaled_value")

  names(predictors_grid_unscaled)[-1] <- paste0("unscaled_", names(predictors_grid_unscaled)[-1])

  # join with geography to have x and y-center

  predictors_grid <- predictors_grid %>%
    left_join(predictors_grid_unscaled, by = "id") %>%
    left_join(geography, by = "id")

  # NOT SO NICE HARD-CODING THIS HERE...WORKAROUND
 predictors_grid$year <- NA
 predictors_grid$x_center <- NA
 predictors_grid$y_center <- NA

  for (predictor_name in predictor_names_add){
    mean_predictor_obs <- mean(predictors_obs[ , paste0("unscaled_", predictor_name)])
    sd_predictor_obs <- sd(predictors_obs[ , paste0("unscaled_", predictor_name)])
    predictors_grid[ ,
                     predictor_name] <- (predictors_grid[ ,
                                         paste0("unscaled_", predictor_name) ] -
                                           mean_predictor_obs) /
                                           sd_predictor_obs
  }

 predictors_grid <- predictors_grid %>%
   dplyr::select(dplyr::all_of(c("id",
                                 predictor_names_for_scaling,
                                 predictor_names_add,
                                 paste0("unscaled_", predictor_names_for_scaling),
                                 paste0("unscaled_", predictor_names_add)))) %>%
   as.data.frame()

  return(predictors_grid)
}
