#' Title
#'
#' @param s1_data 
#' @param rf_model 
#' @param output_dir 
#' @param plots 
#' @param save 
#' @param crs.def 
#'
#' @return
#' @export
#'
#' @examples
#' 
classify_water_surface <- function(s1_data, rf_model, output_dir = "./RESULTS/Test", save = FALSE, plots = c(1, 2), crs.def = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs") {
  
  if(class(rf_model)[2]=="randomForest"){
    cli::cli_bullets(c("v"= "The first column is a {.cls {class(rf_model)[2]}} class"))
    cli::cli_h1("Random Forest Model Information \n The model has been trained with the following parameters:")
    cli::cat_print(rownames(rf_model$importance))
  } else{
    cli::cli_abort(c(
      "!" = "{.arg rf_model} must be a randomForest object.",
      "i" = "You provided an object of class {.cls {class(rf_model)}}."
    ))
  }
  
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) dir.create(output_dir)

  radar_df <- as.data.frame(s1_data, xy = TRUE)
  
  # Store NA positions to restore them later
  na_positions <- which(is.na(radar_df[,5]))
  
  # Remove NA values for prediction
  radar_df_clean <- na.omit(radar_df)
  
  # Rename columns to match training data
  names(radar_df_clean) <- c("X","Y", names(s1_data))
  
  # Predict water surface
  require(randomForest)
  water_prediction <- predict(rf_model, radar_df_clean, type = "prob")
  
  # Create a vector to hold all predictions (including NAs)
  all_predictions <- matrix(NA, nrow(radar_df), 3)
  all_predictions <- as.data.frame(all_predictions)
  all_predictions[-na_positions, 3] <- water_prediction[,2]
  all_predictions[,c(1,2)] <- radar_df[,c(1,2)]
  raster.output <- terra::rast(all_predictions, type="xyz", crs = crs.def)
  
  water_binary <- raster.output > 0.5  # Probability of water > 0.5
  
  # Save water classification
  terra::writeRaster(water_binary, 
              filename = file.path(output_dir, "water_classification_p.tif"), 
              overwrite = TRUE)
  
  # Calculate water surface area
  pixel_area <- raster::res(water_binary)[1] * raster::res(water_binary)[2]  # m²
  water_pixels <- sum(values(water_binary) == 1, na.rm = TRUE)
  water_area_m2 <- water_pixels * pixel_area
  water_area_km2 <- water_area_m2 / 1000000
  
  # Create a simple report
  report <- data.frame(
    total_pixels = length(na.omit(values(water_binary))),
    water_pixels = water_pixels,
    water_percentage = (water_pixels / length(na.omit(values(water_binary)))) * 100,
    water_area_m2 = water_area_m2,
    water_area_km2 = water_area_km2
  )
  
  # Save report
  write.csv(report, file.path(output_dir, "water_surface_report.csv"), row.names = FALSE)
  
  # Visualize results
  mapview(water_binary, col.regions = c("white", "blue"), legend = TRUE)
  
  return(list(classification = water_binary, report = report))
}