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
  
  # Validate s1_data input
  if(missing(s1_data)) {
    cli::cli_abort(c(
      "!" = "{.arg s1_data} is required but missing.",
      "i" = "Please provide Sentinel-1 raster data."
    ))
  }
  
  if(!inherits(s1_data, c("SpatRaster", "RasterStack", "RasterBrick"))) {
    cli::cli_abort(c(
      "!" = "{.arg s1_data} must be a raster object.",
      "i" = "You provided an object of class {.cls {class(s1_data)}}.",
      "i" = "Expected classes: {.cls SpatRaster}, {.cls RasterStack}, or {.cls RasterBrick}."
    ))
  }
  
  # Check if s1_data has enough layers/bands
  n_layers <- terra::nlyr(s1_data)
  if(n_layers < 2) {
    cli::cli_warn(c(
      "!" = "Sentinel-1 data has only {n_layers} layer{?s}.",
      "i" = "Typical Sentinel-1 data should have at least 2 polarizations (VV, VH).",
      "i" = "Classification accuracy may be reduced."
    ))
  }
  
  # Validate rf_model
  if(missing(rf_model)) {
    cli::cli_abort(c(
      "!" = "{.arg rf_model} is required but missing.",
      "i" = "Please provide a trained Random Forest model."
    ))
  }
  
  if(!inherits(rf_model, "randomForest")) {
    cli::cli_abort(c(
      "!" = "{.arg rf_model} must be a randomForest object.",
      "i" = "You provided an object of class {.cls {class(rf_model)}}."
    ))
  }
  
  # Check if model variables match data layers
  model_vars <- rownames(rf_model$importance)
  data_vars <- names(s1_data)
  
  missing_vars <- setdiff(model_vars, data_vars)
  if(length(missing_vars) > 0) {
    cli::cli_abort(c(
      "!" = "Model requires variable{?s} not found in data: {.val {missing_vars}}",
      "i" = "Available variables in data: {.val {data_vars}}",
      "i" = "Please ensure data contains all required model variables."
    ))
  }
  
  # Check if output directory is writable
  parent_dir <- dirname(output_dir)
  if(!dir.exists(parent_dir)) {
    cli::cli_warn(c(
      "!" = "Parent directory {.path {parent_dir}} does not exist.",
      "i" = "Attempting to create directory structure."
    ))
  }
  
  # Validate CRS
  if(!is.character(crs.def) || length(crs.def) != 1) {
    cli::cli_abort(c(
      "!" = "{.arg crs.def} must be a single character string.",
      "i" = "You provided: {.val {crs.def}}"
    ))
  }
  
  # Check CRS validity
  tryCatch({
    terra::crs(crs.def)
  }, error = function(e) {
    cli::cli_abort(c(
      "!" = "Invalid CRS definition: {.val {crs.def}}",
      "i" = "Error: {e$message}"
    ))
  })
  
  
  # Success message for model validation
  if(inherits(rf_model, "randomForest")) {
    cli::cli_bullets(c("v" = "Random Forest model validated successfully"))
    cli::cli_h1("Random Forest Model Information")
    cli::cli_text("The model has been trained with the following variables:")
    cli::cli_ul(rownames(rf_model$importance))
    
    # Model performance info
    if(!is.null(rf_model$err.rate)) {
      oob_error <- rf_model$err.rate[rf_model$ntree, "OOB"]
      cli::cli_bullets(c("i" = "Out-of-bag error rate: {.val {round(oob_error * 100, 2)}%}"))
    }
  }
  
  # Create output directory with error handling
  tryCatch({
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cli::cli_bullets(c("v" = "Created output directory: {.path {output_dir}}"))
    }
  }, error = function(e) {
    cli::cli_abort(c(
      "!" = "Failed to create output directory: {.path {output_dir}}",
      "i" = "Error: {e$message}"
    ))
  })
  
  
  
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