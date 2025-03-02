jh_calculate_distance_between_2_points_function <- function(point_1, point_2){
  sqrt((point_1[1] - point_2[1])^2 + 
         (point_1[2] - point_2[2])^2)
}


jh_calculate_vpa_from_xray_data_function <- function(fem_head_center = c(0,0), 
                                                   vertebral_centroid = c(0.15, 0.3), 
                                                   spine_facing = "left",
                                                   pelvic_tilt = 15){
  
  
  fem_head_to_centroid_length <- jh_calculate_distance_between_2_points_function(point_1 = vertebral_centroid, 
                                                                                 point_2 = fem_head_center)
  
  fem_head_to_centroid_x_length <- jh_calculate_distance_between_2_points_function(point_1 = vertebral_centroid, 
                                                                                   point_2 = c(fem_head_center[1],
                                                                                               vertebral_centroid[2]))
  
  tilt_orientation_modifier <- case_when(
    spine_facing == "left" & fem_head_center[[1]] > vertebral_centroid[[1]] ~ 1,
    spine_facing == "left" & fem_head_center[[1]] < vertebral_centroid[[1]] ~ -1,
    spine_facing == "right" & fem_head_center[[1]] > vertebral_centroid[[1]] ~ -1,
    spine_facing == "right" & fem_head_center[[1]] < vertebral_centroid[[1]] ~ 1
  )
  
  vertebral_tilt <- asin(fem_head_to_centroid_x_length/fem_head_to_centroid_length)*180/pi*tilt_orientation_modifier
  
  vpa <- pelvic_tilt + vertebral_tilt 
  
  return(vpa)
  
}

jh_calculate_segment_angle_between_vertebrae_from_coordinates_function <- function(sp_upper, sa_upper, sp_lower, sa_lower) {
  vec1 <- sa_upper - sp_upper  # Vector for upper vertebra
  vec2 <- sa_lower - sp_lower  # Vector for lower vertebra
  
  # Dot product of the two vectors
  dot_product <- sum(vec1 * vec2)
  
  # Magnitudes of the vectors
  mag_vec1 <- sqrt(sum(vec1^2))
  mag_vec2 <- sqrt(sum(vec2^2))
  
  # Calculate the angle using the dot product formula
  cos_theta <- dot_product / (mag_vec1 * mag_vec2)
  
  # Ensure the value is within the valid range for acos (to avoid numerical errors)
  cos_theta <- min(1, max(-1, cos_theta))
  
  # Compute the angle in radians
  theta_radians <- acos(cos_theta)
  
  # Calculate the cross product to determine if it's lordotic or kyphotic
  cross_product <- (vec1[1] * vec2[2]) - (vec1[2] * vec2[1])
  
  # If the cross product is negative, it is kyphotic (angle opens posteriorly), otherwise lordotic
  
  if(sp_upper[1] < sa_upper[1]){
    sign <- ifelse(cross_product < 0, 1, -1)
  }else{
    sign <- ifelse(cross_product < 0, -1, 1)
  }
  
  # Convert the angle to degrees and apply the sign for lordosis/kyphosis
  theta_degrees <- theta_radians * (180 / pi) * sign
  
  return(theta_degrees)
}


jh_calculate_perpendicular_angle <- function(coord1_x, coord1_y, coord2_x, coord2_y) {
  # Step 1: Calculate the slope of the line connecting the two points
  
  coord1_x <- if_else(is.na(coord1_x), 0, coord1_x)
  coord1_y <- if_else(is.na(coord1_y), 0, coord1_y)
  coord2_x <- if_else(is.na(coord2_x), 0, coord2_x)
  coord2_y <- if_else(is.na(coord2_y), 0, coord2_y)
  
  slope <- (coord2_y - coord1_y) / (coord2_x - coord1_x)
  
  # Step 2: Get the negative reciprocal for the perpendicular slope
  if (slope != 0) {
    perpendicular_slope <- -1 / slope
  } else {
    # If the original slope is zero (horizontal line), perpendicular line is vertical
    perpendicular_slope <- Inf
  }
  
  # Step 3: Calculate the angle in radians
  # If the slope is infinite (vertical line), the angle is 90 degrees
  if (is.infinite(perpendicular_slope)) {
    angle_radians <- pi / 2
  } else {
    angle_radians <- atan(perpendicular_slope)
  }
  
  # Convert the angle to degrees
  angle_degrees <- angle_radians * 180 / pi
  
  return(angle_degrees)
}

jh_calculate_3_point_angle_function <- function(ventral_point_coord = c(0,0), 
                                                vertex_point_coord = c(2, 2),
                                                posterior_point_coord = c(2, 0),
                                                spine_orientation = "left") {
  # Compute vectors
  v1 <- ventral_point_coord - vertex_point_coord
  v2 <- posterior_point_coord - vertex_point_coord
  
  # Compute dot product and magnitudes
  dot_product <- sum(v1 * v2)
  mag_v1 <- sqrt(sum(v1^2))
  mag_v2 <- sqrt(sum(v2^2))
  
  # Compute angle in radians
  angle_rad <- acos(dot_product / (mag_v1 * mag_v2))
  
  # Convert to degrees
  angle_deg <- angle_rad * (180 / pi)
  
  return(angle_deg)
}


jh_rotate_polygon_around_centroid <- function(polygon, angle_degrees) {
  # Step 1: Get the centroid of the polygon
  centroid <- st_centroid(polygon)
  centroid_coords <- st_coordinates(centroid)[1, 1:2]  # Ensure this is a 2D vector (x, y)
  
  # Step 2: Translate the polygon to have the centroid at the origin
  coords <- st_coordinates(polygon)[, 1:2]  # Extract the polygon's coordinates
  translated_coords <- sweep(coords, 2, centroid_coords)  # Subtract centroid from coordinates
  
  # Step 3: Create a rotation matrix
  angle_radians <- angle_degrees * pi / 180  # Convert angle to radians
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                              sin(angle_radians),  cos(angle_radians)),
                            ncol = 2, byrow = TRUE)
  
  # Step 4: Apply the rotation matrix to the translated coordinates
  rotated_coords <- translated_coords %*% rotation_matrix
  
  # Step 5: Translate the polygon back to its original position
  final_coords <- sweep(rotated_coords, 2, centroid_coords, "+")
  
  # Step 6: Rebuild the rotated polygon
  rotated_polygon <- st_polygon(list(final_coords))
  
  return(st_sfc(rotated_polygon, crs = st_crs(polygon)))
}


jh_reformat_vert_tibble <- function(tibble_input) {
  # Convert the tibble to a named list with each point as a vector
  named_list <- list(
    sp = c(tibble_input$x[tibble_input$vert_point == "sp"], tibble_input$y[tibble_input$vert_point == "sp"]),
    sa = c(tibble_input$x[tibble_input$vert_point == "sa"], tibble_input$y[tibble_input$vert_point == "sa"]),
    ia = c(tibble_input$x[tibble_input$vert_point == "ia"], tibble_input$y[tibble_input$vert_point == "ia"]),
    ip = c(tibble_input$x[tibble_input$vert_point == "ip"], tibble_input$y[tibble_input$vert_point == "ip"]),
    centroid = c(tibble_input$x[tibble_input$vert_point == "centroid"], tibble_input$y[tibble_input$vert_point == "centroid"])
  )
  return(named_list)
}



rotate_spine_function <- function(spine_df, angle_degrees, point_of_rotation = c(0, 0)) {
  # Convert angle to radians
  angle_rad <- angle_degrees * pi / 180
  
  # Extract rotation center
  x_center <- point_of_rotation[1]
  y_center <- point_of_rotation[2]
  
  # Rotation matrix components
  cos_theta <- cos(angle_rad)
  sin_theta <- sin(angle_rad)
  
  # Apply rotation to x, y coordinates relative to the rotation point
  spine_df <- spine_df %>%
    mutate(
      x_shifted = x - x_center,
      y_shifted = y - y_center,
      x_rot = x_shifted * cos_theta - y_shifted * sin_theta + x_center,
      y_rot = x_shifted * sin_theta + y_shifted * cos_theta + y_center
    ) %>%
    select(spine_level, vert_point, x_rot, y_rot) %>%
    rename(x = x_rot, y = y_rot)
  
  return(spine_df)
}


jh_construct_geoms_after_planning_function <- function(vertebral_level_tibble, buffer_amount = 0){
  coord_list_df <- vertebral_level_tibble %>%
    mutate(vert_point_coord = map2(.x = x, .y = y, .f = ~ c(.x, .y))) %>%
    select(vert_point, vert_point_coord)
  
  vert_coord_point_list <- coord_list_df$vert_point_coord
  
  names(vert_coord_point_list) <- coord_list_df$vert_point 
  
  if(unique(vertebral_level_tibble$spine_level) == "sacrum"){
    vertebral_geom <- st_polygon(list(rbind(vert_coord_point_list$s1_anterior_superior,
                                            vert_coord_point_list$sac_inf_1, 
                                            vert_coord_point_list$s1_posterior_superior, 
                                            vert_coord_point_list$s1_anterior_superior)))
  }else{
    vertebral_geom <- st_polygon(list(rbind(vert_coord_point_list$sp, 
                                            vert_coord_point_list$sa,
                                            vert_coord_point_list$ia, 
                                            vert_coord_point_list$ip, 
                                            vert_coord_point_list$sp)))
  }
  
  if(buffer_amount > 0){
    geom <- jh_safely_buffer_vert_function(vert_geom = vertebral_geom, 
                                           buffer_amount = buffer_amount)
    return(geom)
  }else{
    return(vertebral_geom)
  }
  
}
