
jh_calculate_distance_between_2_points_function <- function(point_1, point_2){
  sqrt((point_1[1] - point_2[1])^2 + 
         (point_1[2] - point_2[2])^2)
}

# Function to create vertebra as an sf polygon
create_vertebra <- function(centroid_x = 0, centroid_y = 0, width = 5, height = 4) {
  half_height <- height / 2
  half_width <- width / 2
  
  # Define the four corners of the vertebra before rotation
  sp <- c(centroid_x - half_width, centroid_y + half_height)
  sa <- c(centroid_x + half_width, centroid_y + half_height)
  ia <- c(centroid_x + half_width, centroid_y - half_height)
  ip <- c(centroid_x - half_width, centroid_y - half_height)
  
  vert_body <- rbind(sp, sa, ia, ip, sp) ## binds the corners to make a square
  
  return(st_polygon(list(vert_body)))
}

compute_perpendicular_angle <- function(coord1_x, coord1_y, coord2_x, coord2_y) {
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



rotate_polygon_around_centroid <- function(polygon, angle_degrees) {
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


find_sacrum_inf_point_function <- function(s1_posterior_sup = c(0,0), 
                                           s1_midpoint = c(1, 1), 
                                           femoral_heads_center =  c(-1, 0), 
                                           spine_facing = "left") {
  
  length_s1_to_hips <- jh_calculate_distance_between_2_points_function(s1_midpoint, femoral_heads_center)
  
  d_AB <- jh_calculate_distance_between_2_points_function(s1_posterior_sup, s1_midpoint)
  
  # Extract coordinates for points A and B
  x1 <- s1_posterior_sup[1]
  y1 <- s1_posterior_sup[2]
  x2 <- s1_midpoint[1]
  y2 <- s1_midpoint[2]
  
  # Calculate the distance between A and B (sanity check)
  dist_AB <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
  # Check if provided AB distance matches the actual one
  if (abs(dist_AB - d_AB) > 1e-6) {
    stop("Provided distance between A and B does not match the actual distance.")
  }
  
  # Calculate the slope of line AB
  dx <- x2 - x1
  dy <- y2 - y1
  
  # If AB is vertical (dx = 0), then BC is horizontal
  if (dx == 0) {
    C1 <- c(x2 + length_s1_to_hips, y2)  # move right
    C2 <- c(x2 - length_s1_to_hips, y2)  # move left
  } 
  # If AB is horizontal (dy = 0), then BC is vertical
  else if (dy == 0) {
    C1 <- c(x2, y2 + length_s1_to_hips)  # move up
    C2 <- c(x2, y2 - length_s1_to_hips)  # move down
  } 
  else {
    # Slope of the perpendicular line (negative reciprocal)
    slope_perpendicular <- -dx / dy
    
    # Calculate the possible coordinates for C using the perpendicular direction
    angle <- atan(slope_perpendicular)
    
    # Move in both directions along the perpendicular at a distance length_s1_to_hips
    C1_x <- x2 + length_s1_to_hips * cos(angle)
    C1_y <- y2 + length_s1_to_hips * sin(angle)
    
    C2_x <- x2 - length_s1_to_hips * cos(angle)
    C2_y <- y2 - length_s1_to_hips * sin(angle)
    
    C1 <- c(C1_x, C1_y)
    C2 <- c(C2_x, C2_y)
  }
  
  # Return both possible coordinates for C
  # tibble(C1 = list(C1), C2 = list(C2))
  if(str_to_lower(spine_facing) == "left"){
   
    sacrum_inf_x <- if_else(C1_x < femoral_heads_center[1], C2_x, C1_x)
    sacrum_inf_y <- if_else(C1_y < s1_midpoint[2], C1_y, C2_y) 
    
  }else{
    sacrum_inf_x <- if_else(C1_x < s1_midpoint[1], C1_x, C2_x)
    sacrum_inf_y <- if_else(C1_y < s1_midpoint[2], C1_y, C2_y) 
  }

  c(sacrum_inf_x, sacrum_inf_y)
}


compute_inferior_corner <- function(x, y, inferior_x, inferior_y, return_x_or_y = "x") {
  # Calculate the new x and y coordinates for the point 20% along the line
  new_x <- inferior_x + 0.2 * (x - inferior_x)
  new_y <- inferior_y + 0.2 * (y - inferior_y)
  
  # Return the new coordinates as a vector
  # return(c(new_x, new_y))
  if(return_x_or_y == "x"){
    return(new_x)
  }else{
    return(new_y)
  }
}

build_spine_from_coordinates_function <- function(femoral_head_center = c(0,0),
                                                  s1_anterior_superior = c(1, 1),
                                                  s1_posterior_superior = c(1, 2),
                                                  centroid_df = tibble(), 
                                                  spine_facing = "right"){
  
  # spine_facing_modifier <- if_else(str_to_lower(spine_facing) == "right", -1, 1)
  
  s1_endplate_width <- jh_calculate_distance_between_2_points_function(point_1 = s1_anterior_superior, point_2 = s1_posterior_superior)
  
  
  vertebral_heights_df <- centroid_df %>%
    filter(str_detect(spine_point, "centroid|center")) %>%
    mutate(distance_to_cranial_vertebra = pmap(.l = list(..1 = x, ..2 = y, ..3 = lead(x), ..4 = lead(y)), 
                                               .f = ~ jh_calculate_distance_between_2_points_function(point_1 = c(..1, ..2),
                                                                                                      point_2 = c(..3, ..4)))) %>%
    unnest() %>%
    mutate(distance_to_caudal_vertebra = pmap(.l = list(..1 = x, ..2 = y, ..3 = lag(x), ..4 = lag(y)), 
                                              .f = ~ jh_calculate_distance_between_2_points_function(point_1 = c(..1, ..2),
                                                                                                     point_2 = c(..3, ..4)))) %>%
    unnest() %>%
    mutate(vertebral_height = 0.4*distance_to_cranial_vertebra + if_else(spine_point == "l5_centroid", 
                                                                         0.7*distance_to_caudal_vertebra, 0.4*distance_to_caudal_vertebra)) %>%
    # filter(str_detect(spine_point, "centroid")) %>%
    mutate(vertebral_height = if_else(is.na(vertebral_height), lag(vertebral_height), vertebral_height)) %>%
    select(spine_point, x, y, vertebral_height)
  
  
  vert_angles_df <- vertebral_heights_df %>%
    mutate(vert_angle = pmap(.l = list(..1 = lead(x),
                                       ..2 = lead(y),
                                       ..3 = x,
                                       ..4 = y
                                       # ..3 = lag(x),
                                       # ..4 = lag(y)
                                       ), 
                             .f = ~ compute_perpendicular_angle(coord1_x = ..1,
                                                                coord1_y = ..2, 
                                                                coord2_x = ..3, 
                                                                coord2_y = ..4))) %>%
    unnest(vert_angle) %>%
    # mutate(vert_angle = if_else(spine_facing == "left", vert_angle+180, vert_angle)) %>%
    filter(str_detect(spine_point, "centroid")) %>%
    mutate(vert_angle = if_else(spine_point == "c2_centroid", lag(vert_angle), vert_angle)) %>%
    mutate(vert_angle = vert_angle*-1)

  # if(str_to_lower(spine_facing) == "right"){
  #   vert_angles_df <- vert_angles_df %>%
  #     mutate(vert_angle = vert_angle*-1)
  # }
  
  
  buffer_amount <- median(vertebral_heights_df$vertebral_height, na.rm = TRUE)*0.2
  

  # vertebral_coord_dim_sf_df <- vertebral_heights_df %>%
  #   filter(str_detect(spine_point, "centroid")) %>%
  #   mutate(vertebral_width = seq(from = s1_endplate_width, to = 0.6*s1_endplate_width, length = nrow(.))) %>%
  #   mutate(vert_sf = pmap(.l = list(..1 = x, ..2 = y, ..3 = vertebral_width, ..4 = vertebral_height),
  #                         .f = ~ create_vertebra(centroid_x = ..1,
  #                                                centroid_y = ..2,
  #                                                width = ..3,
  #                                                height = ..4))) %>%
  #   left_join(vert_angles_df)%>%
  #   mutate(vert_sf = map(.x = vert_sf, .f = ~ st_geometry(.x))) %>%
  #   mutate(vert_sf_rot = map2(.x = vert_sf, .y = vert_angle,
  #                             .f = ~ rotate_polygon_around_centroid(polygon = .x, angle_degrees = .y))) %>%
  #   mutate(vert_sf_rot = map(.x = vert_sf_rot, .f = ~ st_geometry(st_zm(st_buffer(x = st_buffer(x = .x, dist = -buffer_amount, endCapStyle = "ROUND"),
  #                                                                                 dist = buffer_amount, endCapStyle = "ROUND")))
  #   )
  #   ) %>%
  #   rowwise() %>%
  #   mutate(geometry = st_sfc(vert_sf_rot)) %>%
  #   ungroup()
  
  vert_coord_for_sup_endplates_df <- vertebral_heights_df %>%
    filter(str_detect(spine_point, "centroid")) %>%
    mutate(vertebral_width = seq(from = s1_endplate_width, to = 0.6*s1_endplate_width, length = nrow(.))) %>%
    mutate(vert_sf = pmap(.l = list(..1 = x, ..2 = y, ..3 = vertebral_width, ..4 = vertebral_height),
                          .f = ~ create_vertebra(centroid_x = ..1,
                                                 centroid_y = ..2,
                                                 width = ..3,
                                                 height = ..4))) %>%
    left_join(vert_angles_df)%>%
    mutate(vert_sf = map(.x = vert_sf, .f = ~ st_geometry(.x))) %>%
    mutate(vert_sf_rot = map2(.x = vert_sf, .y = vert_angle,
                              .f = ~ rotate_polygon_around_centroid(polygon = .x, angle_degrees = .y))) %>%
    rowwise() %>%
    mutate(geometry = st_sfc(vert_sf_rot)) %>%
    ungroup()
  
  
  
  
  
  fem_head_sf <- st_buffer(x = st_point(femoral_head_center),
                           dist = buffer_amount * 4,
                           endCapStyle = "ROUND")
  
  
  s1_mid <- c((centroid_df %>% filter(spine_point == "s1_center"))$x, (centroid_df %>% filter(spine_point == "s1_center"))$y)
  
  sac_inf_1 <- find_sacrum_inf_point_function(s1_posterior_sup = s1_posterior_superior, 
                                              s1_midpoint = s1_mid, 
                                              femoral_heads_center = femoral_head_center, 
                                              spine_facing = spine_facing)
  
  sacrum_sf <- st_polygon(list(rbind(s1_anterior_superior, 
                                     sac_inf_1,
                                     s1_posterior_superior, 
                                     s1_anterior_superior)))
  
  sacrum_sf <-  st_buffer(x = st_buffer(x = sacrum_sf, dist = -buffer_amount, endCapStyle = "ROUND"), dist = buffer_amount, endCapStyle = "ROUND")
  
  
  # head_center <- st_point(c(c1_list$sa[[1]] + spine_orientation*1, 
  #                           c1_list$sa[[2]] +4))
  # 
  # skull <- st_linestring(x = rbind(c(head_center[[1]] -spine_orientation*1, head_center[[2]]-.25),
  #                                  c(head_center[[1]] +spine_orientation*1, head_center[[2]])))
  # 
  # skull_sf <-  st_buffer(skull, dist = 5.5, endCapStyle = "ROUND")
  # 
  # 
  # ## triangle:
  # jaw <- st_polygon(list(rbind(c(head_center[[1]], head_center[[2]] -4), ## Most Right point
  #                              c(head_center[[1]] - spine_orientation*5.5, head_center[[2]] -1), ## top point
  #                              c(head_center[[1]] - spine_orientation*5.5, head_center[[2]]-8), ## bottom
  #                              c(head_center[[1]], head_center[[2]] -4)))) ## most right
  # 
  # jaw_sf <-  st_buffer(jaw, dist = 1.2, endCapStyle = "ROUND")
  # 
  # head_sf <- st_union(skull_sf, jaw_sf)
  

  ## sup endplates ##
  
  
  spine_geom_list_endplates <- vert_coord_for_sup_endplates_df$geometry
  
  names(spine_geom_list_endplates) <- str_replace_all(vert_coord_for_sup_endplates_df$spine_point, "centroid", "geom")
  
  sup_endplates_list_df <- map(.x = spine_geom_list_endplates, .f = ~ clean_names(as_tibble(st_coordinates(.x)[1:2,1:2])) %>% mutate(endplate_points = rownames(st_coordinates(.x)[1:2,1:2])))

  names(sup_endplates_list_df) <- str_replace_all(names(spine_geom_list_endplates), pattern = "geom", "superior_endplate")

  sup_endplates_list_df <- map2(.x = sup_endplates_list_df, .y = names(sup_endplates_list_df), .f = ~ .x %>% mutate(spine_level = .y))

  sup_endplates_df <- tibble(spine_level = "s1_superior_endplate",
                             sp_x = s1_posterior_superior[1],
                             sp_y = s1_posterior_superior[2],
                             sa_x = s1_anterior_superior[1],
                             sa_y = s1_anterior_superior[2]) %>%
    union_all(bind_rows(sup_endplates_list_df) %>%
                select(spine_level, endplate_points, x, y) %>%
                filter(endplate_points == "sp") %>%
                select(spine_level, sp_x = x, sp_y = y) %>%
                left_join(bind_rows(sup_endplates_list_df) %>%
                            select(spine_level, endplate_points, x, y) %>%
                            filter(endplate_points == "sa") %>%
                            select(spine_level, sa_x = x, sa_y = y))
  )
  
  revised_corners_df <- sup_endplates_df  %>%
    mutate(caudal_level = lag(spine_level)) %>%
    mutate(inferior_sp_x = lag(sp_x),
           inferior_sp_y = lag(sp_y),
           inferior_sa_x = lag(sa_x),
           inferior_sa_y = lag(sa_y)
    ) %>%
    filter(!is.na(caudal_level)) %>%
    mutate(ia_x = pmap(.l = list(
      ..1 = sa_x,
      ..2 = sa_y,
      ..3 = inferior_sa_x,
      ..4 = inferior_sa_y
    ), .f = ~ compute_inferior_corner(x = ..1, y = ..2, inferior_x = ..3, inferior_y = ..4, return_x_or_y = "x"))) %>%
    unnest(ia_x) %>%
    mutate(ia_y = pmap(.l = list(
      ..1 = sa_x,
      ..2 = sa_y,
      ..3 = inferior_sa_x,
      ..4 = inferior_sa_y
    ), .f = ~ compute_inferior_corner(x = ..1, y = ..2, inferior_x = ..3, inferior_y = ..4, return_x_or_y = "y"))) %>%
    unnest(ia_y)%>%
    mutate(ip_x = pmap(.l = list(
      ..1 = sp_x,
      ..2 = sp_y,
      ..3 = inferior_sp_x,
      ..4 = inferior_sp_y
    ), .f = ~ compute_inferior_corner(x = ..1, y = ..2, inferior_x = ..3, inferior_y = ..4, return_x_or_y = "x"))) %>%
    unnest(ip_x) %>%
    mutate(ip_y = pmap(.l = list(
      ..1 = sp_x,
      ..2 = sp_y,
      ..3 = inferior_sp_x,
      ..4 = inferior_sp_y
    ), .f = ~ compute_inferior_corner(x = ..1, y = ..2, inferior_x = ..3, inferior_y = ..4, return_x_or_y = "y"))) %>%
    unnest(ip_y)
  
  smoothed_vert_df <- revised_corners_df %>%
    select(spine_level, sp_x, sp_y, sa_x, sa_y, ia_x, ia_y, ip_x, ip_y) %>%
    pivot_longer(cols = c(sp_x, sa_x, ia_x, ip_x, sp_y, sa_y, ia_y, ip_y), names_to = "vert_point", values_to = "value") %>%
    mutate(x_or_y = if_else(str_detect(vert_point, "_x"), "x", "y")) %>%
    mutate(vert_point = str_remove_all(vert_point, "_x|_y")) %>%
    pivot_wider(names_from = x_or_y, values_from = value) %>%
    group_by(spine_level) %>%
    nest() %>%
    ungroup() %>%
    mutate(vert_coord_df = map(.x = data, .f = ~ .x %>% union_all(head(.x, n = 1)))) %>%
    select(spine_level, vert_coord_df) %>%
    mutate(geometry = map(.x = vert_coord_df, .f = ~ st_polygon(list(as.matrix(.x %>% select(x, y))))))%>%

    mutate(geometry = map(.x = geometry, .f = ~ st_geometry(st_zm(st_buffer(x = st_buffer(x = .x, dist = -buffer_amount, endCapStyle = "ROUND"),
                                                                            dist = buffer_amount, endCapStyle = "ROUND")))
    )
    ) %>%
    rowwise() %>%
    mutate(geometry = st_sfc(geometry)) %>%
    ungroup()%>%
    mutate(spine_level = str_remove_all(spine_level, "_superior_endplate"))
  
  vertebral_coord_dim_sf_df <- vert_coord_for_sup_endplates_df  %>%
    mutate(spine_level = str_remove_all(spine_point, "_centroid")) %>%
    select(spine_level, spine_point, x, y, vertebral_height, vertebral_width, vert_angle, original_vert = geometry) %>%
    left_join(smoothed_vert_df)
  # 
  spine_geom_list <- vertebral_coord_dim_sf_df$geometry

  names(spine_geom_list) <- str_replace_all(vertebral_coord_dim_sf_df$spine_point, "centroid", "geom")
  
  
  ########### lines ################

  spine_geom_list_for_lines <- vertebral_coord_dim_sf_df$geometry

  names(spine_geom_list_for_lines) <- str_replace_all(vertebral_coord_dim_sf_df$spine_point, "centroid", "geom")

  lines_list <- list()
  l1pa_line <- st_linestring(rbind(st_centroid(spine_geom_list_for_lines$l1_geom),
                                   femoral_head_center,
                                   s1_mid))

  lines_list$l1pa <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
                                                                  line_st_geometry = l1pa_line,
                                                                  distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list_for_lines$l2_geom),
                                                                                                                    femoral_head_center)))/3)
  t9pa_line <- st_linestring(rbind(st_centroid(spine_geom_list_for_lines$t9_geom),
                                   femoral_head_center,
                                   s1_mid))
  
  lines_list$t9pa <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
                                                                  line_st_geometry = t9pa_line,
                                                                  distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list_for_lines$l1_geom),
                                                                                                                    femoral_head_center)))/3)
  
  
  t4pa_line <- st_linestring(rbind(st_centroid(spine_geom_list_for_lines$t4_geom),
                                   femoral_head_center,
                                   s1_mid))

  lines_list$t4pa <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
                                                                  line_st_geometry = t4pa_line,
                                                                  distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list_for_lines$t12_geom),
                                                                                                                    femoral_head_center)))/3)

  
  c2pa_line <- st_linestring(rbind(st_centroid(spine_geom_list_for_lines$c2_geom),
                                   femoral_head_center,
                                   s1_mid))
  
  lines_list$c2pa <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
                                                                  line_st_geometry = c2pa_line,
                                                                  distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list_for_lines$t11_geom),
                                                                                                                    femoral_head_center)))/3)
  ######### computing segment_angles
  
  segment_angle_modifier <- if_else(str_to_lower(spine_facing) == "right", 1, -1)

  
  
  
  xr_segment_angle_df <- vertebral_coord_dim_sf_df %>%
    mutate(segment = str_replace_all(spine_point, "centroid", "segment_angle")) %>%
    select(segment, vert_angle) %>%
    mutate(vert_angle = vert_angle*segment_angle_modifier) %>%
    mutate(segment_angle = vert_angle - lag(vert_angle)) %>%
    mutate(segment_angle = if_else(is.na(segment_angle), vert_angle, segment_angle)) %>%
    select(segment, segment_angle)
  
  segment_angles_list <- as.list(xr_segment_angle_df$segment_angle)
  
  names(segment_angles_list) <- xr_segment_angle_df$segment
  
  segment_angles_list$c1_segment_angle <- segment_angles_list$c2_segment_angle
  

  return(list(fem_head_sf = fem_head_sf,
              sacrum_sf = sacrum_sf,
              vert_geom_df = vertebral_coord_dim_sf_df,
              spine_geom_list = spine_geom_list,
              lines_list = lines_list,
              sup_endplates_df = sup_endplates_df, 
              segment_angles_list = segment_angles_list
              ))
  
  # return(list(
  #   sup_endplates_df = sup_endplates_df, 
  #   revised_corners_df = revised_corners_df,
  #   smoothed_vert_df = smoothed_vert_df,
  #   vert_coord_for_sup_endplates_df = vert_coord_for_sup_endplates_df
  # ))
  
}





# 
# 
# 
# build_spine_xray_outline_from_coordinates_function <- function(femoral_head_center = c(0,0),
#                                                   s1_anterior_superior = c(1, 1),
#                                                   s1_posterior_superior = c(1, 2),
#                                                   centroid_df = tibble(), 
#                                                   spine_facing = "right"){
#   
#   spine_facing_modifier <- if_else(str_to_lower(spine_facing) == "right", -1, 1)
#   
#   s1_endplate_width <- jh_calculate_distance_between_2_points_function(point_1 = s1_anterior_superior, point_2 = s1_posterior_superior)
#   
#   
#   
#   vertebral_heights_df <- centroid_df %>%
#     filter(str_detect(spine_point, "centroid|center")) %>%
#     mutate(distance_to_cranial_vertebra = pmap(.l = list(..1 = x, ..2 = y, ..3 = lead(x), ..4 = lead(y)), 
#                                                .f = ~ jh_calculate_distance_between_2_points_function(point_1 = c(..1, ..2),
#                                                                                                       point_2 = c(..3, ..4)))) %>%
#     unnest() %>%
#     mutate(distance_to_caudal_vertebra = pmap(.l = list(..1 = x, ..2 = y, ..3 = lag(x), ..4 = lag(y)), 
#                                               .f = ~ jh_calculate_distance_between_2_points_function(point_1 = c(..1, ..2),
#                                                                                                      point_2 = c(..3, ..4)))) %>%
#     unnest() %>%
#     mutate(vertebral_height = 0.4*distance_to_cranial_vertebra + if_else(spine_point == "l5_centroid", 
#                                                                          0.7*distance_to_caudal_vertebra, 0.4*distance_to_caudal_vertebra)) %>%
#     # filter(str_detect(spine_point, "centroid")) %>%
#     mutate(vertebral_height = if_else(is.na(vertebral_height), lag(vertebral_height), vertebral_height)) %>%
#     select(spine_point, x, y, vertebral_height)
#   
#   
#   vert_angles_df <- vertebral_heights_df %>%
#     mutate(vert_angle = pmap(.l = list(..1 = lead(x),
#                                        ..2 = lead(y),
#                                        ..3 = x,
#                                        ..4 = y
#                                        # ..3 = lag(x),
#                                        # ..4 = lag(y)
#                                        ), 
#                              .f = ~ compute_perpendicular_angle(coord1_x = ..1,
#                                                                 coord1_y = ..2, 
#                                                                 coord2_x = ..3, 
#                                                                 coord2_y = ..4))) %>%
#     unnest(vert_angle) %>%
#     filter(str_detect(spine_point, "centroid"))
#   
#   
#   if(str_to_lower(spine_facing) == "left"){
#     vert_angles_df <- vert_angles_df %>%
#       mutate(vert_angle = vert_angle*-1)
#   }
#   # 
#   
#   buffer_amount <- median(vertebral_heights_df$vertebral_height, na.rm = TRUE)*0.2
#   
#   # buffer_amount <- mean(vertebral_heights_df$vertebral_height, )*0.2
#   # buffer_amount <- 10
#   
#   vertebral_coord_dim_sf_df <- vertebral_heights_df %>%
#     filter(str_detect(spine_point, "centroid")) %>%
#     mutate(vertebral_width = seq(from = s1_endplate_width, to = 0.6*s1_endplate_width, length = nrow(.))) %>%
#     mutate(vert_sf = pmap(.l = list(..1 = x, ..2 = y, ..3 = vertebral_width, ..4 = vertebral_height),
#                           .f = ~ create_vertebra(centroid_x = ..1,
#                                                  centroid_y = ..2,
#                                                  width = ..3,
#                                                  height = ..4))) %>%
#     left_join(vert_angles_df)%>%
#     mutate(vert_sf = map(.x = vert_sf, .f = ~ st_geometry(.x))) %>%
#     mutate(vert_sf_rot = map2(.x = vert_sf, .y = vert_angle,
#                               .f = ~ rotate_polygon_around_centroid(polygon = .x, angle_degrees = .y))) %>%
#     rowwise() %>%
#     mutate(geometry = st_sfc(vert_sf_rot)) %>%
#     ungroup()
#   
# 
#   
#   fem_head_sf <- st_buffer(x = st_point(femoral_head_center),
#                            dist = buffer_amount * 4,
#                            endCapStyle = "ROUND")
#   
#   
#   s1_mid <- c((centroid_df %>% filter(spine_point == "s1_center"))$x, (centroid_df %>% filter(spine_point == "s1_center"))$y)
#   
#   sac_inf_1 <- find_sacrum_inf_point_function(s1_posterior_sup = s1_posterior_superior, 
#                                               s1_midpoint = s1_mid, 
#                                               femoral_heads_center = femoral_head_center, 
#                                               spine_facing = spine_facing)
#   
#   sacrum_sf <- st_polygon(list(rbind(s1_anterior_superior, 
#                                      sac_inf_1,
#                                      s1_posterior_superior, 
#                                      s1_anterior_superior)))
#   
#   sacrum_sf <-  st_buffer(x = st_buffer(x = sacrum_sf, dist = -buffer_amount, endCapStyle = "ROUND"), dist = buffer_amount, endCapStyle = "ROUND")
#   
#   
#   # head_center <- st_point(c(c1_list$sa[[1]] + spine_orientation*1, 
#   #                           c1_list$sa[[2]] +4))
#   # 
#   # skull <- st_linestring(x = rbind(c(head_center[[1]] -spine_orientation*1, head_center[[2]]-.25),
#   #                                  c(head_center[[1]] +spine_orientation*1, head_center[[2]])))
#   # 
#   # skull_sf <-  st_buffer(skull, dist = 5.5, endCapStyle = "ROUND")
#   # 
#   # 
#   # ## triangle:
#   # jaw <- st_polygon(list(rbind(c(head_center[[1]], head_center[[2]] -4), ## Most Right point
#   #                              c(head_center[[1]] - spine_orientation*5.5, head_center[[2]] -1), ## top point
#   #                              c(head_center[[1]] - spine_orientation*5.5, head_center[[2]]-8), ## bottom
#   #                              c(head_center[[1]], head_center[[2]] -4)))) ## most right
#   # 
#   # jaw_sf <-  st_buffer(jaw, dist = 1.2, endCapStyle = "ROUND")
#   # 
#   # head_sf <- st_union(skull_sf, jaw_sf)
#   
#   spine_geom_list <- vertebral_coord_dim_sf_df$geometry
#   
#   names(spine_geom_list) <- str_replace_all(vertebral_coord_dim_sf_df$spine_point, "centroid", "geom")
#   
#   lines_list <- list()
#   l1pa_line <- st_linestring(rbind(st_centroid(spine_geom_list$l1_geom),
#                                    femoral_head_center,
#                                    s1_mid))
#   
#   lines_list$l1pa_line_curve_sf <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
#                                                                   line_st_geometry = l1pa_line,
#                                                                   distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list$l2_geom),
#                                                                                                                     femoral_head_center)))/3)
#   t4pa_line <- st_linestring(rbind(st_centroid(spine_geom_list$t4_geom),
#                                    femoral_head_center,
#                                    s1_mid))
#   
#   lines_list$t4pa_line_curve_sf <-   jh_plot_angle_curve_function(vertex_vector = femoral_head_center,
#                                                                   line_st_geometry = t4pa_line,
#                                                                   distance_of_curve = st_length(st_linestring(rbind(st_centroid(spine_geom_list$t12_geom),
#                                                                                                                     femoral_head_center)))/3)
#   
#   
#   
#   
#   ## sup endplates ##
#   sup_endplates_list_df <- map(.x = spine_geom_list, .f = ~ clean_names(as_tibble(st_coordinates(.x)[1:2,1:2])) %>% mutate(endplate_points = rownames(st_coordinates(.x)[1:2,1:2])))
#   
#   names(sup_endplates_list_df) <- str_replace_all(names(spine_geom_list), pattern = "geom", "superior_endplate")
#   
#   sup_endplates_list_df <- map2(.x = sup_endplates_list_df, .y = names(sup_endplates_list_df), .f = ~ .x %>% mutate(spine_level = .y))
#   
#   
#   sup_endplates_df <- tibble(spine_level = "s1_superior_endplate",
#                              sp_x = s1_posterior_superior[1],
#                              sp_y = s1_posterior_superior[2],
#                              sa_x = s1_anterior_superior[1],
#                              sa_y = s1_anterior_superior[2]) %>%
#     union_all(bind_rows(sup_endplates_list_df) %>%
#                 select(spine_level, endplate_points, x, y) %>%
#                 filter(endplate_points == "sp") %>%
#                 select(spine_level, sp_x = x, sp_y = y) %>%
#                 left_join(bind_rows(sup_endplates_list_df) %>%
#                             select(spine_level, endplate_points, x, y) %>%
#                             filter(endplate_points == "sa") %>%
#                             select(spine_level, sa_x = x, sa_y = y))
#     )
#   
#   
#   
#   return(list(fem_head_sf = fem_head_sf, 
#               sacrum_sf = sacrum_sf,
#               vert_geom_df = vertebral_coord_dim_sf_df, 
#               spine_geom_list = spine_geom_list, 
#               lines_list = lines_list, 
#               sup_endplates_df = sup_endplates_df))
#   
# }

