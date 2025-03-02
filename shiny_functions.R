generate_spine_level_controls <- function(spine_level, 
                                          spine_level_font_size = 11, 
                                          return_as_full_table = TRUE) {
  
  label_style <- glue("font-size:{paste(spine_level_font_size)}px; color:darkblue; font-weight:bold; text-align:center; margin-top:0; margin-bottom:0")
  
  spine_level_id <- str_to_lower(str_replace_all(spine_level, "-", "_"))
  
  label_percent_width <- 30
  button_percent_width <- (100-label_percent_width)/4
  
  row <- tags$tr(width = "100%",
                 tags$td(width = paste0(button_percent_width, "%"),
                         actionBttn(
                           inputId = paste0(spine_level_id, "_lordosis_down_5"),
                           label = "-5",
                           style = "material-circle",
                           size = "xs"
                         )
                 ),
                 tags$td(width = paste0(button_percent_width, "%"),
                         actionBttn(
                           inputId = paste0(spine_level_id, "_lordosis_down_1"),
                           label = "-1",
                           style = "material-circle",
                           size = "xs"
                         )
                 ),
                 tags$td(width = paste0(label_percent_width, "%"), 
                         tags$div(style = label_style, paste(spine_level))),
                 tags$td(width = paste0(button_percent_width, "%"),
                         actionBttn(
                           inputId = paste0(spine_level_id, "_lordosis_up_1"),
                           label = "+1",
                           style = "material-circle",
                           size = "xs"
                         )
                 ),
                 tags$td(width = paste0(button_percent_width, "%"),
                         actionBttn(
                           inputId = paste0(spine_level_id, "_lordosis_up_5"),
                           label = "+5",
                           style = "material-circle",
                           size = "xs"
                         )
                 ),
  )
  
  if(return_as_full_table == TRUE){
    return(tags$table(width = "100%",
                      row))
  }else{
    return(row)
  }
  
}

# Function to update the dataframe

# update_spine_segmental_planning_df_function <- function(spine_segmental_planning_df = tibble(),
#                                                           spine_interspace_input, change) {
#   spine_segmental_planning_df %>%
#     mutate(adjustment = if_else(spine_interspace == spine_interspace_input, adjustment + change, adjustment))
# }

update_spine_segmental_planning_df_function <- function(spine_segmental_planning_df, spine_interspace_input, change) {
  spine_segmental_planning_df$df <- spine_segmental_planning_df$df %>%
    mutate(adjustment = if_else(spine_interspace == spine_interspace_input, adjustment + change, adjustment))
}


update_spine_segmental_planning_table_observe_button_function <- function(spine_segmental_planning_df, 
                                                                          spine_interspace, session) {
  spine_interspace_id <- gsub("-", "_", tolower(spine_interspace))  # Ensure ID consistency
  
  observeEvent(session$input[[paste0(spine_interspace_id, "_lordosis_down_5")]], {
    spine_segmental_planning_df$df <- update_spine_segmental_planning_df_function(spine_segmental_planning_df = spine_segmental_planning_df,
                                                                                  spine_interspace, -5)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(session$input[[paste0(spine_interspace_id, "_lordosis_down_1")]], {
    spine_segmental_planning_df$df <- update_spine_segmental_planning_df_function(spine_segmental_planning_df = spine_segmental_planning_df,
                                                spine_interspace, -1)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(session$input[[paste0(spine_interspace_id, "_lordosis_up_1")]], {
    spine_segmental_planning_df$df <- update_spine_segmental_planning_df_function(spine_segmental_planning_df = spine_segmental_planning_df,
                                                spine_interspace, 1)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(session$input[[paste0(spine_interspace_id, "_lordosis_up_5")]], {
    spine_segmental_planning_df$df <-  update_spine_segmental_planning_df_function(spine_segmental_planning_df = spine_segmental_planning_df,
                                                spine_interspace, 5)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
}



construct_rod_coordinates_function <- function(final_spine_coordinates_df, uiv = "t4"){
  
  l5_inf_df <- final_spine_coordinates_df %>%
    filter(spine_level == str_to_lower("L5"))%>% 
    filter(vert_point == "sp")
  
  sacrum_sup_df <- final_spine_coordinates_df %>%
    filter(spine_level == "sacrum")%>% 
    filter(vert_point == "s1_posterior_superior")
  
  # inferior_rod_point <- jh_get_point_along_line_function(coord_a = c(l5_inf_df$x, l5_inf_df$y), 
  #                                                        coord_b =  c(sacrum_sup_df$x, sacrum_sup_df$y), 
  #                                                        percent_a_to_b = 1.5)
  
  inferior_rod_point <- c(sacrum_sup_df$x, sacrum_sup_df$y)
  
  uiv_sup_df <- final_spine_coordinates_df %>%
    filter(spine_level == str_to_lower(uiv))%>% 
    filter(vert_point == "sp")
  
  uiv_inf_df <- final_spine_coordinates_df %>%
    filter(spine_level == str_to_lower(uiv))%>% 
    filter(vert_point == "ip")
  
  sup_rod_point <- jh_get_point_along_line_function(coord_a = c(uiv_inf_df$x, uiv_inf_df$y), coord_b = c(uiv_sup_df$x, uiv_sup_df$y), percent_a_to_b = 1.65)
  
  rod_coordinate_points_df <- final_spine_coordinates_df %>%
    filter(y < sup_rod_point[[2]])%>%
    filter(y > inferior_rod_point[[2]]) %>%
    filter(vert_point %in% c("s1_posterior_superior", "sp", "ip"))%>%
    select(spine_level, x, y) %>%
    add_row(spine_level = "superior_rod", x = sup_rod_point[[1]], y = sup_rod_point[[2]]) %>%
    add_row(spine_level = "inferior_rod", x = inferior_rod_point[[1]], y = inferior_rod_point[[2]])%>%
    mutate(x = x + 10)%>%
    arrange(y)
  
  return(rod_coordinate_points_df)
  
}
