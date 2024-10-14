#library(colourpicker)
library(shiny)
# library(reactlog)
library(sf)
library(tidyverse)
library(ggplot2)
library(shinyWidgets)
library(shinyBS)
# library(kableExtra)
# library(rms)
library(svglite)
library(glue)
library(cowplot)
library(janitor)
# library(rms)
# library(Hmisc)
library(cowplot)
# library(assertr)
library(lubridate)
library(shinydashboard)
library(magick)
# library(ggforce)
# library(plotly)

options(shiny.maxRequestSize = 25*1024^2)

source("jh_functions.R", local = TRUE)
source("compute_segment_angles_function_from_sim_data_2024.R", local = TRUE)
source("function_segment_angles_separated.R", local = TRUE)

source("jh_spine_build_NEW_function.R", local = TRUE)

source("jh_prescribing_alignment_functions.R", local = TRUE)

source("spinal_regional_alignment_analysis_by_vpa.R", local = TRUE)

# source("jh_build_spine_by_vertebral_pelvic_angles_only.R", local = TRUE)
source("jh_build_spine_by_vertebral_pelvic_angles_cleaned.R", local = TRUE)
source("prescribing_alignment_by_matching_unfused.R", local = TRUE)
source("xray_segment_angles_model_functions.R", local = TRUE)
source("build_spine_from_coordinates_functions.R", local = TRUE)


jh_calculate_distance_between_2_points_function <- function(point_1, point_2){
  sqrt((point_1[1] - point_2[1])^2 + 
         (point_1[2] - point_2[2])^2)
}

jh_compute_vpa_from_xray_data_function <- function(fem_head_center = c(0,0), 
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

# Define the labels for both modes
get_spine_labels <- function(all_centroids = FALSE) {
  if (all_centroids) {
    return(c("fem_head_center", 
             "s1_anterior_superior", "s1_posterior_superior", 
             "l5_centroid", "l4_centroid", "l3_centroid", "l2_centroid", "l1_centroid", 
             "t12_centroid", "t11_centroid", "t10_centroid", "t9_centroid", 
             "t8_centroid", "t7_centroid", "t6_centroid", "t5_centroid", 
             "t4_centroid", "t3_centroid", "t2_centroid", "t1_centroid", 
             "c7_centroid", "c6_centroid", "c5_centroid", "c4_centroid", "c3_centroid", 
             "c2_centroid"))
  } else {
    return(c("fem_head_center", 
             "s1_anterior_superior", "s1_posterior_superior", 
             "l4_centroid", "l1_centroid", "t9_centroid", "t4_centroid", 
             "t1_centroid", "c2_centroid"))
  }
}

compute_perpendicular_points <- function(x1, y1, x2, y2, distance = 0.01) {
  # Midpoint
  midpoint_x <- (x1 + x2) / 2
  midpoint_y <- (y1 + y2) / 2
  
  # Slope of the line (tangent)
  slope <- (y2 - y1) / (x2 - x1)
  
  # Perpendicular slope (-1 / slope)
  perpendicular_slope <- -1 / slope
  
  # Calculate the change in x and y for perpendicular points
  delta_x <- distance / sqrt(1 + perpendicular_slope^2)
  delta_y <- perpendicular_slope * delta_x
  
  # Two perpendicular points
  point_1 <- c(midpoint_x + delta_x, midpoint_y + delta_y)
  point_2 <- c(midpoint_x - delta_x, midpoint_y - delta_y)
  
  return(tibble(x1 = point_1[1], y1 = point_1[2], x2 = point_2[1], y2 = point_2[2]))
}

calculate_pelvic_incidence_line_coordinates <- function(fem_head_center = c(0,0), 
                                                        s1_anterior, 
                                                        s1_posterior, 
                                                        spine_facing = "left",
                                                        pelvic_tilt = 10, 
                                                        pelvic_incidence_value = 50) {
  
  # Step 1: Calculate the center (midpoint)
  center_x <- (s1_anterior[1] + s1_posterior[1]) / 2
  center_y <- (s1_anterior[2] + s1_posterior[2]) / 2
  center <- c(center_x, center_y)
  
  # Step 2: Calculate the length of the line between s1_anterior and s1_posterior
  line_length <- sqrt((s1_anterior[1] - s1_posterior[1])^2 + (s1_anterior[2] - s1_posterior[2])^2)
  
  # Step 4: Calculate the length of the perpendicular line (5 times the original length)
  extended_length <- 3 * line_length
  
  if (s1_anterior[1] == s1_posterior[1]) {
    # For vertical lines, the perpendicular is horizontal
    dx <- extended_length
    dy <- 0
  } else if (s1_anterior[2] == s1_posterior[2]) {
    # For horizontal lines, the perpendicular is vertical
    dx <- 0
    dy <- extended_length
  } else {
    pt_pi_diff <- pelvic_incidence_value - pelvic_tilt
    
    pt_pi_diff_rad <- pt_pi_diff*(pi/180)
    
    orientation_modifier <- if_else(spine_facing == "left", -1, 1)
    
    dx <- sin(pt_pi_diff_rad)*extended_length*orientation_modifier
    dy <- cos(pt_pi_diff_rad)*extended_length
  }
  
  # Inferior point is displaced from the center by (dx, dy)
  inferior_x <- center[1] - dx
  inferior_y <- center[2] - dy
  inferior <- c(inferior_x, inferior_y)
  
  pi_line_coordinates_df <- tibble(spine_point = c("fem_head_center", "s1_center", "s1_inferior"), 
                                   x = c(fem_head_center[1], 
                                         center[1],
                                         inferior_x),
                                   y = c(fem_head_center[2],
                                         center[2],
                                         inferior_y)
  )
  # Return the center and inferior points as a list
  # return(list(center = center, inferior = inferior))
  pi_line_coordinates_df
}

# all_possible_lumbar_segments_angles_with_lpa_df <- read_csv("all_possible_lumbar_segment_angles_for_lpa.csv")

# reactlog_enable()

spinal_segments_labels_vector <- c('L5-S1', 'L4-L5', 'L3-L4', 'L2-L3', 'L1-L2', 
                                   'T12-L1', 'T11-T12', 'T10-T11', 'T9-T10', 'T8-T9', 'T7-T8', 'T6-T7', 'T5-T6', 'T4-T5', 'T3-T4', 'T2-T3', 'T1-T2',
                                   'C7-T1', 'C6-C7', 'C5-C6', 'C4-C5', 'C3-C4', 'C2-C3', 'C1-C2')

create_spine_rigid_level_input_function <- function(segment_input_label){
  segment_id <- paste0("preop_", str_to_lower(str_replace_all(segment_input_label, pattern = "-", "_")), "_segment")
  rigid_segment_id <- str_replace_all(segment_id, "_segment", "_rigid_xray")
  segment_label <- segment_input_label
  # initial_value <- segment_value_input
  
  # div(
    # class = "segment-input",
    # span(segment_label,
    #      class = "segment-label"),
    div(
      class = "segment-input",
      prettyCheckbox(
        inputId = rigid_segment_id,
        label = segment_label,
        value = FALSE,
        bigger = TRUE,
        status = "danger",
        shape = "curve"
      )
    )
  # )
}



ui <- dashboardPage(
  dashboardHeader(title = "SolaSpine"),
  dashboardSidebar(
    tags$style(HTML("
  .segment-input {
    display: flex;
    align-items: center; /* Align items to the center vertically */
    justify-content: end; /* Ensure space between label and input */
    margin-bottom: 0px; /* Adjusts the spacing between the inputs */
    font-size: 14px;
    color: black;
  }
  .segment-label {
    margin-right: 5px; /* Slightly increase space for the label */
    white-space: nowrap; /* Prevent labels from wrapping */
    font-size: 12px;
    color: black;
  }
  .segment-input .form-group {
    margin-bottom: 1px; /* Reduce default margin-bottom of form-group */
    color: black;
  }
  .custom-numeric-input {
    padding: 0; /* Remove padding from the numeric input container */
    margin: 0; /* Remove margin from the numeric input container */
    text-align: -webkit-left;
  }
  .custom-numeric-input .form-control {
    padding: 2px 5px; /* Adjust padding inside the numeric input */
    margin-bottom: 0px; /* Ensure no extra margin below the input */
    text-align: -webkit-left;
    width: 50px;
  }
")),
    # fileInput("image", "Upload an Image", accept = c('image/png', 'image/jpeg', 'image/jpg')), 
    fileInput("image", "Upload an Image", accept = 'image/'), 
    br(), 
    conditionalPanel(
      condition = "input.xray_file_uploaded == true",
      h4("Xray Orientation:"),
    actionBttn(
      inputId = "spine_orientation_button",
      label = "Facing LEFT", 
      style = "material-flat",
      color = "primary",
      icon = icon("arrow-left")
    )
    ),
    br(),
    conditionalPanel(
      condition = "input.xray_file_uploaded == true",
      h4(strong("Patient Factors:")),
      sliderInput(
        "preop_age",
        "Patient Age:",
        min = 18,
        max = 90,
        value = 60
      ),
      radioGroupButtons(
        inputId = "preop_sex",
        label = NULL,
        individual = TRUE,
        choices = c("Male", "Female"),
        selected = "Female",
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square", 
                       style = "color: steelblue"),
          no = tags$i(class = "fa fa-square-o", 
                      style = "color: steelblue"))
      ),
      hr(),
      
    uiOutput(outputId = "preop_xray_rigid_segments_ui")
    ),
    div(
      style = "display: none;",  # Hide the entire div, including the switch
      switchInput(
        inputId = "xray_file_uploaded",
        size = "mini", label = NULL,
        value = FALSE, 
        onLabel = "Y", 
        offLabel = "N",
      ),
      switchInput(
        inputId = "all_points_recorded",
        size = "mini", label = NULL,
        value = FALSE, 
        onLabel = "Y", 
        offLabel = "N",
      )
    )
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      box(width = 3,
          conditionalPanel(
            condition = "input.xray_file_uploaded == true",
            fluidRow(
            class = "d-flex justify-content-center",  # Center content horizontally
                   tags$div(
                     style = "font-size: 20px; 
                 font-weight: bold; 
                 color: yellow; 
                 font-family: arial; 
                 font-style: italic; 
                 text-align: center; 
                 background-color: black; 
                 padding: 3px; 
                 border-radius: 12px;  /* Rounded corners */
                 display: block;
                 margin-left: 10px;
                 margin-right: 10px;
                     box-sizing: border-box;  /* Include padding and border in the element's width */",
                     htmlOutput(outputId = "xray_click_instructions")
                   ),
                   br()
          )
            ),
          conditionalPanel(
            condition = "input.xray_file_uploaded == true & input.all_points_recorded == false",
            
          fluidRow(
            column(width = 12, 
                   tags$div(
                 id = "image-container",
                 style = "position: relative; width: 350px; height: 700px; overflow: hidden; border: 0px solid #ccc;",
                 tags$img(
                   id = "uploadedImage",
                   src = "",
                   style = "position: absolute; top: 0; left: 0; cursor: crosshair;"
                 )
               ),
               # tags$script(src = "https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"),
               tags$script(HTML("
       $(document).ready(function() {
  let scale = 1;
  let panX = 0, panY = 0;
  let isPanning = false;
  let startX, startY;

  function updateImageTransform() {
    $('#uploadedImage').css({
      'transform-origin': 'top left',
      'transform': `translate(${panX}px, ${panY}px) scale(${scale})`
    });

    // Update positions of all dots to match the transformation
    $('.dot').each(function() {
      const originalX = $(this).data('orig-x');
      const originalY = $(this).data('orig-y');
      const adjustedX = (originalX * scale) + panX;
      const adjustedY = ((imageHeight - originalY) * scale) + panY;

      $(this).css({
        left: adjustedX + 'px',
        top: adjustedY + 'px'
      });
    });
  }

  let imageHeight = null; // We'll determine the height once the image is loaded

  Shiny.addCustomMessageHandler('load-image', function(data) {
    var img = document.getElementById('uploadedImage');
    img.src = data.src;

    // Once the image loads, set the natural height
    img.onload = function() {
      imageHeight = img.naturalHeight;
      // Reset scaling and position when new image is loaded
      scale = 1;
      panX = 0;
      panY = 0;
      updateImageTransform();

      // Remove existing dots when a new image is loaded
      $('.dot').remove();
    };
  });

  Shiny.addCustomMessageHandler('plot-coordinates', function(data) {
    // Remove existing dots
    $('.dot').remove();

    if (!imageHeight) {
      console.error('Image height not set yet.');
      return;
    }

    // Plot all coordinates
    data.coords.forEach(function(coord, index) {
      console.log(`Plotting point ${index + 1}:`, coord);  // Debugging log for each coordinate

      // Create a new dot element
      const dot = $('<div class=\"dot\"></div>');

      // Adjust Y-coordinate for the Cartesian system
      const adjustedX = (coord.x * scale) + panX;
      const adjustedY = ((imageHeight - coord.y) * scale) + panY;  // Adjusting Y-coordinate

      // Debugging log for adjusted positions
      console.log('Adjusted position for dot:', { adjustedX, adjustedY });
      
      // Adjust the dot's CSS to center it on the point clicked
    const dotSize = 10; // This is the width and height of the dot (in pixels)
    const correctionOffset = 0.5; // Small adjustment if there’s still an offset issue

    dot.css({
      position: 'absolute',
      top: (adjustedY - (dotSize / 2)) + correctionOffset + 'px',
      left: (adjustedX - (dotSize / 2)) + correctionOffset + 'px',
      width: dotSize + 'px',
      height: dotSize + 'px',
      'background-color': 'red',
      'border-radius': '50%',
      'pointer-events': 'none', // Ensures dots don't interfere with panning/zooming
      'z-index': 10 // Ensures the dots are layered above the image
    });

      // Store original coordinates for reference during zoom and pan
      dot.data('orig-x', coord.x);
      dot.data('orig-y', coord.y);

      // Append dot to the image container
      $('#image-container').append(dot);
    });

    // Immediately update dot positions to reflect the current zoom and pan state
    updateImageTransform();
  });

  // Handle zoom with the mouse wheel
  $('#image-container').on('wheel', function(e) {
    e.preventDefault();
    const zoomIntensity = 0.1;
    const delta = e.originalEvent.deltaY > 0 ? -1 : 1;
    const previousScale = scale;

    // Update scale
    scale *= (1 + delta * zoomIntensity);
    scale = Math.min(Math.max(0.5, scale), 5);

    // Calculate new pan to keep the zoom centered at mouse position
    const mouseX = e.pageX - $(this).offset().left;
    const mouseY = e.pageY - $(this).offset().top;

    panX = mouseX - (mouseX - panX) * (scale / previousScale);
    panY = mouseY - (mouseY - panY) * (scale / previousScale);

    updateImageTransform();
  });

  // Handle panning with right-click only
  $('#image-container').on('mousedown', function(e) {
    if (e.which === 3) { // Right-click
      isPanning = true;
      startX = e.pageX - panX;
      startY = e.pageY - panY;
      $(this).css('cursor', 'grabbing');
      return false; // Prevent context menu
    }
  });

  $(document).on('mouseup', function() {
    isPanning = false;
    $('#image-container').css('cursor', 'crosshair');
  });

  $(document).on('mousemove', function(e) {
    if (!isPanning) return;
    panX = e.pageX - startX;
    panY = e.pageY - startY;

    updateImageTransform();
  });

  // Prevent the default context menu from appearing on right-click
  $('#image-container').on('contextmenu', function(e) {
    return false;
  });

  // Record click coordinates on left-click
  $('#image-container').on('click', function(e) {
    if (e.which === 1) { // Left-click
      var img = document.getElementById('uploadedImage');
      const rect = img.getBoundingClientRect(); // Get the image's bounding box relative to the viewport

      // Get the click coordinates relative to the image
      const clickX = e.clientX - rect.left;
      const clickY = e.clientY - rect.top;

      // Adjust the coordinates for the current pan and zoom level to get the original image reference frame
      const adjustedX = clickX / scale;
      const adjustedY = clickY / scale;

      // Correcting Y-coordinate (flipping the y-axis based on the height of the image)
      const correctedY = imageHeight - adjustedY;

      // Debugging log for adjusted click positions
      console.log('Adjusted click position:', { adjustedX, correctedY });

      // Send the corrected click coordinates to the Shiny server
      Shiny.setInputValue('xray_click', {x: adjustedX - 1, y: correctedY + 1}, {priority: 'event'});
    }
  });
});
      "))
        )
          )
          ),  
        conditionalPanel(
          condition = "input.all_points_recorded == true",
          tags$head(
            # Include custom CSS to style the container for zoom/pan
            tags$style(HTML("
      #plot-container {
        height: 700px;
        width: 100%;
        overflow: hidden;
        position: relative;
      }
      #xray_plot_image {
        cursor: grab;
        position: relative;
      }
    "))
          ),
          tags$div(
            id = "plot-container",
            tags$img(id = "xray_plot_image", src = "", style = "width: 100%;")
          ),
          tags$script(HTML("
    $(document).ready(function() {
  let plotScale = 1;
  let plotPanX = 0, plotPanY = 0;
  let isPlotPanning = false;
  let plotStartX, plotStartY;

  function updatePlotTransform() {
    $('#xray_plot_image').css({
      'transform-origin': 'top left',
      'transform': `translate(${plotPanX}px, ${plotPanY}px) scale(${plotScale})`
    });
  }
  
    Shiny.addCustomMessageHandler('load-plot-image', function(data) {
    var img = document.getElementById('xray_plot_image');
    img.src = data.src;

    // Once the image loads, set the natural height
    img.onload = function() {
      imageHeight = img.naturalHeight;
      // Reset scaling and position when new image is loaded
      scale = 1;
      panX = 0;
      panY = 0;
      updateImageTransform();

      // Remove existing dots when a new image is loaded
      $('.dot').remove();
    };
  });

  Shiny.addCustomMessageHandler('load-plot-image', function(data) {
    var plotImg = document.getElementById('xray_plot_image');
    plotImg.src = data.src;

    plotImg.onload = function() {
      // Reset scaling and position when new image is loaded
      imageHeight = img.naturalHeight;
      plotScale = 1;
      plotPanX = 0;
      plotPanY = 0;
      updatePlotTransform();
    };
  });

  // Handle zoom with the mouse wheel
  $('#plot-container').on('wheel', function(e) {
    e.preventDefault();
    const zoomIntensity = 0.1;
    const delta = e.originalEvent.deltaY > 0 ? -1 : 1;
    const previousScale = plotScale;

    // Update plotScale
    plotScale *= (1 + delta * zoomIntensity);
    plotScale = Math.min(Math.max(0.5, plotScale), 5);

    // Calculate new pan to keep the zoom centered at mouse position
    const mouseX = e.pageX - $(this).offset().left;
    const mouseY = e.pageY - $(this).offset().top;

    plotPanX = mouseX - (mouseX - plotPanX) * (plotScale / previousScale);
    plotPanY = mouseY - (mouseY - plotPanY) * (plotScale / previousScale);

    updatePlotTransform();
  });
  
  // Handle panning with right-click only
  $('#plot-container').on('mousedown', function(e) {
    if (e.which === 3) { // Right-click
      isPlotPanning = true;
      plotStartX = e.pageX - plotPanX;
      plotStartY = e.pageY - plotPanY;
      $(this).css('cursor', 'grabbing');
      return false; // Prevent context menu
    }
  });

  
    $(document).on('mouseup', function() {
    isPlotPanning = false;
    $('#plot-container').css('cursor', 'crosshair');
  });

  $(document).on('mousemove', function(e) {
    if (!isPlotPanning) return;
    plotPanX = e.pageX - plotStartX;
    plotPanY = e.pageY - plotStartY;

    updatePlotTransform();
  });


  $('#plot-container').on('contextmenu', function(e) {
    return false;
  });
});
  "))
          # plotOutput(outputId = "xray_plot", height = "750px")
        ),
        conditionalPanel(
          condition = "input.xray_file_uploaded == true",
          fluidRow(
            tags$div(
              "zoom with scroll wheel; pan with right click",
              style = "font-size: 8pt; font-style: italic; color: #555; text-align: center;"
            )
          ),
        fluidRow(
          column(
            width = 6,
            actionBttn(
              inputId = "xray_delete_last_point",
              block = TRUE,
              size = "md",
              label = "Delete Last",
              style = "jelly",
              color = "success",
              icon = icon("delete-left")
            )
          ),
          column(
            width = 6,
            actionBttn(
              size = "md",
              inputId = "xray_reset_points",
              block = TRUE,
              label = "Reset",
              style = "unite",
              color = "danger",
              icon = icon("trash-can")
            )
          ),
        )
        )
      ),
      column(width = 3, 
             conditionalPanel(
        condition = "input.all_points_recorded == true",
        box(title = "Preop Alignment:", 
            width = 12,
            plotOutput(outputId = "preop_spine_simulation_plot",
                       height = "750px"),
            tableOutput("alignment_parameters_df")
        ), 
        box(collapsible = TRUE, collapsed = TRUE,
            tableOutput(outputId = "click_coordinates_df"),
            br(),
            textOutput(outputId = "click_coordinates_text"),
            br(),
            hr(),
            h4("Centroid Coordinates:"),
            tableOutput(outputId = "centroid_coordinates_df"), 
            br(), 
            textOutput("xray_centroid_tibble_text")
            )
      )
      ),
      conditionalPanel(
        condition = "input.all_points_recorded == true",
      box(width = 6,
        title = "Plans",
        actionBttn(
          inputId = "compute_plan_xray",
          label = "Compute Plan",
          style = "unite", 
          color = "danger"
        ),
        hr(),
        fluidRow(
          plotOutput(outputId = "spine_plan_lower_t_xray", height = 650),
        ),
        hr(),
        fluidRow(
          plotOutput(outputId = "spine_plan_upper_t_xray", height = 650),
        )
      )
      )
    )
  )
)


############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################
############################################################# SERVER ##########################################################

# Server logic
server <- function(input, output, session) {

  observeEvent(input$image, {
    req(input$image)
    updateSwitchInput(session = session, inputId = "xray_file_uploaded", value = TRUE)
  })
  
  spine_orientation <- reactiveVal("left")
  
  observeEvent(input$spine_orientation_button, ignoreInit = TRUE, {
    # Update the spine orientation value
    if (spine_orientation() == "left") {
      spine_orientation("right")
    } else {
      spine_orientation("left")
    }
  })
  
  observeEvent(spine_orientation(), {
    # Update the spine orientation value
    if (spine_orientation() == "left") {
      updateActionButton(session, "spine_orientation_button", label = "Facing LEFT", icon = icon("arrow-left"))
    } else {
      updateActionButton(session, "spine_orientation_button", label = "Facing RIGHT", icon = icon("arrow-right"))
    }
  })
  
  observeEvent(click_coord_reactive_list$coords, {
    if(length(click_coord_reactive_list$coords)==3){
      fem_head_x <- click_coord_reactive_list$coords$fem_head_center$x
      s1_anterior_superior_x <- click_coord_reactive_list$coords$s1_anterior_superior$x
      s1_posterior_superior_x <- click_coord_reactive_list$coords$s1_posterior_superior$x

      if(s1_anterior_superior_x < s1_posterior_superior_x){
        spine_orientation("left")
      }else{
        spine_orientation("right")
      }
    }

    ## extreme scenarios
  # if(length(click_coord_reactive_list$coords)== 4){
  # 
  #   fem_head_x <- click_coord_reactive_list$coords$fem_head_center$x
  #   s1_anterior_superior_x <- click_coord_reactive_list$coords$s1_anterior_superior$x
  #   s1_posterior_superior_x <- click_coord_reactive_list$coords$s1_posterior_superior$x
  # 
  #   if((s1_anterior_superior_x < s1_posterior_superior_x) & click_coord_reactive_list$coords$l4_centroid$x >  s1_anterior_superior_x){
  #     spine_orientation("right")
  #   }
  #   
  #   if((s1_anterior_superior_x > s1_posterior_superior_x) & click_coord_reactive_list$coords$l4_centroid$x <  s1_anterior_superior_x){
  #     spine_orientation("left")
  #   }
  #   
  # }
    
    
  }
  )
  
  
  # observeEvent(input$spine_orientation_button, {
  #   # Update the spine orientation value
  #   if (spine_orientation() == "left") {
  #     spine_orientation("right")
  #     updateActionButton(session, "spine_orientation_button", label = "Facing RIGHT", icon = icon("arrow-right"))
  #   } else {
  #     spine_orientation("left")
  #     updateActionButton(session, "spine_orientation_button", label = "Facing LEFT", icon = icon("arrow-left"))
  #   }
  # })
  # 
  # observeEvent(input$xray_click, {
  #   if(length(click_coord_reactive_list$coords)==3){
  #     
  #     fem_head_x <- click_coord_reactive_list$coords$fem_head_center$x
  #     s1_anterior_superior_x <- click_coord_reactive_list$coords$s1_anterior_superior$x
  #     s1_posterior_superior_x <- click_coord_reactive_list$coords$s1_posterior_superior$x
  # 
  #     if(s1_anterior_superior_x < s1_posterior_superior_x){
  #       # xray_orientation <- "left"
  #       spine_orientation("left")
  #       updateActionButton(session, "spine_orientation_button", label = "Facing LEFT", icon = icon("arrow-left"))
  #     }else{
  #       spine_orientation("right")
  #       updateActionButton(session, "spine_orientation_button", label = "Facing RIGHT", icon = icon("arrow-right"))
  #       # xray_orientation <- "right"
  #     }
  #   }
    # 
    # if(length(click_coord_reactive_list$coords)== 4){
    #   
    #   fem_head_x <- click_coord_reactive_list$coords$fem_head_center$x
    #   s1_anterior_superior_x <- click_coord_reactive_list$coords$s1_anterior_superior$x
    #   s1_posterior_superior_x <- click_coord_reactive_list$coords$s1_posterior_superior$x
    #   
    #   if((s1_anterior_superior_x < s1_posterior_superior_x) & click_coord_reactive_list$coords$l4_centroid$x >  s1_anterior_superior_x){
    #     spine_orientation("right")
    #     updateActionButton(session, "spine_orientation_button", label = "Facing RIGHT", icon = icon("arrow-right"))
    #     
    #   }else{
    #     spine_orientation("left")
    #     updateActionButton(session, "spine_orientation_button", label = "Facing LEFT", icon = icon("arrow-left"))
    #     
    #   }
    # }
  # }
  # )
  
  observeEvent(input$image, {
    req(input$image)
    
    # img <- image_read(input$image$datapath)
    img_scaled <- image_scale(image_read(input$image$datapath), "400x")  # Scale to 400px width
    
    # Write the scaled image to a temporary file
    temp_file <- tempfile(fileext = ".jpg")
    image_write(img_scaled, path = temp_file, format = "jpeg")
    
    # Encode the scaled image to base64
    image_src <- base64enc::dataURI(file = temp_file, mime = "image/jpeg")
    
    # Create a base64-encoded URI for the uploaded image
    # image_path <- input$image$datapath
    # image_src <- base64enc::dataURI(file = image_path, mime = input$image$type)
    
    # Send the image URI to the UI
    session$sendCustomMessage('load-image', list(src = image_src))
  })
  
  click_coord_reactive_list <- reactiveValues(coords = list(), index = 1)
  
  
  
  # Reset button to clear all points
  observeEvent(input$xray_reset_points, ignoreInit = TRUE, {
    click_coord_reactive_list$coords <- list()
    click_coord_reactive_list$index <- 1
  })
  
  # Button to remove the last recorded point
  observeEvent(input$xray_delete_last_point, ignoreInit = TRUE, {
    if (click_coord_reactive_list$index > 1) {
      click_coord_reactive_list$coords[[click_coord_reactive_list$index - 1]] <- NULL
      click_coord_reactive_list$index <- click_coord_reactive_list$index - 1
    }
    
    if(length(plot_points_coordinates_reactiveval())>0){
      plot_points_list <- plot_points_coordinates_reactiveval()
      
      plot_points_list <- plot_points_list[-length(plot_points_list)]
      
      plot_points_coordinates_reactiveval(plot_points_list)
      
    }
  })
  
  
  plot_points_coordinates_reactiveval <- reactiveVal()
  
  # Store clicks and assign them to the correct label
  observeEvent(input$xray_click, {
    spine_input_labels <- get_spine_labels(FALSE)
    
    # Only proceed if there's a label available for the current index
    if (click_coord_reactive_list$index <= length(spine_input_labels)) {
      target_name <- spine_input_labels[click_coord_reactive_list$index]
      # Create a named list with the click coordinates
      new_click <- list(x = input$xray_click$x, y = input$xray_click$y)
      
      # Append the new click as a named list with the label
      click_coord_reactive_list$coords[[target_name]] <- new_click
      click_coord_reactive_list$index <- click_coord_reactive_list$index + 1
    }
    
    # Debugging step to print the structure of click_coord_reactive_list$coords
    # print(str(click_coord_reactive_list$coords))
    
    # Convert the named list to a list of lists
    coords_for_js <- unname(lapply(click_coord_reactive_list$coords, function(coord) {
      list(x = coord$x, y = coord$y)
    }))
    
    # Debugging step to print coords_for_js
    # print(coords_for_js)
    
    plot_points_coordinates_reactiveval(coords_for_js)
    
    # Send the coordinates to JavaScript
    # session$sendCustomMessage('plot-coordinates', list(coords = coords_for_js))
  })
  
  observeEvent(list(input$xray_click, input$xray_delete_last_point, input$xray_reset_points), {
    
    session$sendCustomMessage('plot-coordinates', list(coords = plot_points_coordinates_reactiveval()))
  })
  
  
  
  xray_instructions_reactiveval <- reactiveVal("x")
  
  # Render instructions dynamically based on the number of recorded clicks
  output$xray_click_instructions <- renderText({
    spine_input_labels <- get_spine_labels(FALSE)
    click_count <- length(click_coord_reactive_list$coords)
    
    if (click_count < length(spine_input_labels)) {
      instruction <- spine_input_labels[click_count + 1]
      
      instruction <- str_replace_all(instruction, "fem_head_center", "Center of Hips")
      instruction <- str_replace_all(instruction, "_superior", "_superior Corner")
      
      instruction <- str_to_title(str_replace_all(instruction, "_", " "))
      
      instruction <- glue("Click:<br>{instruction}")
      
      xray_instructions_reactiveval("x")
      
    } else {
      instruction <- "All points recorded."
      xray_instructions_reactiveval("Completed")
    }
    HTML("<div>", instruction, "</div>")
  })
  
  
  # Create a reactive table to display coordinates
  click_coordinates_df_reactive <- reactive({
    if (length(click_coord_reactive_list$coords) > 0) {
      # Convert the list to a tibble
      tibble(
        spine_point = names(click_coord_reactive_list$coords),
        x = map_dbl(click_coord_reactive_list$coords, "x"),
        y = map_dbl(click_coord_reactive_list$coords, "y")
      )
    } else {
      tibble(spine_point = character(), x = double(), y = double())
    }
  })
  
  output$click_coordinates_df <- renderTable({
    # click_coordinates_df_reactive()
    
    plotting_coord_list <- plot_points_coordinates_reactiveval()
    
    if (length(plotting_coord_list) > 0) {
      # Convert the list to a tibble
      plotting_coord_df <- tibble(
        spine_point = names(click_coord_reactive_list$coords),
        x = map_dbl(plotting_coord_list, "x"),
        y = map_dbl(plotting_coord_list, "y")
      ) 
      
      plotting_coord_df 
        # mutate(y = abs(y - max(plotting_coord_df$y)))
      
    } else {
      tibble(spine_point = character(), x = double(), y = double())
    }
  })
  
  output$click_coordinates_text <- renderText({
    # click_coordinates_df_reactive()
    
    if (length(plot_points_coordinates_reactiveval()) > 0) {
      # Convert the list to a tibble
      s1_center <- jh_get_point_along_line_function(coord_a = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y), 
                                       coord_b = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y), 
                                       percent_a_to_b = 0.5)
      
      click_df <- tibble(spine_point = "s1_center", x = s1_center[1], y = s1_center[2]) %>%
        union_all(click_coordinates_df_reactive())

        spine_point_labels <- glue_collapse(click_df$spine_point, sep = "', '")
        
        x_values <- glue_collapse(click_df$x, sep = ", ")
        y_values <- glue_collapse(click_df$y, sep = ", ")
        
        glue("click_df <- tibble(spine_point = c('{spine_point_labels}'), x = c({x_values}), y = c({y_values}))")
      
      
    } else {
      glue("click_df <- tibble(spine_point)")
    }
  })
  

  ############ COMPUTE CENTROIDS #################
  # Reactive tibble for centroid coordinates
  xray_centroid_coordinates_reactive_df <- reactive({
    # Get clicked coordinates dataframe
    xray_click_coordinates_df <- click_coordinates_df_reactive()
    
    # Labels for spine points that we need to have in the final dataframe
    spine_coordinate_labels <- tibble(
      spine_point = c("fem_head_center", "s1_center", "l5_centroid", "l4_centroid", "l3_centroid",
                      "l2_centroid", "l1_centroid", "t12_centroid", "t11_centroid", "t10_centroid",
                      "t9_centroid", "t8_centroid", "t7_centroid", "t6_centroid", "t5_centroid",
                      "t4_centroid", "t3_centroid", "t2_centroid", "t1_centroid", "c7_centroid",
                      "c6_centroid", "c5_centroid", "c4_centroid", "c3_centroid", "c2_centroid")
    )
    
    # Calculate S1 center based on anterior and posterior clicks
    s1_center <- if ("s1_anterior_superior" %in% names(click_coord_reactive_list$coords) && 
                     "s1_posterior_superior" %in% names(click_coord_reactive_list$coords)) {
      s1_anterior <- click_coord_reactive_list$coords$s1_anterior_superior
      s1_posterior <- click_coord_reactive_list$coords$s1_posterior_superior
      c((s1_anterior[[1]] + s1_posterior[[1]]) / 2, (s1_anterior[[2]] + s1_posterior[[2]]) / 2)
    } else {
      c(NA, NA)
    }
    
    # Build the tibble including the S1 center
    current_coords_df <-  tibble(spine_point = "s1_center", x = s1_center[1], y = s1_center[2]) %>%
      bind_rows(xray_click_coordinates_df %>%
                  filter(spine_point != "s1_anterior_superior", 
                         spine_point != "s1_posterior_superior",
                         spine_point != "fem_head_center"))
    
    if(any(xray_click_coordinates_df$spine_point == "c2_centroid")){
      
      spine_coordinate_labels_df <- tibble(
        spine_point = c("s1_center", "l5_centroid",
                        "l4_centroid", "l3_centroid", "l2_centroid", "l1_centroid", "t12_centroid",
                        "t11_centroid", "t10_centroid", "t9_centroid", "t8_centroid", "t7_centroid",
                        "t6_centroid", "t5_centroid", "t4_centroid", "t3_centroid", "t2_centroid",
                        "t1_centroid", "c7_centroid", "c6_centroid", "c5_centroid", "c4_centroid",
                        "c3_centroid", "c2_centroid")
      )%>%
        mutate(index_count = row_number())
      
      spine_coordinates_short_df <- spine_coordinate_labels_df %>%
        left_join(current_coords_df) %>%
        filter(!is.na(x))
      
      l5_centroid_y_adjustment <- (spine_coordinates_short_df %>% filter(spine_point == "s1_center"))$y + ((spine_coordinates_short_df %>% filter(spine_point == "l4_centroid"))$y - (spine_coordinates_short_df %>% filter(spine_point == "s1_center"))$y)*0.35
      
      head_df <- spine_coordinates_short_df %>%
        filter(spine_point == "c2_centroid") %>%
        mutate(y = y*1.1) %>%
        mutate(spine_point = "head") %>%
        mutate(index_count = index_count + 1)
      
      # spine_y_filled_df <- spine_coordinate_labels_df %>%
      #   left_join(spine_coordinates_short_df) %>%
      #   mutate(y = zoo::na.spline(y)) %>%
      #   # mutate(y = zoo::na.approx(y, rule = 2)) %>%
      #   mutate(y = round(y, 3)) %>%
      #   mutate(y = if_else(spine_point == "l5_centroid", l5_centroid_y_adjustment, y))
      # 
      # 
      # final_coords_df <- spine_y_filled_df %>%
      #   # mutate(x = spline(spine_coordinates_short_df$y, spine_coordinates_short_df$x, xout = spine_y_filled_df$y)$y) %>%
      #   mutate(x = zoo::na.spline(x))%>%
      #   select(spine_point, x, y) %>%
      #   mutate(x = round(x, 3))
      
      final_coords_df <- spine_coordinate_labels_df %>%
        left_join(spine_coordinates_short_df) %>%
        union_all(head_df) %>%
        mutate(x = zoo::na.spline(x))%>%
        filter(spine_point != "head") %>%
        mutate(y = zoo::na.spline(y)) %>%
        mutate(y = round(y, 3)) %>%
        mutate(y = if_else(spine_point == "l5_centroid", l5_centroid_y_adjustment, y)) %>%
        select(spine_point, x, y) %>%
        mutate(x = round(x, 3))
      
    }else{
      final_coords_df <- current_coords_df  %>%
        filter(!is.na(x))  # Return what has been clicked so far
      
    }
    
    
    final_coords_df%>%
      filter(!is.na(x)) 
      # mutate(y = abs(y - max(final_coords_df$y)))
  })
  
  output$centroid_coordinates_df <- renderTable({
    xray_centroid_coordinates_reactive_df()
  })
  
  output$xray_centroid_tibble_text <- renderText({
    spine_point_labels <- glue_collapse(xray_centroid_coordinates_reactive_df()$spine_point, sep = "', '")
    
    x_values <- glue_collapse(xray_centroid_coordinates_reactive_df()$x, sep = ", ")
    y_values <- glue_collapse(xray_centroid_coordinates_reactive_df()$y, sep = ", ")
    
    fem_head_center <- c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y)
    # 
    s1_anterior_superior <- c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y)
    s1_posterior_superior <- c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y)
    
    glue("fem_head_center <- c({fem_head_center[1]}, {fem_head_center[2]})\n  s1_anterior_superior <- c({s1_anterior_superior[1]}, {s1_anterior_superior[2]})\n s1_posterior_superior <- c({s1_posterior_superior[1]}, {s1_posterior_superior[2]})\n \n centroid_df <- tibble(spine_point = c('{spine_point_labels}'), x = c({x_values}), y = c({y_values}))")
  })
  
  
  alignment_parameters_reactivevalues_list <- reactiveValues()
  
  observeEvent(list(input$xray_click,
                    spine_orientation()), {

                      alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
                      
                      if((any(names(alignment_parameters_list) == "pelvic_incidence") == FALSE) & any(xray_centroid_coordinates_reactive_df()$spine_point == "s1_center")){
                        fem_head_center <- c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y)
                        
                        s1_center <- c((click_coord_reactive_list$coords$s1_anterior_superior$x + click_coord_reactive_list$coords$s1_posterior_superior$x)/2,
                                                       (click_coord_reactive_list$coords$s1_anterior_superior$y + click_coord_reactive_list$coords$s1_posterior_superior$y)/2)
                        
                        #   ### COMPUTE PT ###
                          fem_head_to_s1_length <- jh_calculate_distance_between_2_points_function(point_1 = fem_head_center,
                                                                                                   point_2 = s1_center) ## hypotenuse

                          fem_head_to_s1_x_length <- jh_calculate_distance_between_2_points_function(point_1 = fem_head_center,
                                                                                                     point_2 = c(s1_center[[1]], fem_head_center[[2]])) ## opposite
                        
                            pt_orientation_modifier <- case_when(
                              spine_orientation() == "left" & fem_head_center[[1]] < s1_center[[1]] ~ 1,
                              spine_orientation() == "left" & fem_head_center[[1]] > s1_center[[1]] ~ -1,
                              spine_orientation() == "right" & fem_head_center[[1]] > s1_center[[1]] ~ 1,
                              spine_orientation() == "right" & fem_head_center[[1]] < s1_center[[1]] ~ -1
                            )

                            alignment_parameters_reactivevalues_list$pelvic_tilt <- asin(fem_head_to_s1_x_length/fem_head_to_s1_length)*180/pi*pt_orientation_modifier
                          
                              ### COMPUTE SS ###
                              s1_length <- jh_calculate_distance_between_2_points_function(point_1 = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y),
                                                                                           point_2 = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y))

                              s1_x_length <- jh_calculate_distance_between_2_points_function(point_1 = c(click_coord_reactive_list$coords$s1_anterior_superior$x,
                                                                                                         click_coord_reactive_list$coords$s1_posterior_superior$y),
                                                                                             point_2 = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y))

                              alignment_parameters_reactivevalues_list$sacral_slope <- acos(s1_x_length/s1_length)*180/pi

                              ### COMPUTE PI ###
                              alignment_parameters_reactivevalues_list$pelvic_incidence <- alignment_parameters_reactivevalues_list$pelvic_tilt + alignment_parameters_reactivevalues_list$sacral_slope
                            
                      }
                      
                      ## COMPUTE ALL VPAs ##
                      if(any(xray_centroid_coordinates_reactive_df()$spine_point == "c2_centroid")){

                        vpa_df <- xray_centroid_coordinates_reactive_df() %>%
                          filter(spine_point != "s1_center") %>%
                          mutate(vpa = map2(.x = x, .y = y, .f = ~ jh_compute_vpa_from_xray_data_function(fem_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y),
                                                                                                          vertebral_centroid = c(.x, .y),
                                                                                                          spine_facing = spine_orientation(),
                                                                                                          pelvic_tilt = alignment_parameters_reactivevalues_list$pelvic_tilt
                          )
                          )
                          ) %>%
                          unnest() %>%
                          mutate(vpa_label = str_replace_all(spine_point, "_centroid", "pa")) %>%
                          select(vpa_label, vpa)

                        vpa_list <- as.list(vpa_df$vpa)
                        names(vpa_list) <- vpa_df$vpa_label

                        for (name in names(vpa_list)) {
                          alignment_parameters_reactivevalues_list[[name]] <- vpa_list[[name]]
                        }
                      }

                    }
  )

  
  observeEvent(spine_orientation(), {
                      
                      alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
                      
                      if((any(names(alignment_parameters_list) == "pelvic_incidence") == TRUE)){
                        fem_head_center <- c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y)
                        
                        s1_center <- c((click_coord_reactive_list$coords$s1_anterior_superior$x + click_coord_reactive_list$coords$s1_posterior_superior$x)/2,
                                       (click_coord_reactive_list$coords$s1_anterior_superior$y + click_coord_reactive_list$coords$s1_posterior_superior$y)/2)
                        
                        #   ### COMPUTE PT ###
                        fem_head_to_s1_length <- jh_calculate_distance_between_2_points_function(point_1 = fem_head_center,
                                                                                                 point_2 = s1_center) ## hypotenuse
                        
                        fem_head_to_s1_x_length <- jh_calculate_distance_between_2_points_function(point_1 = fem_head_center,
                                                                                                   point_2 = c(s1_center[[1]], fem_head_center[[2]])) ## opposite
                        
                        pt_orientation_modifier <- case_when(
                          spine_orientation() == "left" & fem_head_center[[1]] < s1_center[[1]] ~ 1,
                          spine_orientation() == "left" & fem_head_center[[1]] > s1_center[[1]] ~ -1,
                          spine_orientation() == "right" & fem_head_center[[1]] > s1_center[[1]] ~ 1,
                          spine_orientation() == "right" & fem_head_center[[1]] < s1_center[[1]] ~ -1
                        )
                        
                        alignment_parameters_reactivevalues_list$pelvic_tilt <- asin(fem_head_to_s1_x_length/fem_head_to_s1_length)*180/pi*pt_orientation_modifier
                        
                        ### COMPUTE SS ###
                        s1_length <- jh_calculate_distance_between_2_points_function(point_1 = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y),
                                                                                     point_2 = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y))
                        
                        s1_x_length <- jh_calculate_distance_between_2_points_function(point_1 = c(click_coord_reactive_list$coords$s1_anterior_superior$x,
                                                                                                   click_coord_reactive_list$coords$s1_posterior_superior$y),
                                                                                       point_2 = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y))
                        
                        alignment_parameters_reactivevalues_list$sacral_slope <- acos(s1_x_length/s1_length)*180/pi
                        
                        ### COMPUTE PI ###
                        alignment_parameters_reactivevalues_list$pelvic_incidence <- alignment_parameters_reactivevalues_list$pelvic_tilt + alignment_parameters_reactivevalues_list$sacral_slope
                        
                      }
                      
                      ## COMPUTE ALL VPAs ##
                      if(any(xray_centroid_coordinates_reactive_df()$spine_point == "c2_centroid")){
                        
                        vpa_df <- xray_centroid_coordinates_reactive_df() %>%
                          filter(spine_point != "s1_center") %>%
                          mutate(vpa = map2(.x = x, .y = y, .f = ~ jh_compute_vpa_from_xray_data_function(fem_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y),
                                                                                                          vertebral_centroid = c(.x, .y),
                                                                                                          spine_facing = spine_orientation(),
                                                                                                          pelvic_tilt = alignment_parameters_reactivevalues_list$pelvic_tilt
                          )
                          )
                          ) %>%
                          unnest() %>%
                          mutate(vpa_label = str_replace_all(spine_point, "_centroid", "pa")) %>%
                          select(vpa_label, vpa)
                        
                        vpa_list <- as.list(vpa_df$vpa)
                        names(vpa_list) <- vpa_df$vpa_label
                        
                        for (name in names(vpa_list)) {
                          alignment_parameters_reactivevalues_list[[name]] <- vpa_list[[name]]
                        }
                      }
                      
                    }
  )
  
  # output$alignment_parameters_df_text <- renderText({
  #   
  #   glue_collapse(reactiveValuesToList(alignment_parameters_reactivevalues_list), sep = " \n")
  #   
  # })
  
  output$alignment_parameters_df <- renderTable({

    enframe(reactiveValuesToList(alignment_parameters_reactivevalues_list)) %>%
      mutate(name = str_replace_all(name, "pelvic_tilt", "PT")) %>%
      mutate(name = str_replace_all(name, "pelvic_incidence", "PI")) %>%
      mutate(name = str_replace_all(name, "sacral_slope", "SS")) %>%
      mutate(name = str_to_upper(name))

  })
  
  observeEvent(xray_instructions_reactiveval(), ignoreInit = TRUE, {
    if(xray_instructions_reactiveval() == "Completed"){
    updateSwitchInput(session = session, inputId = "all_points_recorded", value = TRUE)
    }else{
      updateSwitchInput(session = session, inputId = "all_points_recorded", value = FALSE)
      }
  })

  xray_reactive_plot <- reactive({
    if(xray_instructions_reactiveval() == "Completed"){

      xray <- image_scale(image_read(path = input$image$datapath), "400x")
      
      xray_height <- image_info(xray)$height
      xray_width <- image_info(xray)$width
      
      xlim_left <-0.5 - (xray_width/xray_height)/2
      xlim_right <-0.5 + (xray_width/xray_height)/2
      
      spine_colors_df <- xray_centroid_coordinates_reactive_df() %>%
        mutate(spine_point = str_remove_all(spine_point, "_centroid|_center")) %>%
        filter(spine_point != "s1") %>%
        mutate(spine_color = case_when(
          str_detect(spine_point, "c") ~ "lightblue",
          str_detect(spine_point, "t") ~ "lightgreen",
          str_detect(spine_point, "l") ~ "darkblue"
        )
        )
      
      alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
      
      s1_center <- c((click_coord_reactive_list$coords$s1_anterior_superior$x + click_coord_reactive_list$coords$s1_posterior_superior$x)/2,
                     (click_coord_reactive_list$coords$s1_anterior_superior$y + click_coord_reactive_list$coords$s1_posterior_superior$y)/2)
      
      pi_df <- calculate_pelvic_incidence_line_coordinates(fem_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y),
                                                           s1_anterior = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y),
                                                           s1_posterior = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y),
                                                           spine_facing = spine_orientation(),
                                                           pelvic_tilt = alignment_parameters_list$pelvic_tilt,
                                                           pelvic_incidence_value = alignment_parameters_list$pelvic_incidence
      )
      
      spine_coordinates_df <- tibble(spine_point = "s1_center",
                                     x = s1_center[[1]],
                                     y = s1_center[[2]]) %>%
        union_all(click_coordinates_df_reactive()) %>%
        filter(spine_point %in% c("fem_head_center", "s1_anterior_superior", "s1_posterior_superior") == FALSE)
      
      l1pa_df <- tibble(x = c(s1_center[[1]],
                              click_coord_reactive_list$coords$fem_head_center$x,
                              click_coord_reactive_list$coords$l1_centroid$x),
                        y = c(s1_center[[2]],
                              click_coord_reactive_list$coords$fem_head_center$y,
                              click_coord_reactive_list$coords$l1_centroid$y))
      
      t4pa_df <- tibble(x = c(s1_center[[1]],
                              click_coord_reactive_list$coords$fem_head_center$x,
                              click_coord_reactive_list$coords$t4_centroid$x),
                        y = c(s1_center[[2]],
                              click_coord_reactive_list$coords$fem_head_center$y,
                              click_coord_reactive_list$coords$t4_centroid$y))
      
      # xray_plot <- ggdraw()
      
      spine_xr_build_list <- build_spine_from_coordinates_function(femoral_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y), 
                                                                                s1_anterior_superior = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y), 
                                                                                s1_posterior_superior = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y), 
                                                                                centroid_df = xray_centroid_coordinates_reactive_df(), 
                                                                                spine_facing = spine_orientation())
      

      # xray_plot <-  ggdraw(xlim = c(0, xray_width), ylim = c(0, xray_height)) +
      xray_plot <- ggdraw() +
        draw_image(
          scale = xray_height,
                   xray,
                   x = 0,
                   halign = 0, 
                   valign = 0,
                   y = 0,
                   width = 1,
                   # height = 1,
                   clip = FALSE
        ) +
        geom_path(data = xray_centroid_coordinates_reactive_df(),
                  aes(x = x, y = y), color = "lightblue", size = 0.2) +
        geom_point(data = spine_colors_df,
                   aes(x = x, y = y, color = spine_color, fill = spine_color),size = 0.2) +
        scale_fill_identity() +
        scale_color_identity()+
        geom_path(data = pi_df,
                  aes(x = x, y = y), color = "darkgreen", size = 0.25)+
        geom_path(data = l1pa_df,
                  aes(x = x, y = y),
                  color = "darkblue", size = 0.25)+
        geom_path(data = t4pa_df,
                  aes(x = x, y = y),
                  color = "purple", size = 0.25) +
        geom_segment(data = spine_xr_build_list$sup_endplates_df, aes(x = sp_x, y = sp_y, xend = sa_x, yend = sa_y), color = "red", size = 0.2, lineend = "round") + 
        coord_fixed(xlim = c(0, xray_width), ylim = c(0, xray_height))
      
     
      
      xray_plot
      
    }
    
   
    
    
  })
  
  # output$xray_plot <- renderPlot({
  #   xray_reactive_plot()
  # }
  # )
  
  
    observe({
      req(input$image)  # Ensure there's an image uploaded
      
      # Generate the plot
      plot <- xray_reactive_plot()
      
      # Save the plot as a temporary file
      temp_file <- tempfile(fileext = ".png")
      # ggsave(temp_file, plot = plot, width = 8, height = 10, dpi = 150)
      ggsave(temp_file, plot = plot, width = 350, height = 700, units = "px")
      
      # Read the image back and convert it to a base64 string for embedding
      img <- magick::image_read(temp_file)
      img_base64 <- base64enc::dataURI(file = temp_file, mime = "image/png")
      
      # Send the image to the UI
      session$sendCustomMessage('load-plot-image', list(src = img_base64))
    })

  
  
  output$preop_xray_rigid_segments_ui <- renderUI({
    # preop_segment_angles_list <- preop_segment_angles_list_reactive()
    
    segment_angles_input_list <- rev(map(spinal_segments_labels_vector, function(label) {
      create_spine_rigid_level_input_function(segment_input_label = label)
    }))
    # create_spine_rigid_level_input_function
    
    # column(width = 5,
    fluidRow(
      column(width = 12,
      h4("Select Any Fused Levels:"),
      box(width = 12,
          title = "Cervical", 
          collapsible = TRUE, 
          collapsed = TRUE,
          # h5("Check box if Rigid Level"),
          segment_angles_input_list[1:6] %>% tagList()
      ),
      box(width = 12,title = "Thoracic", 
          collapsible = TRUE, 
          collapsed = TRUE,
          # h5("Check box if Rigid Level"),
          segment_angles_input_list[7:19] %>% tagList()
      ),
      box(width = 12,title = "Lumbar", 
          collapsible = TRUE, 
          collapsed = FALSE,
          # h5("Check box if Rigid Level"),
          segment_angles_input_list[20:24]%>% tagList()
      )
    )
    )
  })
  
  preop_rigid_levels_vector_reactive_xray <- reactive({
    # segment_id <- paste0("preop_", str_to_lower(str_replace_all(segment_input_label, pattern = "-", "_")), "_rigid")
    rigid_levels <- c("na")
    if(input$preop_l5_s1_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "L5-S1")} 
    if(input$preop_l4_l5_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "L4-L5")} 
    if(input$preop_l3_l4_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "L3-L4")} 
    if(input$preop_l2_l3_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "L2-L3")} 
    if(input$preop_l1_l2_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "L1-L2")} 
    if(input$preop_t12_l1_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T12-L1")} 
    if(input$preop_t11_t12_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T11-T12")} 
    if(input$preop_t10_t11_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T10-T11")} 
    if(input$preop_t9_t10_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T10-T9")} 
    if(input$preop_t8_t9_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T8-T9")} 
    if(input$preop_t7_t8_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T7-T8")} 
    if(input$preop_t6_t7_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T6-T7")} 
    if(input$preop_t5_t6_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T5-T6")} 
    if(input$preop_t4_t5_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T4-T5")} 
    if(input$preop_t3_t4_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T3-T4")} 
    if(input$preop_t2_t3_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T2-T3")} 
    if(input$preop_t1_t2_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "T1-T2")} 
    if(input$preop_c7_t1_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C7-T1")} 
    if(input$preop_c6_c7_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C6-C7")} 
    if(input$preop_c5_c6_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C5-C6")} 
    if(input$preop_c4_c5_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C4-C5")} 
    if(input$preop_c3_c4_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C3-C4")} 
    if(input$preop_c2_c3_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C2-C3")} 
    if(input$preop_c1_c2_rigid_xray == TRUE){rigid_levels <- append(rigid_levels, "C1-C2")} 
    # set_names(map_chr(spinal_segments_labels_vector, ~ paste0(str_to_lower(strsplit(.x, "-")[[1]][1]), "_segment_angle")))
    if(length(rigid_levels)>0){
      rigid_levels <- map_chr(rigid_levels, ~ paste0(str_to_lower(strsplit(.x, "-")[[1]][1]), "_segment"))
    }
    rigid_levels
  })
  
  spine_build_from_coordinates_reactive <- reactive({
    if(any(names(click_coord_reactive_list$coords) == "c2_centroid")){
    spine_build_list <- build_spine_from_coordinates_function(femoral_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y), 
                                                              s1_anterior_superior = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y), 
                                                              s1_posterior_superior = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y), 
                                                              centroid_df = xray_centroid_coordinates_reactive_df(), 
                                                              spine_facing = spine_orientation())
    }
    
  })
  
  output$preop_spine_simulation_plot <- renderPlot({
    if(xray_instructions_reactiveval() == "Completed"){
      alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
      spine_build_list <- spine_build_from_coordinates_reactive()
      
      ggplot() +
        geom_sf(data = spine_build_list$vert_geom_df,
                color = "black",
                aes(geometry = geometry), 
                fill = "grey90", 
                alpha = 0.9) +
        geom_sf(data = spine_build_list$spine_geom_list$c1_geom, 
                # fill = "grey90", 
                fill = "grey90",
                alpha = 0.5) +
        geom_path(data =  spine_build_list$vert_geom_df, aes(x = x, y = y)) +
        geom_sf(data = spine_build_list$fem_head_sf, 
                fill = "grey90") +
        geom_sf(data = spine_build_list$sacrum_sf, 
                fill = "grey90", 
                alpha = 0.8) +
        geom_sf(data = spine_build_list$lines_list$l1pa, 
                color = "blue")+
        geom_sf(data = spine_build_list$lines_list$t9pa, 
                color = "orange") +
        geom_sf(data = spine_build_list$lines_list$t4pa, 
                color = "purple") +
        geom_sf(data = spine_build_list$lines_list$c2pa, 
                color = "darkgreen") +
        theme_void() +
        labs(title = "Preop Alignment") +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(
            size = 16,
            hjust = 0.5,
            vjust = -0.5,
            face = "bold.italic"
          ),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA)
        ) 
        # geom_segment(data = spine_build_list$sup_endplates_df, aes(x = sp_x, y = sp_y, xend = sa_x, yend = sa_y), color = "red", linewidth = 1, lineend = "round") 
      
      
      
    }
    
    # if(any(names(alignment_parameters_list) == "c2pa")){
      
      # if(any(names(click_coord_reactive_list$coords) == "c2_centroid")){
        
        # xray_centroid_coordinates_reactive_df()
        
        # spine_build_list <- build_spine_from_coordinates_function(femoral_head_center = c(click_coord_reactive_list$coords$fem_head_center$x, click_coord_reactive_list$coords$fem_head_center$y), 
        #                                                  s1_anterior_superior = c(click_coord_reactive_list$coords$s1_anterior_superior$x, click_coord_reactive_list$coords$s1_anterior_superior$y), 
        #                                                  s1_posterior_superior = c(click_coord_reactive_list$coords$s1_posterior_superior$x, click_coord_reactive_list$coords$s1_posterior_superior$y), 
        #                                                  centroid_df = xray_centroid_coordinates_reactive_df(), 
        #                                                  spine_facing = spine_orientation())
        
        
        # spine_simulation_list <- build_full_spine_from_vertebral_pelvic_angles_function(pelv_inc_value = alignment_parameters_list$pelvic_incidence,
        #                                                                                 pt_value = alignment_parameters_list$pelvic_tilt,
        #                                                                                 l1pa_value_input = alignment_parameters_list$l1pa,
        #                                                                                 # l1s1_value_input = alignment_parameters_list$preop_l1s1,
        #                                                                                 # t10_l2_value_input = alignment_parameters_list$preop_t10_l2,
        #                                                                                 t9pa_value_input = alignment_parameters_list$t9pa,
        #                                                                                 t4pa_value_input = alignment_parameters_list$t4pa,
        #                                                                                 c2pa_value_input = alignment_parameters_list$c2pa,
        #                                                                                 # c2_c7_value_input = alignment_parameters_list$preop_c2c7,
        #                                                                                 # input_segment_angles = "yes",
        #                                                                                 # segment_angles_input = segment_angles_list,
        #                                                                                 spine_faces = spine_orientation()
        # )
        # 
        # spine_geoms_df <- spine_simulation_list$spine_df %>%
        #   select(object, geom, geom_alpha)
        
        # measurements_df <- spine_simulation_list$measurements_df
        
        # ggplot() +
        #   geom_sf(data = spine_geoms_df,
        #           color = "black",
        #           aes(geometry = geom,
        #               alpha = geom_alpha,
        #               fill = "grey90")) +
        #   # ylim(-20, 110) +
        #   theme_void() +
        #   labs(title = "Preop Alignment") +
        #   theme(
        #     axis.text = element_blank(),
        #     axis.title = element_blank(),
        #     plot.title = element_text(
        #       size = 16,
        #       hjust = 0.5,
        #       vjust = -0.5,
        #       face = "bold.italic"
        #     ),
        #     plot.background = element_rect(fill = "transparent", colour = NA),
        #     panel.background = element_rect(fill = "transparent", colour = NA)
        #   ) +
        #   draw_text(text = spine_simulation_list$measurements_df$label,
        #             x = spine_simulation_list$measurements_df$x,
        #             y = spine_simulation_list$measurements_df$y, size = 11) +
        #   scale_fill_identity() +
        #   scale_alpha_identity()+
        #   xlim(-30, 30)
      # } 
      # coord_fixed()
      
    # }
    
    
  })
  
  
  
  spine_plan_uiv_t11_option_1_xray <- eventReactive(input$compute_plan_xray, {
    alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
    
    spine_build_list <- spine_build_from_coordinates_reactive()
    
    build_t11_spine_plot_function(pso_option_number = 1, 
                                  preop_age = input$preop_age,
                                  preop_sex = input$preop_sex,
                                  preop_pelvic_incidence = alignment_parameters_list$pelvic_incidence,
                                  preop_pt = alignment_parameters_list$pelvic_tilt,
                                  preop_l1pa = alignment_parameters_list$l1pa,
                                  preop_t9pa = alignment_parameters_list$t9pa,
                                  preop_t4pa = alignment_parameters_list$t4pa,
                                  preop_c2pa = alignment_parameters_list$c2pa,
                                  # l1pa_line_color = input$l1pa_line_color,
                                  # t4pa_line_color = input$t4pa_line_color,
                                  # c2pa_line_color = input$c2pa_line_color, 
                                  # preop_segment_angles_input_list_reactive = spine_build_list$segment_angles_list,
                                  preop_rigid_levels_vector_reactive = preop_rigid_levels_vector_reactive_xray(),
                                  return_list_or_plot = "list"
    )
  })
  
  
  spine_plan_uiv_t11_option_2_xray <- eventReactive(input$compute_plan_xray, {
    alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
    
    spine_build_list <- spine_build_from_coordinates_reactive()
    
    build_t11_spine_plot_function(pso_option_number = 2, 
                                  preop_age = input$preop_age,
                                  preop_sex = input$preop_sex,
                                  preop_pelvic_incidence = alignment_parameters_list$pelvic_incidence,
                                  preop_pt = alignment_parameters_list$pelvic_tilt,
                                  preop_l1pa = alignment_parameters_list$l1pa,
                                  preop_t9pa = alignment_parameters_list$t9pa,
                                  preop_t4pa = alignment_parameters_list$t4pa,
                                  preop_c2pa = alignment_parameters_list$c2pa,
                                  # l1pa_line_color = input$l1pa_line_color,
                                  # t4pa_line_color = input$t4pa_line_color,
                                  # c2pa_line_color = input$c2pa_line_color, 
                                  # preop_segment_angles_input_list_reactive = spine_build_list$segment_angles_list,
                                  preop_rigid_levels_vector_reactive = preop_rigid_levels_vector_reactive_xray(), 
                                  return_list_or_plot = "list"
    )
  })
  
  
  output$spine_plan_lower_t_xray <- renderPlot({
    
    if(spine_plan_uiv_t11_option_1_xray()$pso_level == "none"){
      
      plot_grid(NULL,
                spine_plan_uiv_t11_option_1_xray()$prescribed_plot, 
                nrow = 2, rel_heights = c(0.075, 0.925)
      ) + 
        draw_text(text = "Prescribed Alignment: Lower Thoracic UIV", x = 0.5, y = 0.975, size = 18, fontfact = "bold")
      
    }else{
      plot_grid(NULL, NULL, NULL,
                spine_plan_uiv_t11_option_1_xray()$prescribed_plot, NULL, spine_plan_uiv_t11_option_2_xray()$prescribed_plot, 
                nrow = 2, rel_heights = c(0.075, 0.925), 
                rel_widths = c(1, -0.15, 1)
      ) + 
        draw_text(text = "Prescribed Alignment: Lower Thoracic UIV", x = 0.5, y = 0.975, size = 18, fontfact = "bold")
      
    }
    
    
    
  })
  
  spine_plan_uiv_t4_option_1_xray <- eventReactive(input$compute_plan_xray, {
    alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
    spine_build_list <- spine_build_from_coordinates_reactive()
    
    build_upper_t_uiv_spine_plot_function(pso_option_number = 1, 
                                          preop_age = input$preop_age,
                                          preop_sex = input$preop_sex,
                                          preop_pelvic_incidence = alignment_parameters_list$pelvic_incidence,
                                          preop_pt = alignment_parameters_list$pelvic_tilt,
                                          preop_l1pa = alignment_parameters_list$l1pa,
                                          preop_t9pa = alignment_parameters_list$t9pa,
                                          preop_t4pa = alignment_parameters_list$t4pa,
                                          preop_c2pa = alignment_parameters_list$c2pa,
                                          # l1pa_line_color = input$l1pa_line_color,
                                          # t4pa_line_color = input$t4pa_line_color,
                                          # c2pa_line_color = input$c2pa_line_color, 
                                          # preop_segment_angles_input_list_reactive = spine_build_list$segment_angles_list,
                                          preop_rigid_levels_vector_reactive = preop_rigid_levels_vector_reactive_xray(),
                                          return_list_or_plot = "list"
    )
  })
  
  
  ############## PSO OPTION 2 UIV T4
  spine_plan_uiv_t4_option_2_xray <- eventReactive(input$compute_plan_xray, {
    alignment_parameters_list <- reactiveValuesToList(alignment_parameters_reactivevalues_list)
    
    spine_build_list <- spine_build_from_coordinates_reactive()
    
    build_upper_t_uiv_spine_plot_function(pso_option_number = 2, 
                                          preop_age = input$preop_age,
                                          preop_sex = input$preop_sex,
                                          preop_pelvic_incidence = alignment_parameters_list$pelvic_incidence,
                                          preop_pt = alignment_parameters_list$pelvic_tilt,
                                          preop_l1pa = alignment_parameters_list$l1pa,
                                          preop_t9pa = alignment_parameters_list$t9pa,
                                          preop_t4pa = alignment_parameters_list$t4pa,
                                          preop_c2pa = alignment_parameters_list$c2pa,
                                          # l1pa_line_color = input$l1pa_line_color,
                                          # t4pa_line_color = input$t4pa_line_color,
                                          # c2pa_line_color = input$c2pa_line_color, 
                                          # preop_segment_angles_input_list_reactive = spine_build_list$segment_angles_list,
                                          preop_rigid_levels_vector_reactive = preop_rigid_levels_vector_reactive_xray(), #preop_rigid_levels_vector_reactive()
                                          return_list_or_plot = "list"
    )
  })
  
  output$spine_plan_upper_t_xray <- renderPlot({
    
    if(spine_plan_uiv_t4_option_1_xray()$pso_level == "none"){
      
      plot_grid(NULL,
                spine_plan_uiv_t4_option_1_xray()$prescribed_plot, 
                nrow = 2, rel_heights = c(0.075, 0.925)
      ) + 
        draw_text(text = "Prescribed Alignment: Upper Thoracic UIV", x = 0.5, y = 0.975, size = 18, fontfact = "bold")
      
    }else{
      plot_grid(NULL, NULL, NULL,
                spine_plan_uiv_t4_option_1_xray()$prescribed_plot, NULL, spine_plan_uiv_t4_option_2_xray()$prescribed_plot, 
                nrow = 2, rel_heights = c(0.075, 0.925), 
                rel_widths = c(1, -0.15, 1)
      ) + 
        draw_text(text = "Prescribed Alignment: Upper Thoracic UIV", x = 0.5, y = 0.975, size = 18, fontfact = "bold")
      
    }
    
  })
  
  
  
}

# Run the app
shinyApp(ui = ui, server = server)
