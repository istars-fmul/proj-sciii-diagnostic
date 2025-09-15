library(yaml)
library(plotly)
library(here)
cfg <- yaml.load_file("../config.yml")
font_family <- cfg$style$family
source(here("scripts", "utils.R"))


# plot_gpa function: Plots the alignment of landmarks before and after Generalized Procrustes Analysis (GPA)
# Inputs: 
#   - proc_gpa: GPA object
#   - coordinates_data: Coordinates data (matrix[,,]) - > before the GPA
#   - joinlines: List of lists containing indices for connecting landmarks
#   - output_path: Path to save the plot
# Output:
#   - Combined plot of the alignment of landmarks before and after GPA
#   - Individual plots of the alignment of landmarks before and after GPA
plot_gpa <- function(proc_gpa, coordinates_data, joinlines, output_path) {

    # Generate colors for markers
    marker_colors <- rainbow(dim(coordinates_data)[1])
    # Plot1: Plot before GPA
    plot1 <- plot_ly()


    for(i in 1:dim(coordinates_data)[3]){
        plot1 <- add_trace(plot1,
                        type = "scatter",
                        mode = "markers",
                        x = coordinates_data[, 1, i],
                        y = coordinates_data[, 2, i],
                        marker = list(color = marker_colors, symbol = "circle", size = 5),
                        showlegend = FALSE)
    }

    # Set layout without legend
    plot1 <- layout(plot1,
                title = "",
                xaxis = list(title = "X", autorange = TRUE),
                yaxis = list(title = "Y", autorange = TRUE),
                legend = TRUE)



    # Plot2: Plot after GPA

    plot2 <- plot_ly()

    # Adding point by point of a single patient, to be used as legend
    sample <- proc_gpa$rotated[,,1]
    for (c in 1:dim(sample)[1]){
        label <- c("Sella","Basion","PNS","A Point","B Point","Pogonion","Menton","Gonion","Ramus Point","Distal Aspect of Condyle","Condylion","Nasion")[c]
        plot2 <- add_trace(plot2, 
                        name = label,
                        type = "scatter", 
                        mode = "markers", 
                        x = sample[c,1],
                        y = sample[c,2],
                        marker = list(color = marker_colors[c], symbol = "circle", size = 5),
                        showlegend = TRUE)
    }


    # Adding scatter plot for points after GPA
    for (i in 2:dim(proc_gpa$rotated)[3]){
        plot2 <- add_trace(plot2, 
                        type = "scatter", 
                        mode = "markers", 
                        x = proc_gpa$rotated[,1,i], 
                        y = proc_gpa$rotated[,2,i],
                        marker = list(color = marker_colors, symbol = "circle", size = 5),
                        showlegend = FALSE)
    }

    # Adding lines connecting landmarks
    for (l in seq_along(joinlines)){
        plot2 <- add_trace(plot2,
                        type = "scatter",
                        mode = "lines",
                        x = proc_gpa$mshape[joinlines[[l]], 1],
                        y = proc_gpa$mshape[joinlines[[l]], 2],
                        line = list(color = "black", width = 2),
                        name = "Mean Shape",
                        showlegend = ifelse(l == 1, TRUE, FALSE))
    }

    # Set layout for the plot after GPA
    plot2 <- layout(plot2,
                title = "",
                xaxis = list(title = "X", autorange = TRUE, titlefont = list(family = font_family)),
                yaxis = list(title = "Y", autorange = TRUE, titlefont = list(family = font_family)),
                legend = TRUE)


    # Combine the plots
    combined_plot <- subplot(plot1, plot2) %>%
    layout(
        title = "",
        annotations = list(
        list(
            x = 0.2,
            y = 1.0,
            text = "Before GPA",
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE
        ),
        list( 
            x = 0.8,
            y = 1.0,
            text = "After GPA",
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE
        )
        )
    )

    # Save the combined plot
    save_image_svg(combined_plot, output_path = paste0(output_path, "plot_gpa.svg"), width = 1400, height = 800)

    # Save separe plots
    # save plot 1: before GPA
    save_image_svg(plot1, output_path = paste0(output_path, "plot_gpa_1.svg"), width = 600, height = 600)
    # save plot 2: after GPA
    save_image_svg(plot2, output_path = paste0(output_path, "plot_gpa_2.svg"), width = 700, height = 600)
    
    return(combined_plot)
}