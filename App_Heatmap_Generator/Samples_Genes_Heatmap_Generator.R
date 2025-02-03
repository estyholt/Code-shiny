# Rscript Samples_Genes_Heatmap_Generator.R --expression normalised_expr.txt --metadata metadata.txt --name xxxx


library(tidyverse)
library(shiny)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(shinycssloaders)
library(optparse)
library(shinyjs)
library(grid)



# Command-line options using optparse
option_list <- list(
  make_option(c("--expression"), type = "character", default = "normalised_expr.txt", 
              help = "Path to the expression data file"),
  make_option(c("--metadata"), type = "character", default = "metadata.txt", 
              help = "Path to the metadata file"),
  make_option(c("--name"), type = "character", default = "Heatmap",
              help = "Heatmap name to display in the heatmap title")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Assign command-line arguments to variables
expr_file_path <- opts$expression
meta_file_path <- opts$metadata
heatmap_name <- opts$name

#####################################################################
####################   UI - User Interface   ########################
#####################################################################

ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs
  titlePanel("Interactive Heatmap Generator"),
  
  
  sidebarLayout(
    sidebarPanel(
      # Dropdown for selecting metadata columns to group by (default: based on logic below)
      selectInput("group_by", "Categories to group/cluster by:", 
                  choices = NULL, multiple = TRUE, selectize = TRUE),
      
      
      # Gene selection/pasting dropdown:
      fluidRow(
        column(6,
               selectInput("geneSelect", "Choose Genes to Display:", 
                           choices = NULL, multiple = TRUE, selectize = TRUE)
        ),
        column(6,
               textAreaInput("geneListInput", "Or Paste A List of Genes:", value = "", 
                             placeholder = "Maximum 180 genes",
                             rows = 4),
               # Button to update gene list
               actionButton("updateGenes", "Update Genes")
        )
      ),
      
      textInput("title_input", "Change Heatmap Title", value = heatmap_name),
      
      # Input field to allow the user to change the color scale name
      textInput("scale_name_input", "Change Color Scale Name", value = "Color Scale"),
      
      numericInput("fontsize_row_input", "Gene Names Font Size", value = 10, min = 4, max = 50),
      
      # Checkbox to select top 10 genes with highest STD
      checkboxInput("selectTop", "Select Top 10 Genes (highest std)", value = TRUE),
      
      # Checkbox to remove all genes
      checkboxInput("removeAll", "Remove All Genes", value = FALSE),
      
      # Checkbox input for displaying column names
      checkboxInput("col_names", "Show Sample Names", value = FALSE),
      
      # Checkbox input for clustering columns
      checkboxInput("cluster_cols", "Cluster Columns (hierarchical clustering)", value = FALSE),
      
      # Checkbox input for clustering rows
      checkboxInput("cluster_rows", "Cluster Rows (hierarchical clustering)", value = FALSE),
      
      
      # Radio buttons for scaling options
      radioButtons("scale", "Scale Values:",
                   choices = c("None" = 0, "By Row" = 1, "By Column" = 2),
                   selected = 1)
    ),
    
    mainPanel(
      #A loading spinner (animation) is shown. 
      #Great for cases when it takes some time to load the big data.
      withSpinner(plotOutput("heatmapPlot", height = "800px"))
    )
  )
)

#####################################################################
####################    Server logic   ##############################
#####################################################################

server <- function(input, output, session) {
  
  MAX_GENES = 180
  MAX_CATEGORIES = 8
  
  # Load normalized expression data
  norm_expr_data <- read.table(expr_file_path, sep = "\t", header = TRUE, 
                               row.names = 1, check.names = FALSE)
  norm_expr_data <- as.matrix(norm_expr_data)
  
  # Load metadata
  metadata_data <- read.table(meta_file_path, sep = "\t", header = TRUE, 
                              row.names = 1, check.names = FALSE)
  
  # Top 10 genes with highest standard deviation across the rows:
  top_genes <- head(rownames(norm_expr_data)[order(apply(norm_expr_data, 1, sd), decreasing = TRUE)], 10)
  
  # Populate the metadata dropdown dynamically and set default selection based on "Disease" and "Treatment"
  observe({
    metadata_cols <- colnames(metadata_data)
    
    # Define default group columns
    default_group_by <- NULL
    
    # If both "Disease" and "Treatment" are in the columns, group by both
    if ("Disease" %in% metadata_cols & "Treatment" %in% metadata_cols) {
      default_group_by <- c("Disease", "Treatment")
    }
    # If only "Disease" is in the columns, group by it
    else if ("Disease" %in% metadata_cols) {
      default_group_by <- c("Disease")
    }
    # If only "Treatment" is in the columns, group by it
    else if ("Treatment" %in% metadata_cols) {
      default_group_by <- c("Treatment")
    }
    # If neither "Disease" nor "Treatment" are found, group by the first two columns
    else {
      default_group_by <- metadata_cols[1:2]
    }
    
    # Update the selectInput for metadata grouping
    updateSelectInput(session, "group_by", 
                      choices = metadata_cols, 
                      selected = default_group_by)  # Default: group by "Disease" and/or "Treatment" if available, else first two columns
  })
  # Observe the group_by input and restrict it to a maximum of 8 selections
  observe({
    if (length(input$group_by) > MAX_CATEGORIES) {
      updateSelectInput(session, "group_by", selected = head(input$group_by, MAX_CATEGORIES))
      showNotification(paste("You can select a maximum of ", MAX_CATEGORIES, " categories.", sep = ""), type = "warning")
    }
  })
  
  # Process metadata based on user selection
  process_metadata <- reactive({
    req(input$group_by)  # Ensure user has selected columns
    metadata <- metadata_data
    
    # Convert selected columns to factors
    for (col in input$group_by) {
      if (is.character(metadata[[col]]) || is.factor(metadata[[col]])) {
        metadata[[col]] <- factor(metadata[[col]])
      }
    }
    
    # Order metadata based on selected grouping
    ordered_indices <- do.call(order, metadata[input$group_by])
    
    list(ordered_metadata = metadata[ordered_indices, , drop = FALSE], 
         ordered_indices = ordered_indices)
  })
  
  # Generate annotation colors for metadata
  generate_annotation_colors <- function(metadata) {
    ann_colors <- list()
    color_palettes <- brewer.pal.info[brewer.pal.info$category == "qual", ]
    palette_names <- rownames(color_palettes)
    
    for (i in seq_along(colnames(metadata))) {
      col_name <- colnames(metadata)[i]
      if (is.factor(metadata[[col_name]])) {
        levels <- levels(metadata[[col_name]])
        
        if (length(levels) == 2) {
          colors <- brewer.pal(3, "Set1")[1:2]  
        } else {
          palette <- palette_names[(i - 1) %% length(palette_names) + 1]
          colors <- brewer.pal(n = min(length(levels), color_palettes[palette, "maxcolors"]), name = palette)
        }
        
        ann_colors[[col_name]] <- setNames(colors, levels)
      }
    }
    ann_colors
  }
  
  # Populate gene selection dropdown dynamically
  updateSelectInput(session, "geneSelect", 
                    choices = rownames(norm_expr_data), 
                    selected = top_genes)
  
  # Observe "Select Top Genes" checkbox
  observeEvent(input$selectTop, {
    if (input$selectTop) {
      updateSelectInput(session, "geneSelect", selected = top_genes)  # Select top genes
      updateSelectInput(session, "geneListInput", selected = character(0))  # Clear all selections
      updateCheckboxInput(session, "removeAll", value = FALSE)  # Ensure "Remove All" is unchecked
    }
  })
  
  # Observe "Remove All Genes" checkbox
  observeEvent(input$removeAll, {
    if (input$removeAll) {
      updateSelectInput(session, "geneSelect", selected = character(0))  # Clear all selections
      updateSelectInput(session, "geneListInput", selected = character(0))  # Clear all selections
      updateCheckboxInput(session, "selectTop", value = FALSE)  # Ensure "Select Top" is unchecked
    }
  })
  
  
  
  # What happens when more than MAX_GENES are selected:
  observeEvent(input$geneSelect, {
    # Clear the geneListInput when a gene is selected
    updateTextAreaInput(session, "geneListInput", value = "")
    
    selected_count <- length(input$geneSelect)
    # Prevent selection of more than MAX_GENES genes and show message
    if (selected_count > MAX_GENES) {
      updateSelectInput(session, "geneSelect", selected = head(input$geneSelect, MAX_GENES))  # Limit to MAX_GENES genes
      output$geneLimitMessage <- renderText("You can select a maximum of 180 genes.")
      
    } else {
      output$geneLimitMessage <- renderText("")  # Clear message if within limit
    }
    # What happens when there is a manual selection of genes:
    if (setequal(input$geneSelect, top_genes)) {
      updateCheckboxInput(session, "selectTop", value = TRUE)   # Check "Select Top Genes"
      updateCheckboxInput(session, "removeAll", value = FALSE)  # Uncheck "Remove All Genes"
    } else if (selected_count == 0) {
      updateCheckboxInput(session, "selectTop", value = FALSE)  # Uncheck "Select Top Genes"
      updateCheckboxInput(session, "removeAll", value = TRUE)   # Check "Remove All Genes"
    } else {
      updateCheckboxInput(session, "selectTop", value = FALSE)  # Uncheck "Select Top Genes"
      updateCheckboxInput(session, "removeAll", value = FALSE)  # Uncheck "Remove All Genes"
    }
  })
  
  observeEvent(input$updateGenes, {
    # Split the input by newlines, commas, spaces, quotes (any type), or tabs:
    genes <- unlist(strsplit(input$geneListInput, '[\n,\'\" \t]+'))
    
    # Remove empty strings (in case of extra newlines or spaces)
    genes <- genes[genes != ""]
    
    genes_upper <- toupper(genes)  # Convert input genes to uppercase for matching
    rownames_upper <- toupper(rownames(norm_expr_data))  # Convert row names to uppercase for matching
    
    # Match genes in uppercase
    valid_genes_upper <- genes_upper[genes_upper %in% rownames_upper]
    
    # Find the corresponding original gene names from your dataset
    valid_genes <- rownames(norm_expr_data)[match(valid_genes_upper, rownames_upper)]
    
    valid_genes <- unique(valid_genes)  # Ensure unique genes
    # Identify genes that don't exist in the data
    invalid_genes <- genes_upper[!(genes_upper %in% rownames_upper)]
    invalid_genes <- invalid_genes[invalid_genes != ""]  # Remove empty strings from invalid_genes
    
    if (length(invalid_genes) > 0) {
      showNotification(paste("Warning: The following genes do not exist in the data:", 
                             paste(invalid_genes, collapse = ", ")), type = "warning")
    }
    
    if (length(valid_genes) > MAX_GENES) {
      valid_genes <- valid_genes[1:MAX_GENES]  # Limit to 180 genes
    }
    
    if (length(valid_genes) > 0) {
      updateSelectInput(session, "geneSelect", selected = valid_genes)
    }
  })
  
  
  # Reactive expression for selected genes
  selected_data <- reactive({
    req(input$geneSelect)
    norm_expr_data[input$geneSelect, , drop = FALSE]
  })
  
  
  
  
  # Render the heatmap
  output$heatmapPlot <- renderPlot({
    req(selected_data(), input$group_by)
    
    metadata_info <- process_metadata()
    data <- selected_data()
    
    ordered_indices <- metadata_info$ordered_indices
    ordered_metadata <- metadata_info$ordered_metadata
    ordered_metadata <- as.data.frame(ordered_metadata)
    
    # Order expression data columns
    data_ordered <- data[, ordered_indices, drop = FALSE]
    
    # Determine scale option
    scale_option <- if (input$scale == 1) "row" else if (input$scale == 2) "column" else "none"
    
    # Generate annotation colors
    ann_colors <- generate_annotation_colors(ordered_metadata)
    
    # Generate the heatmap:
    ComplexHeatmap::pheatmap(data_ordered,
                             main = input$title_input,
                             name = input$scale_name_input,
                             border_color = "grey60",
                             show_rownames = TRUE,
                             show_colnames = input$col_names,
                             scale = scale_option,
                             annotation_col = ordered_metadata,
                             annotation_colors = ann_colors,
                             cluster_cols = input$cluster_cols,
                             cluster_rows = input$cluster_rows,
                             fontsize_number = 8,
                             fontsize_row = input$fontsize_row_input, 
                             na_col = "white")
    
    
  })
}

# Run the application
shinyApp(ui = ui, server = server)
