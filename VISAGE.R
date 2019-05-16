# Code to run VISAGE from viability data to .rnk file for GSEA.  Put this script in a directory that contains:
# VISAGE_functions.R
# VISAGE_plots.R
# visage_params.txt
# And a subdirectory containing:
#   Expression data for cell lines
#   CL_{cell_line_name}.txt viability data for every cell line.
# See the user's guide available at https://github.com/yaffelab/visage/
# 
# If you use VISAGE, please cite: 
# VISAGE Reveals a Targetable Mitotic Spindle Vulnerability in Cancer Cells.
# Patterson JC, Joughin BA, Prota AE, Mühlethaler T, Jonas OH, Whitman MA, Varmeh S, Chen S, 
# Balk SP, Steinmetz MO, Lauffenburger DA, Yaffe MB.
# Cell Syst. 2019 
# PMID: XXXXXXXX


# Install and load necessary packages.
#https://cran.r-project.org/web/packages/MASS/index.html
if(!require(MASS)) {
  install.packages("MASS")
  library(MASS)
}

# visage_params.txt is R code to set parameters.
source("visage_params.txt") 

# Create the output directory if it does not already exist.
if (!file.exists(outdir)) {
  dir.create(outdir)
}

# Load custom functions and classes.
source("VISAGE_functions.R")
source("VISAGE_plots.R")

# Log transform dose data.
logdoseA = log10(doseA)
logdoseB = log10(doseB)

# Read in cell line_data and calculate synergy.  These must be files in a subdirectory
# called "cell_line_data" and each file must be of the form "CL_{anything}.txt".  
# Each contains on line 1 the name of the cell line that will be used by the software 
# and matched against expression datasets.  On line 2 is the number of replicates contained.
# Then each line is a row of dose reponse data.
cell_line_files = Sys.glob(paste(inputdir, "/CL_*.txt", sep="")) # A list of all cell line files, only identifies those files that start with "CL_" and end with ".txt"
cell_line_data = vector("list", length(cell_line_files)) # Cell line data structure that holds the results of each cell line. See VISAGE_functions.R for list of attributes.

for (cell_line_index in 1:length(cell_line_files)) { 
    # Read each file in (in VISAGE_functions.R).
    cell_line_data[[cell_line_index]] <- read_cell_line_file(cell_line_files[cell_line_index], length(doseA), length(doseB), logdoseA, logdoseB)
}

# Compile cell line sensitivity and synergy metrics: AUC for each drug, unscaled synergy and scaled/shifted synergy.
# Create & initialize the data structure with data from the first cell line.
compiled_metrics <- data.frame(a = cell_line_data[[1]]@name, b = cell_line_data[[1]]@aucA, c = cell_line_data[[1]]@aucB, d = cell_line_data[[1]]@synergy, e = cell_line_data[[1]]@ss_syn) 
for (metric_index in 2:length(cell_line_data)) {
  # Add data from each further cell line.
  compiled_metrics <- rbind(compiled_metrics, data.frame(a=cell_line_data[[metric_index]]@name, b=cell_line_data[[metric_index]]@aucA, c=cell_line_data[[metric_index]]@aucB, d=cell_line_data[[metric_index]]@synergy, e=cell_line_data[[metric_index]]@ss_syn))
}
# Give rows and columns of the compiled metrics data structure names.
colnames(compiled_metrics) <- c("Cell Line", paste(drugA, "AUC",sep="_"), paste(drugB, "AUC",sep="_"), "untransfomed_synergy", "scaled_shifted_synergy")
rownames(compiled_metrics) <- compiled_metrics[,1]

# Write out cell line sensitivity and synergy for each cell line.
write.matrix(compiled_metrics, file = paste(outdir, "compiled_metrics.txt", sep="/"), sep = "\t") 

# Get a list of cell line names.
cell_line_list <- rownames(compiled_metrics)

# Read in expression data for the cell lines of interest (in VISAGE_functions.R).
expression <- read_expression_data(paste(inputdir, gene_expr_file, sep="/"), cell_line_list)

# Warn the user about missing cell lines.
missing_cell_lines <- setdiff(cell_line_list, colnames(expression))
if(length(missing_cell_lines)>0) {
  warning("Warning: The following cell lines were not found in expression data:\n", paste(missing_cell_lines, "\n"), call.=FALSE)
}

# For each sensitivity and synergy metric, calculate correlation with each gene.
for(metric_index in 2:5) {
  metric_vector <- compiled_metrics[colnames(expression)[2:length(colnames(expression))],metric_index]
  correlation <- cor(t(as.matrix(expression[,2:30])), metric_vector)
  compiled_correlations <- data.frame(Gene = expression$Description, Correlation = correlation)
  compiled_correlations <- compiled_correlations[order(compiled_correlations$Correlation, decreasing=TRUE),]
  
  # Write out sorted correlation vectors.
  write.table(compiled_correlations, file=paste(outdir, "/genes_ranked_by_viability_", colnames(compiled_metrics)[metric_index], "_correlation.rnk", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}

# If make_plots is set to TRUE in visage_params.txt, make a PDF on sensitivity & synergy for each cell line.
# If make_plots is set to FALSE in visage_params.txt, do not.
if(make_plots) {

  # Make an images subdirectory in the output directory.  
  if (!file.exists(paste(outdir, "/images", sep=""))) {
    dir.create(paste(outdir, "/images", sep=""))
  }
  
  # Create a data structure conaining all plots for all cell lines (in VISAGE_plots.R).
  cell_line_plots <- make_all_plots(cell_line_data, logdoseA, logdoseB)

  # Generate PDFs.
  for (cell_line_index in 1:length(cell_line_plots)) {
    
    # For the purposes of spacing, add blank plots after drug A slice plots to make drug B slices start on a new page.
    for (whitespace in 1:(8-(length(cell_line_plots[[cell_line_index]]@drugA_slices) %% 8))) {
      cell_line_plots[[cell_line_index]]@drugA_slices[[length(cell_line_plots[[cell_line_index]]@drugA_slices)+1]] <- ggplot() + theme_void()
    }

    # For the purposes of spacing, add blank plots after drug B slice plots to end at the bottom of a page.
    for (whitespace in 1:(8-(length(cell_line_plots[[cell_line_index]]@drugB_slices) %% 8))) {
      cell_line_plots[[cell_line_index]]@drugB_slices[[length(cell_line_plots[[cell_line_index]]@drugB_slices)+1]] <- ggplot() + theme_void()
    }

    # For the purposes of spacing, intersperse blank plots among surface plots to fill 1 column and start drug A slice plots on a new page.    
    cell_line_plots[[cell_line_index]]@expected[[2]] <- ggplot() + theme_void()
    cell_line_plots[[cell_line_index]]@observed[[2]] <- ggplot() + theme_void()
    cell_line_plots[[cell_line_index]]@diff[[2]] <- ggplot() + theme_void()
    cell_line_plots[[cell_line_index]]@diff[[3]] <- ggplot() + theme_void()
    cell_line_plots[[cell_line_index]]@diff[[4]] <- ggplot() + theme_void()
    
    # Generate a PDF file of: surfaces, drug A slices, and drug B slices.
    ggexport(plotlist=c(cell_line_plots[[cell_line_index]]@expected,
                        cell_line_plots[[cell_line_index]]@observed,
                        cell_line_plots[[cell_line_index]]@diff,
                        cell_line_plots[[cell_line_index]]@drugA_slices, cell_line_plots[[cell_line_index]]@drugB_slices), 
             ncol=2, nrow=4, filename=paste(outdir, "/images/", cell_line_data[[cell_line_index]]@name, ".pdf", sep=""))
  }
}