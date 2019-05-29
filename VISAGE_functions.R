# Define a class for containing cell line data.
cell_line <- setClass("cell_line", slots = c(name="character", num_reps="numeric", raw_data="list", norm_data="list", 
                                             avg="matrix", sterr="matrix", max_kill="numeric",
                                             scale_shift_data="matrix", aucA="numeric", aucB="numeric",
                                             expected="matrix", ss_expected="matrix", synergy="numeric", ss_syn="numeric"))

# Function for reading a single cell line file, doing the necessary data transformations and calculating the relevent statistics
read_cell_line_file <- function(cell_line_filename, numdoses_A, numdoses_B, logdoses_A, logdoses_B) {

  # Create an "empty" cell line data structure.
  result <- cell_line() 
  
  # Read the contents of the cell line file for parsing.
  cl_file <- readLines(cell_line_filename)
  
  # Read the cell line name.
  result@name <- cl_file[1] 
  cl_file <- cl_file[-1]

  # Read the number of replicates.  
  result@num_reps <- as.numeric(cl_file[1]) 
  cl_file <- cl_file[-1]
  
  # Read the viability data itself, expecting repeated replicates of drugA doses in columns
  # and drugB doses in rows.
  for (replicate_index in 1:result@num_reps) { 
    result@raw_data[[replicate_index]] = matrix(nrow=numdoses_B, ncol=numdoses_A)
    for (doseB_index in 1:numdoses_B) {
        result@raw_data[[replicate_index]][doseB_index,] <- as.numeric(c(strsplit(cl_file[1], '\t')[[1]]))
        cl_file<-cl_file[-1]
    }

    # Normalize each replicate to its own zero drug control.
    result@norm_data[[replicate_index]] = result@raw_data[[replicate_index]] / result@raw_data[[replicate_index]][1,1] # Normalize raw data
  }

  # Average viability across replicates
  result@avg <- apply(simplify2array(result@norm_data), 1:2, mean)
  
  # Standard error of viability across replicates.
  result@sterr <- apply(simplify2array(result@norm_data), 1:2, sd) / sqrt(result@num_reps)

  # Max kill, any dose combo, in this cell line.
  result@max_kill <- min(result@avg)
  
  # Transformation of viability data such that no drug = 1, and most potent dose combo = 0.
  result@scale_shift_data <- result@avg * (1/(1-result@max_kill)) + (1-(1/(1-result@max_kill)))

  # Calculate drug AUC values as a fraction of total possible killing.

  result@aucA <- 0
  result@aucB <- 0
  drugA_auc <- 0
  drugA_totalarea <- 0
  # Sum up raw AUC of drug A and raw total area.
  for(doseA_index in 2:(numdoses_A-1)) {
    drugA_auc <- drugA_auc + (logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*mean(result@avg[1,doseA_index:(doseA_index+1)])
    drugA_totalarea <- drugA_totalarea + (logdoses_A[doseA_index+1]-logdoses_A[doseA_index])
  }
  # Calculate AUC of drug A as raw AUC over total area.
  result@aucA <- drugA_auc/drugA_totalarea
  
  drugB_auc <- 0
  drugB_totalarea <- 0
  # Sum up raw AUC of drug B and raw total area.
  for(doseB_index in 2:(numdoses_B-1)) {
    drugB_auc <- drugB_auc + (logdoses_B[doseB_index+1]-logdoses_B[doseB_index])*mean(result@avg[doseB_index:(doseB_index+1),1])
    drugB_totalarea <- drugB_totalarea + (logdoses_B[doseB_index+1]-logdoses_B[doseB_index])
  }
  # Calculate AUC of drug B as raw AUC over total area.
  result@aucB <- drugB_auc/drugB_totalarea
  
  # Expected surfaces - unscaled and scaled&shifted.
  result@expected <- result@avg
  result@ss_expected <- result@scale_shift_data
  # Expected normalized viability is product of observed single drug doses.
  for (doseB_index in 2:numdoses_B) {
    result@expected[doseB_index,] <- result@expected[1,] * result@expected[doseB_index,1]
    result@ss_expected[doseB_index,] <- result@ss_expected[1,] * result@ss_expected[doseB_index,1]
  }
  
  # Synergy calculations.
  max_vol <- 0
  e_vol <- 0
  o_vol <- 0
  # Sum up total volume, and volumes under expected and observed surfaces for the unscaled data.
  for (doseA_index in 2:(numdoses_A-1)) {
    for (doseB_index in 2:(numdoses_B-1)) {
      max_vol <- max_vol + ((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index]))
      e_vol <- e_vol + (((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index])) * (mean(result@expected[doseB_index:(doseB_index+1),doseA_index:(doseA_index+1)])))
      o_vol <- o_vol + (((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index])) * (mean(result@avg[doseB_index:(doseB_index+1),doseA_index:(doseA_index+1)])))
    }
  }
  # Synergy is difference between expected and observed volumes over total volume.
  result@synergy <- (e_vol-o_vol)/max_vol

  max_ss_vol <- 0
  e_ss_vol <- 0
  o_ss_vol <- 0
  # Sum up total volume, and volumes under expected and observed surfaces for the scaled/shifted data.
  for (doseA_index in 2:(numdoses_A-1)) {
    for (doseB_index in 2:(numdoses_B-1)) {
      max_ss_vol <- max_ss_vol + ((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index]))
      e_ss_vol <- e_ss_vol + (((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index])) * (mean(result@ss_expected[doseB_index:(doseB_index+1),doseA_index:(doseA_index+1)])))
      o_ss_vol <- o_ss_vol + (((logdoses_A[doseA_index+1]-logdoses_A[doseA_index])*(logdoses_B[doseB_index+1]-logdoses_B[doseB_index])) * (mean(result@scale_shift_data[doseB_index:(doseB_index+1),doseA_index:(doseA_index+1)])))
    }
  }
  # Synergy is difference between expected and observed volumes over total volume.
  result@ss_syn <- (e_ss_vol-o_ss_vol)/max_ss_vol
  
  return(result)
}

# Function for reading a gene expression file and keeping relevant columns.
read_expression_data <- function(exp_filename, cl_names) {
  # Install/load a necessary package.
  # https://cran.r-project.org/web/packages/readr/index.html
  if(!require(readr)) {
    install.packages("readr")
    library("readr")
  }
  # Read in the raw expression data
  expression <- read_delim(exp_filename, "\t", escape_double = FALSE, trim_ws = TRUE)

  # Remove each column that is not titled exactly as a cell line under investigation.
  for (column_index in length(colnames(expression)):1){
    if(column_index != 1 & !(colnames(expression)[column_index] %in% cl_names)){
      expression<-expression[,-column_index]
    }
  }
  
  # Remove unnamed rows and rows that have a number in the genename column with no other data.
  for (row_index in nrow(expression):1) {
    # If there is no gene name,            or the gene name is a number (suppress warnings here              and there's at least one following blank.
    if ( is.na(expression[row_index,1]) | ((!is.na(suppressWarnings(as.numeric(expression[row_index,1])))) & sum(is.na(expression[row_index,])) >=1)) {
      expression<-expression[-row_index,]
    }
  }

    return(expression)
}


