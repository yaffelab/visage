# Functions for plotting:
# A slice of expected/observed viability
# Expected/observed viability surfaces.

# Install/load necessary packages.
# https://cran.r-project.org/web/packages/ggplot2/index.html
# https://cran.r-project.org/web/packages/ggpubr/index.html
# https://cran.r-project.org/web/packages/reshape2/index.html
# https://cran.r-project.org/web/packages/akima/index.html

if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(ggpubr)) {
  install.packages("ggpubr")
  library(ggpubr)
}
if(!require(reshape2)) {
  install.packages("reshape2")
  library(reshape2)
}
if(!require(akima)){
  install.packages("akima")
  library(akima)
}

# Data structure for holding all plots.
cell_line_plot_set <- setClass("cell_line_plot_set", slots = c(drugA_slices = "list", drugB_slices = "list", 
                                                       expected="list", observed="list", diff="list"))

# Function to generate plots and store in a data structure for printing to PDF by main script.
make_all_plots <- function(data, logdoseA, logdoseB){
  
  # Find maximum and minimums for synergy and maximum viability to set common plot axes.
  max_viab_per_cl <- c()
  max_diff_per_cl <- c()
  min_diff_per_cl <- c()
  # Maximum and minimum synergy and maximum viability for each cell line.
  for (cl in 1:length(data)) {
    max_viab_per_cl[cl] <- max(max(data[[cl]]@avg), max(data[[cl]]@expected))
    max_diff_per_cl[cl] <- max(data[[cl]]@expected-data[[cl]]@avg)
    min_diff_per_cl[cl] <- min(data[[cl]]@expected-data[[cl]]@avg)
  }
  # Global max/min synergy and max viablity.
  max_viab <- max(max_viab_per_cl)
  max_diff <- max(max_diff_per_cl)
  min_diff <- min(min_diff_per_cl)
  
  # Calculate the width of the dose spaces (for getting errorbar widths right).
  rangeA <- logdoseA[length(logdoseA)]-logdoseA[2]
  rangeB <- logdoseB[length(logdoseB)]-logdoseB[2]
  cell_line_plots = vector("list", length(data))

  # Procedure to follow for each cell line.  
  for (i in 1:length(data)) {
    
    # Create a new empty plotset for this cell line, and read in its data.
    this_plotset <- cell_line_plot_set()
    this_cell_line <- data[[i]]
    
    # Assemble a dataset "viability", holding expected, observed, se, and e-o
    # for each dose pair.
    viability <- this_cell_line@avg
    viability <- melt(viability) # Melt turns rectangular data to a list of intersections
    colnames(viability) <- c("drugB", "drugA", "observed")
    viability$drugB <- logdoseB[viability$drugB]
    viability$drugA <- logdoseA[viability$drugA]
    expected <- this_cell_line@expected
    expected <- melt(expected) # Melt turns rectangular data to a list of intersections
    colnames(expected) <- c("drugB", "drugA", "expected")
    expected$drugB <- logdoseB[expected$drugB]
    expected$drugA <- logdoseA[expected$drugA]
    viability <- merge(viability, expected, sort=FALSE) # Merge combines two datasets using common columns.
    se <- this_cell_line@sterr
    se <- melt(se) # Melt turns rectangular data to a list of intersections
    colnames(se) <- c("drugB", "drugA", "se")
    se$drugB <- logdoseB[se$drugB]
    se$drugA <- logdoseA[se$drugA]
    viability <- merge(viability, se, sort=FALSE) # Merge combines two datasets using common columns.
    viability <- cbind(viability, viability$expected-viability$observed) # Add a columnd for difference.
    colnames(viability)[length(colnames(viability))] <- "difference"
    viability$daf <- as.factor(viability$drugA) # Create "factor" versions of the dose columns.
    viability$dbf <- as.factor(viability$drugB)
    
    ##########################################
    # A single plot of drug A effects (observed) at all drug B doses.***
    this_plot <- ggplot()
       
    for (b_dose in 1:length(logdoseB)) {
      this_plot <- this_plot + 
        geom_line(data=viability[viability$drugB==logdoseB[b_dose]&viability$drugA>logdoseA[1],], aes(x=drugA,y=observed, color=dbf),size=.5) +
        geom_point(data=viability[viability$drugB==logdoseB[b_dose]&viability$drugA>logdoseA[1],], aes(x=drugA, y=observed, color=dbf), size=1) +
        geom_errorbar(data=viability[viability$drugB==logdoseB[b_dose]&viability$drugA>logdoseA[1],],aes(x=drugA, ymin=observed-se, ymax=observed+se, color=dbf), size=.5, width=.03*rangeA)
    }
    this_plot <- this_plot + scale_color_brewer(palette = "Spectral", labels=(10^logdoseB[length(logdoseB):1])) +
      xlab(paste("Log [",drugA,"] (M)",sep="")) +
      ylab("Relative Viability") + 
      labs(color=paste(drugB, "Dose"), title=paste(drugA, "Effect by", drugB, "Dose"), subtitle=data[[i]]@name) + 
      guides(color=guide_legend(reverse=TRUE)) +
      theme_light()+
      theme(
        plot.title = element_text(size=8,margin=margin(c(0,0,0,0))),
        legend.text = element_text(size = 7),
        legend.margin=margin(c(0,0,0,0)),
        legend.key.height =  unit(.04, "in"),
        legend.box.margin=margin(c(0,0,0,0)),
        legend.spacing.y=unit(.04, "in"),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=7, face="bold"),
        axis.title.y = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
    coord_cartesian(ylim=c(0,max_viab_per_cl[i]))
    
    this_plotset@drugA_slices[[1]] <- this_plot

    ##########################################
    # A series of plots depicting drug A effect (observed & expected) at each drug B dose.***
    for (b_dose in 2:length(logdoseB)) {
      numtext <- strsplit(as.character(10^logdoseB[b_dose]), "e")
      numtext[[1]][2] <- as.character(as.numeric(numtext[[1]][2]))
      
      #Each drug A dose plot for Drug B effect.
      this_plot <- ggplot(viability[viability$drugB==logdoseB[b_dose]&viability$drugA>logdoseA[1],]) +
        geom_line(aes(x=drugA, y=expected, color="c"), size=.5) +
        geom_point(aes(x=drugA, y=expected, color="c"), size=1) +
        geom_line(data=viability[viability$drugB==logdoseB[1]&viability$drugA>logdoseA[1],], aes(x=drugA, y=observed, color="a"), size=.5) +
        geom_point(data=viability[viability$drugB==logdoseB[1]&viability$drugA>logdoseA[1],], aes(x=drugA, y=observed, color="a"), size=1) +
        geom_errorbar(data=viability[viability$drugB==logdoseB[1]&viability$drugA>logdoseA[1],],aes(x=drugA, ymin=observed-se, ymax=observed+se, color="a"), size=.5, width=.03*rangeA) +
        geom_line(aes(x=drugA, y=observed, color="b"), size=.5) + 
        geom_point(aes(x=drugA, y=observed, color="b"), size=1) + 
        geom_errorbar(aes(x=drugA, ymin=observed-se, ymax=observed+se, color="b"), size=.5, width=.03*rangeA) +
        xlab(paste("Log [",drugA,"] (M)",sep="")) +
        ylab("Relative Viability") + 
        theme_light() +
        theme(
          plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
          legend.text = element_text(size = 7),
          legend.margin=margin(c(0,0,0,-10)),
          legend.box.margin=margin(c(0,0,0,0)),
          legend.spacing.x=unit(0.05, "in"),
          axis.title.x = element_text(size=7, face="bold"),
          axis.title.y = element_text(size=7, face="bold"),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6), 
          plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
        coord_cartesian(ylim=c(0,max_viab_per_cl[i])) +
        labs(title=bquote(.(paste(drugA, "Effect at")) ~ .(numtext[[1]][1]) ~ x ~ 10^.(numtext[[1]][2]) ~ M ~ .(drugB)), subtitle=data[[i]]@name) +
        scale_color_manual(name = NULL, 
                           values =c('a'='black','b'='red', 'c'='gray50'), 
                           labels = c(paste("No", drugB),
                                      bquote(.(numtext[[1]][1]) ~ x ~ 10^.(numtext[[1]][2]) ~ M ~ .(drugB)),
                                      'Expected Combination'))
      
      this_plotset@drugA_slices[[a=b_dose]] <- this_plot
      
    }
    ##########################################
    # A series of plots depicting drug B effect (observed & expected) at each drug A dose.***
    this_plot <- ggplot()
    
    for (a_dose in 1:length(logdoseA)) {
      this_plot <- this_plot + 
        geom_line(data=viability[viability$drugA==logdoseA[a_dose]&viability$drugB>logdoseB[1],], aes(x=drugB,y=observed, color=daf),size=.5) +
        geom_point(data=viability[viability$drugA==logdoseA[a_dose]&viability$drugB>logdoseB[1],], aes(x=drugB, y=observed, color=daf), size=1) +
        geom_errorbar(data=viability[viability$drugA==logdoseA[a_dose]&viability$drugB>logdoseB[1],],aes(x=drugB, ymin=observed-se, ymax=observed+se, color=daf), size=.5, width=.03*rangeB)
    }
    this_plot <- this_plot + scale_color_brewer(palette = "Spectral", labels=(10^logdoseA[length(logdoseA):1])) +
      xlab(paste("Log [",drugB,"] (M)",sep="")) +
      ylab("Relative Viability") + 
      labs(color=paste(drugA, "Dose"), title=paste(drugB, "Effect by", drugA, "Dose"), subtitle=data[[i]]@name) + 
      guides(color=guide_legend(reverse=TRUE)) +
      theme_light()+
      theme(
        plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
        legend.text = element_text(size = 7),
        legend.margin=margin(c(0,0,0,0)),
        legend.key.height =  unit(.04, "in"),
        legend.box.margin=margin(c(0,0,0,0)),
        legend.spacing.y=unit(.04, "in"),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=7, face="bold"),
        axis.title.y = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) + 
      coord_cartesian(ylim=c(0,max_viab_per_cl[i]))
      
    
    this_plotset@drugB_slices[[1]] <- this_plot

    ##########################################
    # A series of plots depicting drug B effect (observed & expected) at each drug A dose.***
    for (a_dose in 2:length(logdoseA)) {
      numtext <- strsplit(as.character(10^logdoseA[a_dose]), "e")
      numtext[[1]][2] <- as.character(as.numeric(numtext[[1]][2]))

      this_plot <- ggplot(viability[viability$drugA==logdoseA[a_dose]&viability$drugB>logdoseB[1],]) +
        geom_line(aes(x=drugB, y=expected, color="c"), size=.5) +
        geom_point(aes(x=drugB, y=expected, color="c"), size=1) +
        geom_line(data=viability[viability$drugA==logdoseA[1]&viability$drugB>logdoseB[1],], aes(x=drugB, y=observed, color="a"), size=.5) +
        geom_point(data=viability[viability$drugA==logdoseA[1]&viability$drugB>logdoseB[1],], aes(x=drugB, y=observed, color="a"), size=1) +
        geom_errorbar(data=viability[viability$drugA==logdoseA[1]&viability$drugB>logdoseB[1],],aes(x=drugB, ymin=observed-se, ymax=observed+se, color="a"), size=.5, width=.03*rangeB) +
        geom_line(aes(x=drugB, y=observed, color="b"), size=.5) + 
        geom_point(aes(x=drugB, y=observed, color="b"), size=1) + 
        geom_errorbar(aes(x=drugB, ymin=observed-se, ymax=observed+se, color="b"), size=.5, width=.03*rangeB) +
        xlab(paste("Log [",drugB,"] (M)",sep="")) +
        ylab("Relative Viability") + 
        theme_light() +
        theme(
          plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
          legend.text = element_text(size = 7),
          legend.margin=margin(c(0,0,0,-10)),
          legend.box.margin=margin(c(0,0,0,0)),
          legend.spacing.x=unit(0.05, "in"),
          axis.title.x = element_text(size=7, face="bold"),
          axis.title.y = element_text(size=7, face="bold"),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
        coord_cartesian(ylim=c(0,max_viab_per_cl[i])) +
        labs(title=bquote(.(paste(drugB, "Effect at")) ~ .(numtext[[1]][1]) ~ x ~ 10^.(numtext[[1]][2]) ~ M ~ .(drugA)), subtitle=data[[i]]@name) +
        scale_color_manual(name = NULL, 
                           values =c('a'='black','b'='red', 'c'='gray50'), 
                           labels = c(paste("No", drugA),
                                      bquote(.(numtext[[1]][1]) ~ x ~ 10^.(numtext[[1]][2]) ~ M ~ .(drugA)),
                                      'Expected Combination'))
      
      this_plotset@drugB_slices[[a_dose]] <- this_plot
      
    }
    
    ##########################################
    # Make surfaces for observed & expected viability and the difference.
    # This requires linear interpolation of the dose-response matrices onto a regular grid.

    # Remove the -Inf data.
    vtrunc <- viability[viability$drugA!=-Inf & viability$drugB !=-Inf,]

    # # # # # # # # # # #
    # Expected Viability
    exp.interp <- interp(x = vtrunc$drugA, y = vtrunc$drugB, z = vtrunc$expected, nx = 50, ny = 50)
    exp.interp.xyz <- as.data.frame(interp2xyz(exp.interp))
    exp_plot <- ggplot(exp.interp.xyz, aes(x = x, y = y, fill = z)) + 
      geom_raster(interpolate=TRUE) + 
      #scale_fill_distiller(palette="Spectral", limits=c(0,1), direction=1, name = "Relative\nViability") +
      scale_fill_distiller(palette="Spectral", limits=c(0,max_viab), direction=1, name = "Relative\nViability") +
      theme_classic() +
      theme(plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
            legend.margin=margin(c(0,0,0,0)),
            legend.box.margin=margin(c(0,0,0,0)),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            axis.title.x = element_text(size=7, face="bold"),
            axis.title.y = element_text(size=7, face="bold"),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
      labs(title="Expected Relative Viability", subtitle=data[[i]]@name) +
      xlab(paste("Log [",drugA,"] (M)", sep="")) + 
      ylab(paste("Log [",drugB,"] (M)", sep=""))

    # # # # # # # # # # #
    # Observed Viability
    obs.interp <- interp(x = vtrunc$drugA, y = vtrunc$drugB, z = vtrunc$observed, nx = 50, ny = 50)
    obs.interp.xyz <- as.data.frame(interp2xyz(obs.interp))
    obs_plot <- ggplot(obs.interp.xyz, aes(x = x, y = y, fill = z)) + 
      geom_raster(interpolate=TRUE) + 
      #scale_fill_distiller(palette="Spectral", limits=c(0,1), direction=1, name = "Relative\nViability") +
      scale_fill_distiller(palette="Spectral", limits=c(0,max_viab), direction=1, name = "Relative\nViability") +
      theme_classic() +
      theme(plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
            legend.margin=margin(c(0,0,0,0)),
            legend.box.margin=margin(c(0,0,0,0)),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            axis.title.x = element_text(size=7, face="bold"),
            axis.title.y = element_text(size=7, face="bold"),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
      labs(title="Observed Relative Viability", subtitle=data[[i]]@name) +
      xlab(paste("Log [",drugA,"] (M)", sep="")) + 
      ylab(paste("Log [",drugB,"] (M)", sep=""))

    # # # # # # # # # # #
    # Difference, Expected minus observed.
    diff.interp <- interp(x = vtrunc$drugA, y = vtrunc$drugB, z = vtrunc$difference, nx = 50, ny = 50)
    diff.interp.xyz <- as.data.frame(interp2xyz(diff.interp))
    diff_plot <- ggplot(diff.interp.xyz, aes(x = x, y = y, fill = z)) + 
      geom_raster(interpolate=TRUE) + 
      scale_fill_gradient2(high="red", mid="white", low="blue", limits=c(min_diff,max_diff), name = "Relative\nViability\nDifference", midpoint=0) +
      theme_classic() +
      theme(plot.title = element_text(size=8, margin=margin(c(0,0,0,0))),
            legend.margin=margin(c(0,0,0,0)),
            legend.box.margin=margin(c(0,0,0,0)),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            axis.title.x = element_text(size=7, face="bold"),
            axis.title.y = element_text(size=7, face="bold"),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            plot.subtitle = element_text(size=5, margin=margin(c(0,0,0,0)))) +
      labs(title="Difference (Exp-Obs) in Relative Viability", subtitle=data[[i]]@name) +
      xlab(paste("Log [",drugA,"] (M)", sep="")) + 
      ylab(paste("Log [",drugB,"] (M)", sep=""))
    
    this_plotset@expected[[1]] <- exp_plot
    this_plotset@observed[[1]] <- obs_plot
    this_plotset@diff[[1]] <- diff_plot

    # Store the whole set of plots for this cell line.    
    cell_line_plots[[i]] <- this_plotset
  }

  return(cell_line_plots)
}




