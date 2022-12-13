############################################################################
#   Title:  Epicapture hybridization capture stats
#   Description: Script used to plot target coverage stats (R)
#   Author: Miljana Tanic (m.tanic@ucl.ac.uk) 
#   Created: September 2019
#   Last edited: 18-04-2020 # barplots on coverage stats
###############################################################################

library(tidyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(beanplot)
library(RColorBrewer)
library("devtools")
# install_github("ndphillips/yarrr")
library("yarrr")
library(ggplot2)
library(wesanderson) # works in my R session only 





#===================================
# Set working directory
#==================================

  Epicapture <- "~/RDS_C2c/EpiCapture/"
  setwd(Epicapture)

  RESULTS <- paste0(Epicapture,"RESULTS/")
  list.files(paste0(RESULTS, "3_TargetCaptureMetrics/"))


#========================================================
# Import hybcapture data output from Picard
#========================================================


  # Import metric on target capture efficiencies for 3 hybridization platforms: Target enrichment stats sumarry_Final dataset.xlsx, sheet4
    # hyb_capture.df <- read.delim("clipboard", header = TRUE, stringsAsFactors = TRUE )
    hyb_capture.df <- read.delim(paste0(RESULTS, "3_TargetCaptureMetrics/Target enricment stats_4R.txt"), header=TRUE, stringsAsFactors = TRUE)


    dim(hyb_capture.df)
    str(hyb_capture.df)

    save(hyb_capture.df, file=paste0(RESULTS,"3_TargetCaptureMetrics/hyb_capture.RData"))

    load(paste0(RESULTS,"3_TargetCaptureMetrics/hyb_capture.RData"))

#========================================================
# Plot coverage by sequencing depth per sample for each platform 
#========================================================

    # Arrange for 3 plots - 1 per row/ 3 per row
      par(mfrow=c(1,1))

    # plot for each Platform
      names(hyb_capture.df)
      Agilent="Agilent \n SureSelect MethylSeq"
      Roche= "Roche \n SeqCap EpiGiant"
      Illumina="Illumina \n TruSeq Methyl Cature EPIC"

      head(hyb_capture.df[,c(1:4,21:28)])
      pct_cov.df <- hyb_capture.df[,c("PLATFORM_SAMPLE","Platform", "SAMPLE","PF_READS", "PCT_TARGET_BASES_1X", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X", "PCT_TARGET_BASES_50X", "PCT_TARGET_BASES_100X")]

      head(pct_cov.df)
      library(dplyr)

      x <- c(1,2,10,20,30,40,50,100)


    #-------------------
    # Agilent:
    #-------------------

      pct_cov.df %>% filter(Platform=="Agilent") -> y

      rownames(y) <- y[,3]
      colnames(y)
      y <- y[,5:12]

      # set color pallete:

      names(wes_palettes)
      col = wes_palette("Zissou1", 16, type = "continuous")

      # plot the first line:
      plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of capture targets >= Depth", xlab = "Depth of coverege", main = Agilent)

      # plot x-axis
      depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
      axis(side = 1, at = x, labels = depth, tcl = -0.2)
      abline(v=x, col="gray", lty=2)
      abline(h=0.8, col="red", lty=1)
      abline(h=0.6, col="blue", lty=1)

      # Add a second line
      for (i in 2:16){
        lines(x, y[i,], type = "l", lty=1, col= col[i])
      }

      # Add a legend to the plot
      legend("topright", legend=rownames(y),
            col=col, lty = 1, cex=0.7)

    #----------------
    # Roche:
    #----------------

      pct_cov.df %>% filter(Platform=="Roche") -> y

      rownames(y) <- y[,3]
      colnames(y)
      y <- y[,5:12]

      # set color pallete:
      colors <- RColorBrewer::brewer.pal(16,)

      names(wes_palettes)
      col = wes_palette("Zissou1", 16, type = "continuous")

      # plot the first line:
      plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of targets", xlab = "Depth of coverege", main = Roche)

      # plot x-axis
      depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
      axis(side = 1, at = x, labels = depth, tcl = -0.2)
      abline(v=x, col="gray", lty=2)
      abline(h=0.8, col="red", lty=1)
      abline(h=0.6, col="blue", lty=1)

      # Add a second line
      for (i in 2:16){
        lines(x, y[i,], type = "l", lty=1, col= col[i])
      }

      # Add a legend to the plot
      legend("topright", legend=rownames(y),
            col=col, lty = 1, cex=0.7)


    #------------
    # Illumina:
    #------------

      pct_cov.df %>% filter(Platform=="Illumina") -> y

      rownames(y) <- y[,3]
      colnames(y)
      y <- y[,5:12]

      # set color pallete:
      colors <- RColorBrewer::brewer.pal(16,)

      names(wes_palettes)
      col = wes_palette("Zissou1", 16, type = "continuous")

      # plot the first line:
      plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of targets", xlab = "Depth of coverege", main = Illumina)

      # plot x-axis
      depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
      axis(side = 1, at = x, labels = depth, tcl = -0.2)
      abline(v=x, col="gray", lty=2)
      abline(h=0.8, col="red", lty=1)
      abline(h=0.6, col="blue", lty=1)

      # Add a second line
      for (i in 2:16){
        lines(x, y[i,], type = "l", lty=1, col= col[i])
      }

      # Add a legend to the plot
      legend("topright", legend=rownames(y),
            col=col, lty = 1, cex=0.7)

#========================================================
# Plot agregated hybridization stats for each platform 
#========================================================

  #------------------------------------------
  # make new dataset and calculate summaries
  #------------------------------------------

    # hyb_summary <- data.frame(total_unq_mapped=numeric(), pct_on_target=numeric(), pct_off_target=numeric(), pct_on_near_bait_vs_selected=numeric(), mean_bait_cov=numeric(), median_target_cov=numeric(), median_target_cov=numeric(), pct_bases_onbait=numeric(), pct_bases_ontarget=numeric(), fold_enrichment=numeric(), unifirmity=numeric())

    # conveert to tibble
    hyb_capture <- as_tibble(hyb_capture.df)
    names(hyb_capture)

    # group by platform
    by_platform <- hyb_capture %>% 
                    group_by(Platform)

    groups(by_platform)
    group_vars(by_platform)

    # calculate mean and median for one columns
    hyb_capture %>% 
                group_by(Platform) %>%
                summarise(total_unq_mapped.mean=mean(PF_UQ_READS_ALIGNED),
                          total_unq_mapped.sd=sd(PF_UQ_READS_ALIGNED))

    # calculate mean and median for all numeric columns
    hyb_summary <- hyb_capture %>% 
                    group_by(Platform) %>%
                    summarise_if(is.numeric, list(~mean(.), ~sd(.)), na.rm = TRUE)


  #-----------------
  # Set colors
  #----------------


    # Selecting colors using yarr (pirateplot)
    piratepal(palette= "all")
    piratepal("google") 
    # blue         red      yellow       green 
    # "#3D79F3FF" "#E6352FFF" "#F9B90AFF" "#34A74BFF" 
    col.platforms <- c("#3D79F3FF",  "#E6352FFF", "#34A74BFF", "#7570b3" , "#F9B90AFF") # GOOD LOKING DIVERGING PALLETE

    # make color plaette transparent uing yarr transparent() function:
    col.platforms <- transparent(orig.col = col.platforms, trans.val = 0.3) # BEAUTIFUL :)

    # check
    barplot(x, col= col.platforms2)
    legend(x = 'topleft', legend=col.platforms, fill=col.platforms, cex = 0.8)


  #------------------------------------------
  # plot hybridization efficiency 
  #------------------------------------------

    # make dataset for plotting and set levels for Platform
    data <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), pct_on_target.mean=hyb_summary$PCT_SELECTED_BASES_mean, pct_on_target.sd=hyb_summary$PCT_SELECTED_BASES_sd)

    # plot barchart with error bars - Standard deviation
      # single color
      ggplot(data) +
        geom_bar( aes(x=platform, y=pct_on_target.mean), stat="identity", fill="forestgreen", alpha=0.5) +
        geom_errorbar( aes(x=platform, ymin=pct_on_target.mean-pct_on_target.sd, ymax=pct_on_target.mean+pct_on_target.sd), width=0.4, colour="orange", alpha=0.9, size=1) + # add error bars
        ggtitle("hybridization efficiency") +
        xlab("platform") +
        ylab("Fraction of reads on/near bait") +
        theme_classic()

      # by platform colors
      ggplot(data) +
        geom_bar( aes(x=platform, y=pct_on_target.mean, fill=platform), stat="identity", alpha=0.5) +
        geom_errorbar( aes(x=platform, ymin=pct_on_target.mean-pct_on_target.sd, ymax=pct_on_target.mean+pct_on_target.sd), width=0.4, colour="gray56", alpha=0.9, size=0.8) + # add error bars
        scale_fill_manual(values=col.platforms[1:3]) + 
        scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1), labels=c("20", "40", "60", "80", "100")) +
        ggtitle("hybridization efficiency") +
        xlab("platform") +
        ylab("Reads (percent)") +
        theme_classic()

  #------------------------------------------
  # plot hybridization efficiency 
  #------------------------------------------

    # make dataset for plotting and set levels for Platform
    data.mean <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
          pct_on_target.mean=hyb_summary$PCT_SELECTED_BASES_mean, 
          pct_off_target.mean=hyb_summary$PCT_OFF_BAIT_mean)

    data.sd <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
          pct_on_target.sd=hyb_summary$PCT_SELECTED_BASES_sd, 
          pct_off_target.sd=hyb_summary$PCT_OFF_BAIT_sd)


    # transform to long format
    data.mean_long <- gather(data.mean, key = "mean", value = "mean.pct", c("pct_on_target.mean", "pct_off_target.mean"), na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data.sd_long <- gather(data.sd, key = "sd", value = "sd.pct", c("pct_on_target.sd", "pct_off_target.sd"), na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data_long <- cbind(data.mean_long, data.sd_long[,2:3])
    data_long

    # plot barchart with error bars - Standard deviation
      cols <- c("forestgreen","tomato")

      ggplot(data_long, aes(x=platform, y=mean.pct, fill=mean)) +
        geom_bar(stat="identity", position=position_dodge(),  alpha=0.5) +
        geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
        scale_fill_manual(values=cols, name="", labels=c("on target", "off target")) + 
        scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=c("10", "20","30", "40", "50", "60", "70", "80","90", "100")) +
        ggtitle("hybridization efficiency") +
        xlab("") +
        ylab("Reads (percent)") +
        theme_classic()

  #--------------------------------------------------------
  # plot mean target coverage, Fold Enrichment, Uniformity 
  #--------------------------------------------------------

    # make dataset for plotting and set levels for Platform
    data.mean <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
          mean_coverage.mean=hyb_summary$MEAN_TARGET_COVERAGE_mean, 
          fold_enrichment.mean=hyb_summary$FOLD_ENRICHMENT_mean,
          uniformity.mean=hyb_summary$FOLD_80_BASE_PENALTY_mean)

    data.sd <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
          mean_coverage.sd=hyb_summary$MEAN_TARGET_COVERAGE_sd, 
          fold_enrichment.sd=hyb_summary$FOLD_ENRICHMENT_sd,
          uniformity.sd=hyb_summary$FOLD_80_BASE_PENALTY_sd)

    # transform to long format
    data.mean_long <- gather(data.mean, key = "mean", value = "mean.pct", c("mean_coverage.mean", "fold_enrichment.mean", "uniformity.mean"), na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data.sd_long <- gather(data.sd, key = "sd", value = "sd.pct", c("mean_coverage.sd", "fold_enrichment.sd", "uniformity.sd"), na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data_long <- cbind(data.mean_long, data.sd_long[,2:3])
    data_long


    # plot barchart with error bars - Standard deviation
      cols <- c("lightsteelblue3","thistle3", "snow3")

      # ggplot(data_long, aes(x=platform, y=mean.pct, fill=mean)) +
      #   geom_bar(stat="identity", position=position_dodge(),  alpha=0.5) +
      #   geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
      #   # scale_fill_manual(values=cols, name="", labels=c("on target", "off target")) + 
      #   # scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=c("10", "20","30", "40", "50", "60", "70", "80","90", "100")) +
      #   # ggtitle("hybridization efficiency") +
      #   xlab("") +
      #   ylab("Reads (percent)") +
      #   theme_classic()

      ggplot(data_long, aes(x=platform, y=mean.pct, fill=mean)) +
        facet_grid( . ~ mean) +
        geom_bar(stat="identity", position=position_dodge(),  alpha=0.9) +
        geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
        scale_fill_manual(values=cols, name="", labels=c("Mean Coverage", "Fold Enrichment", "Uniformity")) + 
        # scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=c("10", "20","30", "40", "50", "60", "70", "80","90", "100")) +
        # ggtitle("hybridization efficiency") +
        xlab("") +
        ylab("") +
        theme_classic()


  #--------------------------------------------------------
  # plot mean target coverage
  #--------------------------------------------------------

    data_long %>%
          filter(mean=="mean_coverage.mean") %>%
          ggplot(aes(x=platform, y=mean.pct, fill=mean)) +
            geom_bar( stat="identity", fill="lightsteelblue3", alpha=0.9) +
            geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
            scale_y_continuous(breaks=c(10, 20, 30, 40, 50), labels=c("10x", "20x", "30x", "40x", "50x")) +
            ggtitle("Mean coverage") +
            xlab("") +
            ylab("Coverage") +
            theme_classic()

  #--------------------------------------------------------
  # plot  fold enrichment
  #--------------------------------------------------------

    data_long %>%
          filter(mean=="fold_enrichment.mean") %>%
          ggplot(aes(x=platform, y=mean.pct, fill=mean)) +
            geom_bar( stat="identity", fill="thistle3", alpha=0.9) +
            geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
            scale_y_continuous(breaks=c(5, 10, 15, 20, 25, 30)) +
            ggtitle("Fold Enrichment") +
            xlab("") +
            ylab("Fold") +
            theme_classic()

  #--------------------------------------------------------
  # plot uniformity
  #--------------------------------------------------------

    data_long %>%
          filter(mean=="uniformity.mean") %>%
          ggplot(aes(x=platform, y=mean.pct, fill=mean)) +
            geom_bar( stat="identity", fill="snow3", alpha=0.9) +
            geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
            scale_y_continuous(breaks=c(1, 2, 3, 4)) +
            ggtitle("Uniformity") +
            xlab("") +
            ylab("") +
            theme_classic()

  #------------------------------------------
  # plot average coverage per platform
  #------------------------------------------

    # make dataset for plotting and set levels for Platform
    data <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
            cov_1x.mean=hyb_summary$PCT_TARGET_BASES_1X_mean, cov_1x.sd=hyb_summary$PCT_TARGET_BASES_1X_sd,
            cov_2x.mean=hyb_summary$PCT_TARGET_BASES_2X_mean, cov_2x.sd=hyb_summary$PCT_TARGET_BASES_2X_sd,
            cov_10x.mean=hyb_summary$PCT_TARGET_BASES_10X_mean, cov_10x.sd=hyb_summary$PCT_TARGET_BASES_10X_sd,
            cov_20x.mean=hyb_summary$PCT_TARGET_BASES_20X_mean, cov_20x.sd=hyb_summary$PCT_TARGET_BASES_20X_sd,
            cov_30x.mean=hyb_summary$PCT_TARGET_BASES_30X_mean, cov_30x.sd=hyb_summary$PCT_TARGET_BASES_30X_sd,
            cov_40x.mean=hyb_summary$PCT_TARGET_BASES_40X_mean, cov_40x.sd=hyb_summary$PCT_TARGET_BASES_40X_sd,
            cov_50x.mean=hyb_summary$PCT_TARGET_BASES_50X_mean, cov_50x.sd=hyb_summary$PCT_TARGET_BASES_50X_sd,
            cov_100x.mean=hyb_summary$PCT_TARGET_BASES_100X_mean, cov_100x.sd=hyb_summary$PCT_TARGET_BASES_100X_sd)
    

    # make dataset for plotting and set levels for Platform
    data.mean <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
            cov_1x.mean=hyb_summary$PCT_TARGET_BASES_1X_mean, 
            cov_2x.mean=hyb_summary$PCT_TARGET_BASES_2X_mean, 
            cov_10x.mean=hyb_summary$PCT_TARGET_BASES_10X_mean, 
            cov_20x.mean=hyb_summary$PCT_TARGET_BASES_20X_mean, 
            cov_30x.mean=hyb_summary$PCT_TARGET_BASES_30X_mean, 
            cov_40x.mean=hyb_summary$PCT_TARGET_BASES_40X_mean, 
            cov_50x.mean=hyb_summary$PCT_TARGET_BASES_50X_mean,
            cov_100x.mean=hyb_summary$PCT_TARGET_BASES_100X_mean)

    data.sd <- data.frame(platform=factor(hyb_summary$Platform, levels=c("Agilent", "Illumina", "Roche")), 
            cov_1x.sd=hyb_summary$PCT_TARGET_BASES_1X_sd,
            cov_2x.sd=hyb_summary$PCT_TARGET_BASES_2X_sd,
            cov_10x.sd=hyb_summary$PCT_TARGET_BASES_10X_sd,
            cov_20x.sd=hyb_summary$PCT_TARGET_BASES_20X_sd,
            cov_30x.sd=hyb_summary$PCT_TARGET_BASES_30X_sd,
            cov_40x.sd=hyb_summary$PCT_TARGET_BASES_40X_sd,
            cov_50x.sd=hyb_summary$PCT_TARGET_BASES_50X_sd,
            cov_100x.sd=hyb_summary$PCT_TARGET_BASES_100X_sd)

    # transform to long format
    data.mean_long <- gather(data.mean, key = "mean", value = "mean.pct", 2:9, na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data.sd_long <- gather(data.sd, key = "sd", value = "sd.pct", 2:9, na.rm = FALSE,  convert = FALSE, factor_key = TRUE)

    data_long <- cbind(data.mean_long, data.sd_long[,2:3])
    data_long

    # plot by platform colors


    # New facet label names for supp variable
    mean.labs <- c( "1x", "2x", "10x", "20x", "30x", "40x", "50x", "100x")
    names(mean.labs) <- c("cov_1x.mean", "cov_2x.mean", "cov_10x.mean", "cov_20x.mean","cov_30x.mean","cov_40x.mean", "cov_50x.mean", "cov_100x.mean")

    # Create the plot
    library(scales)
          ggplot(data_long, aes(x=platform, y=mean.pct, fill=platform)) +
            facet_grid(. ~ mean, labeller =labeller(mean=mean.labs)) +
            geom_bar(stat="identity", position=position_dodge(),  alpha=0.5) +
            geom_errorbar( aes(x=platform, ymin=mean.pct-sd.pct, ymax=mean.pct+sd.pct), position=position_dodge(0.9), width=0.4, colour="gray21", alpha=0.9, size=0.5) + # add error bars
            scale_fill_manual(values=col.platforms[1:3]) + 
            # scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=c("10", "20","30", "40", "50", "60", "70", "80","90", "100")) +
            scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=percent) +
            ggtitle("Percent of Targets covered at specific depth") +
            xlab("") +
            ylab("percent") +
            theme_classic() +        
            theme(axis.text.x = element_blank())



        geom_pirate(aes(colour = Platform, fill=Platform), points = FALSE, bars = TRUE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                    points_params = list(shape = 19, alpha = 0.2),
                    lines_params = list(size = 0.8))
