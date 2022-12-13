#################################################################
#                                                               #
# Project: Epicapture - Targeted BS-seq platform comparison     #
# Description: Number of unique CpGs covered by each platform   #
# Author: Miljana Tanic                                         #
# Date last edited: 05-03-2019                                  #
#                                                               #
#################################################################

getwd()

library(readr)
library(data.table)

#===================================
# Working directory
#==================================

  Epicapture <- "~/RDS_C2c/EpiCapture/"
  setwd(Epicapture)

  RESULTS <- paste0(Epicapture,"RESULTS/")
  list.files(RESULTS)


#--------------------------------------
# Path to each platform
#--------------------------------------

  path2Agilent <- paste0(Epicapture, "Agilent_SureSelect/bedGraphs_Coverage")
  path2Illumina <- paste0(Epicapture, "Illumina_EPIC/bedGraphs_Coverage")
  path2Roche <-  paste0(Epicapture, "Roche_Nimblegen/bedGraphs_Coverage")
  path2Diagenode <-  paste0(Epicapture, "Diagenode_RRBS/bedGraphs_Coverage")
  path2Nugen <-  paste0(Epicapture, "NuGen_RRBS/bedGraphs_Coverage")


#--------------------------------------
# Sample list
#--------------------------------------

  # samples <- list.files(pattern="bismark.cov") # Total number of CpG sites for both strands independently

  # Total number of CpG sites - output from coverage2cytosine with information from top and bottom strand merged into one
  samples.agilent <- list.files(path=path2Agilent, pattern="CpG_evidence.cov") 
  samples.roche <- list.files(path=path2Roche, pattern="CpG_evidence.cov") 
  samples.illumina <- list.files(path=path2Illumina, pattern="CpG_evidence.cov") 
  samples.diagenode <- list.files(path=path2Diagenode, pattern="CpG_evidence.cov") 
  samples.nugen <- list.files(path=path2Nugen, pattern="CpG_evidence.cov") 

#---------------------------------------
# # Test sample reading of coverage file 
#---------------------------------------
  # sample <- fread(samples[1])
  # class(sample)
  # head(sample)
  # names(sample) <- c("chrom", "start", "end", "%_methylation", "count_M", "count_UM")
  # head(sample)
  # # test$col1 <- paste0("chr", test$col1)
  # sample$chrom <- paste0("chr", sample$chrom)
  # head(sample)
  # # Total number of reads per CpG site:
  # sample$sum <- sample$count_M + sample$count_UM
  # head(sample)
  # dim(sample)
  # dim(sample[sum>=5])


#--------------------------------------------------------------------
# Make function to read all sample CpG coverage files to a data.frame:
#--------------------------------------------------------------------

countCpGs <- function(samples, df){
  
  df<- data.frame(Total_CpGs=integer(),
                  MeanDepthPerCpG=double(),
                  NoCpGs_5x=integer(),
                  PercentCpG_5x=double(),
                  NoCpGs_10x=integer(),
                  PercentCpG_10x=double(),
                  NoCpGs_30x=integer(),
                  PercentCpG_30x=double())
  
  for (i in samples) {
    # Read coverage file:
    sample <- fread(i)
    names(sample) <- c("chrom", "start", "end", "%_methylation", "count_M", "count_UM")
    
    # test$col1 <- paste0("chr", test$col1)
    sample$chrom <- paste0("chr", sample$chrom)
    
    # Total number of reads per CpG site:
    sample$sum <- sample$count_M + sample$count_UM
    
    # Total number of CpGs (covered at least once):
    df[i,1] <- dim(sample)[1]
    
    # Mean readDepth per CpG
    df[i,2] <- mean(sample$sum)
    
    # No of CpGs with Reads >=5x
    
    sample.5x <- sample[sample$sum>=5,]
    df[i,3] <-   dim(sample.5x)[1]
    # Percent CpGs with >5x read depth
    df[i,4] <- dim(sample.5x)[1]/dim(sample)[1]
    
    # No of CpGs with Reads >=10x  
    sample.10x <- sample[sample$sum>=10,]
    df[i,5] <- dim(sample.10x)[1]
    # Percent CpGs with >=10x read depth
    df[i,6] <- dim(sample.10x)[1]/dim(sample)[1]
    
    # No of CpGs with Reads >=30x
    sample.30x <-  sample[sample$sum>=30,]
    df[i,7] <- dim(sample.30x)[1]
    # Percent CpGs with >=30x read depth
    df[i,8] <-  dim(sample.30x)[1]/dim(sample)[1]
  }  
  
  head(df)
  return(df)
}

#---------------------------------
# Read samples from each platform:
#---------------------------------

    setwd(path2Agilent)
    df.agilent <- countCpGs(samples.agilent)

    setwd(path2Roche)
    df.roche <- countCpGs(samples.roche)

    setwd(path2Illumina)
    df.illumina <- countCpGs(samples.illumina)

    setwd(path2Diagenode)
    df.diagenode <- countCpGs(samples.diagenode)

    setwd(path2Nugen)
    df.nugen <- countCpGs(samples.nugen)


  # Save RData objects:
    setwd(Epicapture)
    save(df.agilent, df.roche, df.illumina, df.diagenode, df.nugen, file=paste0(RESULTS,"2_SequencingSummaryMetrics/CpG count.RData"))

  # Merge dataframes to one:

    df <- rbind(df.agilent, df.roche, df.illumina, df.diagenode, df.nugen)

  ## Rename samples - remove file extention:
    rownames(df) <-gsub("_R1_val_1_bismark_bt2_pe.bismark.CpG_report.merged_CpG_evidence.cov", "",rownames(df) )
    rownames(df) <-gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.CpG_report.merged_CpG_evidence.cov", "",rownames(df) )
    rownames(df) <-gsub("_R1_val_1_trimmed_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.CpG_report.merged_CpG_evidence.cov", "",rownames(df) )
    rownames(df) <-gsub("_R1_val_1_trimmed.fq.gz_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.CpG_report.merged_CpG_evidence.cov", "",rownames(df) )
    rownames(df) <-gsub("_100M", "",rownames(df) )


    rownames(df) # check if ok

  # Split sample names and platform name
    df$sample <- rownames(df)

    library(dplyr)
    library(tidyr)

  df <- df %>% separate(sample, c("Platform", "Sample"), "_") #%>% Passes object on left hand side as first argument (or .argument) of function on righthand side
  head(df)

  # Define platforms order:
    Platforms <- c("Agilent", "Roche", "Illumina", "Diagenode", "Nugen")

  # Define factors and level order:
    df$Sample <- as.factor(df$Sample)
    df$Platform <- factor(df$Platform, levels = Platforms)
    df

    Platforms%in%levels(df$Platform) # platforms present in dataframe

#----------------------------------------------------
# Save csv file and Rdata objects for future analysis
#----------------------------------------------------

  setwd(Epicapture)
  #write.table(df, "Coverage_summary_byStrand.txt", sep="\t")
  write.table(df,file=paste0(RESULTS, "/2_SequencingSummaryMetrics/Coverage_summary_mrgCpG.txt", sep="\t"))

  # Save all objects
  save.image(file= "CpG_coverage.RData")

  # # Save some, but not all objects:
  # save(df.illumina, data, file="Illumina_CpGcoverage_stats.RData")
  # 
  # # these can be reloaded into another R session using the load() function (w/ SAME NAMES):
  # load(file="Diagenode_CpGcoverage_stats.RData") 
  # 
  # #  Save a single object you can use the saveRDS() function:
  # 
  # saveRDS(df, file="Diagenode_CpG_coverage_df.rds")
  # # You can load these into your R session using the readRDS() function, but you will need to assign the result into a the desired variable:
  # 
  # df.Roche <- readRDS("Roche_CpG_coverage_df.rds")
  # data.Roche <- readRDS("Roche_CpG_coverage_data.rds")


#===========================================
# Statistics and plotting on merged dataset:
#===========================================

  dplyr::tbl_df(df) # Converts data to tbl class - easier to examine than data frames

  glimpse(df) # Information dense summary of tbl data
  View(df) # View data set in spreadsheet-like display 
  levels(df$Platform)

  # Summarise data  into single row of values - mean of all platforms
    summarise_all(df[,1:8], funs(mean)) 


  # Caluclate mean number of CpGs and stanadrd deviation BY PLATFORM:
    mean <- tapply(df[,1], df$Platform, mean)
    mean <- mean[c(1,3,4,2)] # reorder
    sd <- tapply(df[,1], df$Platform, sd)
    sd <- sd[c(1,3,4,2)]

  # Calculate mean for each column by platform: 
    apply(df[,1:8], 2, tapply, df$Platform, mean) # combine tapply with apply to allow us to calculate clustered values from several different columns of a data frame

  # Calculate standard deviation for each column by platform: 
    apply(df[,1:8], 2, tapply, df$Platform, sd) # combine tapply with apply to allow us to calculate clustered values from several different columns of a data frame



  # Subset CpG counts to dataframe:
    CpG.counts <- df[,c(1,3,5,7,9:10)]
    # reorder by Platform

  # Calculate mean & sd for each platform
    CpG.counts.mean <- apply(CpG.counts[1:4], 2, tapply, CpG.counts$Platform, mean)
    CpG.counts.sd <- apply(CpG.counts[1:4], 2, tapply, CpG.counts$Platform, sd)
  # reorder:
    CpG.counts.mean <- CpG.counts.mean[c(1,4,3,2),]
    CpG.counts.sd <- CpG.counts.sd[c(1,4,3,2),]

#---------------------------
# Reshape data for plotting:
#---------------------------
  df.long <- gather(CpG.counts, 1:4, key="Coverage" , value = "CpG_count" )
  str(df.long)
  tbl_df(df.long)

  # change names
    df.long$Coverage <- gsub("Total_CpGs", "1x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_5x", "5x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_10x", "10x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_30x" , "30x", df.long$Coverage)

  # Order factors:
    View(df.long)
    str(df.long)
    df.long
    df.long$Coverage <- factor(df.long$Coverage, levels = c("1x", "5x", "10x", "30x"))
    df.long



#-----------------------------
# Set colors for each platform
#------------------------------

  # make dataframe to test color palletes
  x <- c(1,2,3,4,5)
  as.data.frame(x)
  names(x) <- Platforms

  # Determine the order of colors for Platforms
  # # [1] "Agilent"   "Roche"     "Illumina"  "Diagenode" "Nugen" 
  # col.platforms <- c( "blue", "red", "purple","green3",  "yellow")
  # barplot(x, col= col.platforms)
  # legend(x = 'topleft', legend=col.platforms, fill=col.platforms, cex = 0.8)



  # show baseR colors
  #colors()

  # choose palette:
  # library(colorspace)
  # hcl_palettes(plot=TRUE)
  #col.platforms <- sequential_hcl(5,"Plasma") # select 5 colors from sequentian palette


  # # Selected form http://colorbrewer2.org/#type=qualitative&scheme=Accent&n=5
  # col.platforms <- c("#80b1d3",  "#fb8072", "#bebada", "#8dd3c7", "#ffffb3") # looks good
  # barplot(x, col= col.platforms)
  # legend(x = 'topleft', legend=col.platforms, fill=col.platforms, cex = 0.8)


  # Selecting colors using yarr (pirateplot)
    library(yarrr)
    piratepal(palette= "all")
    piratepal("google") 
    # blue         red      yellow       green 
    # "#3D79F3FF" "#E6352FFF" "#F9B90AFF" "#34A74BFF" 
    col.platforms <- c("#3D79F3FF",  "#E6352FFF", "#34A74BFF", "#7570b3" , "#F9B90AFF") # GOOD LOKING DIVERGING PALLETE
    barplot(x, col= col.platforms)
    legend(x = 'topleft', legend=col.platforms, fill=col.platforms, cex = 0.8)

  # make color plaette transparent uing yarr transparent() function:
    col.platforms <- transparent(orig.col = col.platforms, trans.val = 0.3) # BEAUTIFUL :)
    barplot(x, col= col.platforms)

#----------------------

  library(RColorBrewer)
  display.brewer.all() # The palettes are composed of 8-12 distinct colors, but if you have more than 12 categories to graph, you can use the colorRampPalette() function with any of the sequential or diverging palettes. This ramps the color at the necessary interval to create as many hues as your data calls for.
  # col.platforms <- brewer.pal(length(Platforms), "Set1")

  # Select colors for a specified number of variables from any color pallete:
    YlGnBu <- brewer.pal(9, "YlGnBu")
    col.samples <- colorRampPalette(YlGnBu) (20)
    col.samples <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61",  "#fee08b",  "#ffffbf","#c7e9b4", "#7fcdbb", "#2171b5",  "#08519c",  "#08306b", "#41b6c4",  "#1d91c0",  "#225ea8",  "#253494",  "#081d58")


#-------------------------
# Set plotting parameters:
#-------------------------
  ?par


#---------------------
# Ploting with BaseR
#--------------------

  # Gropuped barchart by platform
    barplot(CpG.counts.mean, col=col.platforms[1:4], beside=TRUE, xlab="Number of CpG sites covered at specific depth", font.lab=2 ) 
    legend(x = 'topright', legend=rownames(CpG.counts.mean), fill=col.platforms[1:4], cex = 0.8)
    title('Number of CpG sites covered at specific depth')

  #--------------------------------------------------
  # DOESN'T WORK == Bbecause of using average values 
  # # Get the stacked barplot for percentages:
  # df$percentCpG_under5x <- (df$Total_CpGs - df$NoCpGs_5x)/df$Total_CpGs
  # head(df)
  # colnames(df)
  # CpG.pct <-  df[,c(11,4,6,8:10)]
  # head(CpG.pct)
  # CpG.pct.mean <- apply(CpG.pct[1:4], 2, tapply, CpG.pct$Platform, mean)
  # CpG.pct.mean <- t(CpG.pct.mean[c(1,4,3,2),])
  # barplot(CpG.pct.mean, col=col.platforms[1:4], beside=FALSE, xlab="Percent of CpG sites covered at specific depth", font.lab=2 ) 
  #--------------------------------------------------



#---------------------
# Pirate plot
#--------------------


  library("devtools")
  install_github("ndphillips/yarrr")
  library("yarrr")

  defoult.par <- par()
  par(mar=c(7,8,5,5), las=1)


  # For one coverage 
  # pirateplot(formula = Total_CpGs ~ Platform, 
  #            data = CpG.counts, 
  #            xlab = "Total Number of CpGs covered at 1x", 
  #            ylab = "", 
  #            main = col.platforms,
  #            pal=col.platforms,
  #            gl.col = "gray",
  #            cex.lab = 1, #size of labels
  #            cex.axis = 0.8, # size of axes
  #            cex.names = 0.8 , # size of bean names
  #            theme = 2)

  # Plotting number of Cpgs covered at each depth treshold stratified by Platform:
    pirateplot(formula = CpG_count ~ Platform  + Coverage, 
              data = df.long, 
              xlab = "Coverage", 
              ylab = "", 
              main = "Number of CpGs covered  by each platform at specific depth", # Title
              pal = col.platforms[1:5], # select colors for variables
              cex.lab = 1, #size of labels
              cex.axis = 0.8, # size of axes
              cex.names = 0.7 , # size of bean names,
              bean.lty = 1, # type of line for the bean,
              #xaxt ="n", # don't plot xaxis labels
              theme = 2 # set te theme of the plot
              )
    #legend(x= 21, y=max(df.long$CpG_count), legend=levels(df.long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")
    legend(x= "topright", legend=levels(df.long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")
    mtext(side = 2, " # CpG ", line = 5)
    dev.print(pdf, file= "No CpGs covered per platform.pdf")


#--------------------
# Plot with ggplot2
#----------------------

#--------------------------------------------  
  # Plot
  library(ggplot2)
  devtools::install_github("mikabr/ggpirate")
  library(ggpirate)

  p <- ggplot(df.long, aes(Platform,CpG_count, colour=Platform))  

  # Option-1
  theme_set(theme_light())
  p  + facet_grid(~Coverage) +
    geom_point(position = "jitter") +
    xlab("Coverage") +
    ylab("Count") +
    ggtitle(" Number of CpGs covered at specific depth per platform") +
    #geom_boxplot(alpha=0, colour="black") +
    geom_pirate(aes(colour = Platform, fill=Platform), points = FALSE, bars = TRUE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                points_params = list(shape = 19, alpha = 0.2),
                lines_params = list(size = 0.8)) +
    scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
    scale_fill_manual(values=col.platforms[1:5])
    
  dev.print(pdf, file= "No CpGs covered per platform_ggplot-v1.pdf")

  # Option-2:
  theme_set(theme_light())
  p +  facet_grid(~Coverage) +
    geom_point(position = "jitter") +
    xlab("Coverage") +
    ylab("Count") +
    ggtitle(" Number of CpGs covered at specific depth per platform") +
    geom_boxplot(alpha=0, colour="black") +
    #geom_pirate(aes(colour = Platform), points = FALSE, bars = FALSE, violins = FALSE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                # points_params = list(shape = 19, alpha = 0.2),
                # lines_params = list(size = 0.8)) +
    scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
    scale_fill_manual(values=col.platforms[1:5])

  dev.print(pdf, file= "No CpGs covered per platform_ggplot-v2.pdf")


  # Option-3:
  theme_set(theme_light())
  p  + facet_grid(~Coverage) +
    geom_point(position = "jitter") +
    xlab("Coverage") +
    ylab("Count") +
    ggtitle(" Number of CpGs covered at specific depth per platform") +
    #geom_boxplot(alpha=0, colour="black") +
    geom_pirate(aes(colour = Platform), points = FALSE, bars = FALSE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                points_params = list(shape = 19, alpha = 0.2),
                lines_params = list(size = 0.8)) +
    scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
    scale_fill_manual(values=col.platforms[1:5])

  dev.print(pdf, file= "No CpGs covered per platform_ggplot-v3.pdf")



#=========================================
# Save workspace immage to continue later:
#=========================================

  #  The save.image() function will save all objects currently in your R session:

  save.image(file="Epicapture_Count CpGs and Read depth.RData") 

  # These objects can then be loaded back into a new R session using the load() function:
  load(file="Epicapture_Count CpGs and Read depth.RData")
