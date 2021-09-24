#################################################################
#                                                               #
# Project: Epicapture - Targeted BS-seq platform comparison     #
# Description: Number of unique CpGs covered by each platform   #
# Author: Miljana Tanic                                         #
# Date last edited: 05-03-2019                                  #
#                                                               #
#################################################################

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

  # Total number of CpG sites - output from coverage2cytosine with information from top and bottom strand merged into one
  samples.agilent <- list.files(path=path2Agilent, pattern="CpG_evidence.cov") 
  samples.roche <- list.files(path=path2Roche, pattern="CpG_evidence.cov") 
  samples.illumina <- list.files(path=path2Illumina, pattern="CpG_evidence.cov") 
  samples.diagenode <- list.files(path=path2Diagenode, pattern="CpG_evidence.cov") 
  samples.nugen <- list.files(path=path2Nugen, pattern="CpG_evidence.cov") 

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
  tbl_df(df.long)

  # change names
    df.long$Coverage <- gsub("Total_CpGs", "1x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_5x", "5x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_10x", "10x", df.long$Coverage)
    df.long$Coverage <- gsub( "NoCpGs_30x" , "30x", df.long$Coverage)

  # Order factors:
    df.long$Coverage <- factor(df.long$Coverage, levels = c("1x", "5x", "10x", "30x"))
    
  # Selecting colors using yarr (pirateplot)
    library(yarrr)
    piratepal(palette= "all")
    piratepal("google") 
    col.platforms <- c("#3D79F3FF",  "#E6352FFF", "#34A74BFF", "#7570b3" , "#F9B90AFF") 
  # make color plaette transparent uing yarr transparent() function:
    col.platforms <- transparent(orig.col = col.platforms, trans.val = 0.3)


  # Select colors for a specified number of variables from any color pallete:
    library(RColorBrewer)
    YlGnBu <- brewer.pal(9, "YlGnBu")
    col.samples <- colorRampPalette(YlGnBu) (20)
    col.samples <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61",  "#fee08b",  "#ffffbf","#c7e9b4", "#7fcdbb", "#2171b5",  "#08519c",  "#08306b", "#41b6c4",  "#1d91c0",  "#225ea8",  "#253494",  "#081d58")
    
#---------------------
# Pirate plot
#--------------------

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
    legend(x= "topright", legend=levels(df.long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")
    mtext(side = 2, " # CpG ", line = 5)
    dev.print(pdf, file= "No CpGs covered per platform.pdf")

