


library(reshape2)
library(tidyr)


# Make annotation of experimentally covered regions as pirate plot and as jaccard similarity index correlogram plot

 #load GrangesList objects made fro  CpG_merged coverage files (some ranges have consecutive CpGs!)
load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/mrg_CpGcoverage_grl.PLATFORM.RData"))

#==================================================
# define Functions:
#==================================================

  # find overlaping pairs intersection size - code form HelloRanges package
  findOverlapSize <- function(platform, annotation){
    pairs <- findOverlapPairs(platform, annotation, ignore.strand = TRUE) # find overlaping pairs
    ans <- pintersect(pairs, ignore.strand=TRUE) # get actual intersectiong regions
    sum(width(reduce(ans)))  # Calculate overlap between targeted regions in kilobases
    }

  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  # Get sequences from GRangesList object:  
  grl.seq <- getSeq(genome, platforms_design.grl) # The return values are DNAStringSetList objects.

    # count number of CpGs in overlaping region
  findOverlapCpGCount <- function(platform, annotation){
    pairs <- findOverlapPairs(platform, annotation, ignore.strand = TRUE)
    ans <- pintersect(pairs, ignore.strand=TRUE)
    x <- getSeq(genome, reduce(ans)) # # The return values are DNAStringSetList objects
    sum(vcountPattern("CG", x)) # count number of CG occurancesin the sequence
  }

  countOverlaps()

  # jaccard statistics
  intersects <- intersect(gr_a, gr_b, ignore.strand = TRUE)
  intersection <- sum(width(intersects))
  union <- sum(width(union(gr_a, gr_b, ignore.strand = TRUE)))
  ans <- DataFrame(intersection, union, jaccard = intersection/union, n_intersections = length(intersects))


  # Note that both of these are strand-specific, although findOverlaps has an ignore.strand option.

#=============================================
# 1. Jaccard Statistic
#============================================

    # jaccard statistic, a measure of similarity between two tracks. It is defined as the total width of their intersection over the total width of their union

 # Define a function for the Jaccard statistic
    jaccardIndex <- function(gr_a,gr_b) {
        intersects <- intersect(gr_a, gr_b, ignore.strand = TRUE)
        intersection <- sum(width(intersects))
        union <- sum(width(union(gr_a, gr_b, ignore.strand = TRUE)))
        ans <- DataFrame(intersection, union, jaccard = intersection/union, n_intersections = length(intersects))
        ans
    }
    # Be careful with the strand in your GRanges object!

    # Compute the statistics over all pairs of samples in parallel
    jaccard_matrix <- outer(files, files, function(a, b) mcmapply(jaccard, a, b)) ## wrapper function for parallel computing

    # Make the plot
    library(gplots)
    library(RColorBrewer)
    heatmap.2(jaccard_matrix, col = brewer.pal(9, "Blues"))

#-------------------------------------------------------


    # make one huge GRanges List for all platforms & samples
    grl.all_platforms <- c(grl.Agilent, grl.Roche, grl.Illumina, grl.Diagenode, grl.Nugen)
    class(grl.all_platforms)

    # convert to Simle list object
    files <- as(grl.all_platforms, "SimpleList")
    
    # run jaccardIndex function to get pairwise similarity distance matrix
    jaccard_matrix <- outer(files, files, function(a, b) mcmapply(jaccardIndex, a, b)) # was doing somethingbut doesn't work with mc.cores
# this outputs a matrix of quastion marks characters [?]

    # this outputs a matrix of numbers
    jaccard <- apply(jaccard_matrix, 1:2, function(x) x[[1]]$jaccard) 

    write.table(jaccard,sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/jaccard_matrix.txt"))

    heatmap.2(jaccard, col = brewer.pal(9, "Blues"))

  # I used Intervene Intervene Shiny App  to generate the plots
  https://asntech.shinyapps.io/intervene/




#================================================================================================
# 2. make a dataframe object for each platform containing # CpGs in each feature for each sample
#================================================================================================



# make a nested for loop to calculate #CpGs annotated in each sample per each platform

  findOverlapCpGCount <- function(platform, annotation){
    pairs <- findOverlapPairs(platform, annotation, ignore.strand = TRUE)
    ans <- pintersect(pairs, ignore.strand=TRUE)
    x <- getSeq(genome, reduce(ans)) # # The return values are DNAStringSetList objects
    sum(vcountPattern("CG", x)) # count number of CG occurancesin the sequence
  }

 pairs <- findOverlapPairs(unlist(grl.Agilent[1]), unlist(annotations_genes.grl[1]), ignore.strand = TRUE)

#======================
# Genes + enhancers
#======================
  
    updateObject(annotations_genes.grl,  verbose=TRUE)
    annotations_genes_enhanc.grl <- annotations_genes.grl 
    annotations_genes_enhanc.grl[["enhancers"]] <- annotations_enhancers_fantom.gr
    updateObject(annotations_genes_enhanc.grl,  verbose=TRUE)
    names(annotations_genes_enhanc.grl)

    # save objects
    save(annotations_genes_enhanc.grl,  file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/annotations_genes_enhanc.grl.RData"))

    load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/annotations_genes_enhanc.grl"))


    # find #  CpGs:
      df_Agilent <- data.frame(matrix(vector(),ncol=length(grl.Agilent)))
      names(df_Agilent) <- names(grl.Agilent)
      df_Agilent
      for (n in names(grl.Agilent)) {
        for (i in (names(annotations_genes_enhanc.grl))) {
            df_Agilent[i, n] <- findOverlapCpGCount(unlist(grl.Agilent[n]), unlist(annotations_genes_enhanc.grl[i]))
        }
      }

      df_Roche <- data.frame(matrix(vector(),ncol=length(grl.Roche)))
      names(df_Roche) <- names(grl.Roche)
      df_Roche
      for (n in names(grl.Roche)) {
        for (i in (names(annotations_genes_enhanc.grl))) {
            df_Roche[i, n] <- findOverlapCpGCount(unlist(grl.Roche[n]), unlist(annotations_genes_enhanc.grl[i]))
        }
      }

      df_Illumina <- data.frame(matrix(vector(),ncol=length(grl.Illumina)))
      names(df_Illumina) <- names(grl.Illumina)
      df_Illumina

      for (n in names(grl.Illumina)) {
        for (i in (names(annotations_genes_enhanc.grl))) {
            df_Illumina[i, n] <- findOverlapCpGCount(unlist(grl.Illumina[n]), unlist(annotations_genes_enhanc.grl[i]))
        }
      }

      df_Diagenode <- data.frame(matrix(vector(),ncol=length(grl.Diagenode)))
      names(df_Diagenode) <- names(grl.Diagenode)
      df_Diagenode
      for (n in names(grl.Diagenode)) {
        for (i in (names(annotations_genes_enhanc.grl))) {
            df_Diagenode[i, n] <- findOverlapCpGCount(unlist(grl.Diagenode[n]), unlist(annotations_genes_enhanc.grl[i]))
        }
      }

      df_Nugen <- data.frame(matrix(vector(),ncol=length(grl.Nugen)))
      names(df_Nugen) <- names(grl.Nugen)
      df_Nugen
      for (n in names(grl.Nugen)) {
        for (i in (names(annotations_genes_enhanc.grl))) {
            df_Nugen[i, n] <- findOverlapCpGCount(unlist(grl.Nugen[n]), unlist(annotations_genes_enhanc.grl[i]))
        }
      }

      df_Total <- data.frame(Total=integer())
      for (i in seq_along(names(annotations_genes_enhanc.grl))){ 
        x <- getSeq(genome, unlist(annotations_genes_enhanc.grl[i]))
        df_Total[i,"Total"]  <-sum(vcountPattern("CG", x))} 

    # save all files
      save(df_Agilent,df_Roche,df_Illumina,df_Diagenode,df_Nugen, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_genes_enhanc_platform.RData"))

      load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_genes_enhanc_platform.RData"))

    # Combine platform dataframes to one  
        df_all <- cbind(df_Agilent, df_Roche,df_Illumina,df_Diagenode,df_Nugen)
        names(df_all)
        rownames(df_all)
        df_all$Features <- rownames(df_all) # add features column
        df_all
        dim(df_all)

    # melt to long format
        df_all_long <- melt(df_all)
        names(df_all_long) <- c("Feature","Platform_sample", "CountCpG"  )
        head(df_all_long)
        dim(df_all_long)
        # same thing just using tidyr function
        #  df_all_long <- tidyr::gather(df_all, "Platform_sample", "CountCpG", 1:(ncol(df_all)-1))

    # split 
        df_all_long <- df_all_long %>% separate(Platform_sample, c("Platform", "Sample"), "_") #%>% Passes object on left hand side as first argument (or .argument) of function on righthand side

    # ready for ploting
        head(df_all_long)

    # save file
      save(df_all_long, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_long_genes_enhanc.RData"))

    # Define platforms order:
        Platforms <- c("Agilent", "Roche", "Illumina", "Diagenode", "Nugen")

    # Define factors and level order:
      df_all_long$Sample <- as.factor(df_all_long$Sample)
      df_all_long$Platform <- factor(df_all_long$Platform, levels = Platforms)
      df_all_long
      names(df_all_long)
      head(df_all_long)

      Platforms%in%levels(df_all_long$Platform) # platforms present in dataframe

      df_all_long$Feature <- as.factor(df_all_long$Feature)
      # reorder levels
      levels(df_all_long$Feature) <- levels(df_all_long$Feature)[c(1,2,4,7,5:6,3)]  
      # change level names
      # library(plyr)
      df_all_long$Feature <- revalue(df_all_long$Feature, c("enhancers"="FANTOM5 \n enhancers","hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
          "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"))


    # Save csv file and Rdata objects for future analysis
    write.table(df_all_long, file=paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/AllSamples_CpGCount.txt", sep="\t"))
    save(df_all_long, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_long.PLATFORM.RData"))

    load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_long.PLATFORM.RData"))
  

#---------------------
# Pirate plot
#--------------------


  library("devtools")
  install_github("ndphillips/yarrr")
  library("yarrr")

    # Selecting colors using yarr (pirateplot)
    library(yarrr)
    piratepal(palette= "all")
    piratepal("google") 
    # blue         red      yellow       green 
    # "#3D79F3FF" "#E6352FFF" "#F9B90AFF" "#34A74BFF" 
    col.platforms <- c("#3D79F3FF",  "#E6352FFF", "#34A74BFF", "#7570b3" , "#F9B90AFF") # GOOD LOKING DIVERGING PALLETE


  # Plotting number of Cpgs  in each feature type stratified by Platform:

  # option-1

  pirateplot(formula = CountCpG ~ Platform  + Feature, 
                data = df_all_long, 
                ylab = "Number of CpGs", 
                main = "Number of CpGs covered per platform", # Title
                pal = col.platforms[1:5], # select colors for variables
                theme = 2,   # set theme as 0 to fully customize
                bar.f.col = gray(.8), # bar filling color
                bar.f.o = 0.3, # bars fill opacity
                cex.axis = 0.8, # size of axes
                cex.names = 0.3 , # size of bean names,
                cex.lab = 0.8 # size of labels
              ) 

    legend(x= "topright", 
          # inset=c(-0.25,0),   
          legend=levels(df_all_chrom_long$Platform), fill=col.platforms, cex = 0.7, pt.cex=1, xpd = NA,bty="n") 

    dev.print(pdf, width=12, height=5, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No exp covered CpGs per feature type stratified by Platform_genic.pdf"))

  # option-2
    pirateplot(formula = CountCpG ~ Platform  + Feature, 
              data = df_all_long, 
             # xlab = "Feature", 
              ylab = "Number of CpGs", 
              # ylim = c(1e7,1.05e8),
              main = "Number of CpGs covered per  platform by feature type", # Title
              pal = col.platforms[1:5], # select colors for variables
              theme = 0,   # set theme as 0 to fully customize
              # bean.lty = 1, # type of line for the bean,
              inf.f.o = 0.3, # Inference fill opacity. 0=Turn off inf fill
              inf.b.o = 0.5, # Inference border opacity 0=Turn off inf border
              inf.f.col = "white", # Inf fill col
              inf.b.col = "black", # Inf border col
              point.o = .4,   #  point opacity
              point.pch = 21,
              point.cex = .7,
              point.bg = "white",
              point.col = "black", # Black points
              bar.f.col = gray(.8), # bar filling color
              bar.f.o = 0.3, # bars fill opacity
              bar.b.col="white", # bar border color
              bar.b.o=1, # bars border opacity
              bean.f.o = .6, # Bean fill opacity - Light bean filling
              # bean.b.col ="black",  # Bean border color 
              bean.b.o = .8, # Bean border opacity - Light bean border
              avg.line.col = "black", # avg line col
              avg.line.o = 0.5, #  average line opacitu 0=Turn off 
              cex.lab = 0.8, # size of labels
              yaxt = "n", # add custom axis labels later
              # xaxt ="n", # don't plot xaxis labels
              cex.axis = 0.8, # size of axes
              cex.names = 0.7 , # size of bean names,
              ) 
        # add legend
        legend(x= "topright", 
        # inset=c(-0.25,0), 
        legend=levels(df_all_long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")  

        # add custom y-axis
        axis(2,cex.axis=0.8,  las=1, at=c(2e5,4e5,6e5,8e5, 1e6, 1.2e6,1.4e6,1.6e6, 1.8e6, 2e6 ), 
                labels=c("200 K", "400 K",  "600 K", "800 K", "1 M", "1.2 M","1.4 M", "1.6 M", "1.8 M", "2 M" )) 
        # # Add custom x-axis
        # axis(1,cex.axis=0.9, at=(5*c(1,2,3,4,5,6,7)), 
        #         labels=c("FANTOM5 \n enhancers", "1 to 5kb",  "promoters", "5'UTRs", "exons", "introns","3'UTRs")) 
        

    dev.print(pdf, width=7, height=5, file= paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_genic_v2.pdf"))





    #-------------------
    # Plot with ggplot2
    #-------------------
    library(ggplot2)
    devtools::install_github("mikabr/ggpirate")
    library(ggpirate)

    p <- ggplot(df_all_long, aes(x=Platform, y=CountCpG )) + facet_wrap(~Feature)
    
    p +
      geom_pirate(aes(colour=Platform,fill=Platform), show.legend=TRUE) +
      theme_light() 

    p  +
      geom_point(position = "jitter") +
      xlab("") +
      ylab("No CpGs") +
      ggtitle(" Number of CpGs  per platform") +
      #geom_boxplot(alpha=0, colour="black") +
      geom_pirate(aes(colour = Platform, fill = Platform), points = FALSE, bars = TRUE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                  points_params = list(shape = 19, alpha = 0.2),
                  lines_params = list(size = 0.8)) +
      scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
      scale_fill_manual(values=col.platforms[1:5]) +
      theme_minimal() 

    ggsave(width=12, height=6, paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_genic_ggplot.pdf"))



#======================
# CpG Islands
#======================
  
  updateObject(annotations_CpGs.grl,  verbose=TRUE)
  head(annotations_CpGs.grl)
  names(annotations_CpGs.grl)



    # find #  CpGs:
      df_Agilent_CGI <- data.frame(matrix(vector(),ncol=length(grl.Agilent)))
      names(df_Agilent_CGI) <- names(grl.Agilent)
      df_Agilent_CGI
      for (n in names(grl.Agilent)) {
        for (i in (names(annotations_CpGs.grl))) {
            df_Agilent_CGI[i, n] <- findOverlapCpGCount(unlist(grl.Agilent[n]), unlist(annotations_CpGs.grl[i]))
        }
      }

      df_Roche_CGI <- data.frame(matrix(vector(),ncol=length(grl.Roche)))
      names(df_Roche_CGI) <- names(grl.Roche)
      df_Roche_CGI
      for (n in names(grl.Roche)) {
        for (i in (names(annotations_CpGs.grl))) {
            df_Roche_CGI[i, n] <- findOverlapCpGCount(unlist(grl.Roche[n]), unlist(annotations_CpGs.grl[i]))
        }
      }

      df_Illumina_CGI <- data.frame(matrix(vector(),ncol=length(grl.Illumina)))
      names(df_Illumina_CGI) <- names(grl.Illumina)
      df_Illumina_CGI
      for (n in names(grl.Illumina)) {
        for (i in (names(annotations_CpGs.grl))) {
            df_Illumina_CGI[i, n] <- findOverlapCpGCount(unlist(grl.Illumina[n]), unlist(annotations_CpGs.grl[i]))
        }
      }

      df_Diagenode_CGI <- data.frame(matrix(vector(),ncol=length(grl.Diagenode)))
      names(df_Diagenode_CGI) <- names(grl.Diagenode)
      df_Diagenode_CGI
      for (n in names(grl.Diagenode)) {
        for (i in (names(annotations_CpGs.grl))) {
            df_Diagenode_CGI[i, n] <- findOverlapCpGCount(unlist(grl.Diagenode[n]), unlist(annotations_CpGs.grl[i]))
        }
      }

      df_Nugen_CGI <- data.frame(matrix(vector(),ncol=length(grl.Nugen)))
      names(df_Nugen_CGI) <- names(grl.Nugen)
      df_Nugen_CGI
      for (n in names(grl.Nugen)) {
        for (i in (names(annotations_CpGs.grl))) {
            df_Nugen_CGI[i, n] <- findOverlapCpGCount(unlist(grl.Nugen[n]), unlist(annotations_CpGs.grl[i]))
        }
      }

      df_Total_CGI <- data.frame(Total=integer())
      for (i in seq_along(names(annotations_CpGs.grl))){ 
        x <- getSeq(genome, unlist(annotations_CpGs.grl[i]))
        df_Total_CGI[i,"Total"]  <-sum(vcountPattern("CG", x))} 

    # save all files
      save(df_Agilent,df_Roche,df_Illumina,df_Diagenode,df_Nugen, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_CGI_platform.RData"))

      load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_CGI_platform.RData"))

    # Combine platform dataframes to one  
        df_all_CGI <- cbind(df_Agilent_CGI, df_Roche_CGI,df_Illumina_CGI,df_Diagenode_CGI,df_Nugen_CGI)
        names(df_all_CGI)
        rownames(df_all_CGI)
        df_all_CGI$Features <- rownames(df_all_CGI) # add features column
        df_all_CGI
        dim(df_all_CGI)

    # melt to long format
        df_all_CGI_long <- melt(df_all_CGI)
        names(df_all_CGI_long) <- c("Feature","Platform_sample", "CountCpG"  )
        head(df_all_CGI_long)
        dim(df_all_CGI_long)

    # split 
        df_all_CGI_long <- df_all_CGI_long %>% separate(Platform_sample, c("Platform", "Sample"), "_") #%>% Passes object on left hand side as first argument (or .argument) of function on righthand side

    # ready for ploting
        head(df_all_CGI_long)

    # Define platforms order:
        Platforms <- c("Agilent", "Roche", "Illumina", "Diagenode", "Nugen")

  # Define factors and level order:
    df_all_CGI_long$Sample <- as.factor(df_all_CGI_long$Sample)
    df_all_CGI_long$Platform <- factor(df_all_CGI_long$Platform, levels = Platforms)
    df_all_CGI_long
    names(df_all_CGI_long)

    Platforms%in%levels(df_all_CGI_long$Platform) # platforms present in dataframe

      df_all_CGI_long$Feature <- as.factor(df_all_CGI_long$Feature)
      # reorder levels
      levels(df_all_CGI_long$Feature) <- levels(df_all_CGI_long$Feature)[c(2,4,3,1)]  
      # change level names
      # library(plyr)
      df_all_CGI_long$Feature <- revalue(df_all_CGI_long$Feature, c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shores"="shores", "hg19_cpg_shelves"="shelves","hg19_cpg_inter"="inter CGI"))

    # Save csv file and Rdata objects for future analysis
    save(df_all_CGI_long, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_CGI_long.PLATFORM.RData"))

    load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_CGI_long.PLATFORM.RData"))

  #-----------------
  # pirate lot
  #---------------

    pirateplot(formula = CountCpG ~ Platform  + Feature, 
              data = df_all_CGI_long, 
              ylab = "Number of CpGs", 
                main = "Number of CpGs covered per platform", # Title
                pal = col.platforms[1:5], # select colors for variables
                theme = 2,   # set theme as 0 to fully customize
                bar.f.col = gray(.8), # bar filling color
                bar.f.o = 0.3, # bars fill opacity
                cex.axis = 0.8, # size of axes
                cex.names = 0.3 , # size of bean names,
                cex.lab = 0.8 # size of labels
              ) 

    legend(x= "topright", 
          # inset=c(-0.25,0),   
          legend=levels(df_all_chrom_long$Platform), fill=col.platforms, cex = 0.7, pt.cex=1, xpd = NA,bty="n") 

        # # add custom y-axis
        # axis(2,cex.axis=0.8,  las=1, at=c(2e5,4e5,6e5,8e5, 1e6, 1.2e6,1.4e6,1.6e6, 1.8e6, 2e6 ), 
        #         labels=c("200 K", "400 K",  "600 K", "800 K", "1 M", "1.2 M","1.4 M", "1.6 M", "1.8 M", "2 M" )) 
        

    dev.print(pdf, width=6, height=5, file= paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No exp covered CpGs per feature type stratified by Platform_CGI.pdf"))


  # option-2
    pirateplot(formula = CountCpG ~ Platform  + Feature, 
              data = df_all_CGI_long, 
              # xlab = "Feature", 
              ylab = "Number of CpGs", 
              # ylim = c(1e7,1.05e8),
              main = "Number of CpGs covered per  platform by feature type", # Title
              pal = col.platforms[1:5], # select colors for variables
              theme = 0,   # set theme as 0 to fully customize
              # bean.lty = 1, # type of line for the bean,
              inf.f.o = 0.3, # Inference fill opacity. 0=Turn off inf fill
              inf.b.o = 0.5, # Inference border opacity 0=Turn off inf border
              inf.f.col = "white", # Inf fill col
              inf.b.col = "black", # Inf border col
              point.o = .4,   #  point opacity
              point.pch = 21,
              point.cex = .7,
              point.bg = "white",
              point.col = "black", # Black points
              bar.f.col = gray(.8), # bar filling color
              bar.f.o = 0.3, # bars fill opacity
              bar.b.col="white", # bar border color
              bar.b.o=1, # bars border opacity
              bean.f.o = .6, # Bean fill opacity - Light bean filling
              # bean.b.col ="black",  # Bean border color 
              bean.b.o = .8, # Bean border opacity - Light bean border
              avg.line.col = "black", # avg line col
              avg.line.o = 0.5, #  average line opacitu 0=Turn off 
              cex.lab = 0.8, # size of labels
              # yaxt = "n", # add custom axis labels later
              # xaxt ="n", # don't plot xaxis labels
              cex.axis = 0.8, # size of axes
              cex.names = 0.7 , # size of bean names,
              ) 
        # add legend
        legend(x= "topright", 
        # inset=c(-0.25,0), 
        legend=levels(df_all_long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")  

        # add custom y-axis
        # axis(2,cex.axis=0.8,  las=1, at=c(2e5,4e5,6e5,8e5, 1e6, 1.2e6,1.4e6,1.6e6, 1.8e6, 2e6 ), 
        #         labels=c("200 K", "400 K",  "600 K", "800 K", "1 M", "1.2 M","1.4 M", "1.6 M", "1.8 M", "2 M" )) 
        

    dev.print(pdf, width=8, height=5, file= paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_CGI_v2.pdf"))

    #-------------------
    # Plot with ggplot2
    #-------------------
    library(ggplot2)
    devtools::install_github("mikabr/ggpirate")
    library(ggpirate)

    p <- ggplot(df_all_CGI_long, aes(x=Platform, y=CountCpG )) + facet_wrap(~Feature)
    
    p +
      geom_pirate(aes(colour=Platform,fill=Platform), show.legend=TRUE) +
      theme_light() 

    p  +
      geom_point(position = "jitter") +
      xlab("") +
      ylab("No CpGs") +
      ggtitle(" Number of CpGs  per platform") +
      #geom_boxplot(alpha=0, colour="black") +
      geom_pirate(aes(colour = Platform, fill = Platform), points = FALSE, bars = TRUE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                  points_params = list(shape = 19, alpha = 0.2),
                  lines_params = list(size = 0.8)) +
      scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
      scale_fill_manual(values=col.platforms[1:5]) +
      theme_minimal() 

    ggsave(width=12, height=6, paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_CGI_ggplot.pdf"))


#========================
# Annotate chromHMM marks:
#========================

 annotations_chromHMM.grl <- GRangesList( "ActivePromoter"=annotations_ActivePromoter.gr, "Heterochrom"=annotations_Heterochrom.gr, "Insulator"=annotations_Insulator.gr, "PoisedPromoter"=annotations_PoisedPromoter.gr, "Repressed"=annotations_Repressed.gr, "StrongEnhancer"=annotations_StrongEnhancer.gr,  "TxnElongation"=annotations_TxnElongation.gr,  "TxnTransition"=annotations_TxnTransition.gr,  "WeakEnhancer"=annotations_WeakEnhancer.gr,  "WeakPromoter"=annotations_WeakPromoter.gr,  "WeakTxn"=annotations_WeakTxn.gr)

  names(annotations_chromHMM.grl)
 annotations_chromHMM.grl <-  annotations_chromHMM.grl[seqnames(annotations_chromHMM.grl) == c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM") ]

  annotations_chromHMM.grl <-  endoapply(annotations_chromHMM.grl, reduce) 


    # find #  CpGs:
      df_Agilent_chrom <- data.frame(matrix(vector(),ncol=length(grl.Agilent)))
      names(df_Agilent_chrom) <- names(grl.Agilent)
      df_Agilent_chrom
      for (n in names(grl.Agilent)) {
        for (i in (names(annotations_chromHMM.grl))) {
            df_Agilent_chrom[i, n] <- findOverlapCpGCount(unlist(grl.Agilent[n]), unlist(annotations_chromHMM.grl[i]))
        }
      }

      df_Roche_chrom <- data.frame(matrix(vector(),ncol=length(grl.Roche)))
      names(df_Roche_chrom) <- names(grl.Roche)
      df_Roche_chrom
      for (n in names(grl.Roche)) {
        for (i in (names(annotations_chromHMM.grl))) {
            df_Roche_chrom[i, n] <- findOverlapCpGCount(unlist(grl.Roche[n]), unlist(annotations_chromHMM.grl[i]))
        }
      }

      df_Illumina_chrom <- data.frame(matrix(vector(),ncol=length(grl.Illumina)))
      names(df_Illumina_chrom) <- names(grl.Illumina)
      df_Illumina_chrom

      for (n in names(grl.Illumina)) {
        for (i in (names(annotations_chromHMM.grl))) {
            df_Illumina_chrom[i, n] <- findOverlapCpGCount(unlist(grl.Illumina[n]), unlist(annotations_chromHMM.grl[i]))
        }
      }

      df_Diagenode_chrom <- data.frame(matrix(vector(),ncol=length(grl.Diagenode)))
      names(df_Diagenode_chrom) <- names(grl.Diagenode)
      df_Diagenode_chrom
      for (n in names(grl.Diagenode)) {
        for (i in (names(annotations_chromHMM.grl))) {
            df_Diagenode_chrom[i, n] <- findOverlapCpGCount(unlist(grl.Diagenode[n]), unlist(annotations_chromHMM.grl[i]))
        }
      }

      df_Nugen_chrom <- data.frame(matrix(vector(),ncol=length(grl.Nugen)))
      names(df_Nugen_chrom) <- names(grl.Nugen)
      df_Nugen_chrom
      for (n in names(grl.Nugen)) {
        for (i in (names(annotations_chromHMM.grl))) {
            df_Nugen_chrom[i, n] <- findOverlapCpGCount(unlist(grl.Nugen[n]), unlist(annotations_chromHMM.grl[i]))
        }
      }

      df_Total_chrom <- data.frame(Total=integer())
      for (i in seq_along(names(annotations_chromHMM.grl))){ 
        x <- getSeq(genome, unlist(annotations_chromHMM.grl[i]))
        df_Total_chrom[i,"Total"]  <-sum(vcountPattern("CG", x))} 

    # save all files
      save(df_Agilent_chrom,df_Roche_chrom,df_Illumina_chrom,df_Diagenode_chrom,df_Nugen_chrom, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_chrom_platform.RData"))

      load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_chrom_platform.RData"))

    # Combine platform dataframes to one  
        df_all_chrom <- cbind(df_Agilent_chrom, df_Roche_chrom,df_Illumina_chrom,df_Diagenode_chrom,df_Nugen_chrom)
        names(df_all_chrom)
        rownames(df_all_chrom)
        df_all_chrom$Features <- rownames(df_all_chrom) # add features column
        df_all_chrom
        dim(df_all_chrom)

    # melt to long format
        df_all_chrom_long <- melt(df_all_chrom)
        names(df_all_chrom_long) <- c("Feature","Platform_sample", "CountCpG"  )
        head(df_all_chrom_long)
        dim(df_all_chrom_long)

    # split 
        df_all_chrom_long <- df_all_chrom_long %>% separate(Platform_sample, c("Platform", "Sample"), "_") #%>% Passes object on left hand side as first argument (or .argument) of function on righthand side

    # ready for ploting
        head(df_all_chrom_long)

    # Define platforms order:
        Platforms <- c("Agilent", "Roche", "Illumina", "Diagenode", "Nugen")

  # Define factors and level order:
    df_all_chrom_long$Sample <- as.factor(df_all_chrom_long$Sample)
    df_all_chrom_long$Platform <- factor(df_all_chrom_long$Platform, levels = Platforms)
    df_all_chrom_long
    names(df_all_chrom_long)

    Platforms%in%levels(df_all_chrom_long$Platform) # platforms present in dataframe

    df_all_chrom_long$Feature <- as.factor(df_all_chrom_long$Feature)
      # reorder levels
      levels(df_all_chrom_long$Feature) <- levels(df_all_chrom_long$Feature)[c(1,10,4,6,9,3,8,7,11,5,2)]  


      # change level names
    # library(plyr)
    df_all_chrom_long$Feature <- revalue(df_all_chrom_long$Feature, c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn"))

    

    # Save csv file and Rdata objects for future analysis
    save(df_all_chrom_long, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_chrom_long.PLATFORM.RData"))

    load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/df_all_chrom_long.PLATFORM.RData"))


#-----------------
  # pirate lot
  #---------------
 # option-1

  pirateplot(formula = CountCpG ~ Platform  + Feature, 
                data = df_all_chrom_long, 
                ylab = "Number of CpGs", 
                main = "Number of CpGs covered per platform", # Title
                pal = col.platforms[1:5], # select colors for variables
                theme = 2,   # set theme as 0 to fully customize
                bar.f.col = gray(.8), # bar filling color
                bar.f.o = 0.3, # bars fill opacity
                cex.axis = 0.8, # size of axes
                cex.names = 0.3 , # size of bean names,
                cex.lab = 0.8 # size of labels
              ) 

    legend(x= "topright", 
          # inset=c(-0.25,0),   
          legend=levels(df_all_chrom_long$Platform), fill=col.platforms, cex = 0.7, pt.cex=1, xpd = NA,bty="n")  

    dev.print(pdf, width=18, height=5, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No exp covered CpGs per feature type stratified by Platform_chrom.pdf"))

  # option-2
    pirateplot(formula = CountCpG ~ Platform  + Feature, 
              data = df_all_chrom_long, 
              # xlab = "Feature", 
              ylab = "Number of CpGs", 
              # ylim = c(1e7,1.05e8),
              main = "Number of CpGs covered per  platform by feature type", # Title
              pal = col.platforms[1:5], # select colors for variables
              theme = 0,   # set theme as 0 to fully customize
              # bean.lty = 1, # type of line for the bean,
              inf.f.o = 0.3, # Inference fill opacity. 0=Turn off inf fill
              inf.b.o = 0.5, # Inference border opacity 0=Turn off inf border
              inf.f.col = "white", # Inf fill col
              inf.b.col = "black", # Inf border col
              point.o = .4,   #  point opacity
              point.pch = 21,
              point.cex = .7,
              point.bg = "white",
              point.col = "black", # Black points
              bar.f.col = gray(.8), # bar filling color
              bar.f.o = 0.3, # bars fill opacity
              bar.b.col="white", # bar border color
              bar.b.o=1, # bars border opacity
              bean.f.o = .6, # Bean fill opacity - Light bean filling
              # bean.b.col ="black",  # Bean border color 
              bean.b.o = .8, # Bean border opacity - Light bean border
              avg.line.col = "black", # avg line col
              avg.line.o = 0.5, #  average line opacitu 0=Turn off 
              cex.lab = 0.8, # size of labels
              # yaxt = "n", # add custom axis labels later
              # xaxt ="n", # don't plot xaxis labels
              cex.axis = 0.8, # size of axes
              cex.names = 0.7 , # size of bean names,
              ) 
        # add legend
        legend(x= "topright", 
        # inset=c(-0.25,0), 
        legend=levels(df_all_chrom_long$Platform), fill=col.platforms, cex = 0.8, xpd = NA,bty="n")  

        # add custom y-axis
        # axis(2,cex.axis=0.8,  las=1, at=c(2e5,4e5,6e5,8e5, 1e6, 1.2e6,1.4e6,1.6e6, 1.8e6, 2e6 ), 
        #         labels=c("200 K", "400 K",  "600 K", "800 K", "1 M", "1.2 M","1.4 M", "1.6 M", "1.8 M", "2 M" )) 
        # # Add custom x-axis
        # axis(1,cex.axis=0.9, at=(5*c(1,2,3,4,5,6,7)), 
        #         labels=c("FANTOM5 \n enhancers", "1 to 5kb",  "promoters", "5'UTRs", "exons", "introns","3'UTRs")) 
        

    dev.print(pdf, width=11, height=5, file= paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_chrom_v2.pdf"))





    #-------------------
    # Plot with ggplot2
    #-------------------
    library(ggplot2)
    devtools::install_github("mikabr/ggpirate")
    library(ggpirate)

    p <- ggplot(df_all_long, aes(x=Platform, y=CountCpG )) + facet_wrap(~Feature)
    
    p +
      geom_pirate(aes(colour=Platform,fill=Platform), show.legend=TRUE) +
      theme_light() 

    p  +
      geom_point(position = "jitter") +
      xlab("") +
      ylab("No CpGs") +
      ggtitle(" Number of CpGs  per platform") +
      #geom_boxplot(alpha=0, colour="black") +
      geom_pirate(aes(colour = Platform, fill = Platform), points = FALSE, bars = TRUE, violins = TRUE, # Each of the layers can be turned off, e.g. for just means and confidence intervals:
                  points_params = list(shape = 19, alpha = 0.2),
                  lines_params = list(size = 0.8)) +
      scale_color_manual(values=col.platforms[1:5]) + # select colors for dots
      scale_fill_manual(values=col.platforms[1:5]) +
      theme_minimal() 

    ggsave(width=12, height=6, paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_genic_ggplot.pdf"))

        

    dev.print(pdf, width=16, height=5, file= paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/No CpGs covered per platform_chrom.pdf"))






