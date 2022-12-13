############################################################################
#   Title: Epicapture Target Regions - 3.step 
#   Description: Script used to annotate  Target capture design files across 5 platforms
#   Author: Miljana Tanic (m.tanic@ucl.ac.uk) 
#   Created: November 2018
#   Edited: Mart 2019
#   Edited: Dec 2019 - 
####################################################################################

# !!! If packages, commands or something doesn't work - disconect from the  server and reconect, it will usually "miracolously" start working!

#=========================================
# Load libraries and set working directory
#=========================================


  source("https://bioconductor.org/biocLite.R")
  # biocLite("rtracklayer")
  #biocLite("annotatr") 

  library(rtracklayer) # not necessary for importing BED files, can do it manualy crating GRanges object from dataframe 
  library(dplyr)
  library(readr)
  library(annotatr)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19) # sequence of the hg19 genome
  library(ggplot2)
  library(HelloRanges)
  library(reshape2)
  library(yarrr)



#------------------------------------------
# Set Path to files
#------------------------------------------

     # set work directly on R server for increased speed
    #setwd("/home/regmani@ad.ucl.ac.uk/data/EpiCapture/")
    setwd("/mnt/254b78b9-76b4-422d-84b1-cc632bff60f7/regmani/Epicapture")
    list.files()

    getwd() 

    # set resource and RESULTS directory
    Epicapture <- "/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/"
    list.files(Epicapture)

    BEDfiles <- paste0(Epicapture,"BEDfiles/")
    list.files(BEDfiles)


    path2Agilent_BED <- paste0(BEDfiles, "AgilentBEDfiles/S03770311/")
    path2Roche_BED <- paste0(BEDfiles,"RocheBEDfiles/")
    path2Illumina_BED <- paste0(BEDfiles,"IlluminaBEDfiles/")
    path2RRBS_BED <- paste0(BEDfiles,"RRBS_in_silico_BEDfile/")

    list.files(path2Agilent_BED)

    RESULTS <- Epicapture <- paste0(Epicapture,"RESULTS/")

#------------------------------------------
# Load annotation files and Platform files:
#------------------------------------------

    # Load design files:
    load(file=paste0(BEDfiles,"Platform_Design_GRanges.RData")) 

    # Load feature annotations:
    load(file=paste0(BEDfiles,"FeatureAnnotations/Genes_annotations_GRangesList.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/CpG_annottations_GRangesList.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/Enhancers_annotations_GRanges.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/lncRNA_annotations_GRanges.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/ChromHMM-byState_GRanges.RData")) 
    # load(file=paste0(BEDfiles, "FeatureAnnotations/CpG_locations_hg19_GRanges.RData"))
    # load(file=paste0(BEDfiles, "FeatureAnnotations/all_cpgs.RData"))
    # load(file=paste0(BEDfiles,"FeatureAnnotations/hg19_CpGsites_GRanges.RData"))
    # load(file=paste0(BEDfiles,"FeatureAnnotations/UCSC_CGI_GRanges.RData"))

#==================================================
# find overlaping pairs intersection size Function:
#==================================================

  # code form HelloRanges package
  findOverlapSize <- function(platform, annotation){
    pairs <- findOverlapPairs(platform, annotation, ignore.strand = TRUE) # find overlaping pairs
    ans <- pintersect(pairs, ignore.strand=TRUE) # get actual intersectiong regions
    sum(width(reduce(ans)))  # Calculate overlap between targeted regions in kilobases
    }

  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  # Get sequences from GRangesList object:  
  grl.seq <- getSeq(genome, platforms_design.grl) # The return values are DNAStringSetList objects.


  findOverlapCpGCount <- function(platform, annotation){
    pairs <- findOverlapPairs(platform, annotation, ignore.strand = TRUE)
    ans <- pintersect(pairs, ignore.strand=TRUE)
    x <- getSeq(genome,reduce(ans)) # # The return values are DNAStringSetList objects
    sum(vcountPattern("CG", x)) # count number of CG occurancesin the sequence
  }


#==========================
# Set colors for Plotting
#==========================

    # Selecting colors using yarr (pirateplot)
      library(yarrr)
      piratepal(palette= "all")
      piratepal("google") 
      # blue         red      yellow       green 
      # "#3D79F3FF" "#E6352FFF" "#F9B90AFF" "#34A74BFF" 

      # platforms <- c( "Agilent", "Roche", "Illumina", "Diagenode", "Nugen")
      col.platforms <- c("#E6352FFF", "#3D79F3FF", "#34A74BFF", "#7570b3" , "#F9B90AFF") # GOOD LOKING DIVERGING PALLETE
      barplot(rep(1,length(col.platforms)), col= col.platforms)
      legend(x = 'topleft', legend=col.platforms, fill=col.platforms, cex = 0.8)
    
    # make color plaette transparent uing yarr transparent() function:
      col.platforms <- transparent(orig.col = col.platforms, trans.val = 0.2) # BEAUTIFUL :)

    # change default colors in ggplot2
      opts <- options()  # save old options
      library(RColorBrewer)
      feature.cols <- brewer.pal(n = 8, name = "Dark2")
      barplot(rep(1,length(feature.cols)), col= feature.cols)


      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 
      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

#------------------------
# Plot data with ggplot2:
#------------------------


#==========================================
# Annotate platforms for all Feature types:
#==========================================


#======================
# Genes + enhancers
#======================

  # Make GRangesList for all annotation files:
    grl <- GRangesList(enhancers=annotations_enhancers_fantom.gr, lncRNA=annotations_lncRNA.gr)
    grl <- annotations_genes.grl 
    grl[["enhancers"]] <- annotations_enhancers_fantom.gr
    updateObject(grl,  verbose=TRUE)
    names(grl)


    # names(grl) <- c("1to5kb", "3'UTR", "5'UTR","exons","introns", "promoters", "CpG intergenic", "CpG islands", "CpG shelves", "CpG shores", "enhancers")

  # Make data frame to store size values:
    df <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())

  # Calculate overlap - intersection between platform design file and anntation file:
    for (i in (names(grl))) {
      df[i,"Agilent"] <- findOverlapSize(Agilent_SureSelect.gr, unlist(grl[i]))}

    for (i in (names(grl))) {
      df[i,"Roche"] <- findOverlapSize(Roche_EpiGiant.gr, unlist(grl[i]))}

    for (i in (names(grl))) {
      df[i,"Illumina"] <- findOverlapSize(Illumina_EPIC.gr, unlist(grl[i]))}

    for (i in (names(grl))) {
      df[i,"RRBS"] <- findOverlapSize(RRBS.gr.sub, unlist(grl[i]))}

  # Calculate total size of each feature in annotation file:
    df$Total <- sum(width(reduce(grl)))

  # Inspect summary dataframe
    head(df)
    rownames(df)
  # Save dataframe with summary
    save(df, file = "design_features_size.RData")
    write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_features_size.txt"))

  # Transform dataframe to long format for plotting:
    df_long <- melt(t(df[,-1]))
    names(df_long) <- c("Platform", "Feature", "Size_bp" )
    head(df_long)
    df_long$Size_Mb <- round(df_long$Size_bp/1000000,2)


    #------------------------------------------
    # Plot Feature size by Platform - barplots:
    #------------------------------------------

      gene_annots <- ggplot(data=df_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
  
        scale_fill_brewer(palette="Set1") +
        # scale_color_manual(values=piratepal("appletv")[1:5]) + # select colors for dots
        # scale_fill_manual(values=col.platforms[1:4])+
        ylab("Genomic size (Mb)") + 
        ggtitle("Breadth of coverage of genomic features by platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                          labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters", "enhancers"="FANTOM5 \n enhancers"
                                )) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=7, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_GenomicFeatures_ByFeature_size.pdf"))

    #-------------------------------------
    # Features size by Platform - stacked
    #-------------------------------------

      # Stacked Plot by Platform
        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=df_long, aes(x=Platform, y=Size_Mb, fill=Feature)) +
          geom_bar(stat = "identity", colour="black", size=0.25)+
          ylab("Genomic size (Mb)") + 
          ggtitle("Breadth of coverage of genomic features by platform") + 
          scale_fill_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",
                                  "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("enhancers"="FANTOM5 enhancers","hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                  "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters")) + # Change the name of items
          theme_minimal()
        gene_annots
        ggsave(width=4, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_GenomicFeatures_byPlatform_size.pdf"))

    #--------------------------------------
    # Percent features covered by Platform
    #--------------------------------------


      # Calculate percentages:
      head(df)
      df_pct <- apply(df, 2, function(x) x/df$Total)
      head(df_pct)
      df_pct <- abs(df_pct)

      # Plot percentages:
      df_pct_long <- melt(t(df_pct[,-1]))
      names(df_pct_long) <-  c("Platform", "Feature", "Percent_covered" )

      annots_pct <- ggplot(data=df_pct_long, aes(x=Feature, y=Percent_covered, fill=Platform))  +
            geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
            # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
            # scale_y_continuous("",formatter="percent") + 
            scale_fill_brewer(palette="Set1") +
            # scale_fill_hue(c = 40) +
            # scale_color_manual(values=piratepal("appletv")[1:5]) + # select colors for dots
            # scale_fill_manual(values=col.platforms[1:4])+
            ylab("% covered") + 
            xlab("") + 
            ggtitle("Percent genomic features covered by each platform") + 
            theme(axis.title.x=element_blank()) +
            # coord_flip()+ # make horizontal chart
            scale_x_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                    "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                              labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                    "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters", "enhancers"="FANTOM5 \n enhancers")) + # Change the name of items
            theme_classic()
      annots_pct
     
    # Save figure
      ggsave(width=7, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_GenomicFeatures_ByFeature_size_pct.pdf"))

    #-------------------------------------
    # Piechart % features by Platform
    #-------------------------------------

      # Compute percentages of each features covered by each platform
        df_platform_pct <- apply(df, 2, function(x) pct=x/sum(x))

      # calculate percentage:
          df_platform_long <- melt(t(df_platform_pct[,-1]))
          names(df_platform_long) <- c("Platform", "Feature", "pct" )
          head(df_platform_long)

      gene_annots <- ggplot(data=df_platform_long, aes(x="", y=pct, fill=Feature)) +
              facet_grid(. ~ Platform) + 
              geom_bar(width = 1, stat = "identity",  color="white") +
              coord_polar("y", start=0) +
              ylab("Genomic size (bp)") + 
              scale_fill_discrete(name="Genomic features",
                        limits=c("enhancers","hg19_genes_1to5kb","hg19_genes_5UTRs",  
                        "hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                        labels=c("enhancers"="FANTOM5 enhancers","hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                        "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"),
                        guide = guide_legend(label = TRUE, label.position = "right",
                                              legend.theme	= element_text(size = 8)
                        )) + # Change the name of items
              theme(legend.position="bottom") +
              theme_void() # remove background, grid, numeric labels
        gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_GenomicFeatures_byPlatform_pct.pdf"))

  #------------------------------------------------------
  # Number of CpGs covered in each feature  by Platform
  #------------------------------------------------------

    # find #  CpGs:
      df_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


      for (i in (names(grl))) {
        df_cpg[i,"Agilent"] <- findOverlapCpGCount(Agilent_SureSelect.gr, unlist(grl[i]))}

      for (i in (names(grl))) {
        df_cpg[i,"Roche"] <- findOverlapCpGCount(Roche_EpiGiant.gr, unlist(grl[i]))}

      for (i in (names(grl))) {
        df_cpg[i,"Illumina"] <- findOverlapCpGCount(Illumina_EPIC.gr, unlist(grl[i]))}

      for (i in (names(grl))) {
        df_cpg[i,"RRBS"] <- findOverlapCpGCount(RRBS.gr.sub, unlist(grl[i]))}

      for (i in seq_along(names(grl))){ 
        x <- getSeq(genome, unlist(grl[i]))
        df_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}


      head(df_cpg)
      save(df_cpg, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_genomic_cpg.RData"))
      write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_GenomicFeatures_CpGs.txt"))

    # Transform dataframe to long format for plotting:
      df_cpg_long <- melt(t(df_cpg[,-1]))
      names(df_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
      head(df_cpg_long)
      df_cpg_long$No_CpGs_M <- round(df_cpg_long$No_CpGs/1000000,2)

    #----------------------------------------
    # Plot Features by Platform - barplots:
    #----------------------------------------

      gene_annots <- ggplot(data=df_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        # scale_fill_manual(values=col.platforms[1:4]) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per genomic feature by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                          labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters", "enhancers"="FANTOM5 \n enhancers")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=7, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_GenomicFeatures_ByFeature_NoCpGs.pdf"))
    
    #---------------------------------------------
    # Stacked Plot NoCpGs per Feature by Platform
    #----------------------------------------------

        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=df_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
          geom_bar(stat = "identity", colour="black", size=0.25)+
          ylab("# CpGs") + 
          ggtitle("Number of CpGs per genomic feature by platform") + 
          scale_fill_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",
                                  "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("enhancers"="FANTOM5 enhancers","hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                  "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters")) + # Change the name of items
          theme_minimal()
        gene_annots
        ggsave(width=5, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_GenomicFeatures_byPlatform_NoCpGs.pdf"))

    #----------------------------------------------
    # Percent CpGs covered per Feature by Platform
    #----------------------------------------------

      # Calculate percentages of total features covered by each platform
        head(df_cpg)
        df_cpg_pct <- apply(df_cpg, 2, function(x) x/df_cpg$Total)
        head(df_cpg_pct)

      # Plot percentages:
      df_cpg_pct_long <- melt(t(df_cpg_pct[,-1]))
      names(df_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

      gene_annots_pct <- ggplot(data=df_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
        #          position = position_dodge(0.9), size=3.5)+
        scale_fill_brewer(palette="Set1")+
        ylab("% CpGs covered") + 
        xlab("") +
        ggtitle("Percent CpGs covered per features by platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers", "hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                          labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters", "enhancers"="FANTOM5 \n enhancers")) + # Change the name of items
        theme_classic()
      gene_annots_pct
      ggsave(width=7, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_GenomicFeatures_ByFeature_pctCpGs.pdf"))
  
    #----------------------------------------
    # Piechart % CpGs in features by Platform
    #----------------------------------------
 
      # Compute percentages of each features covered by each platform
        df_platform_pct_CpGs <- apply(df_cpg, 2, function(x) pct=x/sum(x))

      # calculate percentage:
          df_platform_pct_CpGs_long <- melt(t(df_platform_pct_CpGs[,-1]))
          names(df_platform_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
          head(df_platform_pct_CpGs_long)

      gene_annots <- ggplot(data=df_platform_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
              facet_grid(. ~ Platform) + 
              geom_bar(width = 1, stat = "identity",  color="white") +
              coord_polar("y", start=0) +
              ylab("% CpGs covered") + 
              ggtitle("Percent CpGs covered per features by platform") + 
              scale_fill_discrete(name="Genomic features",
                        limits=c("enhancers","hg19_genes_1to5kb","hg19_genes_5UTRs",  
                        "hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                        labels=c("enhancers"="FANTOM5 \n enhancers","hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                        "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"),
                        guide = guide_legend(label = TRUE, label.position = "right",
                                              legend.theme	= element_text(size = 8)
                        )) + # Change the name of items
              theme(legend.position="bottom") +
              theme_void() # remove background, grid, numeric labels
        gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_GenomicFeatures_byPlatform_pctCpGs.pdf"))


#===========================
# Annotate Genic Features:
#===========================

  # Inspect loaded object:
  # head(annotations_genes.gr)
  # table(annotations_genes.gr$type)

  head(annotations_genes.grl)
  names(annotations_genes.grl)
  #df_genes <- data.frame(Total=sum(width(reduce(annotations_genes.grl))))
  

  # Convert to GRanges object
  # annotations_genes.gr <- unlist(annotations_genes.grl) # Downstream coercing to dataframe doesn't work with unlisted GRabgesList object

  # # check levels in GRanges objects metadata column type:
  # levels(as.factor(mcols(annotations_genes.gr)$type)) 

  # # can be used for subseting by type of feature:
  # annotations_genes.gr[mcols(annotations_genes.gr)$type=="hg19_genes_5UTRs",]

  # Coerce to dataframe and calculate size of each feature:
  # df_genes <- data.frame(annotations_genes.gr)
  # tapply(as.numeric(df_genes$width), df_genes$type, sum)

  #---------------------------
  # Features size by Platform
  #---------------------------

    # find overlaping pairs:
      df_genes <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


      for (i in (names(annotations_genes.grl))) {
        df_genes[i,"Agilent"] <- findOverlapSize(Agilent_SureSelect.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_genes[i,"Roche"] <- findOverlapSize(Roche_EpiGiant.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_genes[i,"Illumina"] <- findOverlapSize(Illumina_EPIC.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_genes[i,"RRBS"] <- findOverlapSize(RRBS.gr.sub, unlist(annotations_genes.grl[i]))}

      df_genes$Total <- sum(width(reduce(annotations_genes.grl)))

      head(df_genes)
      save(df_genes, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_genes_size.RData"))
      write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_GenomicFeatures_size.txt"))
      load(file=paste0(RESULTS,"1_PlatformDesignDiferences/df_genes_size.RData"))


    # Transform dataframe to long format for plotting:
      df_genes_long <- melt(t(df_genes[,-1]))
      names(df_genes_long) <- c("Platform", "Feature", "Size_bp" )
      head(df_genes_long)
      df_genes_long$Size_Mb <- round(df_genes_long$Size_bp/1000000,2) 

    #------------------------------------------
    # Plot Feature size by Platform - barplots:
    #------------------------------------------

      gene_annots <- ggplot(data=df_genes_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
  
        scale_fill_brewer(palette="Set1") +
        # scale_color_manual(values=piratepal("appletv")[1:5]) + # select colors for dots
        # scale_fill_manual(values=col.platforms[1:4])+
        ylab("Genomic size (Mb)") + 
        ggtitle("Breadth of coverage of genomic features by platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                          labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                                )) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Genic_ByFeature_size.pdf"))

    #-------------------------------------
    # Features size by Platform - stacked
    #-------------------------------------

      # Stacked Plot by Platform
        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=df_genes_long, aes(x=Platform, y=Size_Mb, fill=Feature)) +
          geom_bar(stat = "identity", colour="black", size=0.25)+
          ylab("Genomic size (Mb)") + 
          ggtitle("Breadth of coverage of genomic features by platform") + 
          scale_fill_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",
                                  "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                  "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                                  )) + # Change the name of items
          theme_minimal()
        gene_annots
        ggsave(width=5, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Genic_byPlatform_size.pdf"))

    #--------------------------------------
    # Percent features covered by Platform
    #--------------------------------------

      # Calculate percentages of total features covered by each platform
        head(df_genes)
        df_genes_pct <- apply(df_genes, 2, function(x) x/df_genes$Total)
        head(df_genes_pct)

      # Plot percentages:
      df_genes_pct_long <- melt(t(df_genes_pct[,-1]))
      names(df_genes_pct_long) <-  c("Platform", "Feature", "Percent_covered" )

      gene_annots_pct <- ggplot(data=df_genes_pct_long, aes(x=Feature, y=Percent_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
        #          position = position_dodge(0.9), size=3.5)+
        scale_fill_brewer(palette="Set1") +
        ylab("% covered") + 
        ggtitle("Percent genomic features covered by each platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                            "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                      labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                            "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                              )) + # Change the name of items

        theme_classic()
      gene_annots_pct
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Genic_ByFeature_sizepct.pdf"))

    #-------------------------------------
    # Piechart % features by Platform
    #-------------------------------------

            # Compute percentages of each features covered by each platform
              df_platform_pct <- apply(df_genes, 2, function(x) pct=x/sum(x))

            # calculate percentage:
                df_platform_long <- melt(t(df_platform_pct[,-1]))
                names(df_platform_long) <- c("Platform", "Feature", "pct" )
                head(df_platform_long)

            gene_annots <- ggplot(data=df_platform_long, aes(x="", y=pct, fill=Feature)) +
                    facet_grid(. ~ Platform) + 
                    geom_bar(width = 1, stat = "identity",  color="white") +
                    coord_polar("y", start=0) +
                    ylab("Genomic size (bp)") + 
                    ggtitle("Breadth of coverage of genomic features by platform") + 
                    scale_fill_discrete(name="Genic features",
                              limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs",  
                              "hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                              labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                              "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"),
                              guide = guide_legend(label = TRUE, label.position = "right",
                                                    legend.theme	= element_text(size = 8)
                              )) + # Change the name of items
                    theme(legend.position="bottom") +
                    theme_void() # remove background, grid, numeric labels
              gene_annots
            ggsave(width=6, height=5,paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_Genic_byPlatform_pct.pdf"))

  #------------------------------------------------------
  # Number of CpGs covered in each feature  by Platform
  #------------------------------------------------------

    # find #  CpGs:
      df_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


      for (i in (names(annotations_genes.grl))) {
        df_cpg[i,"Agilent"] <- findOverlapCpGCount(Agilent_SureSelect.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_cpg[i,"Roche"] <- findOverlapCpGCount(Roche_EpiGiant.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_cpg[i,"Illumina"] <- findOverlapCpGCount(Illumina_EPIC.gr, unlist(annotations_genes.grl[i]))}

      for (i in (names(annotations_genes.grl))) {
        df_cpg[i,"RRBS"] <- findOverlapCpGCount(RRBS.gr.sub, unlist(annotations_genes.grl[i]))}

      for (i in seq_along(names(annotations_genes.grl))){ 
        x <- getSeq(genome, unlist(annotations_genes.grl[i]))
        df_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}


      head(df_cpg)
      save(df_cpg, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_genes_cpg.RData"))
      write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_Genic_CpGs.txt"))

    # Transform dataframe to long format for plotting:
      df_cpg_long <- melt(t(df_cpg[,-1]))
      names(df_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
      head(df_cpg_long)
      df_cpg_long$No_CpGs_M <- round(df_cpg_long$No_CpGs/1000000,2)

    #----------------------------------------
    # Plot Features by Platform - barplots:
    #----------------------------------------

      gene_annots <- ggplot(data=df_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        # scale_fill_manual(values=col.platforms[1:4]) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per genomic feature by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters","hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                          labels=c("hg19_genes_1to5kb"="1 to 5kb","hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns","hg19_genes_promoters"="promoters"
                                )) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Genic_ByFeature_NoCpGs.pdf"))
    
    #----------------------------------------
    # Stacked Plot NoCpGs per Feature by Platform
    #----------------------------------------

        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=df_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
          geom_bar(stat = "identity", colour="black", size=0.25)+
          ylab("# CpGs") + 
          ggtitle("Number of CpGs per genomic feature by platform") + 
          scale_fill_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",
                                  "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                  "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                                  )) + # Change the name of items
          theme_minimal()
        gene_annots
        ggsave(width=5, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Genic_byPlatform_NoCpGs.pdf"))

    #----------------------------------------------
    # Percent CpGs covered per Feature by Platform
    #----------------------------------------------

      # Calculate percentages of total features covered by each platform
        head(df_cpg)
        df_cpg_pct <- apply(df_cpg, 2, function(x) x/df_cpg$Total)
        head(df_cpg_pct)

      # Plot percentages:
      df_cpg_pct_long <- melt(t(df_cpg_pct[,-1]))
      names(df_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

      gene_annots_pct <- ggplot(data=df_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
        #          position = position_dodge(0.9), size=3.5)+
        scale_fill_brewer(palette="Set1")+
        ylab("% CpGs covered") + 
        ggtitle("Percent CpGs covered per features by platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                            "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                      labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                            "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                              )) + # Change the name of items

        theme_classic()
      gene_annots_pct
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Genic_ByFeature_pctCpGs.pdf"))
  
    #----------------------------------------
    # Piechart % CpGs in features by Platform
    #----------------------------------------
 
      # Compute percentages of each features covered by each platform
        df_platform_pct_CpGs <- apply(df_cpg, 2, function(x) pct=x/sum(x))

      # calculate percentage:
          df_platform_pct_CpGs_long <- melt(t(df_platform_pct_CpGs[,-1]))
          names(df_platform_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
          head(df_platform_pct_CpGs_long)

      gene_annots <- ggplot(data=df_platform_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
              facet_grid(. ~ Platform) + 
              geom_bar(width = 1, stat = "identity",  color="white") +
              coord_polar("y", start=0) +
              ylab("% CpGs covered") + 
              ggtitle("Percent CpGs covered per features by platform") + 
              scale_fill_discrete(name="Genic features",
                        limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs",  
                        "hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                        labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                        "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"),
                        guide = guide_legend(label = TRUE, label.position = "right",
                                              legend.theme	= element_text(size = 8)
                        )) + # Change the name of items
              theme(legend.position="bottom") +
              theme_void() # remove background, grid, numeric labels
        gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_byPlatform_pctCpGs.pdf"))

#================
# Annotate CpGs:
#================

  # Inspect loaded object:
  head(annotations_CpGs.gr)
  head(annotations_CpGs.grl)
  names(annotations_CpGs.grl)


  # find overlaping pairs:

    df_CGI <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


    for (i in (names(annotations_CpGs.grl))) {
        df_CGI[i,"Agilent"] <- findOverlapSize(Agilent_SureSelect.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
        df_CGI[i,"Roche"] <- findOverlapSize(Roche_EpiGiant.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
        df_CGI[i,"Illumina"] <- findOverlapSize(Illumina_EPIC.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
        df_CGI[i,"RRBS"] <- findOverlapSize(RRBS.gr.sub, unlist(annotations_CpGs.grl[i]))}

    # df_CGI$Total <- sum(width(reduce(annotations_CpGs.grl, ignore.strand=T))) # this gives negative values for inter-CGI, while  values for other regions are the same
    df_CGI$Total <- sapply(annotations_CpGs.grl, function(x) sum(width(reduce(x))))

    head(df_CGI)
    save(df_CGI, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_CGI_size.RData"))
    write.table(df_CGI, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_CGI_size.txt"))
    load(file=paste0(RESULTS,"1_PlatformDesignDiferences/df_CGI_size.RData"))


  # Transform dataframe to long format for plotting:
    df_CGI_long <- melt(t(df_CGI[,-1]))
    names(df_CGI_long) <- c("Platform", "Feature", "Size_bp" )
    head(df_CGI_long)
    df_CGI_long$Size_Mb <- round(df_CGI_long$Size_bp/1000000)

  #------------------------------------------
  # Plot Feature size by Platform - barplots:
  #------------------------------------------

    gene_annots <- ggplot(data=df_CGI_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
      # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
      scale_fill_brewer(palette="Set1") +
      ylab("Genomic size (Mb)") + 
      ggtitle("Breadth of coverage of CGI by each platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
      scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shores", "hg19_cpg_shelves",  
                                "hg19_cpg_inter"), # Change the order of items
                      labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shores"="shores", "hg19_cpg_shelves"="shelves","hg19_cpg_inter"="inter CGI")) + # Change the name of items
      theme_classic()
    gene_annots
    ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_CGI_ByFeature_size.pdf"))

  #----------------------------------------------
  # Features size by Platform - stacked  barplot
  #----------------------------------------------

    # Stacked Plot by Platform
      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=df_CGI_long, aes(x=Platform, y=Size_Mb, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25) +
        ylab("Genomic size (Mb)") + 
        ggtitle("Breadth of coverage of CGI by each platform") + 
        scale_fill_discrete(labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shores"="shores","hg19_cpg_shelves"="shelves",  "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(width=5, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_CGI_byPlatform_size.pdf"))

  #--------------------------------------
  # Percent features covered by Platform
  #--------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_CGI)
      df_CGI_pct <- apply(df_CGI, 2, function(x) x/df_CGI$Total)
      head(df_CGI_pct)

    # Transform to long format
      df_CGI_pct_long <- melt(t(df_CGI_pct[,-1]))
      names(df_CGI_pct_long) <-  c("Platform", "Feature", "Percent_covered" )

    # Plot percentages:
      gene_annots_pct <- ggplot(data=df_CGI_pct_long, aes(x=Feature, y=Percent_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        scale_fill_brewer(palette="Set1") +
        ylab("% covered") + 
        ggtitle("Percent genomic features covered by each platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shelves","hg19_cpg_shores",  
                                    "hg19_cpg_inter"), # Change the order of items
                        labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                                    "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots_pct
    ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_CGI_ByFeature_sizepct.pdf"))

  #-------------------------------------
  # Piechart % features by Platform
  #-------------------------------------

    # Compute percentages of each features covered by each platform
    df_platform_pct <- apply(df_CGI, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_long <- melt(t(df_platform_pct[,-1]))
        names(df_platform_long) <- c("Platform", "Feature", "pct" )
        head(df_platform_long)

    gene_annots <- ggplot(data=df_platform_long, aes(x="", y=pct, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("Genomic size (bp)") + 
            ggtitle("Breadth of coverage of genomic features by platform") + 
            scale_fill_discrete(name="CGI features",
                      labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                            "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI") ,
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
    gene_annots
    ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_CGI_byPlatform_pct.pdf"))

#------------------------------------------------------
# Number of CpGs covered in each feature  by Platform
#------------------------------------------------------

  # find #  CpGs:
    df_CGI_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


    for (i in (names(annotations_CpGs.grl))) {
      df_CGI_cpg[i,"Agilent"] <- findOverlapCpGCount(Agilent_SureSelect.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
      df_CGI_cpg[i,"Roche"] <- findOverlapCpGCount(Roche_EpiGiant.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
      df_CGI_cpg[i,"Illumina"] <- findOverlapCpGCount(Illumina_EPIC.gr, unlist(annotations_CpGs.grl[i]))}

    for (i in (names(annotations_CpGs.grl))) {
      df_CGI_cpg[i,"RRBS"] <- findOverlapCpGCount(RRBS.gr.sub, unlist(annotations_CpGs.grl[i]))}

    for (i in seq_along(names(annotations_CpGs.grl))){ 
      x <- getSeq(genome, unlist(annotations_CpGs.grl[i]))
      df_CGI_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

    head(df_CGI_cpg)
    save(df_CGI_cpg, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_CGI_cpg.RData"))
    write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_CGI_CpGs.txt"))

  # Transform dataframe to long format for plotting:
    df_CGI_cpg_long <- melt(t(df_CGI_cpg[,-1]))
    names(df_CGI_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
    head(df_CGI_cpg_long)
    df_CGI_cpg_long$No_CpGs_M <- round(df_CGI_cpg_long$No_CpGs/1000000,2)

  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=df_CGI_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shelves","hg19_cpg_shores",  
                                    "hg19_cpg_inter"), # Change the order of items
                        labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                                    "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_CGI_ByFeature_NoCpGs.pdf"))
  
  #----------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #----------------------------------------

      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=df_CGI_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGsper feature category by platform") + 
        scale_fill_discrete(labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves", "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(width=5, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_CGI_byPlatform_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_CGI_cpg)
      df_CGI_cpg_pct <- apply(df_CGI_cpg, 2, function(x) x/df_CGI_cpg$Total)
      head(df_CGI_cpg_pct)

    # Plot percentages:
    df_CGI_cpg_pct_long <- melt(t(df_CGI_cpg_pct[,-1]))
    names(df_CGI_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=df_CGI_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_brewer(palette="Set1")+
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shelves","hg19_cpg_shores",  
                                    "hg19_cpg_inter"), # Change the order of items
                        labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                                    "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
      theme_classic()
    gene_annots_pct
    ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_CGI_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      df_platform_CGI_pct_CpGs <- apply(df_CGI_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_CGI_pct_CpGs_long <- melt(t(df_platform_CGI_pct_CpGs[,-1]))
        names(df_platform_CGI_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(df_platform_CGI_pct_CpGs_long)

    gene_annots <- ggplot(data=df_platform_CGI_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("% CpGs covered") + 
            ggtitle("Percent CpGs covered per features by platform") + 
            scale_fill_discrete(name="CGI features",
                      labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                            "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI") ,
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
      gene_annots
    ggsave(width=6, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_CGI_byPlatform_pctCpGs.pdf"))

#===============================
# Annotate Enhancers and lncRNA:
#===============================

 # Inspect loaded object:
        head(annotations_enhancers_fantom.gr)
        head(annotations_lncRNA.gr)
        annotations_regul.grl <- GRangesList(enhancers=annotations_enhancers_fantom.gr,lncRNA=annotations_lncRNA.gr)

        df_regul <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())
        df_regul <- data.frame(Total=sum(width(reduce(annotations_regul.grl))))

        for (i in (names(annotations_regul.grl))) {
        df_regul[i,"Agilent"] <- findOverlapSize(Agilent_SureSelect.gr, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        df_regul[i,"Roche"] <- findOverlapSize(Roche_EpiGiant.gr, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        df_regul[i,"Illumina"] <- findOverlapSize(Illumina_EPIC.gr, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        df_regul[i,"RRBS"] <- findOverlapSize(RRBS.gr.sub, unlist(annotations_regul.grl[i]))}

        df_regul$Total <-  sapply(annotations_regul.grl, function(x) sum(width(reduce(x))))

        head(df_regul)
        save(df_regul, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_regul_size.RData"))
        write.table(df_regul, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_regul_size.txt"))
        load( file=paste0(RESULTS,"1_PlatformDesignDiferences/df_regul_size.RData"))


    # Calculate percentages
        head(df_regul)
        df_regul_pct <- apply(df_regul, 2, function(x) x/df_regul$Total)
        head(df_regul_pct)

    # Transform dataframe to long format for plotting:
        df_regul_long <- melt(t(df_regul[,-1]))
        names(df_regul_long) <- c("Platform", "Feature", "Size_bp" )
        head(df_regul_long)
        df_regul_long$Size_Mb <- round(df_regul_long$Size_bp/1000000)

    #------------------------------------------
    # Plot Feature size by Platform - barplots:
    #------------------------------------------

      gene_annots <- ggplot(data=df_regul_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
        #facet_grid(. ~Feature)
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        # scale_color_manual(values=piratepal("appletv")[1:5]) + # select colors for dots
        # scale_fill_manual(values=col.platforms[1:5])+
        ylab("Genomic size (Mb)") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Regul_ByFeature_size.pdf"))

    #-------------------------------------
    # Features size by Platform - stacked
    #-------------------------------------
      # NOT INFORMATIVE

      # Stacked Plot by Platform
        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=df_regul_long, aes(x=Platform, y=Size_Mb, fill=Feature)) +
          geom_bar(stat = "identity", colour="black", size=0.25)+
          ylab("Genomic size (Mb)") + 
          ggtitle("Breadth of coverage of genomic features by platform") + 
          scale_fill_discrete( labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
          theme_minimal()
        gene_annots
        ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Regul_byPlatform_size.pdf"))

  #--------------------------------------
  # Percent features covered by Platform
  #--------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_regul)
      df_regul_pct <- apply(df_regul, 2, function(x) x/df_regul$Total)
      head(df_regul_pct)

    # Transform to long format
      df_regul_pct_long <- melt(t(df_regul_pct[,-1]))
      names(df_regul_pct_long) <-  c("Platform", "Feature", "Percent_covered" )

    # Plot percentages:
      gene_annots_pct <- ggplot(data=df_regul_pct_long, aes(x=Feature, y=Percent_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        scale_fill_brewer(palette="Set1") +
        ylab("% covered") + 
        ggtitle("Percent genomic features covered by each platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_minimal()
      gene_annots_pct
    ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Regul_ByFeature_sizepct.pdf"))


  #-------------------------------------
  # Piechart % features by Platform
  #-------------------------------------
    # NOT THAT INFORMATIVE  

    # Compute percentages of each features covered by each platform
    df_platform_pct <- apply(df_regul, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_long <- melt(t(df_platform_pct[,-1]))
        names(df_platform_long) <- c("Platform", "Feature", "pct" )
        head(df_platform_long)

    gene_annots <- ggplot(data=df_platform_long, aes(x="", y=pct, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("Genomic size (bp)") + 
            ggtitle("Breadth of coverage of genomic features by platform") + 
            scale_fill_discrete( labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA"), # Change the name of items
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
    gene_annots
    ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_Regul_byPlatform_pct.pdf"))

#------------------------------------------------------
# Number of CpGs covered in each feature  by Platform
#------------------------------------------------------

  # find #  CpGs:
    df_regul_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


    for (i in (names(annotations_regul.grl))) {
      df_regul_cpg[i,"Agilent"] <- findOverlapCpGCount(Agilent_SureSelect.gr, unlist(annotations_regul.grl[i]))}

    for (i in (names(annotations_regul.grl))) {
      df_regul_cpg[i,"Roche"] <- findOverlapCpGCount(Roche_EpiGiant.gr, unlist(annotations_regul.grl[i]))}

    for (i in (names(annotations_regul.grl))) {
      df_regul_cpg[i,"Illumina"] <- findOverlapCpGCount(Illumina_EPIC.gr, unlist(annotations_regul.grl[i]))}

    for (i in (names(annotations_regul.grl))) {
      df_regul_cpg[i,"RRBS"] <- findOverlapCpGCount(RRBS.gr.sub, unlist(annotations_regul.grl[i]))}

    for (i in seq_along(names(annotations_regul.grl))){ 
      x <- getSeq(genome, unlist(annotations_regul.grl[i]))
      df_regul_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

    head(df_regul_cpg)
    save(df_regul_cpg, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_regul_cpg.RData"))
    write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_Regul_CpGs.txt"))

  # Transform dataframe to long format for plotting:
    df_regul_cpg_long <- melt(t(df_regul_cpg[,-1]))
    names(df_regul_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
    head(df_regul_cpg_long)
    df_regul_cpg_long$No_CpGs_M <- round(df_regul_cpg_long$No_CpGs/1000000,2)

  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=df_regul_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Regul_ByFeature_NoCpGs.pdf"))
  
  #----------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #----------------------------------------

      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=df_regul_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        scale_fill_discrete( labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Regul_byPlatform_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_regul_cpg)
      df_regul_cpg_pct <- apply(df_regul_cpg, 2, function(x) x/df_regul_cpg$Total)
      head(df_regul_cpg_pct)

    # Plot percentages:
    df_regul_cpg_pct_long <- melt(t(df_regul_cpg_pct[,-1]))
    names(df_regul_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=df_regul_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_brewer(palette="Set1")+
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
      theme_minimal()
    gene_annots_pct
    ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Regul_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      df_platform_regul_pct_CpGs <- apply(df_regul_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_regul_pct_CpGs_long <- melt(t(df_platform_regul_pct_CpGs[,-1]))
        names(df_platform_regul_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(df_platform_regul_pct_CpGs_long)

    gene_annots <- ggplot(data=df_platform_regul_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("% CpGs covered") + 
            ggtitle("Percent CpGs covered per features by platform") + 
            scale_fill_discrete( labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA"), # Change the name of items
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
      gene_annots
    ggsave(paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_Regul_byPlatform_pctCpGs.pdf"))

#================
# Non coding RNA
#================

ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/9.0/genome_coordinates/Homo_sapiens.GRCh38.bed.gz

#========================
# Annotate chromHMM marks:
#========================


 annotations_chromHMM.grl <- GRangesList( "ActivePromoter"=annotations_ActivePromoter.gr, "Heterochrom"=annotations_Heterochrom.gr, "Insulator"=annotations_Insulator.gr, "PoisedPromoter"=annotations_PoisedPromoter.gr, "Repressed"=annotations_Repressed.gr, "StrongEnhancer"=annotations_StrongEnhancer.gr,  "TxnElongation"=annotations_TxnElongation.gr,  "TxnTransition"=annotations_TxnTransition.gr,  "WeakEnhancer"=annotations_WeakEnhancer.gr,  "WeakPromoter"=annotations_WeakPromoter.gr,  "WeakTxn"=annotations_WeakTxn.gr)

  names(annotations_chromHMM.grl)
 annotations_chromHMM.grl <-  annotations_chromHMM.grl[seqnames(annotations_chromHMM.grl) == c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM") ]

  sapply(annotations_chromHMM.grl, function(x) length(x))  
  sapply(annotations_chromHMM.grl, function(x) length(reduce(x)))  
  sapply(annotations_chromHMM.grl, function(x) sum(width((x))))  
  sapply(annotations_chromHMM.grl, function(x) sum(width(reduce(x))))  

  annotations_chromHMM.grl <-  endoapply(annotations_chromHMM.grl, reduce) 


  df_chrom <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


  for (i in (names(annotations_chromHMM.grl))) {
    df_chrom[i,"Agilent"] <- findOverlapSize(Agilent_SureSelect.gr, unlist(annotations_chromHMM.grl[i]))}

  for (i in (names(annotations_chromHMM.grl))) {
    df_chrom[i,"Roche"] <- findOverlapSize(Roche_EpiGiant.gr, unlist(annotations_chromHMM.grl[i]))}

  for (i in (names(annotations_chromHMM.grl))) {
    df_chrom[i,"Illumina"] <- findOverlapSize(Illumina_EPIC.gr, unlist(annotations_chromHMM.grl[i]))}

  for (i in (names(annotations_chromHMM.grl))) {
    df_chrom[i,"RRBS"] <- findOverlapSize(RRBS.gr.sub, unlist(annotations_chromHMM.grl[i]))}

  df_chrom$Total <- sapply(annotations_chromHMM.grl, function(x) sum(width(reduce(x))))

  head(df_chrom)
  rownames(df_chrom)
  save(df_chrom, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_chrom.RData"))
  write.table(df_chrom, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_Chrom_Features_size.txt"))
  load(file=paste0(RESULTS,"1_PlatformDesignDiferences/df_chrom.RData"))


  # Transform dataframe to long format for plotting:
    df_chrom_long <- melt(t(df_chrom[,-1]))
    names(df_chrom_long) <- c("Platform", "Feature", "Size_bp" )
    head(df_chrom_long)
    df_chrom_long$Size_Mb <- round(df_chrom_long$Size_bp/1000000,2)

  #------------------------------------------
  # Plot Feature size by Platform - barplots:
  #------------------------------------------


    gene_annots <- ggplot(data=df_chrom_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
      # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
      scale_fill_brewer(palette="Set1") +
      ylab("Genomic size (Mb)") + 
      ggtitle("Breadth of coverage of chromatine state by each platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
      scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
        labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
      theme_classic()
    gene_annots
    ggsave(width=12, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Chrom_ByFeature_size.pdf"))

  #----------------------------------------------
  # Features size by Platform - stacked  barplot
  #----------------------------------------------

    # Stacked Plot by Platform
      feature.cols <- brewer.pal(n=11, "Set3")
      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=df_chrom_long, aes(x=Platform, y=Size_Mb, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25) +
        ylab("Genomic size (Mb)") + 
        ggtitle("Breadth of coverage of chromatine state by each platform") + 
        scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(width=5, height=5,paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Chrom_byPlatform_size.pdf"))

  #--------------------------------------
  # Percent features covered by Platform
  #--------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_chrom)
      df_chrom_pct <- apply(df_chrom, 2, function(x) x/df_chrom$Total)
      head(df_chrom_pct)

    # Transform to long format
      df_chrom_pct_long <- melt(t(df_chrom_pct[,-1]))
      names(df_chrom_pct_long) <-  c("Platform", "Feature", "Percent_covered" )

    # Plot percentages:
      gene_annots_pct <- ggplot(data=df_chrom_pct_long, aes(x=Feature, y=Percent_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        scale_fill_brewer(palette="Set1") +
        ylab("% covered") + 
        ggtitle("Percent genomic features covered by each platform") + 
        theme(axis.title.x=element_blank()) +
        # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
          theme_classic()
      gene_annots_pct
    ggsave(width=12, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Chrom_ByFeature_sizepct.pdf"))

  #-------------------------------------
  # Piechart % features by Platform
  #-------------------------------------

    # Compute percentages of each features covered by each platform
    df_platform_pct <- apply(df_chrom, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_long <- melt(t(df_platform_pct[,-1]))
        names(df_platform_long) <- c("Platform", "Feature", "pct" )
        head(df_platform_long)

    gene_annots <- ggplot(data=df_platform_long, aes(x="", y=pct, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("Genomic size (bp)") + 
            ggtitle("Breadth of coverage of genomic features by platform") + 
            scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn"),
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
    gene_annots
    ggsave(width=6, height=5,paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_Chrom_byPlatform_pct.pdf"))

#------------------------------------------------------
# Number of CpGs covered in each feature  by Platform
#------------------------------------------------------

  # find #  CpGs:
    df_chrom_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), RRBS=integer())


    for (i in (names(annotations_chromHMM.grl))) {
      df_chrom_cpg[i,"Agilent"] <- findOverlapCpGCount(Agilent_SureSelect.gr, unlist(annotations_chromHMM.grl[i]))}

    for (i in (names(annotations_chromHMM.grl))) {
      df_chrom_cpg[i,"Roche"] <- findOverlapCpGCount(Roche_EpiGiant.gr, unlist(annotations_chromHMM.grl[i]))}

    for (i in (names(annotations_chromHMM.grl))) {
      df_chrom_cpg[i,"Illumina"] <- findOverlapCpGCount(Illumina_EPIC.gr, unlist(annotations_chromHMM.grl[i]))}

    for (i in (names(annotations_chromHMM.grl))) {
      df_chrom_cpg[i,"RRBS"] <- findOverlapCpGCount(RRBS.gr.sub, unlist(annotations_chromHMM.grl[i]))}

    for (i in seq_along(names(annotations_chromHMM.grl))){ 
      x <- getSeq(genome, unlist(annotations_chromHMM.grl[i]))
      df_chrom_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

    head(df_chrom_cpg)
    save(df_chrom_cpg, file=paste0(RESULTS,"1_PlatformDesignDiferences/df_chrom_cpg.RData"))
    write.table(df, sep="\t", file=paste0(RESULTS,"1_PlatformDesignDiferences/design_Chrom_CpGs.txt"))

  # Transform dataframe to long format for plotting:
    df_chrom_cpg_long <- melt(t(df_chrom_cpg[,-1]))
    names(df_chrom_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
    head(df_chrom_cpg_long)
    df_chrom_cpg_long$No_CpGs_M <- round(df_chrom_cpg_long$No_CpGs/1000000,2)

  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=df_chrom_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(width=12, height=5, file=paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Chrom_ByFeature_NoCpGs.pdf"))

  #----------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #----------------------------------------
      gene_annots <- ggplot(data=df_chrom_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_minimal()
      gene_annots
      ggsave(width=5, height=5,paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/StakedBarplot_Chrom_byPlatform_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(df_chrom_cpg)
      df_chrom_cpg_pct <- apply(df_chrom_cpg, 2, function(x) x/df_chrom_cpg$Total)
      head(df_chrom_cpg_pct)

    # Plot percentages:
    df_chrom_cpg_pct_long <- melt(t(df_chrom_cpg_pct[,-1]))
    names(df_chrom_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=df_chrom_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="white", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_brewer(palette="Set1")+
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
      theme_classic()
    gene_annots_pct
    ggsave(width=12, height=5, paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Barplot_Chrom_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      df_platform_Chrom_pct_CpGs <- apply(df_chrom_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        df_platform_Chrom_pct_CpGs_long <- melt(t(df_platform_Chrom_pct_CpGs[,-1]))
        names(df_platform_Chrom_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(df_platform_Chrom_pct_CpGs_long)

    gene_annots <- ggplot(data=df_platform_Chrom_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
            facet_grid(. ~ Platform) + 
            geom_bar(width = 1, stat = "identity",  color="white") +
            coord_polar("y", start=0) +
            ylab("% CpGs covered") + 
            ggtitle("Percent CpGs covered per features by platform") + 
            scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn"),
                      guide = guide_legend(label = TRUE, label.position = "right",
                                            legend.theme	= element_text(size = 8)
                      )) + # Change the name of items
            theme(legend.position="bottom") +
            theme_void() # remove background, grid, numeric labels
      gene_annots
    ggsave(width=6, height=5,paste0(RESULTS, "1_PlatformDesignDiferences/1.4_DesignFeatureAnnotations/Piechart_Chrom_byPlatform_pctCpGs.pdf"))

