############################################################################
#   Title: Epicapture - annotate expperimentaly covered regions 
#   Description: Script used to annotate  Target capture coverage files across 5 platforms
#   Author: Miljana Tanic (m.tanic@ucl.ac.uk) 
#   Created: November 2018
#   Last edited: Mart 2019
#   edited: Dec 2019 # plot all annotations
#####################################################################################

#=========================================
# Load libraries and set working directory
#=========================================

  library(GenomicRanges)
  library(AnnotationHub)
  library(Rsamtools)
  library(data.table)
  library(genomation)
  library(rtracklayer) # not necessary for importing BED files, can do it manualy crating GRanges object from dataframe 
  library(dplyr)
  library(readr)
  library(annotatr)
  library(BSgenome.Hsapiens.UCSC.hg19) # sequence of the hg19 genome
  library(ggplot2)
  library(HelloRanges)
  library(reshape2)
  library(yarrr)
  library(GenomeInfoDb) # Contains data and functions that define and allow translation between different chromosome sequence naming conventions (e.g., "chr1" versus "1")

#------------------------------------------
# Set Path to files
#------------------------------------------

  # set work directly on R server for increased speed
    setwd("/home/regmani@ad.ucl.ac.uk/data/EpiCapture/")
    setwd("/mnt/254b78b9-76b4-422d-84b1-cc632bff60f7/regmani/Epicapture")

    list.files()

    getwd() 

  # set resource and RESULTS directory
    Epicapture <- "/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/"
    list.files(Epicapture)

    RESULTS <- paste0(Epicapture,"RESULTS/")
    list.files(RESULTS)

    BEDfiles <- paste0(Epicapture,"BEDfiles/")
    list.files(BEDfiles)

      path2Agilent_BED <- paste0(BEDfiles, "AgilentBEDfiles/S03770311/")
      path2Roche_BED <- paste0(BEDfiles,"RocheBEDfiles/")
      path2Illumina_BED <- paste0(BEDfiles,"IlluminaBEDfiles/")
      path2RRBS_BED <- paste0(BEDfiles,"RRBS_in_silico_BEDfile/")


  #   # get working directory and paths from script "Epicapture_Sample names and Paths.R"
  #     platforms <- paste0(Epicapture,"PLATFORMS/")

  #     path2Agilent <- paste0(platforms, "Agilent_SureSelect/")
  #     path2Illumina <- paste0(platforms, "Illumina_EPIC/")
  #     path2Roche <-  paste0(platforms, "Roche_Nimblegen/")
  #     path2Diagenode <-  paste0(platforms, "Diagenode_RRBS/")
  #     path2Nugen <-  paste0(platforms, "NuGen_RRBS/")

    # # Path to each platform
    # Agilent_bedcov <- paste0(path2Agilent, "bedGraphs_Coverage/")
    # Illumina_bedcov <- paste0(path2Illumina, "bedGraphs_Coverage/")
    # Roche_bedcov <-  paste0(path2Roche, "bedGraphs_Coverage/")
    # Diagenode_bedcov <-  paste0(path2Diagenode, "bedGraphs_Coverage/")
    # Nugen_bedcov <-  paste0(path2Nugen, "bedGraphs_Coverage/")

    # # Platform prefix
    # agilent.prefix <- "Agilent_"
    # illumina.prefix <- "Illumina_"
    # roche.prefix <- "Roche_"
    # diagenode.prefix <- "Diagenode_"
    # nugen.prefix <- "NuGEN_" 

    # # The bedGraph output (optional) looks like this (tab-delimited, 0-based start, 1-based end coords):
    # bedgraph.sufix <- "_R1_val_1_trimmed_bismark_bt2_pe_n6dupsRemoved_NameSorted.bedGraph"


#---------------------------------------------------
# Load annotation files and MethylKit GRanges files:
#---------------------------------------------------

    # Load design files:
    load(file=paste0(BEDfiles,"Platform_Design_GRanges.RData")) 

    # Load feature annotations:
    load(file=paste0(BEDfiles,"FeatureAnnotations/Genes_annotations_GRangesList.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/extended_annotations_GRangesList.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/CpG_annottations_GRangesList.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/CpG_annottations_GRanges.RData"))
    load(file=paste0(BEDfiles,"FeatureAnnotations/Enhancers_annotations_GRanges.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/lncRNA_annotations_GRanges.RData")) 
    load(file=paste0(BEDfiles,"FeatureAnnotations/ChromHMM-byState_GRanges.RData")) 
    # load(file=paste0(BEDfiles, "FeatureAnnotations/CpG_locations_hg19_GRanges.RData"))
    # load(file=paste0(BEDfiles, "FeatureAnnotations/all_cpgs.RData"))
    # load(file=paste0(BEDfiles,"FeatureAnnotations/hg19_CpGsites_GRanges.RData"))
    # load(file=paste0(BEDfiles,"FeatureAnnotations/UCSC_CGI_GRanges.RData"))

#-----------------------------------------------------------------
# Load experimentaly covered GRangesList objects for each platform
#-----------------------------------------------------------------

  #   # Load GRangesList objects 
  #       # grl.PLATFORM was generated from my.DB.PLATFORM object using MethylKit 
  #       # CpGs  covered at min 10x 
  #       # CpGs at both strands included! start=end and separated by strand info!

  #     load(file=paste0(RESULTS,"5_MethylKit/","grl.Agilent.RData"))
  #     load(file=paste0(RESULTS,"5_MethylKit/","grl.Roche.RData"))
  #     load(file=paste0(RESULTS,"5_MethylKit/","grl.Illumina.RData"))
  #     load(file=paste0(RESULTS,"5_MethylKit/","grl.Diagenode.RData"))
  #     load(file=paste0(RESULTS,"5_MethylKit/","grl.Nugen.RData"))

  #   # Roche_Hela_1 has failed so exclud it form analysis
  #   grl.Roche <-  grl.Roche[-6]

  # #-------------------------------------------------------------
  # # Find union of all covered regions - max breadth of coverage:
  # #-------------------------------------------------------------
  #   reduce(unlist(grl.Agilent)) # this is the union

  #   # get number of all total CpGs covered - union
  #   gr.Agilent<- unlist(grl.Agilent)
  #   gr.Roche<- unlist(grl.Roche)
  #   gr.Illumina<- unlist(grl.Illumina)
  #   gr.Diagenode<- unlist(grl.Diagenode)
  #   gr.Nugen<- unlist(grl.Nugen)


  #   sum(width(reduce(gr.Agilent)))
  #   sum(width(reduce(gr.Roche)))
  #   sum(width(reduce(gr.Illumina)))
  #   sum(width(reduce(gr.Diagenode)))
  #   sum(width(reduce(gr.Nugen)))

  # #------------------------------------------------
  # # Find overlap (intersection) between several GRanges objects:
  # #-----------------------------------------------

  # # Find intersection of all CpGs covered by each platform - 
  #   # using MethylKit's functionality to make a merged dataset and then convert to GRanges
  #     meth.Agilent=unite(my.DB.Agilent, destrand=TRUE)
  #     meth.Roche=unite(my.DB.Roche, destrand=TRUE)
  #     meth.Illumina=unite(my.DB.Illumina, destrand=TRUE)
  #     meth.Diagenode=unite(my.DB.Diagenode, destrand=TRUE)
  #     meth.Nugen=unite(my.DB.Nugen, destrand=TRUE)

  #     save(meth.Agilent, meth.Roche, meth.Illumina, meth.Diagenode,meth.Nugen, file="meth.PLATFORM.RData")

  #   # convert to GRanges

  #   gr.common.Agilent <- as(meth.Agilent, "GRanges")
  #   gr.common.Roche <- as(meth.Roche, "GRanges")
  #   gr.common.Illumina <- as(meth.Illumina, "GRanges")
  #   gr.common.Diagenode <- as(meth.Diagenode, "GRanges")
  #   gr.common.Nugen <- as(meth.Nugen, "GRanges")

  #   # get number of all common CpGs covered - intersection
  #   length(gr.common.Agilent)
  #   length(gr.common.Roche)
  #   length(gr.common.Illumina)
  #   length(gr.common.Diagenode)
  #   length(gr.common.Nugen)



    # returns an object of overlappingPeaks, which contains there elements: venn_cnt, peaklist (a list of overlapping peaks or unique peaks), and overlappingPeaks (a list of data frame consists of the annotation of all the overlapping peaks).
    # ol <- findOverlapsOfPeaks(grl.Agilent) # max 5 peak lists!
    # head(ol)


#==============================
# Sample list 
#==============================

  all <- c("Ref-gDNA-500-1", 
                    "Ref-gDNA-500-2",
                    "Ref-gDNA-recommended-1",
                    "Ref-gDNA-recommended-2",
                    "Coriell-NA12878-K12-1",
                    "Coriell-NA12878-K12-2",
                    "Hela-1",
                    "Hela-2",
                    "ZYMO-FM",
                    "ZYMO-UM",
                    "DNAm-5pct",
                    "DNAm-10pct",
                    "T24",
                    "253J",
                    "RT112",
                    "RT112-CP")


  tech_repl <- c("Ref-gDNA-500-1", 
                "Ref-gDNA-500-2",
                "Ref-gDNA-recommended-1",
                "Ref-gDNA-recommended-2",
                "Coriell-NA12878-K12-1",
                "Coriell-NA12878-K12-2",
                "Hela-1",
                "Hela-2")



#--------------------------------------------------------------
# Import bedgraph files from each platform:
#--------------------------------------------------------------

    # to import bedgraph files use rtracklayer:
    # imp = import(f1, format="bedGraph")



    #! bedGraph doesn't have covergae info for filtering! Better use cov.files they also have BED like info but can be filtered:
    # Try imorting with Methyl Kit and use those files for annotation later :

    # # List files:
    # print(samples.agilent <- list.files(path=path2Agilent, pattern="bismark.cov"))
    # print(samples.roche <- list.files(path=path2Roche, pattern="bismark.cov"))
    # print(samples.illumina <- list.files(path=path2Illumina, pattern="bismark.cov"))
    # print(samples.diagenode <- list.files(path=path2Diagenode, pattern="bismark.cov"))  
    # print(samples.nugen <- list.files(path=path2Nugen, pattern="bismark.cov"))

    # import function from rtracklayer, which can load pretty much any kind of genomic data into the appropriate typevof Bioconductor object
    #x <- import.bedGraph(paste0(path2Agilent,"/",samples.agilent[1]))

    # fread()

    #------------
    # Agilent
    #------------

    # # Sample names:
    # samples_names.agilent <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph","",samples.agilent)

    # # Import multiple bedGraph samples:
    # grl <-GRangesList()
    # sample.list.agilent <- list()
    # for (i in length(samples.agilent) ) {
    #   sample.list.agilent[i] <- fread(paste0(path2Agilent,"/",samples.agilent[i]))
    # }

    # # Coerce to GrangesList object
    # grl.Agilent <- GRangesList(file.list)
    # class(grl.Agilent)
    # names(grl.Agilent)

#==============================
# Import merged CpG coverage file
#==============================
  #--------------------------------------
  # Sample list
  #--------------------------------------

  # samples <- list.files(pattern="bismark.cov") # Total number of CpG sites for both strands independently

  # Total number of CpG sites - output from coverage2cytosine with information from top and bottom strand merged into one
  samples.agilent <- list.files(path=Agilent_bedcov, pattern="CpG_evidence.cov") 
  samples.roche <- list.files(path=Roche_bedcov, pattern="CpG_evidence.cov") 
  samples.illumina <- list.files(path=Illumina_bedcov, pattern="CpG_evidence.cov") 
  samples.diagenode <- list.files(path=Diagenode_bedcov, pattern="CpG_evidence.cov") 
  samples.nugen <- list.files(path=Nugen_bedcov, pattern="CpG_evidence.cov") 

  #--------------------------------------
  # read samples 
  #--------------------------------------

    # make a list object
    list.Agilent <- list()
    list.Roche <- list()
    list.Illumina <- list()
    list.Diagenode <- list()
    list.Nugen <- list()

    #test
    # for (i in seq_along(samples.agilent)) {
    #   print(samples.agilent[i])
    # }

    # list.Agilent[[1]] <- fread(paste0(Agilent_bedcov,samples.agilent[1]))
    # names(list.Agilent[[1]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
  
  # read merged cov files
    for (i in seq_along(samples.agilent)) {
        list.Agilent[[i]] <- fread(paste0(Agilent_bedcov,samples.agilent[i]))
        names(list.Agilent[[i]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
      }
    for (i in seq_along(samples.roche)) {
        list.Roche[[i]] <- fread(paste0(Roche_bedcov,samples.roche[i]))
        names(list.Roche[[i]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
      }
    for (i in seq_along(samples.illumina)) {
        list.Illumina[[i]] <- fread(paste0(Illumina_bedcov,samples.illumina[i]))
        names(list.Illumina[[i]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
      }
    for (i in seq_along(samples.diagenode)) {
        list.Diagenode[[i]] <- fread(paste0(Diagenode_bedcov,samples.diagenode[i]))
        names(list.Diagenode[[i]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
      }
    for (i in seq_along(samples.nugen)) {
        list.Nugen[[i]] <- fread(paste0(Nugen_bedcov,samples.nugen[i]))
        names(list.Nugen[[i]]) <- c("chr", "start", "end", "pct_meth", "count_M", "count_UM")
      }
  # make new sample names
        samples.agilent<- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.agilent)
        samples.roche <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.roche)
        samples.illumina <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.illumina)
        samples.diagenode <- gsub("_R1_val_1_bismark_bt2_pe.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.diagenode)
        samples.nugen <- gsub("_R1_val_1_trimmed_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.nugen)
       "R1_val_1_trimmed.fq.gz_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.CpG_report.merged_CpG_evidence.cov"
        samples.nugen <- gsub("_R1_val_1_trimmed.fq.gz_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.CpG_report.merged_CpG_evidence.cov", "", samples.nugen)

        samples.diagenode  <- gsub("Coriell-NA12878", "Coriell-NA12878-K12", samples.diagenode)  
        samples.nugen <- gsub("Nugen_", "", samples.nugen)
        samples.nugen <- gsub("_100M", "", samples.nugen)
        samples.nugen <- gsub("Ref.", "Ref-", samples.nugen)


  names(list.Agilent) <- samples.agilent
  names(list.Roche) <- samples.roche
  names(list.Illumina) <- samples.illumina
  names(list.Diagenode) <- samples.diagenode
  names(list.Nugen) <- samples.nugen


  # add a coverage column
  for (i in seq_along(list.Agilent)){
    list.Agilent[[i]][["cov"]] <- list.Agilent[[i]][["count_M"]] + list.Agilent[[i]][["count_UM"]]
  }
  for (i in seq_along(list.Roche)){
    list.Roche[[i]][["cov"]] <- list.Roche[[i]][["count_M"]] + list.Roche[[i]][["count_UM"]]
  }
  for (i in seq_along(list.Illumina)){
    list.Illumina[[i]][["cov"]] <- list.Illumina[[i]][["count_M"]] + list.Illumina[[i]][["count_UM"]]
  }
  for (i in seq_along(list.Diagenode)){
    list.Diagenode[[i]][["cov"]] <- list.Diagenode[[i]][["count_M"]] + list.Diagenode[[i]][["count_UM"]]
  }
  for (i in seq_along(list.Nugen)){
    list.Nugen[[i]][["cov"]] <- list.Nugen[[i]][["count_M"]] + list.Nugen[[i]][["count_UM"]]
  }

  # save 
  save(list.Agilent, list.Roche, list.Illumina, list.Diagenode, list.Nugen, file="mrg_CpGcoverage_list.PLATFORM.RData")

  load(file="mrg_CpGcoverage_list.PLATFORM.RData")

#===================================================
# filter by coverage >10x and convert to GRangesList
#===================================================

  #-----------------------------
  # Filter by coverage >=10x
  #-----------------------------
      # test
      # list.Nugen[[1]] %>% filter(cov >=10)

    for (i in seq_along(list.Agilent)){
      list.Agilent[[i]]<-  list.Agilent[[i]] %>% filter(cov >=10)}

    for (i in seq_along(list.Roche)){
      list.Roche[[i]]<-  list.Roche[[i]] %>% filter(cov >=10)}

    for (i in seq_along(list.Illumina)){
      list.Illumina[[i]]<-  list.Illumina[[i]] %>% filter(cov >=10)}

    for (i in seq_along(list.Diagenode)){
      list.Diagenode[[i]]<-  list.Diagenode[[i]] %>% filter(cov >=10)}

    for (i in seq_along(list.Nugen)){
      list.Nugen[[i]]<-  list.Nugen[[i]] %>% filter(cov >=10)}

  # save 
  save(list.Agilent, list.Roche, list.Illumina, list.Diagenode, list.Nugen, file="mrg_CpGcoverage_list.PLATFORM.RData")

  load(file="mrg_CpGcoverage_list_10x.PLATFORM.RData")

  #-----------------------------
  # Convert to GRangesList
  #-----------------------------

  makeGRangesFromDataFrame(df,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

  grl.Agilent <- GRangesList()
  for (i in names(list.Agilent)) {
    grl.Agilent[[i]] <- makeGRangesFromDataFrame(list.Agilent[[i]],
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)
  }

  grl.Roche <- GRangesList()
  for (i in names(list.Roche)) {
    grl.Roche[[i]] <- makeGRangesFromDataFrame(list.Roche[[i]],
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)
  }
  grl.Illumina <- GRangesList()
  for (i in names(list.Illumina)) {
    grl.Illumina[[i]] <- makeGRangesFromDataFrame(list.Illumina[[i]],
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)
  }
  grl.Diagenode <- GRangesList()
  for (i in names(list.Diagenode)) {
    grl.Diagenode[[i]] <- makeGRangesFromDataFrame(list.Diagenode[[i]],
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)
  }
  grl.Nugen <- GRangesList()
  for (i in names(list.Nugen)) {
    grl.Nugen[[i]] <- makeGRangesFromDataFrame(list.Nugen[[i]],
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field=c("end", "stop"),
                         starts.in.df.are.0based=FALSE)
  }

    # Roche_Hela_1 has failed so exclud it form analysis
    grl.Roche <-  grl.Roche[-6]

    # add Platform names to Nugen
        names(grl.Nugen) <- paste0("Nugen_", names(grl.Nugen))


  # Get number of CpGs covered - lapply and sapply return a dataframe and a list object respectively!
    sapply(grl.Agilent,length)
    sapply(grl.Roche,length)
    sapply(grl.Illumina,length)
    sapply(grl.Diagenode,length)
    sapply(grl.Nugen,length)

  # add "chr" to seqnames - ALWAYS USE ENDOAPPLY() TO RETURN A GRANGESLIST OBJECT
    BiocManager::install("diffloop")
    library(diffloop)

    grl.Agilent <- endoapply(grl.Agilent,addchr) 
    grl.Roche <- endoapply(grl.Roche,addchr) 
    grl.Illumina <- endoapply(grl.Illumina,addchr) 
    grl.Diagenode <- endoapply(grl.Diagenode,addchr) 
    grl.Nugen <- endoapply(grl.Nugen,addchr) 

  # save 
  save(grl.Agilent, grl.Roche, grl.Illumina, grl.Diagenode, grl.Nugen, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/mrg_CpGcoverage_grl.PLATFORM.RData"))

  load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/mrg_CpGcoverage_grl.PLATFORM.RData"))

  #------------------------------------------------
  # find union - all CpGs covered by each platform:
  #------------------------------------------------

    # using reduce results in consecutive CpGs being merged into one region

    ugr.Agilent <- reduce(unlist(grl.Agilent))
    length(ugr.Agilent)
    head(ugr.Agilent)

    ugr.Roche <- reduce(unlist(grl.Roche))
    length(ugr.Roche)
    head(ugr.Roche)

    ugr.Illumina <- reduce(unlist(grl.Illumina))
    length(ugr.Illumina)
    head(ugr.Illumina)
    
    ugr.Diagenode <- reduce(unlist(grl.Diagenode))
    length(ugr.Diagenode)
    head(ugr.Diagenode)
    
    ugr.Nugen <- reduce(unlist(grl.Nugen))
    length(ugr.Nugen)
    head(ugr.Nugen)

  #------------------------------------------------------------
  # find intersection - commonly covered CpGs  by each platform:
  #-------------------------------------------------------------

  # subsetByOverlaps method simply subsets the first GRanges object to include only those that overlap the second.
      subsetByOverlaps(gr, g)

  # use intersect instead

    for (i in names(grl.Agilent)) {
    cgr.Agilent <- intersect(ugr.Agilent, grl.Agilent[[i]])}
    length(cgr.Agilent)  

    for (i in names(grl.Roche)) {
    cgr.Roche <- intersect(ugr.Roche, grl.Roche[[i]])}
    length(cgr.Roche)  
    
    for (i in names(grl.Illumina)) {
    cgr.Illumina <- intersect(ugr.Illumina, grl.Illumina[[i]])}
    length(cgr.Illumina)  
    
    for (i in names(grl.Diagenode)) {
    cgr.Diagenode <- intersect(ugr.Diagenode, grl.Diagenode[[i]])}
    length(cgr.Diagenode)  
    
    for (i in names(grl.Nugen)) {
    cgr.Nugen <- intersect(ugr.Nugen, grl.Nugen[[i]])}
    length(cgr.Nugen)  

  # # add "chr" to seqnames to make it compatible with annotation files
    # #  addchr takes a loops object or GRanges object and simply adds 'chr' to seqnames

    #   BiocManager::install("diffloop")
    #   library(diffloop)
    #   ## S4 method for signature 'GRanges'
    #   ugr.Agilent <- addchr(ugr.Agilent)
    #   ugr.Roche <- addchr(ugr.Roche)
    #   ugr.Illumina <- addchr(ugr.Illumina)
    #   ugr.Diagenode <- addchr(ugr.Diagenode)
    #   ugr.Nugen <- addchr(ugr.Nugen)

    #   cgr.Agilent <- addchr(cgr.Agilent)
    #   cgr.Roche <- addchr(cgr.Roche)
    #   cgr.Illumina <- addchr(cgr.Illumina)
    #   cgr.Diagenode <- addchr(cgr.Diagenode)
    #   cgr.Nugen <- addchr(cgr.Nugen)

  # save objects
    save(ugr.Agilent,ugr.Roche,ugr.Illumina,ugr.Diagenode, ugr.Nugen, cgr.Agilent,cgr.Roche, cgr.Illumina, cgr.Diagenode, cgr.Nugen,  file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/Total_n_common_CpGs.RData"))

    load(file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/Total_n_common_CpGs.RData"))


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
    x <- getSeq(genome, reduce(ans)) # # The return values are DNAStringSetList objects
    sum(vcountPattern("CG", x)) # count number of CG occurancesin the sequence
  }

  countOverlaps()

  # jaccard statistics
  intersects <- intersect(gr_a, gr_b, ignore.strand = TRUE)
  intersection <- sum(width(intersects))
  union <- sum(width(union(gr_a, gr_b, ignore.strand = TRUE)))
  ans <- DataFrame(intersection, union, jaccard = intersection/union, n_intersections = length(intersects))


  findOverlapPairs #convenience for creating a  Pairs object that matches up the overlapping ranges: 
 
  # Another way of getting at overlap information is to use %over% which returns a logical vector of which ranges in the first argument overlapped any ranges in the second.
  gr1 %over% gr2
  ## [1] FALSE FALSE  TRUE  TRUE FALSE
  gr1[gr1 %over% gr2]

  # Note that both of these are strand-specific, although findOverlaps has an ignore.strand option.


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

    scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 
    scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 


#===========================================================================
# 1. Annotate Union of all CpGs coverd by each platforms for each Feature type:
#===========================================================================


#===========================
# Annotate Genomic Features:
#===========================

  #----------------------------------------------------------------------------------
  # Calculate overlap - intersection between platform design file and annotation file:
  #----------------------------------------------------------------------------------

    u_df_genes <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

    for (i in (names(annotations_genes.grl))) {
      u_df_genes[i,"Agilent"] <- findOverlapSize(ugr.Agilent, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      u_df_genes[i,"Roche"] <- findOverlapSize(ugr.Roche, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      u_df_genes[i,"Illumina"] <- findOverlapSize(ugr.Illumina, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      u_df_genes[i,"Diagenode"] <- findOverlapSize(ugr.Diagenode, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      u_df_genes[i,"Nugen"] <- findOverlapSize(ugr.Nugen, unlist(annotations_genes.grl[i]))}


    u_df_genes$Total <- sum(width(reduce(annotations_genes.grl)))

    head(u_df_genes)
    save(u_df_genes, file=paste0(RESULTS,"1_PlatformDesignDiferences/u_df_genes_size.RData"))
    write.table(u_df_genes, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/union_genic_size.txt"))

    load(file=paste0(RESULTS,"1_PlatformDesignDiferences/u_df_genes_size.RData"))

    # Transform dataframe to long format for plotting:
      u_df_genes_long <- melt(t(u_df_genes[,-1]))
      names(u_df_genes_long) <- c("Platform", "Feature", "Size_bp" )
      head(u_df_genes_long)
      u_df_genes_long$Size_Mb <- round(u_df_genes_long$Size_bp/1000000,2)

    #------------------------------------------
    # Plot Feature size by Platform - barplots:
    #------------------------------------------

      gene_annots <- ggplot(data=u_df_genes_long, aes(x=Feature, y=Size_Mb, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_brewer(palette="Set1") +
        # scale_color_manual(values=piratepal("appletv")[1:5]) + # select colors for dots
        # scale_fill_manual(values=col.platforms[1:5])+
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
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Genic_ByFeature_size.pdf"))

  #--------------------------------------------
  # find # of CpGs in each feature category
  #--------------------------------------------

    # Count # of CpG in each annotation category for UNION

        u_df_genes_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_genes.grl))) {
        u_df_genes_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_genes.grl[i]))}

        for (i in (names(annotations_genes.grl))) {
        u_df_genes_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_genes.grl[i]))}

        for (i in (names(annotations_genes.grl))) {
        u_df_genes_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_genes.grl[i]))}

        for (i in (names(annotations_genes.grl))) {
        u_df_genes_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_genes.grl[i]))}

        for (i in (names(annotations_genes.grl))) {
        u_df_genes_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_genes.grl[i]))}

        for (i in seq_along(names(annotations_genes.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_genes.grl[i])))
            u_df_genes_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(u_df_genes_cpg)
        save(u_df_genes_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/u_df_genes_cpg.RData"))
        write.table(u_df_genes, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/union_genic_cpg.txt"))

    # Transform dataframe to long format for plotting:
        u_df_genes_cpg_long <- melt(t(u_df_genes_cpg[,-1]))
        names(u_df_genes_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(u_df_genes_cpg_long)
        u_df_genes_cpg_long$No_CpGs_M <- round(u_df_genes_cpg_long$No_CpGs/1000000,2)


    #----------------------------------------
    # Plot Features by Platform - barplots:
    #----------------------------------------

        gene_annots <- ggplot(data=u_df_genes_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of total CpGs per genomic feature by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                                )) + # Change the name of items
        theme_classic()
        gene_annots
        ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Genic_ByFeature_NoCpGs.pdf"))

    #----------------------------------------
    # Stacked Plot NoCpGs per Feature by Platform
    #----------------------------------------

        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=u_df_genes_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
            geom_bar(stat = "identity", colour="black", size=0.25)+
            ylab("# CpGs") + 
            ggtitle("Number of CpGs per genomic feature by platform") + 
            scale_fill_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                    "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns","hg19_genes_promoters"="promoters"
                                    )) + # Change the name of items
             theme_classic()
        gene_annots
        ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_StakedBarplot_byPlatform_NoCpGs.pdf"))


    #----------------------------------------------
    # Percent CpGs covered per Feature by Platform
    #----------------------------------------------

      # Calculate percentages of total features covered by each platform
        head(u_df_genes_cpg)
        u_df_genes_cpg_pct <- apply(u_df_genes_cpg, 2, function(x) x/u_df_genes_cpg$Total)
        head(u_df_genes_cpg_pct)

      # Plot percentages:
      u_df_genes_cpg_pct_long <- melt(t(u_df_genes_cpg_pct[,-1]))
      names(u_df_genes_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

      gene_annots_pct <- ggplot(data=u_df_genes_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
        #          position = position_dodge(0.9), size=3.5)+
        scale_fill_manual(values=col.platforms) +
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
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Genic_ByFeature_pctCpGs.pdf"))


    #----------------------------------------
    # Piechart % CpGs in features by Platform
    #----------------------------------------

      # Compute percentages of each features covered by each platform
        u_df_genes_cpg_pct_CpGs <- apply(u_df_genes_cpg, 2, function(x) pct=x/sum(x))
        head(u_df_genes_cpg_pct_CpGs)

      # calculate percentage:
          u_df_genes_cpg_pct_CpGs_long <- melt(t(u_df_genes_cpg_pct_CpGs[,-1]))
          names(u_df_genes_cpg_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
          head(u_df_genes_cpg_pct_CpGs_long)

      gene_annots <- ggplot(data=u_df_genes_cpg_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Piechart_byPlatform_pctCpGs.pdf"))



#================
# Annotate CpGs:
#================

  # Inspect loaded object:
  head(annotations_CpGs.grl)
  updateObject(annotations_CpGs.grl,  verbose=TRUE)
  names(annotations_CpGs.grl)


  #--------------------------------------------
  # find # of CpGs in each feature category
  #--------------------------------------------

    # Count # of CpG in each annotation category for UNION

        u_df_CGI_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_CpGs.grl))) {
        u_df_CGI_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        u_df_CGI_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        u_df_CGI_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        u_df_CGI_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        u_df_CGI_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_CpGs.grl[i]))}

        for (i in seq_along(names(annotations_CpGs.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_CpGs.grl[i])))
            u_df_CGI_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(u_df_CGI_cpg)
        save(u_df_CGI_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/u_df_CGI_cpg.RData"))
        write.table(u_df_CGI_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/union_CGI_cpg.txt"))

    # Transform dataframe to long format for plotting:
        u_df_CGI_cpg_long <- melt(t(u_df_CGI_cpg[,-1]))
        names(u_df_CGI_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(u_df_CGI_cpg_long)
        u_df_CGI_cpg_long$No_CpGs_M <- round(u_df_CGI_cpg_long$No_CpGs/1000000,2)


  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=u_df_CGI_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shelves","hg19_cpg_shores",  
                                    "hg19_cpg_inter"), # Change the order of items
                        labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                                    "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_CGI_ByFeature_NoCpGs.pdf"))

  #---------------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #--------------------------------------------

      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=u_df_CGI_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGsper feature category by platform") + 
        scale_fill_discrete(labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves", "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_StakedBarplot_CGI_byPlatform_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(u_df_CGI_cpg)
      u_df_CGI_cpg_pct <- apply(u_df_CGI_cpg, 2, function(x) x/u_df_CGI_cpg$Total)
      head(u_df_CGI_cpg_pct)

    # Plot percentages:
    u_df_CGI_cpg_pct_long <- melt(t(u_df_CGI_cpg_pct[,-1]))
    names(u_df_CGI_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=u_df_CGI_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_CGI_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      u_df_platform_CGI_pct_CpGs <- apply(u_df_CGI_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        u_df_platform_CGI_pct_CpGs_long <- melt(t(u_df_platform_CGI_pct_CpGs[,-1]))
        names(u_df_platform_CGI_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(u_df_platform_CGI_pct_CpGs_long)

    gene_annots <- ggplot(data=u_df_platform_CGI_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Piechart_CGI_byPlatform_pctCpGs.pdf"))

#===============================
# Annotate Enhancers and lncRNA:
#===============================

 # Inspect loaded object:
        head(annotations_enhancers_fantom.gr)
        head(annotations_lncRNA.gr)
        annotations_regul.grl <- GRangesList(enhancers=annotations_enhancers_fantom.gr,lncRNA=annotations_lncRNA.gr)

    # Count # of CpG in each annotation category for UNION

        u_df_regul_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_regul.grl))) {
        u_df_regul_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        u_df_regul_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        u_df_regul_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        u_df_regul_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        u_df_regul_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_regul.grl[i]))}

        for (i in seq_along(names(annotations_regul.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_regul.grl[i])))
            u_df_regul_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(u_df_regul_cpg)
        save(u_df_regul_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/u_df_regul_cpg.RData"))
        write.table(u_df_regul_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/union_regul_cpg.txt"))

    # Transform dataframe to long format for plotting:
       u_df_regul_cpg_long <- melt(t(u_df_regul_cpg[,-1]))
        names(u_df_regul_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(u_df_regul_cpg_long)
        u_df_regul_cpg_long$No_CpGs_M <- round(u_df_regul_cpg_long$No_CpGs/1000000,2)

  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=u_df_regul_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Regul_ByFeature_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(u_df_regul_cpg)
      u_df_regul_cpg_pct <- apply(u_df_regul_cpg, 2, function(x) x/u_df_regul_cpg$Total)
      head(u_df_regul_cpg_pct)

    # Plot percentages:
    u_df_regul_cpg_pct_long <- melt(t(u_df_regul_cpg_pct[,-1]))
    names(u_df_regul_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=u_df_regul_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
      theme_classic()
    gene_annots_pct
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Regul_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      df_platform_regul_pct_CpGs <- apply(u_df_regul_cpg, 2, function(x) pct=x/sum(x))

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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Piechart_Regul_byPlatform_pctCpGs.pdf"))

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

    # Count # of CpG in each annotation category for UNION

        u_df_chrom_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_chromHMM.grl))) {
        u_df_chrom_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        u_df_chrom_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        u_df_chrom_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        u_df_chrom_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        u_df_chrom_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_chromHMM.grl[i]))}

        for (i in seq_along(names(annotations_chromHMM.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_chromHMM.grl[i])))
            u_df_chrom_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(u_df_chrom_cpg)
        save(u_df_chrom_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/u_df_chrom_cpg.RData"))
        write.table(u_df_chrom_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/union_chrom_cpg.txt"))

    # Transform dataframe to long format for plotting:
       u_df_chrom_cpg_long <- melt(t(u_df_chrom_cpg[,-1]))
        names(u_df_chrom_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(u_df_chrom_cpg_long)
        u_df_chrom_cpg_long$No_CpGs_M <- round(u_df_chrom_cpg_long$No_CpGs/1000000,2)


  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=u_df_chrom_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Chrom_ByFeature_NoCpGs.pdf"))

  #----------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #----------------------------------------
        feature.cols <- brewer.pal(n=11, "Set3")
      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=u_df_chrom_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_StakedBarplot_Chrom_byPlatform_NoCpGs.pdf"))


  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(u_df_chrom_cpg)
      u_df_chrom_cpg_pct <- apply(u_df_chrom_cpg, 2, function(x) x/u_df_chrom_cpg$Total)
      head(u_df_chrom_cpg_pct)

    # Plot percentages:
    u_df_chrom_cpg_pct_long <- melt(t(u_df_chrom_cpg_pct[,-1]))
    names(u_df_chrom_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=u_df_chrom_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
    gene_annots_pct
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Barplot_Chrom_ByFeature_pctCpGs.pdf"))


  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      u_df_platform_Chrom_pct_CpGs <- apply(u_df_chrom_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        u_df_platform_Chrom_pct_CpGs_long <- melt(t(u_df_platform_Chrom_pct_CpGs[,-1]))
        names(u_df_platform_Chrom_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(u_df_platform_Chrom_pct_CpGs_long)

    gene_annots <- ggplot(data=u_df_platform_Chrom_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/u_Piechart_Chrom_byPlatform_pctCpGs.pdf"))












#=======================================================================================
# 2. Annotate Intersection (common) CpGs coverd by each platforms for each Feature type:
#=======================================================================================

 
  # Count # of CpG in each annotation category for INTERSECTION (common CpGs)

    c_df_genes_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

    for (i in (names(annotations_genes.grl))) {
      c_df_genes_cpg[i,"Agilent"] <- findOverlapCpGCount(cgr.Agilent, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      c_df_genes_cpg[i,"Roche"] <- findOverlapCpGCount(cgr.Roche, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      c_df_genes_cpg[i,"Illumina"] <- findOverlapCpGCount(cgr.Illumina, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      c_df_genes_cpg[i,"Diagenode"] <- findOverlapCpGCount(cgr.Diagenode, unlist(annotations_genes.grl[i]))}

    for (i in (names(annotations_genes.grl))) {
      c_df_genes_cpg[i,"Nugen"] <- findOverlapCpGCount(cgr.Nugen, unlist(annotations_genes.grl[i]))}

    for (i in seq_along(names(annotations_genes.grl))){ 
        x <- getSeq(genome, reduce(unlist(annotations_genes.grl[i])))
        c_df_genes_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

    head(c_df_genes_cpg)
    save(c_df_genes_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/c_df_genes_cpg.RData"))
    write.table(c_df_genes_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/common_genic_cpg.txt"))

  # Transform dataframe to long format for plotting:
    c_df_genes_cpg_long <- melt(t(c_df_genes_cpg[,-1]))
    names(c_df_genes_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
    head(c_df_genes_cpg_long)
    c_df_genes_cpg_long$No_CpGs_M <- round(c_df_genes_cpg_long$No_CpGs/1000000,2)

    #----------------------------------------
    # Plot Features by Platform - barplots:
    #----------------------------------------

        gene_annots <- ggplot(data=c_df_genes_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of total CpGs per genomic feature by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters",  
                                "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns", "hg19_genes_promoters"="promoters"
                                )) + # Change the name of items
        theme_classic()
        gene_annots
        ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Genic_ByFeature_NoCpGs.pdf"))

    #----------------------------------------
    # Stacked Plot NoCpGs per Feature by Platform
    #----------------------------------------
        feature.cols <- brewer.pal(n = 8, name = "Dark2")
        scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

        gene_annots <- ggplot(data=c_df_genes_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
            geom_bar(stat = "identity", colour="black", size=0.25)+
            ylab("# CpGs") + 
            ggtitle("Number of CpGs per genomic feature by platform") + 
            scale_fill_discrete(limits=c("hg19_genes_1to5kb","hg19_genes_5UTRs","hg19_genes_promoters", "hg19_genes_exons", "hg19_genes_introns","hg19_genes_3UTRs" ), # Change the order of items
                            labels=c("hg19_genes_1to5kb"="1 to 5kb", "hg19_genes_3UTRs"="3'UTRs",
                                    "hg19_genes_5UTRs"="5'UTRs", "hg19_genes_exons"="exons", "hg19_genes_introns"="introns","hg19_genes_promoters"="promoters"
                                    )) + # Change the name of items
             theme_classic()
        gene_annots
        ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_StakedBarplot_byPlatform_NoCpGs.pdf"))


    #----------------------------------------------
    # Percent CpGs covered per Feature by Platform
    #----------------------------------------------

      # Calculate percentages of total features covered by each platform
        head(c_df_genes_cpg)
        c_df_genes_cpg_pct <- apply(c_df_genes_cpg, 2, function(x) x/c_df_genes_cpg$Total)
        head(c_df_genes_cpg_pct)

      # Plot percentages:
      c_df_genes_cpg_pct_long <- melt(t(c_df_genes_cpg_pct[,-1]))
      names(c_df_genes_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

      gene_annots_pct <- ggplot(data=c_df_genes_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
        #          position = position_dodge(0.9), size=3.5)+
        scale_fill_manual(values=col.platforms) +
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
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Genic_ByFeature_pctCpGs.pdf"))


    #----------------------------------------
    # Piechart % CpGs in features by Platform
    #----------------------------------------

      # Compute percentages of each features covered by each platform
        c_df_genes_cpg_pct_CpGs <- apply(c_df_genes_cpg, 2, function(x) pct=x/sum(x))
        head(c_df_genes_cpg_pct_CpGs)

      # calculate percentage:
          c_df_genes_cpg_pct_CpGs_long <- melt(t(c_df_genes_cpg_pct_CpGs[,-1]))
          names(c_df_genes_cpg_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
          head(c_df_genes_cpg_pct_CpGs_long)

      gene_annots <- ggplot(data=c_df_genes_cpg_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Piechart_byPlatform_pctCpGs.pdf"))



#================
# Annotate CpGs:
#================

  # Inspect loaded object:
  head(annotations_CpGs.grl)
  updateObject(annotations_CpGs.grl,  verbose=TRUE)
  names(annotations_CpGs.grl)


  #--------------------------------------------
  # find # of CpGs in each feature category
  #--------------------------------------------

    # Count # of CpG in each annotation category for UNION

        c_df_CGI_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_CpGs.grl))) {
        c_df_CGI_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        c_df_CGI_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        c_df_CGI_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        c_df_CGI_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_CpGs.grl[i]))}

        for (i in (names(annotations_CpGs.grl))) {
        c_df_CGI_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_CpGs.grl[i]))}

        for (i in seq_along(names(annotations_CpGs.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_CpGs.grl[i])))
            c_df_CGI_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(c_df_CGI_cpg)
        save(c_df_CGI_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/c_df_CGI_cpg.RData"))
        write.table(c_df_CGI_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/common_CGI_cpg.txt"))

    # Transform dataframe to long format for plotting:
        c_df_CGI_cpg_long <- melt(t(c_df_CGI_cpg[,-1]))
        names(c_df_CGI_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(c_df_CGI_cpg_long)
        c_df_CGI_cpg_long$No_CpGs_M <- round(c_df_CGI_cpg_long$No_CpGs/1000000,2)


  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=c_df_CGI_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("hg19_cpg_islands","hg19_cpg_shelves","hg19_cpg_shores",  
                                    "hg19_cpg_inter"), # Change the order of items
                        labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves",
                                    "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_CGI_ByFeature_NoCpGs.pdf"))

  #---------------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #--------------------------------------------

      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=c_df_CGI_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGsper feature category by platform") + 
        scale_fill_discrete(labels=c("hg19_cpg_islands"="CpG Islands", "hg19_cpg_shelves"="shelves", "hg19_cpg_shores"="shores", "hg19_cpg_inter"="inter CGI")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_StakedBarplot_CGI_byPlatform_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(c_df_CGI_cpg)
      c_df_CGI_cpg_pct <- apply(c_df_CGI_cpg, 2, function(x) x/c_df_CGI_cpg$Total)
      head(c_df_CGI_cpg_pct)

    # Plot percentages:
    c_df_CGI_cpg_pct_long <- melt(t(c_df_CGI_cpg_pct[,-1]))
    names(c_df_CGI_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=c_df_CGI_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_CGI_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      c_df_platform_CGI_pct_CpGs <- apply(c_df_CGI_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        c_df_platform_CGI_pct_CpGs_long <- melt(t(c_df_platform_CGI_pct_CpGs[,-1]))
        names(c_df_platform_CGI_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(c_df_platform_CGI_pct_CpGs_long)

    gene_annots <- ggplot(data=c_df_platform_CGI_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Piechart_CGI_byPlatform_pctCpGs.pdf"))

#===============================
# Annotate Enhancers and lncRNA:
#===============================

 # Inspect loaded object:
        head(annotations_enhancers_fantom.gr)
        head(annotations_lncRNA.gr)
        annotations_regul.grl <- GRangesList(enhancers=annotations_enhancers_fantom.gr,lncRNA=annotations_lncRNA.gr)

    # Count # of CpG in each annotation category for common

        c_df_regul_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_regul.grl))) {
        c_df_regul_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        c_df_regul_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        c_df_regul_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        c_df_regul_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_regul.grl[i]))}

        for (i in (names(annotations_regul.grl))) {
        c_df_regul_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_regul.grl[i]))}

        for (i in seq_along(names(annotations_regul.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_regul.grl[i])))
            c_df_regul_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(c_df_regul_cpg)
        save(c_df_regul_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/c_df_regul_cpg.RData"))
        write.table(c_df_regul_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/common_regul_cpg.txt"))

    # Transform dataframe to long format for plotting:
       c_df_regul_cpg_long <- melt(t(c_df_regul_cpg[,-1]))
        names(c_df_regul_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(c_df_regul_cpg_long)
        c_df_regul_cpg_long$No_CpGs_M <- round(c_df_regul_cpg_long$No_CpGs/1000000,2)

  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=c_df_regul_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Regul_ByFeature_NoCpGs.pdf"))

  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(c_df_regul_cpg)
      c_df_regul_cpg_pct <- apply(c_df_regul_cpg, 2, function(x) x/c_df_regul_cpg$Total)
      head(c_df_regul_cpg_pct)

    # Plot percentages:
    c_df_regul_cpg_pct_long <- melt(t(c_df_regul_cpg_pct[,-1]))
    names(c_df_regul_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=c_df_regul_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("enhancers","lncRNA"), # Change the order of items
                          labels=c("enhancers"="FANTOM5 enhancers", "lncRNA"="GENCODE lncRNA")) + # Change the name of items
      theme_classic()
    gene_annots_pct
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Regul_ByFeature_pctCpGs.pdf"))

  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      df_platform_regul_pct_CpGs <- apply(c_df_regul_cpg, 2, function(x) pct=x/sum(x))

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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Piechart_Regul_byPlatform_pctCpGs.pdf"))

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

    # Count # of CpG in each annotation category for UNION

        c_df_chrom_cpg <- data.frame(Total=integer(), Agilent=integer(), Roche=integer(), Illumina=integer(), Diagenode=integer(), Nugen=integer())

        for (i in (names(annotations_chromHMM.grl))) {
        c_df_chrom_cpg[i,"Agilent"] <- findOverlapCpGCount(ugr.Agilent, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        c_df_chrom_cpg[i,"Roche"] <- findOverlapCpGCount(ugr.Roche, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        c_df_chrom_cpg[i,"Illumina"] <- findOverlapCpGCount(ugr.Illumina, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        c_df_chrom_cpg[i,"Diagenode"] <- findOverlapCpGCount(ugr.Diagenode, unlist(annotations_chromHMM.grl[i]))}

        for (i in (names(annotations_chromHMM.grl))) {
        c_df_chrom_cpg[i,"Nugen"] <- findOverlapCpGCount(ugr.Nugen, unlist(annotations_chromHMM.grl[i]))}

        for (i in seq_along(names(annotations_chromHMM.grl))){ 
            x <- getSeq(genome, reduce(unlist(annotations_chromHMM.grl[i])))
            c_df_chrom_cpg[i,"Total"]  <-sum(vcountPattern("CG", x))}

        head(c_df_chrom_cpg)
        save(c_df_chrom_cpg, file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/c_df_chrom_cpg.RData"))
        write.table(c_df_chrom_cpg, sep="\t", file=paste0(RESULTS,"8_ExpCoveredCpGsAnnotation/common_chrom_cpg.txt"))

    # Transform dataframe to long format for plotting:
       c_df_chrom_cpg_long <- melt(t(c_df_chrom_cpg[,-1]))
        names(c_df_chrom_cpg_long) <- c("Platform", "Feature", "No_CpGs" )
        head(c_df_chrom_cpg_long)
        c_df_chrom_cpg_long$No_CpGs_M <- round(c_df_chrom_cpg_long$No_CpGs/1000000,2)


  #----------------------------------------
  # Plot Features by Platform - barplots:
  #----------------------------------------
    # scale_fill_discrete <- function(...) {scale_fill_manual(..., values = platform.cols)} 

      gene_annots <- ggplot(data=c_df_chrom_cpg_long, aes(x=Feature, y=No_CpGs, fill=Platform))  +
        geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
        # geom_text(aes(label=Size_bp), vjust=1.6, color="white", position = position_dodge(0.9), size=2)+
        scale_fill_manual(values=col.platforms) +
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        theme(axis.title.x=element_blank()) +
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Chrom_ByFeature_NoCpGs.pdf"))

  #----------------------------------------
  # Stacked Plot NoCpGs per Feature by Platform
  #----------------------------------------
        feature.cols <- brewer.pal(n=11, "Set3")
      scale_fill_discrete <- function(...) {scale_fill_manual(..., values = feature.cols)} 

      gene_annots <- ggplot(data=c_df_chrom_cpg_long, aes(x=Platform, y=No_CpGs, fill=Feature)) +
        geom_bar(stat = "identity", colour="black", size=0.25)+
        ylab("# CpGs") + 
        ggtitle("Number of CpGs per feature category by platform") + 
        scale_fill_discrete(labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
      gene_annots
      ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_StakedBarplot_Chrom_byPlatform_NoCpGs.pdf"))


  #----------------------------------------------
  # Percent CpGs covered per Feature by Platform
  #----------------------------------------------

    # Calculate percentages of total features covered by each platform
      head(c_df_chrom_cpg)
      c_df_chrom_cpg_pct <- apply(c_df_chrom_cpg, 2, function(x) x/c_df_chrom_cpg$Total)
      head(c_df_chrom_cpg_pct)

    # Plot percentages:
    c_df_chrom_cpg_pct_long <- melt(t(c_df_chrom_cpg_pct[,-1]))
    names(c_df_chrom_cpg_pct_long) <-  c("Platform", "Feature", "PercentCpGs_covered" )

    gene_annots_pct <- ggplot(data=c_df_chrom_cpg_pct_long, aes(x=Feature, y=PercentCpGs_covered, fill=Platform))  +
      geom_bar(stat = "identity", position = "dodge", colour="black", size=0.25)+
      #geom_text(aes(label=Percent_covered), vjust=1.6, color="white",
      #          position = position_dodge(0.9), size=3.5)+
      scale_fill_manual(values=col.platforms) +
      ylab("% CpGs covered") + 
      ggtitle("Percent CpGs covered per features by platform") + 
      theme(axis.title.x=element_blank()) +
      # coord_flip()+ # make horizontal chart
        scale_x_discrete(limits=c("ActivePromoter","WeakPromoter", "PoisedPromoter","StrongEnhancer","WeakEnhancer", "Insulator", "TxnTransition","TxnElongation", "WeakTxn","Repressed", "Heterochrom"), # Change the order of items
          labels=c("ActivePromoter"="Active Promoter", "Heterochrom"="Heterochromatine", "Insulator"="Insulator","PoisedPromoter"="Poised Promoter", "Repressed"="Repressed", "StrongEnhancer"="Strong Enhancer", "TxnElongation"="Txn Elongation", "TxnTransition"="TxnTransition", "WeakEnhancer"="Weak Enhancer", "WeakPromoter"="Weak Promoter","WeakTxn"="Weak Txn")) + # Change the name of items
        theme_classic()
    gene_annots_pct
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Barplot_Chrom_ByFeature_pctCpGs.pdf"))


  #----------------------------------------
  # Piechart % CpGs in features by Platform
  #----------------------------------------

    # Compute percentages of each features covered by each platform
      c_df_platform_Chrom_pct_CpGs <- apply(c_df_chrom_cpg, 2, function(x) pct=x/sum(x))

    # calculate percentage:
        c_df_platform_Chrom_pct_CpGs_long <- melt(t(c_df_platform_Chrom_pct_CpGs[,-1]))
        names(c_df_platform_Chrom_pct_CpGs_long) <- c("Platform", "Feature", "pctCpGs" )
        head(c_df_platform_Chrom_pct_CpGs_long)

    gene_annots <- ggplot(data=c_df_platform_Chrom_pct_CpGs_long, aes(x="", y=pctCpGs, fill=Feature)) +
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
    ggsave(paste0(RESULTS, "8_ExpCoveredCpGsAnnotation/c_Piechart_Chrom_byPlatform_pctCpGs.pdf"))



#=======================================
# Sumarry statistics and  test
#=======================================

 my.tbl <- df_all_CGI_long %>% 
          dplyr::filter(Feature == "CpG Islands") %>%
          dplyr::group_by(Platform) %>%  
          dplyr::summarize(
                          count = n(),
                          Mean = mean(CountCpG, na.rm=TRUE),
                          SD= sd(CountCpG, na.rm=TRUE))


# A tibble: 5 x 2
  Platform      Mean
  <fct>        <dbl>
1 Agilent    803313.
2 Roche      991107 
3 Illumina  1283460.
4 Diagenode 1820686.
5 Nugen     1573965.

# ANOVA test hypotheses:

# Null hypothesis: the means of the different groups are the same
# Alternative hypothesis: At least one sample mean is not equal to the others.

my_data <- df_all_CGI_long %>% 
          dplyr::filter(Feature == "CpG Islands") %>%
          dplyr::group_by(Platform) 

# Compute the analysis of variance
res.aov <- aov(CountCpG ~ Platform, data = my_data)
# Summary of the analysis
summary(res.aov)


# make function

myFun_tbl <- function(df,Filter) {
  my.tbl <- df %>% 
            dplyr::filter(Feature == Filter) %>%
            dplyr::group_by(Platform) %>%  
            dplyr::summarize(
                          count = n(),
                          Mean = mean(CountCpG, na.rm=TRUE),
                          SD= sd(CountCpG, na.rm=TRUE))

  return(my.tbl)
}

myFun_anova <- function(df, Filter){
  my_data <- df %>% 
          dplyr::filter(Feature == Filter) %>%
          dplyr::group_by(Platform) 

  res.aov <- aov(CountCpG ~ Platform, data = my_data)
  summary(res.aov)
}

# run functions on all and save output in a txt file:
myFun_tbl(df=df_all_CGI_long, Filter="CpG Islands")
myFun_anova(df=df_all_CGI_long, Filter="CpG Islands")


myFun_tbl(df=df_all_CGI_long, Filter="shores")
myFun_anova(df=df_all_CGI_long, Filter="shores")

myFun_tbl(df=df_all_CGI_long, Filter="inter CGI")
myFun_anova(df=df_all_CGI_long, Filter="inter CGI")


myFun_anova(df=df_all_long, Filter="hg19_genes_promoters")

myFun_tbl(df=df_all_long, Filter="enhancers")
myFun_anova(df=df_all_long, Filter="enhancers")

myFun_tbl(df=df_all_long, Filter="enhancers")
myFun_anova(df=df_all_long, Filter="enhancers")

myFun_tbl(df=df_all_chrom_long, Filter="Active Promoter")
myFun_anova(df=df_all_chrom_long, Filter="Active Promoter")

myFun_tbl(df=df_all_chrom_long, Filter="Weak Promoter")
myFun_anova(df=df_all_chrom_long, Filter="Weak Promoter")

myFun_tbl(df=df_all_chrom_long, Filter="Weak Promoter")
myFun_anova(df=df_all_chrom_long, Filter="Weak Promoter")
