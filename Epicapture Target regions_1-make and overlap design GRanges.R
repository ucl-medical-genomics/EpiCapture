############################################################################
#   Title: Epicapture Target Regions - 1.step 
#   Description: Script used to import and calculate & visualize overlap for  Target capture BED designs across 5 platforms
#   Author: Miljana Tanic (m.tanic@ucl.ac.uk) 
#   Created: November 2018
#    edited: Mart 2019
#   Last edited: Dec 2019 - count CpG overlap

#####################################################################################


#=========================================
# Load libraries and set working directory
#=========================================


  source("https://bioconductor.org/biocLite.R")

  # biocLite("rtracklayer")
  library(rtracklayer) # not necessary for importing BED files, can do it manualy crating GRanges object from dataframe library(dplyr)
  library(readr)
  # biocLite("annotatr") 
  library(annotatr)
  library(BSgenome.Hsapiens.UCSC.hg19) # sequence of the hg19 genome
  library(ggplot2)
  library(ChIPpeakAnno) # for multiple pairs overlaps

#------------------------------------------
# Set Path to files
#------------------------------------------

  # setwd("~/home/regmani@ad.ucl.ac.uk/EpiCapture/BEDfiles")
  setwd("~/RDS_C2c/EpiCapture/BEDfiles")
  list.files()


  path2Agilent <- "AgilentBEDfiles/S03770311/"
  path2Roche <- "RocheBEDfiles/"
  path2Illumina <- "IlluminaBEDfiles/"
  path2RRBS <- "RRBS_in_silico_BEDfile/"

  Epicapture <- "/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/"
  list.files(Epicapture)

  RESULTS <- paste0(Epicapture,"RESULTS/")
  list.files(RESULTS)


#=========================================
# Import Target capture BED files:
#=========================================

  #-----------------------------------------
  # Roche design BED file:
  #-----------------------------------------

    list.files(path2Roche, pattern = "bed")

    Roche_EpiGiant_design <- read.table(paste0(path2Roche,"130912_HG19_CpGiant_4M_EPI.bed"), header = FALSE, sep="\t")
    head(Roche_EpiGiant_design)

    # Make GRanges object:
    Roche_EpiGiant.gr <- GRanges(seqnames = Roche_EpiGiant_design$V1, IRanges(start=Roche_EpiGiant_design$V2, end = Roche_EpiGiant_design$V3)) 
    Roche_EpiGiant.gr
    # check what chromosomes are covered:
    seqnames(Roche_EpiGiant.gr)
    seqlevels(Roche_EpiGiant.gr)
    #extract a vector with chromosome names:
    chr.list <- levels(seqnames(Roche_EpiGiant.gr))

    # get number of regions covered:
    length(Roche_EpiGiant.gr)
    # get the total length without double counting anything:
    sum(width(reduce(Roche_EpiGiant.gr))) 


  #-----------------------------------------
  # Agilent BED files:
  #-----------------------------------------

    list.files(path2Agilent, pattern = "bed")

    Agilent_regons <- read.table(paste0(path2Agilent, "S03770311_Regions.bed"), skip = 2, header= FALSE, sep= "\t")
    head(Agilent_regons)
    dim(Agilent_regons)

    # make GRanges object:
    Agilent_SureSelect.gr <- GRanges(seqnames = Agilent_regons$V1, IRanges(start=Agilent_regons$V2, end = Agilent_regons$V3)) 
    Agilent_SureSelect.gr
    genome(Agilent_SureSelect.gr) <- "hg19"
    seqnames(Agilent_SureSelect.gr)
    seqlevels(Agilent_SureSelect.gr)
    #extract a vector with chromosome names:
    chr.list.agilent <- levels(seqnames(Agilent_SureSelect.gr))
    # get number of regions covered
    length(Agilent_SureSelect.gr) 
    # get the total length without double counting anything
    sum(width(reduce(Agilent_SureSelect.gr))) 

  #-----------------------------------------
  # Illumina BED files:
  #-----------------------------------------

    list.files(path2Illumina, pattern = "bed")

    Illumina_EPIC_2016_regions <-  read.table(paste0(path2Illumina,"EPIC_target_regions_2016_02_26.bed"), header = FALSE, sep="\t")
    Illumina_EPIC_regions <- read.table(paste0(path2Illumina, "truseq-methyl-capture-epic-manifest-file.bed"), header=FALSE, sep="\t")
    dim(Illumina_EPIC_regions) # 437792 regions
    dim(Illumina_EPIC_2016_regions) # 435423 regions - this one was used for the beta trial

    Illumina_EPIC.gr <- GRanges(seqnames=Illumina_EPIC_2016_regions$V1, IRanges(start=Illumina_EPIC_2016_regions$V2, end=Illumina_EPIC_2016_regions$V3))
    Illumina_EPIC.gr 
    genome(Illumina_EPIC.gr) <- "hg19"
    seqnames(Illumina_EPIC.gr)
    #extract a vector with chromosome names:
    chr.list.illumina <- levels(seqnames(Illumina_EPIC.gr))

    length(Illumina_EPIC.gr)
    sum(width(reduce(Illumina_EPIC.gr)))


    Illumina_EPIC.new.gr <- GRanges(seqnames=Illumina_EPIC_regions$V1, IRanges(start=Illumina_EPIC_regions$V2, end=Illumina_EPIC_regions$V3))
    length(Illumina_EPIC.new.gr)
    sum(width(reduce(Illumina_EPIC.new.gr)))


  #-----------------------------------------
  # RRBS in silico BED files:
  #-----------------------------------------

    # In silico digestion og hg19 with MspI restriction enzyme:
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(plyr)
    library(ggplot2)
    library(reshape2)
    library(scales)

    # Identify the MspI recongnition sites for each chromosomal entry
    # Generate a dataframe with the length of MspI digested fragments
    mdf=data.frame();
    for (i in seq_along(Hsapiens)){
      #print(paste("Processing ",seqnames(Hsapiens)[i], sep=""))
      m<-matchPattern("CCGG", Hsapiens[[i]])
      starts<-start(gaps(m))
      ends<-end(gaps(m))
      temp_df<-data.frame(start=starts-4,end=ends,chr=seqnames(Hsapiens)[i]) #actually end = ends
      temp_df$start<-replace(temp_df$start, temp_df$start == -3, 0)
      temp_df<-temp_df[c("chr","start","end")]
      mdf<-rbind(mdf,temp_df)
    }

    # Extract the digested fragment length

    #  methylation insensitive restriction enzyme MspI, which recognizes CCGG.
    # As a result of partial fragmentation during bisulfte conversion, PCR, and effciency of cluster generation, only a subset of these fragments, typically under 300 bp in length, is sequenced. 
    # Note how only a small amount of the digested DNA is represented in the under 300 bp fraction. 
    # However these smaller fragments are derived from genomic DNA that has a high frequency of MspI sites and therefore a high frequency of potential CpG methylation sites. 
    # This is how RRBS achieves reduced representation
    # Adaptor artifacts typically appear as a peak of approximately 145 bp.

    # Keep only the fragments in the range of 40-600 bp
    mdf$width=mdf$end-mdf$start
    ml<-mdf[mdf$width>39&mdf$width<601,]
    dim(ml)
    counts<-ddply(ml,.(width), nrow)
    # Create a plot of the frequency of the fragment lengths (y-axis is logarithmic)
    p<-ggplot(counts,aes(x=width, y=V1))+geom_line()
    p+scale_y_continuous(trans=log2_trans())

    #My modifications - labels
    p<-ggplot(counts,aes(x=width, y=V1))+geom_line() + labs(title="hg19 MspI in silico digest") + labs( x= "Fragment length (bp)", y="Frequency")
    p+scale_y_continuous(trans=log2_trans()) + scale_x_continuous(breaks=seq(40,600,50))

    # Save as .BED file
    RRBS.bed <-ml[,-4]
    library(dplyr)
    RRBS.bed_sub <- RRBS.bed[!RRBS.bed$chr=="chrUn_gl000249",]

    result <- RRBS.bed[ grep("chrUn_", RRBS.bed$chr) , ]  # Note need to use `a$` to get at the `x`

    # Filter out all but Chr1-22, chrX and chrY
    RRBS.bed_sub <- subset(RRBS.bed, !grepl("chrUn_", RRBS.bed$chr))
    RRBS.bed_sub <- subset(RRBS.bed_sub, !grepl("*_random", RRBS.bed_sub$chr))
    RRBS.bed_sub <- subset(RRBS.bed_sub, !grepl("*_hap?", RRBS.bed_sub$chr))
    RRBS.bed_sub <- subset(RRBS.bed_sub, !grepl("chrM", RRBS.bed_sub$chr))
    tail(RRBS.bed_sub)

    write.table(RRBS.bed_sub, file="RRBS_in_silico.bed", sep="\t",row.name=FALSE, quote=FALSE, col.names = FALSE)

    #------------------------------------------------------------------------------------------
    # [ Optional - if already generated bed]

    # Imorort from bed file
    RRBS.bed_sub <-  read.table(paste0(path2RRBS,"RRBS_in_silico.bed"), header = FALSE, sep="\t")
    RRBS.ir <- IRanges(start=RRBS.bed_sub$V2, end = RRBS.bed_sub$V3, names = RRBS.bed_sub$V1)
    #------------------------------------------------------------------------------------------

    # make an IRanges object from data frame
    RRBS.ir <- IRanges(start=mdf$start, end = mdf$end, names = mdf$chr)
    RRBS.ir
    names(RRBS.ir)

    # make a GRanges object
    RRBS.gr <- GRanges(seqnames = names(RRBS.ir), ranges = RRBS.ir)
    genome(RRBS.gr) <- "hg19"
    RRBS.gr 
    seqnames(RRBS.gr)
    seqlengths(RRBS.gr) # not defined for any platform gr?
    sum(width(RRBS.gr)) # gets width of all ranges

    # Filter by fragment length!
    RRBS.gr.sub.plot <- RRBS.gr[width(RRBS.gr)>39 & width(RRBS.gr)<601] # select only fragments between 40-600bp for plotting frequency
    RRBS.gr.sub <- RRBS.gr[width(RRBS.gr)>39 & width(RRBS.gr)<301] # select only fragments between 40-300bp for target annotation
    #  104983071 bp

    # Select only first 100bp from each end of the restricted fragment:
    RRBS.gr.sub.short <- RRBS.gr.sub[width(RRBS.gr.sub)<=200]
    RRBS.gr.sub.200 <- RRBS.gr.sub[width(RRBS.gr.sub)>200]
    RRBS.gr.sub.R1.100 <- resize(RRBS.gr.sub.200, 100) # get first 100 bp
    RRBS.gr.sub.R2.100 <- restrict(RRBS.gr.sub.200, start=(end(RRBS.gr.sub.200)-100)) # get last 100 bp

    # merge granges: 
    RRBS.gr.sub <- c(RRBS.gr.sub.short,RRBS.gr.sub.R1.100,RRBS.gr.sub.R2.100)
    #  97848603 bp

    RRBS.gr.sub
    length(RRBS.gr.sub)
    sum(width(reduce(RRBS.gr.sub))) # reduce overlaping fragments to to count each bp only once


    # Plot frequency distribution of RRBS.gr.sub s - Works on GRanges object! 
    table(width(RRBS.gr.sub.plot)) # count occurences of each fragment size (fragment size frequency)
    RRBS.counts <- as.data.frame(table(width(RRBS.gr.sub.plot))) # convert to a data frame
    p <- ggplot(RRBS.counts, 
                aes(x=as.numeric(RRBS.counts$Var1), # x-axis generated by table is a factor so it needs to be converted to numeric for it to work with ggplot
                    y=RRBS.counts$Freq)) + 
      geom_line() + 
      labs(title="hg19 MspI in silico digest") + 
      labs( x= "Fragment length (bp)", y="Frequency") 
    p+scale_y_continuous(trans="log2") + scale_x_continuous(breaks=seq(40,600,50))


#--------------------------------------
# Save several R objects to load later:
#--------------------------------------
  save(RRBS.gr.sub, Illumina_EPIC.gr, Illumina_EPIC.new.gr, Agilent_SureSelect.gr, Roche_EpiGiant.gr, file="Platform_Design_GRanges.RData")

  # these can be reloaded into another R session using the load() function:
  load(file="Platform_Design_GRanges.RData") 


#============================================
# Find overlap between Platform design files:
#============================================

  source("https://bioconductor.org/biocLite.R")
  biocLite("ChIPpeakAnno")
  library(GenomicRanges)
  library(Biostrings)
  library(ChIPpeakAnno)
  library(HelloRanges)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # Load GRanges files for each platform
  load(file="Platform_Design_GRanges.RData")

  # Load the BSgenome package for homo sapiens.
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19


#-----------------------------------
# Make platforms GRangesList object
#----------------------------------

  platforms_design.grl <- GRangesList(Agilent_SureSelect.gr, Roche_EpiGiant.gr, Illumina_EPIC.new.gr, RRBS.gr.sub)
  names(platforms_design.grl) <- c("Agilent_SureSelect", "Roche_EpiGiant", "Illumina_EPIC", "RRBS_sub")

  # Get target design size of each platform
  sum(width(reduce(platforms_design.grl)))
  # lapply(platforms_design.grl, function(x) sum(width(reduce(x))) )

  # Get sequences from GRangesList object:  
  grl.seq <- getSeq(genome, platforms_design.grl) # The return values are DNAStringSetList objects.

  grl.seq2 <- getSeq(genome, reduce(platforms_design.grl)) # The return values are DNAStringSetList objects.

  # Count number of occurances for CG in a given region of list object:
  sapply(grl.seq, function(x) sum(vcountPattern("CG", x)))
  
  # get number of regions (probes/fragments) covered:
  lapply(platforms_design.grl, function(x) length(x))
  
  # get the total length without double counting anything:
  lapply(platforms_design.grl, function(x) sum(width(reduce(x))))

  # get the average length of target regions:
  lapply(platforms_design.grl, function(x) mean(width(reduce(x)))
  lapply(platforms_design.grl, function(x) min(width(reduce(x))))
  lapply(platforms_design.grl, function(x) max(width(reduce(x))))


#------------------------------------------------
# Find overlaps between several GRanges objects:
#-----------------------------------------------

  # returns an object of overlappingPeaks, which contains there elements: venn_cnt, peaklist (a list of overlapping peaks or unique peaks), and overlappingPeaks (a list of data frame consists of the annotation of all the overlapping peaks).
  ol <- findOverlapsOfPeaks(Agilent_SureSelect.gr, Roche_EpiGiant.gr,Illumina_EPIC.new.gr, RRBS.gr.sub)
  head(ol)

#--------------------------------------------------------------------
# Visualize the number of common and specific peaks with Venn diagram
#--------------------------------------------------------------------

  makeVennDiagram(ol, connectedPeaks="keepAll",
                NameOfPeaks=c("Agilent", "Roche", "Illumina", "RRBS"),
                main="Overlap of platform breadth of coverage by design",
                fill=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3" ), # circle fill color
                col=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3" ), #circle border color
                cat.col=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3") )

  # # Or in one shot:
  # res <- makeVennDiagram(Peaks=list(Agilent_SureSelect.gr, Roche_EpiGiant.gr, Illumina_EPIC.gr, RRBS.gr.sub), ignore.strand=TRUE,
  #                       NameOfPeaks=c("Agilent", "Roche", "Illumina", "RRBS"),
  #                       main="Venn Diagram for platform coverage design",
  #                       fill=c(1,2,3,4), scaled=TRUE, euler.d=FALSE)
  # # Get Venn counts:
  # res$vennCounts


  # A pie chart is used to demonstrate the overlap features of the common peaks.
  pie1(table(ol$overlappingPeaks[["Agilent_SureSelect.gr///Roche_EpiGiant.gr"]]$overlapFeature)) 


#-------------------------------------------
# Visualize intersecting regions with UpSetR 
#-------------------------------------------

  # biocLite("UpSetR")
  library(UpSetR)
  #biocLite("ComplexHeatmap")
  library(devtools)
  # install_github("jokergoo/ComplexHeatmap") # latest version
  library(ComplexHeatmap)


  # union mode: 1 means in that set and 0 is not taken into account. When there are multiple 1, the relationship is OR. Then, 1 1 0 means a set of elements in set A or B, and they can also in C or not in C (union(A, B)). Under this mode, the seven combination sets can overlap.

  # When the sets are genomic regions, the size is calculated as the sum of the width of regions in each set (or in other words, the total number of base pairs).

  # List 
  lt <- list("Agilent" =Agilent_SureSelect.gr,
            "Roche" =  Roche_EpiGiant.gr,
            "Illumina"= Illumina_EPIC.new.gr,
            "RRBS"=RRBS.gr.sub)

  m = make_comb_mat(lt)
  # m = m[comb_size(m) > 500000]
  
  # save
  write.table(comb_size(m), file=paste0(RESULTS,"1_PlatformDesignDiferences/UpSet_matrix.txt"))


  # plot  
  UpSet(m, 
        set_order = c("Agilent", "Roche", "Illumina", "RRBS"), 
        top_annotation = upset_top_annotation(
          m,
          axis_param= list(at=c(0,2e7,4e7,6e7,8e7),
                          labels=c("0MB", "20MB","40MB",  "60MB", "80MB")),
          height = unit(4, "cm")
        ),
        right_annotation = upset_right_annotation(
          m,
          axis_param= list(at=c(0, 5e7,10e7), 
                          labels=c("0MB", "50MB", "100MB"),
                          labels_rot=0),
          width=unit(4, "cm")),
        pt_size = unit(5, "mm"), lwd = 3,
        #comb_col = c("olivedrab", "mediumslateblue", "black")[comb_degree(m)] # doesn't show all 4 intersection
  )
  # export as 10' x 5' pdf plot


#----------------------------------------
# Alternative Venn - doesn't work nicely
#----------------------------------------


  # After finding the overlapping peaks, you can use annotatePeakInBatch to annotate the overlapping peaks with the genomic features in the AnnotationData within certain distance away specified by maxgap, which is 5kb in the following example.


  # But if you have them as a single GRanges stratified by type metadata, you could make it GRangesList and plot it like this:

  # peaks1$type<-"TF1"
  # peaks2$type<-"TF2"
  # peaks3$type<-"TF3"
  # gr <- c(peaks1, peaks2, peaks3) # like your data
  # 
  # grl <- splitAsList(gr, gr$type)
  # res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl))


  # # If you want a nice and colorful Venn diagram, then check https://github.com/js229/Vennerable/blob/master/README.md to install Venerable (I could not get it from Rforge) and do it like this:
  # biocLite(c("RBGL","graph"))
  # install.packages("devtools"); library(devtools);
  # install_github("js229/Vennerable"); library(Vennerable);
  # vignette("Venn")
  # 
  # library(Vennerable)
  # venn_cnt2venn <- function(venn_cnt){
  #   n <- which(colnames(venn_cnt)=="Counts") - 1
  #   SetNames=colnames(venn_cnt)[1:n]
  #   Weight=venn_cnt[,"Counts"]
  #   names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  #   Venn(SetNames=SetNames, Weight=Weight)
  # }
  # v <- venn_cnt2venn(res$vennCounts)
  # plot(v)

  # venn_cnt <- res$vennCounts
  # n <- which(colnames(venn_cnt)=="Counts") - 1
  # SetNames=colnames(venn_cnt)[1:n]
  # Weight=venn_cnt[,"Counts"]
  # names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  # v <- Venn(SetNames=SetNames, Weight=Weight)
  # plot(v)

  # 
  # # Following works for getting the exact intersects between all the ranges.
  # Reduce(intersect, list(gr, gr1, gr2))


#=========================================================================
# Find overlap wth known CpG sites and compare coverage between Platforms
#=========================================================================

    Epicapture <- "/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/"
    list.files(Epicapture)

    BEDfiles <- paste0(Epicapture,"BEDfiles/")
    list.files(BEDfiles)

    load(file=paste0(BEDfiles,"FeatureAnnotations/hg19_CpGsites_GRanges.RData"))

# make subsets og CpG sites per platform
  # subsetByOverlaps method simply subsets the first GRanges object to include only those that overlap the second.
      subsetByOverlaps(gr, g)

  platforms_design.grl <- GRangesList(Agilent_SureSelect.gr, Roche_EpiGiant.gr, Illumina_EPIC.new.gr, RRBS.gr.sub)

     cpgr.Agilent_SureSelect <-  subsetByOverlaps(cpgr, Agilent_SureSelect.gr, ignore.strand=TRUE)
     cpgr.Roche_EpiGiant <-  subsetByOverlaps(cpgr, Roche_EpiGiant.gr, ignore.strand=TRUE)
     cpgr.Illumina_EPIC <-  subsetByOverlaps(cpgr, Illumina_EPIC.gr, ignore.strand=TRUE)
     cpgr.RRBS <-  subsetByOverlaps(cpgr,  RRBS.gr.sub, ignore.strand=TRUE)
   
   # length of each cpgr.PLATFORM object should mach count of CpGs per platform

  length(cpgr.Agilent_SureSelect)
  length(cpgr.Roche_EpiGiant)
  length(cpgr.Illumina_EPIC) 
  length(cpgr.RRBS)     

    #  > sapply(grl.seq, function(x) sum(vcountPattern("CG", x)))
    # [1] 3153816 2806466 3346505 4070782

    # > length(cpgr.Agilent_SureSelect)
    # [1] 3172538
    # > length(cpgr.Roche_EpiGiant)
    # [1] 2818345
    # > length(cpgr.Illumina_EPIC) 
    # [1] 3353173
    # > length(cpgr.RRBS)     
    # [1] 4079434

  # reduce merges consecutive CpG sites
    # > length(reduce(cpgr.Agilent_SureSelect))
    # [1] 2949002
    # > length(reduce(cpgr.Roche_EpiGiant))
    # [1] 2620071
    # > length(reduce(cpgr.Illumina_EPIC))
    # [1] 3143186
    # > length(reduce(cpgr.RRBS))
    # [1] 3817962

  #-------------------------------------------
  # Visualize intersecting regions with UpSetR 
  #-------------------------------------------


  # List 
  lt_cpg <- list("Agilent" =cpgr.Agilent_SureSelect,
            "Roche" =  cpgr.Roche_EpiGiant,
            "Illumina"= cpgr.Illumina_EPIC,
            "RRBS"=cpgr.RRBS)

  m_cpg = make_comb_mat(lt)
  # m = m[comb_size(m) > 500000]

  # save
  write.table(m_cpg, file=)

  # plot  
  UpSet(m_cpg, 
        set_order = c("Agilent", "Roche", "Illumina", "RRBS"), 
        top_annotation = upset_top_annotation(
          m,
          axis_param= list(at=c(0,2e7,4e7,6e7),
                          labels=c("0MB", "10 M","20 M",  "30 M")),
          height = unit(4, "cm")
        ),
        right_annotation = upset_right_annotation(
          m,
          axis_param= list(at=c(0,2e7,4e7,6e7,8e7,1e8),
                          labels=c("0MB", "10 M","20 M",  "30 M", "40 M", "50M" ),
                          labels_rot=0),
          width=unit(4, "cm")),
        pt_size = unit(5, "mm"), lwd = 3,
        #comb_col = c("olivedrab", "mediumslateblue", "black")[comb_degree(m)] # doesn't show all 4 intersection
  )

  # intersection size plots the actual size of CpG regions not the count of common CpGs so the output is "double the count" since each CpG has size of 2bp!
  # can just divede by two


#--------------------------------------------------------------------
# Visualize the number of common and specific CpGs with Venn diagram
#--------------------------------------------------------------------

    ol_cpg <- findOverlapsOfPeaks(cpgr.Agilent_SureSelect, cpgr.Roche_EpiGiant,cpgr.Illumina_EPIC, cpgr.RRBS)
    head(ol_cpg)


    makeVennDiagram(ol_cpg, connectedPeaks="keepAll",
                NameOfPeaks=c("Agilent", "Roche", "Illumina", "RRBS"),
                main="Overlap of CpGs covered by design",
                fill=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3" ), # circle fill color
                col=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3" ), #circle border color
                cat.col=c("#3D79F3FF", "#E6352FFF", "#34A74BFF","#7570b3") )
