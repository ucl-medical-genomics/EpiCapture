############################################################################
#   Title: Epicapture Target Regions - 2.step 
#   Description: Script used to make annotation files for genomic features
#   Author: Miljana Tanic (m.tanic@ucl.ac.uk) 
#   Created: November 2018
#   Last edited: Mart 2019
#############################################################################

# !!! If packages, commands or something doesn't work - disconect from the  server and reconect, it will usually "miracolously" start working!

#=========================================
# Load libraries and set working directory
#=========================================

source("https://bioconductor.org/biocLite.R")

# biocLite("rtracklayer")
library(rtracklayer) # not necessary for importing BED files, can do it manualy crating GRanges object from dataframe library(dplyr)
library(readr)
# biocLite("annotatr")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("annotatr")

library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg19) # sequence of the hg19 genome
library(ggplot2)
library(AnnotationHub)

#------------------------------------------
# Set Path to files
#------------------------------------------

#setwd("~/home/regmani@ad.ucl.ac.uk/EpiCapture/BEDfiles")
setwd("~/RDS_C2c/EpiCapture/BEDfiles/FeatureAnnotations")

list.files()

#=========================================
# Make Anotations - GRanges from anotate R 
#==========================================

#------------------------------------------------- 
# Install annotater for anotation of genomic regions:
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("annotatr")
library(annotatr)

# annotate_regions() requires a GRanges object (either the result of  read_regions() or an existing object), a GRanges object of the annotations, and a logical value indicating whether to ignore.strand when calling  GenomicRanges::findOverlaps().
# The positive integer minoverlap is also passed to GenomicRanges::findOverlaps() and specifies the minimum overlap required for a region to be assigned to an annotation.

# Before annotating regions, they must be built with build_annotations() which requires a character vector of desired annotation codes.
# There are three types of annotations available to annotatr:
# 1. Built-in annotations including CpG annotations, genic annotations, enhancers, GENCODE lncRNAs, and chromatin states from chromHMM.
# 2. AnnotationHub annotations include any GRanges resource within the Bioconductor AnnotationHub web resource.
# 3. Custom annotations provided by the user. 


#-----------------------
# CpG Annotations
# The CpG islands are the basis for all CpG annotations, and are given by the  AnnotationHub package for the given organism.
# CpG shores are defined as 2Kb upstream/downstream from the ends of the CpG islands, less the CpG islands.
# CpG shelves are defined as another 2Kb upstream/downstream of the farthest upstream/downstream limits of the CpG shores, less the CpG islands and CpG shores. 
# The remaining genomic regions make up the inter-CGI annotation.

hg19_cpgs  # shortcut annotates regions to CpG islands, CpG shores, CpG shelves, and inter-CGI. 


#-----------------------
# Genic Annotations
# The genic annotations are determined by functions from GenomicFeatures and data from the TxDb.* and org.*.eg.db packages. 
# Genic annotations include 1-5Kb upstream of the TSS, the promoter (< 1Kb upstream of the TSS), 5’UTR, first exons, exons, introns, CDS, 3’UTR, and intergenic regions (the intergenic regions exclude the previous list of annotations).
# Also included in genic annotations are intronexon and exonintron boundaries. These annotations are 200bp up/down stream of any boundary between an exon and intron. Important to note, is that the boundaries are with respect to the strand of the gene.
# Non-intergenic gene annotations include Entrez ID and gene symbol information where it exists. The org.*.eg.db packages for the appropriate organisms are used to provide gene IDs and gene symbols.
# The genic annotations have populated tx_id, gene_id, and symbol columns. Respectively they are, the knownGene transcript name, Entrez Gene ID, and gene symbol.

hg19_basicgenes # shortcut annotates regions to 1-5Kb, promoters, 5’UTRs, exons, introns, and 3’UTRs.

#------------------------------
# FANTOM5 Permissive Enhancers
# FANTOM5 permissive enhancers were determined from bi-directional CAGE transcription as in Andersson et al. (2014), and are downloaded and processed for hg19 and mm9 from the FANTOM5 resource. Using the rtracklayer::liftOver() function, enhancers from hg19 are lifted to hg38, and mm9 to mm10.


#------------------------------
# GENCODE lncRNA transcripts
# The long non-coding RNA (lncRNA) annotations are from GENCODE for hg19, hg38, and mm10. The lncRNA transcripts are used, and we eventually plan to include the lncRNA introns/exons at a later date. The lncRNA annotations have populated  tx_id, gene_id, and symbol columns. Respectively they are, the Ensembl transcript name, Entrez Gene ID, and gene symbol. As per the transcript_type field in the GENCODE anntotations, the biotypes are given in the id column.


#--------------------------------
# Chromatin states from ChromHMM
# Chromatin states determined by chromHMM (Ernst and Kellis (2012)) in hg19 are available for nine cell lines (Gm12878, H1hesc, Hepg2, Hmec, Hsmm, Huvec, K562, Nhek, and Nhlf) via the UCSC Genome Browser tracks. Annotations for all states can be built using a shortcut like hg19_Gm12878-chromatin, or specific chromatin states can be accessed via codes like hg19_chromatin_Gm12878-StrongEnhancer or  hg19_chromatin_Gm12878-Repressed.

#--------------------------------
# AnnotationHub Annotations
# The AnnotationHub web resource provides a central location where genomic files (e.g., VCF, bed, wig) and other resources from standard locations (e.g., UCSC, Ensembl) can be discovered. The resource includes metadata about each resource, e.g., a textual description, tags, and date of modification. The client creates and manages a local cache of files retrieved by the user, helping with quick and reproducible access.

#--------------------------------
# Custom Annotations
# Users may load their own annotations from BED files using the  read_annotations() function, which uses the rtracklayer::import() function. 
# The output is a GRanges with mcols() for id, tx_id, gene_id, symbol, and  type. 
# If a user wants to include tx_id, gene_id, and/or symbol in their custom annotations they can be included as extra columns on a BED input file.

## Custom annotation objects are given names of the form genome_custom_name

gencode_lncRNA.bed <- import("/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/BEDfiles/FeatureAnnotations/gencode.v19.long_noncoding_RNAs.bed")
head(gencode_lncRNA.bed)
read_annotations(con=gencode_lncRNA.bed, genome = "hg19", name = "gencode_lncRNA", format="bed") 

print(annotatr_cache$list_env())


#======================
# Define Annotations
#=====================

    builtin_genomes()
    builtin_annotations()

    # Make list of built-in annotations for hg19
    hg19_annotations <- grepl("hg19", builtin_annotations()) # grep logical --> logical vector
    # index and subset annotations for hg19:
    hg19_annotations <- builtin_annotations()[hg19_annotations]

#-----------------------------
# Build the gene annotations:
#-----------------------------

    #  Built in annotates regions to 1-5Kb, promoters, 5’UTRs, exons, introns, and 3’UTR:
    annotations_genes.gr = build_annotations(genome = 'hg19', annotations = "hg19_basicgenes")

    annotations_genes.gr
    # get names of chromosomes:
    seqnames(annotations_genes.gr) # has all chromosomes as in hg19 including chrUn
    # calculate reduced region size:
    sum(as.numeric(width(reduce(annotations_genes.gr))))  # size of all genic reagions in bp
    # show metadata columns:
    elementMetadata(annotations_genes.gr)

    # subset GRanges according to defined features (CGI, genic, regulatory, etc...)
    # mcols(annotations_genes.gr)$type=="hg19_genes_cds" # 2 subset - logical 
    annotations_genes.gr[mcols(annotations_genes.gr)$type=="hg19_genes_cds" , ]

    # Make GRangesList object by splitting by type of annotation
    annotations_genes.grl <- split(annotations_genes.gr, mcols(annotations_genes.gr)$type)
    names(annotations_genes.grl) # get names of the list

    # You can always unlist back to GRanges object:
    unlist(annotations_genes.grl)


# calculate size in bp of each feature - use reduce not to count overlaping regions twice
    #sum(as.numeric((width(reduce(annotations_genic.gr[mcols(annotations_genic.gr)$type=="hg19_genes_5UTRs" , ])))))
    # better way:
    sum(width(reduce(annotations_genes.grl)))

    # It is counting both strands! (double the size of the genome? it would make sense since genes are direrently defined on different strands - the iformation load is actualy for double stranded DNA - countng 6 Gb of information) - solved by using reduce not to count overlaping regions


    # Note inclusion of custom annotation, and use of shortcuts
    annot_genes= c("hg19_genes_1to5kb", "hg19_genes_promoters", "hg19_genes_cds" , "hg19_genes_5UTRs", "hg19_genes_exons", "hg19_genes_firstexons" , "hg19_genes_introns" , "hg19_genes_intronexonboundaries", "hg19_genes_exonintronboundaries" , "hg19_genes_3UTRs", "hg19_genes_intergenic"   )
    annotations_genic.gr = build_annotations(genome = 'hg19', annotations = annot_genes)

    # Make GRangesList object by splitting by type of annotation
    annotations_genic.grl <- split(annotations_genic.gr, mcols(annotations_genic.gr)$type)
    names(annotations_genic.grl)

    # getSeq() doesnt't work if there is any undefined strand "*" - need to remove intergenic info
    #intergenic.gr <- annotations_genic.grl[8]  


#------------------------------------
# Build the CpG islands annotations:
#------------------------------------

    # annotates regions to CpG islands, CpG shores, CpG shelves, and inter-CGI:
    annotations_CpGs.gr = build_annotations(genome = 'hg19', annotations = "hg19_cpgs")
    # Make GRangesList object by splitting by type of annotation
    annotations_CpGs.grl <- split(annotations_CpGs.gr, mcols(annotations_CpGs.gr)$type)
    names(annotations_CpGs.grl)

    # Custom:
    annot_CpG = c("hg19_cpg_islands", "hg19_cpg_shores", "hg19_cpg_shelves", "hg19_cpg_inter")
    annotations_CGI.gr = build_annotations(genome = 'hg19', annotations = "hg19_cpg_islands")
    annotations_shores.gr = build_annotations(genome = 'hg19', annotations = "hg19_cpg_shores")
    annotations_shelves.gr = build_annotations(genome = 'hg19', annotations = "hg19_cpg_shelves")
    annotations_interCGI.gr = build_annotations(genome = 'hg19', annotations = "hg19_cpg_inter")
    sum(width(reduce(annotations_interCGI.gr)))

    annotations_CpGs.grl <- GRangesList("hg19_cpg_islands"=annotations_CGI.gr,"hg19_cpg_shores"=annotations_shores.gr,"hg19_cpg_shelves"=annotations_shelves.gr, "hg19_cpg_inter"=annotations_interCGI.gr )

    # Size of each CGI feature
    # sum(width(reduce(annotations_CpGs.grl))) # sum of widths of each grl element BUT gives wrong data for inter-CGI!

    # this is the correct way of doing this and gives good numbers!
    lapply(annotations_CpGs.grl, function(x) sum(width(reduce(x))))

#------------------------------------
# Build the Regulatory annotations:
#------------------------------------

    # Enhancers:
    annot_enhancers_fantom <- c("hg19_enhancers_fantom")         
    annotations_enhancers_fantom.gr = build_annotations(genome = 'hg19', annotations = annot_enhancers_fantom)
    sum(as.numeric(width(reduce(annotations_enhancers_fantom.gr))))


    # Long non coding RNA:
    path2lncRNA <- "/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/BEDfiles/FeatureAnnotations"
    file = "gencode.v19.long_noncoding_RNAs.bed"
    # load their own annotations from BED files using the  read_annotations() function
    annotations_lncRNA.gr <- read_annotations(con = paste0(path2lncRNA,"/", file), genome="hg19", name="lncRNA", format = "bed")
    sum(as.numeric(width(reduce(annotations_lncRNA.gr))))

    #Using the build_ah_annots() function, users can turn any resource of class  GRanges into an annotation for use in annotatr eg. from AnnotationHub

#------------------------------------
# Chromatine annotations:
#------------------------------------

    annot_chromHMM <- c("hg19_Gm12878-chromatin",               
                        "hg19_H1hesc-chromatin", "hg19_Hepg2-chromatin" , "hg19_Hmec-chromatin","hg19_Hsmm-chromatin"  , "hg19_Huvec-chromatin" , "hg19_K562-chromatin" ,  "hg19_Nhek-chromatin", "hg19_Nhlf-chromatin")
    #list of all annotations for chromatin - this will list also stratifed by type of chromatin mark
    # hg19_annotations_chromHMM <- grepl("chromatin", hg19_annotations) # grep logical --> logical vector
    annotations_chromHMM.gr  = build_annotations(genome = 'hg19', annotations = annot_chromHMM)
    sum(as.numeric(width(reduce(annotations_chromHMM.gr))))

#-----------------------------------------------------
# subset GRanges according to type of chromatine mark
#-----------------------------------------------------
    grepl("Insulator", mcols(annotations_chromHMM.gr)$type)
    annotations_chromHMM.gr[ grepl("Insulator", mcols(annotations_chromHMM.gr)$type)==TRUE,]

    # Make a GRangesList object from regulatory annotations - split for each cell line a chromatine mark:
    annotations_chromHMM.grl = split(annotations_chromHMM.gr,mcols(annotations_chromHMM.gr)$type)
    names(annotations_chromHMM.grl) # 107 elements

    # grep logical --> logical vector
    print(grepl("Insulator", hg19_annotations))

    #  subset annotations that contain same type of chromatine mark for all analyzed cell lines:  
    annotations_ActivePromoter.gr  = annotations_chromHMM.gr[ grepl("ActivePromoter", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_WeakPromoter.gr  = annotations_chromHMM.gr[ grepl("WeakPromoter", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_PoisedPromoter.gr  =annotations_chromHMM.gr[ grepl("PoisedPromoter", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_StrongEnhancer.gr  = annotations_chromHMM.gr[ grepl("StrongEnhancer", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_WeakEnhancer.gr  = annotations_chromHMM.gr[ grepl("WeakEnhancer", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_Insulator.gr  = annotations_chromHMM.gr[ grepl("Insulator", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_TxnElongation.gr  =annotations_chromHMM.gr[ grepl("TxnElongation", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_TxnTransition.gr  = annotations_chromHMM.gr[ grepl("TxnTransition", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_WeakTxn.gr  = annotations_chromHMM.gr[ grepl("WeakTxn", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_Repressed.gr  = annotations_chromHMM.gr[ grepl("Repressed", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    annotations_Heterochrom.gr  = annotations_chromHMM.gr[ grepl("Heterochrom/lo", mcols(annotations_chromHMM.gr)$type)==TRUE,]
    #annotations_Repetitive-CNV.gr  = annotations_chromHMM.gr[ grepl("Repetitive", mcols(annotations_chromHMM.gr)$type)==TRUE,] # object 'annotations_Repetitive' not found


    # Calculate the size of each region:
    sum(width(reduce(annotations_chromHMM.gr[grepl("Repetitive/CNV", mcols(annotations_chromHMM.gr)$type)==TRUE,])))

    annotations_regul.grl <- GRangesList(annotations_ActivePromoter.gr, annotations_WeakPromoter.gr, annotations_PoisedPromoter.gr, annotations_StrongEnhancer.gr, annotations_WeakEnhancer.gr, annotations_Insulator.gr, annotations_TxnElongation.gr, annotations_TxnTransition.gr, annotations_WeakTxn.gr, annotations_WeakTxn.gr, annotations_Repressed.gr, annotations_Heterochrom.gr)

    # Get sizes of each chromHHM feature in each cell line:
    file <- sum(width(reduce(annotations_regul.grl)))
    write.table(x=file, file="ChromHHM_size.txt", sep="\t")


#-----------------------------------------------
# Save Annotation GRanges objects to load later:
#-----------------------------------------------
    save(annotations_genes.grl, file="Genes_annotations_GRangesList.RData")
    save(annotations_genes.gr, file="Genes_annotations_GRanges.RData")

    save(annotations_genic.grl, file="Genic-extended_annotations_GRangesList.RData")
    save(annotations_CpGs.grl, file="CpG_annottations_GRangesList.RData")
    save(annotations_CpGs.gr, file="CpG_annottations_GRanges.RData")

    save(annotations_enhancers_fantom.gr, file="Ehnahcers_annotations_GRanges.RData")
    save(annotations_lncRNA.gr, file="lncRNA_annotations_GRanges.RData")
    save(annotations_regul.grl, file="ChromHMM-byState_GRangesList.RData")
    save(annotations_chromHMM.grl, file="ChromHMM-byCellLine-byState_GRangesList.RData")
    save(annotations_ActivePromoter.gr, annotations_WeakPromoter.gr, annotations_PoisedPromoter.gr, annotations_StrongEnhancer.gr, annotations_WeakEnhancer.gr, annotations_Insulator.gr, annotations_TxnElongation.gr, annotations_TxnTransition.gr, annotations_WeakTxn.gr, annotations_WeakTxn.gr, annotations_Repressed.gr, annotations_Heterochrom.gr, file="ChromHMM-byState_GRanges.RData")

    # these can be reloaded into another R session using the load() function:
    load(file="Genes_annotations_GRangesList.RData") 
    load(file="extended_annotations_GRangesList.RData") 
    load(file="CpG_annottations_GRangesList.RData") 
    load(file="Enhancers_annotations_GRanges.RData") 
    load(file="lncRNA_annotations_GRanges.RData") 
    load(file="ChromHMM-byState_GRangesList.RData") 
    load(file="ChromHMM-byCellLine-byState_GRangesList.RData") 
    load(file="ChromHMM-byState_GRanges.RData") 


#=====================
# AnnotationHub
#=====================
    library(AnnotationHub)
    ah <- AnnotationHub()

    # data provaiders
    unique(ah$dataprovider)

    # species available
    head(unique(ah$species))

    # type of data availeble
    head(unique(ah$rdataclass))

    # make queries for spcific data type
    dm <- query(ah, c("GRanges", "UCSC", "Homo sapiens"))
    dm["AH5012"]

    #-----------------------------------------------
    # get the CpG islands positions from UCSC, most easily from AnnotationHub
    #-----------------------------------------------

    query(ah, c("cpg","hg19"))
    cpgs_ucsc.gr <- ah[["AH5086"]]

    # save CGI annotations
    save(cpgs_ucsc.gr, file="UCSC_CGI_GRanges.RData")

    #-----------------------------------------------
    # get the Roadmap data from AnnotatioHub
    #-----------------------------------------------

    epiFiles <- query(ah, "EpigenomeRoadMap")

    # get a more descriptive idea of these different files
    sort(table(epiFiles$description), decreasing=TRUE)




#=========================================================
# get  positions of all CpG sites genome postion for hg19
#=========================================================

    library(BSgenome.Hsapiens.UCSC.hg19)  
    chrs <- names(Hsapiens)[1:24]
    
    cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
    
    cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
    
    cpgr

    # save CpG sites annotations
    save(cpgr, file="hg19_CpGsites_GRanges.RData")

#=========================================================
# Get anotations from Enseble Regulatory build
#=========================================================

    # source("http://bioconductor.org/biocLite.R")
    # biocLite("biomaRt")
    # library(biomaRt)
    # listMarts()

    # listEnsembl()   
    # listEnsembl("GRCh37")
    # listEnsemblArchives(name="GRCh37")


    # # GRC38 Ensembl Regulatory Build:
    # #ensembl <- wget(http://www.ensembl.org/biomart/martview/29726bf6c14acd6b0f252bd47f0064c5?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_regulatory_feature.default.regulatory_feature.chromosome_name|hsapiens_regulatory_feature.default.regulatory_feature.chromosome_start|hsapiens_regulatory_feature.default.regulatory_feature.chromosome_end|hsapiens_regulatory_feature.default.regulatory_feature.feature_type_name&FILTERS=&VISIBLEPANEL=resultspanel)

    # ensembl <- useMart("ensembl")

    # ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)


    # The Ensembl Regulatory Build (Zerbino et al. 2015) contains a genome-wide set of regions that are likely to be involved in gene regulation. The classification of these regions are:

    #     Promoters
    #     Promoter flanking regions
    #     Enhancers
    #     CTCF binding sites
    #     Transcription factor binding sites
    #     Open chromatin regions

    # Cell type specific activity levels

        # For each cell type the regulatory features are assigned labels to describe their activity levels. For each cell type and each regulatory feature defined in the Regulatory Build, we check whether there is a segmentation state or experimental evidence (TF binding from ChIP-seq data or open chromatin from DNase-seq data) that could be used to test the activity of this annotation. Based on this then the activity levels can be one of the following:

        # ACTIVE, when the feature displays an active epigenetic signature.
        # POISED, when the feature displays a epigenetic signature with the potential to be activated.
        # REPRESSED, when the feature is epigenetically repressed.
        # INACTIVE, when the region bears no epigenetic modifications from the ones included in the Regulatory Build.
        # NA, when there is no available data in the cell type for this feature.


#=========================================================
# Repetitive elements
#=========================================================

    # download the RepeatMasker results performed on hg19 from the UCSC Table Browser tool. 
    # I’ve downloaded the results as "all fields from selected table" file and named it hg19_rmsk.gz
    $ gunzip hg19_rmsk.gz

    # use shell bedops to convert to BED file and keep lass annotation

    # convert it to sorted BED data 
    $ rmsk2bed < hg19_rmsk

    # Use bedops –merge to merge elements, e.g.: rmsk2bed < foo.out | bedops --merge - > merged_repeatmasker_elements.bed

    #  we strip the header element, but adding the --keep-header option will preserve this header as a BED element that uses _header as a chromosome name:

    # Note the conversion from 1- to 0-based coordinate indexing, in the transition from RepeatMasker annotation output to BED. BEDOPS supports operations on input with any coordinate indexing, but the coordinate change made here is believed to be convenient for most end users.

    