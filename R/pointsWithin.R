#' Helper functionalities to determine if SNPs are within CNVs or other ranges
#'
#' @author: Jean Fan
#'

library(GenomicRanges)
library(BiocParallel)

#' From data frame of ranges to GRanges object
#'
#' @examples
#' cnv <- data.frame(
#'    'chr' = c('chr1', 'chr2'),
#'    'startPos' = c(5, 0),
#'    'endPos' = c(100, 20000)
#'    )
#' range2GRanges(df)
#'
range2GRanges <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df[,1],
        ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
}

#' From data frame of positions to GRanges object
#'
#' @examples
#' snps <- data.frame(
#' 'chr' = c('chr2', 'chr1'),
#'    'pos' = c(20005, 7)
#'    )
#' pos2GRanges(snps)
#'
pos2GRanges <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df[,1],
        ranges=IRanges(start = df[,2], width=1)
    )
    return(gr)
}


#' Determine if a point is within a range
#'
#' @examples
#' snp1 <- data.frame(
#'     'chr' = 'chr1',
#'     'pos' = 7
#'     )
#' cnv <- data.frame(
#'     'chr' = c('chr1', 'chr2', 'chr3'),
#'     'startPos' = c(5, 0, 1),
#'     'endPos' = c(100, 20000, 2)
#'     )
#' pointsWithin(snp1, cnv)
#'
#' # no hits
#' snp2 <- data.frame(
#'     'chr' = 'chrX',
#'     'pos' = 7
#'     )
#' pointsWithin(snp2, cnv)
#'
pointsWithin <- function(pos, ranges) {
    posGRanges <- pos2GRanges(pos)
    rangesGRanges <- range2GRanges(ranges)
    overlap <- GenomicRanges::findOverlaps(rangesGRanges, posGRanges)
    # which of the ranges did the position hit
    hit <- GenomicRanges::queryHits(overlap)
    return(hit)
}

#' Determine whether a set of SNPs are in a CNV region or other ranges
#' Super slow...see insideCNVs instead
#'
#' @param posList data frame ($chr / $pos)
#' @param cnvs : data frame ($chr / $start / $end)
#' @return a T/F vector of whether the positions are in the CNVs
#'
#' @examples
#' snpList <- data.frame(
#'     'chr' = c('chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr3'),
#'     'pos' = c(7, 100, 101, 0, 10, 5)
#'     )
#' cnv <- data.frame(
#'     'chr' = c('chr1', 'chr2', 'chr3'),
#'     'startPos' = c(5, 0, 1),
#'     'endPos' = c(100, 20000, 2)
#'     )
#' insideRanges(snpList, cnv)
#'
insideRange <- function(posList, range) {
    hits <- unlist(bplapply(seq_len(nrow(posList)), function(i) {
        pos <- posList[i,]
        hit <- pointsWithin(pos, ranges)
        # if hit ie. within range, then will have length 1
        # if no hit, then will have length 0
        length(hit)!=0
    }))
    return(hits)
}

#' Returns a vector of T/F depending on whether SNPs are in CNV regions
#'
#' @param posList data frame ($chr / $pos)
#' @param cnvs : data frame ($chr / $start / $end)
#' @return a T/F vector of whether the positions are in the CNVs
#'
#' @examples
#' snpList <- data.frame(
#'     'chr' = c('chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr3'),
#'     'pos' = c(7, 100, 101, 0, 10, 5)
#'     )
#' cnv <- data.frame(
#'     'chr' = c('chr1', 'chr2', 'chr3'),
#'     'startPos' = c(5, 0, 1),
#'     'endPos' = c(100, 20000, 2)
#'     )
#' insideCnvs(snpList, cnv)
#'
insideCnvs <- function(posList, cnvs) {
    posGRanges <- pos2GRanges(posList)
    rangesGRanges <- range2GRanges(cnvs)
    overlap <- GenomicRanges::findOverlaps(rangesGRanges, posGRanges)
    # which of the ranges did the position hit
    hit <- rep(FALSE, nrow(posList))  # initialize
    hit[GenomicRanges::subjectHits(overlap)] <- TRUE
    return(hit)
}

# returns a factor vector across pos (chr:pos vector)
# indicating if the pos is in one of the genes (row number)

#' In this example gtf, the chromosomes do not have the 'chr' prefix
#' posList <- data.frame(
#'     'chr' = c('1', '1', '1', '2', '2', '3'),
#'     'pos' = c(11869+7, 11869+100, 14363+101, 29554+0, 52473+10, 62948+5)
#'     )
#' file <- 'Homo_sapiens.GRCh37.75.gtf'
#' gtf <- read.table(gtfFile, header=F, stringsAsFactors=F, sep='\t')
geneFactors <- function(posList, gtf, fill=TRUE, gene=TRUE) {
    nam <- paste(posList[,1], posList[,2], sep=':') # preserve names
    posGRanges <- pos2GRanges(posList)
    # gene level or exon level
    if(gene) {
        gtf <- gtf[gtf[,3]=="gene", c(1,4,5)]
    } else {
        gtf <- gtf[gtf[,3]=="exon", c(1,4,5)]
    }
    rangesGRanges <- range2GRanges(gtf)
    overlap <- GenomicRanges::findOverlaps(rangesGRanges, posGRanges)
    # initialize to NAs
    hit <- rep(NA, nrow(posList))
    # get SNP positions that are in genes
    posInGenes <- GenomicRanges::subjectHits(overlap)
    # get which row is associated with those genes
    genesHit <- GenomicRanges::queryHits(overlap)
    hit[posInGenes] <- genesHit
    names(hit) <- nam
    # for those that are still NAs (poor annotation?), assing unique ids
    if(fill) {
        print('NAs:')
        print(table(is.na(hit)))
        hit[is.na(hit)] <- nrow(gtf) + seq(1, sum(is.na(hit)))
    }
    else {
        hit <- na.omit(hit) # remove snps that do not fall in genes (due to poorly annotated UTRs or unannoated genes)
    }
    return(hit)
}
