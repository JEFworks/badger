# Sample helpers to parse vcf files to get to SNP info matrices

library(Rsamtools)
library(GenomicFiles)
library(parallel)

# Read VCF
sampleGetSnps <- function() {
    library(VariantAnnotation)
    vcfFile <- "data-raw/MM16.ug.raw.snp.vqsr.vcf.gz"
    scanVcfHeader(vcfFile)

    # scan normal sample for genotype
    param = ScanVcfParam(
        samples=c("MM16_N", "MM16_T"),
        which=GRanges("1", IRanges(40000000, 50000000)))
    vcf = readVcf(vcfFile, "hg19", param=param)

    # get hets only
    hets <- which(geno(vcf)[["GT"]][,"MM16_N"] == "0/1" & geno(vcf)[["GT"]][,"MM16_T"] != ".")
    vcf <- vcf[hets,]

    # filter for high quality
    vi <- elementMetadata(vcf)$FILTER == "PASS"
    vcf <- vcf[vi,]

    # get final info
    snps <- as.data.frame(rowData(vcf))
    alleles <- do.call(rbind, geno(vcf)[["AD"]][,"MM16_T"])
    l <- alleles[,2]
    cov <- rowSums(alleles)

    rownames(snps) <- names(l) <- names(cov) <- paste0(as.character(snps[,1]), ":", as.character(snps[,2]))
}


#' Process MuTect output from WES
#'
#' @param snps.file MuTect .mut file
#' @param out.file RData output file location and name
#'
#' @return outputs positions and allele information for somatic, germline het, germline homo ref and germline homo alt snps
#'  saves output to out.file
#'
process.bulk <- function(snps.file, out.file) {
    # Read in SNPs from exome
    snps <- read.delim(snps.file,header=T,sep="\t",stringsAsFactors=F, comment.char="#")

    # Focus on somatic snps
    somatic.inds <- which(snps$judgement=='KEEP')
    nonsomatic.inds <- which(snps$judgement=='REJECT')
    snps.somatic <- snps[somatic.inds, c('contig', 'position', 'ref_allele', 'alt_allele')]
    print('Number of somatic variants:')
    print(length(somatic.inds))

    # Focus on heterozygous SNPs
    germline.call <- unlist(lapply(snps$normal_best_gt, function(x) length(unique(strsplit(x, '')[[1]]))))
    germline.hets.inds <- which(germline.call==2)  # genotype must have 2 alleles
    germline.hets.inds <- intersect(germline.hets.inds, nonsomatic.inds)  # just in case, make sure non-somatic; shouldn't matter
    snps.het <- snps[germline.hets.inds, c('contig', 'position', 'ref_allele', 'alt_allele')]
    print('Number of heterozygous germline SNPs:')
    print(length(germline.hets.inds))

    # Focus on homozygous SNPs
    germline.homo.inds <- which(germline.call==1)  # genotype must have 1 allele
    # reference
    germline.homo.ref <- intersect(which(unlist(lapply(snps$normal_f, function(x) x==0))), germline.homo.inds)  # alt allele fraction = 0
    germline.homo.ref <- intersect(germline.homo.ref, nonsomatic.inds)  # just in case, make sure non-somatic; shouldn't matter
    snps.homo.ref <- snps[germline.homo.ref, c('contig', 'position', 'ref_allele', 'alt_allele')]
    print('Number of homozygous reference germline SNPs:')
    print(length(germline.homo.ref))
    # alt
    germline.homo.alt <- intersect(which(unlist(lapply(snps$normal_f, function(x) x==1))), germline.homo.inds)  # alt allele fraction = 1
    germline.homo.alt <- intersect(germline.homo.alt, nonsomatic.inds)  # just in case, make sure non-somatic; shouldn't matter
    snps.homo.alt <- snps[germline.homo.alt, c('contig', 'position', 'ref_allele', 'alt_allele')]
    print('Number of homozygous alternate germline SNPs:')
    print(length(germline.homo.alt))

    save(snps.somatic, snps.het, snps.homo.ref, snps.homo.alt, file=out.file)
    print('Done!')
}

#' Get alternative allele count for positions of interest
#'
#' @param posDf data.frame with positions as chr:pos, first column as reference amino acid, second column as alternative amino acid
#' @param bamFile bam file
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @return
#'   altAlleleCount alternative allele count information for each position of interest
#'   refAlleleCount reference allele count information for each position of interest
#'
getAlleleCount <- function(alleleInfo, bamFile, indexFile, verbose=F) {

    # Split posName into components for GRanges
    chrs <- alleleInfo[,1]
    pos <- as.numeric(alleleInfo[,2])
    names <- paste(alleleInfo[,1], alleleInfo[,2], sep=":")

    refNames <- paste(alleleInfo[,1], alleleInfo[,2], alleleInfo[,3], sep=":")
    altNames <- paste(alleleInfo[,1], alleleInfo[,2], alleleInfo[,4], sep=":")

    if (verbose) {
        print("Getting coverage for...")
        print(names)
    }

    # Set pileup options
    # no max depth
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=TRUE,
        max_depth=10000000,
        min_base_quality=0,
        min_mapq=20
    )
    # Positions of interest
    gr <- GRanges(
        seqnames = chrs,
        IRanges(pos, width=1)  # Should all be SNVs so width=1
    )

    if (verbose) {
        print("Getting pileup...")
    }

    # Get pileup
    pu <- pileup(
        file=bamFile,
        index=indexFile,
        scanBamParam=ScanBamParam(which=gr),
        pileupParam=pp
    )
    rownames(pu) <- paste(pu$seqnames, pu$pos, pu$nucleotide, sep=':')  # Create unique identifiers

    if (verbose) {
        print("Getting allele read counts...")
    }

    # Pileup only returns non-zero read counts so fill in those that have no info
    altCount <- pu[altNames, ]$count
    altCount[is.na(altCount)] <- 0
    names(altCount) <- names

    # Pileup only returns non-zero read counts so fill in those that have no info
    refCount <- pu[refNames, ]$count
    refCount[is.na(refCount)] <- 0
    names(refCount) <- names

    print("Done!")
    return(list(refCount, altCount))
}

getCoverage <- function(alleleInfo, bamFile, indexFile, verbose=F) {

    # Split posName into components for GRanges
    chrs <- alleleInfo[,1]
    pos <- as.numeric(alleleInfo[,2])
    names <- paste(alleleInfo[,1], alleleInfo[,2], sep=":")

    if (verbose) {
        print("Getting coverage for...")
        print(names)
    }

    # Set pileup options
    # Do not distinguish between strands or nucleotides
    # no max depth
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        max_depth=10000000,
        min_base_quality=0,
        min_mapq=20
    )
    # Positions of interest
    gr <- GRanges(
        seqnames = chrs,
        IRanges(pos, width=1)  # Should all be SNVs so width=1
    )

    if (verbose) {
        print("Getting pileup...")
    }

    # Get pileup
    pu <- pileup(
        file=bamFile,
        index=indexFile,
        scanBamParam=ScanBamParam(which=gr),
        pileupParam=pp
    )
    rownames(pu) <- paste(pu$seqnames, pu$pos, sep=':')  # Create unique identifiers

    if (verbose) {
        print("Getting coverage counts...")
    }

    # Pileup only returns non-zero read counts so fill in those that have no info
    totCount <- pu[names,]$count
    totCount[is.na(totCount)] <- 0
    names(totCount) <- names

    print("Done!")
    return(totCount)
}

# helper function
getMats <- function(snps, bamFiles, baiFiles, bamNames, out.file) {
    snps.info <- mclapply(seq(along=bamFiles), function(i) {
        bam = bamFiles[i]
        ind = baiFiles[i]
        cov <- getCoverage(snps, bam, ind)
        ref.alt <- getAlleleCount(snps, bam, ind)
        ref <- ref.alt[[1]]
        alt <- ref.alt[[2]]
        return(list(ref, alt, cov))
    }, mc.cores = 10)
    snps.ref <- do.call(cbind, mclapply(seq(along=bamFiles), function(i) {
        snps.info[[i]][[1]]
    }, mc.cores = 10))
    snps.alt <- do.call(cbind, mclapply(seq(along=bamFiles), function(i) {
        snps.info[[i]][[2]]
    }, mc.cores = 10))
    snps.cov <- do.call(cbind, mclapply(seq(along=bamFiles), function(i) {
        snps.info[[i]][[3]]
    }, mc.cores = 10))
    colnames(snps.ref) <- colnames(snps.alt) <- colnames(snps.cov) <- bamNames
    head(snps.cov)
    save(snps.ref, snps.alt, snps.cov, file=out.file)
}
