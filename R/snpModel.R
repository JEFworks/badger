# Functionalities related to snp model


#' Get alternative allele count for positions of interest
#'
#' @param gr GenomicRanges object for positions of interest
#' @param bamFile bam file
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @return
#'   refCount reference allele count information for each position of interest
#'   altCount alternative allele count information for each position of interest
#'
#' @examples
#' \dontrun{
#' # Sites of interest (chr1:4600000, chr2:2000)
#' gr <- GRanges(c('chr1', 'chr2'), IRanges(start=c(4600000, 2000), width=1))
#' # we can get the coverage at these SNP sites from our bams
#' path <- '../data-raw/bams/'
#' files <- list.files(path = path)
#' files <- files[grepl('.bam$', files)]
#' alleleCounts <- lapply(files, function(f) {
#'    bamFile <- paste0(path, f)
#'    indexFile <- paste0(path, paste0(f, '.bai'))
#'    getAlleleCount(gr, bamFile, indexFile)
#' })
#' altCounts <- do.call(cbind, lapply(1:length(gr), function(i) alleleCounts[[i]][[1]]))
#' refCounts <- do.call(cbind, lapply(1:length(gr), function(i) alleleCounts[[i]][[2]]))
#' colnames(altCounts) <- colnames(refCounts) <- files
#' }
#'
#' @export
#'
getAlleleCount <- function (gr, bamFile, indexFile, verbose = FALSE) {
    names <- paste(gr)
    if (verbose) {
        print("Getting allele counts for...")
        print(names)
    }
    
    pp <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = TRUE, max_depth = 1e+07, min_base_quality = 20, min_mapq = 10)
    if (verbose) {
        print("Getting pileup...")
    }
    pu <- pileup(file = bamFile, index = indexFile, scanBamParam = ScanBamParam(which = gr), pileupParam = pp)
    
    if (verbose) {
        print("Getting allele read counts...")
    }
    refCount <- unlist(lapply(seq_along(names), function(i) {
        b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(data.frame(gr$REF)$value[i])
        if(length(b)==0) { return(0) } # neither allele observed
        else if(sum(b)==0) { return(0) } # alt allele observed only
        else { return(pu[pu$which_label==names[i],]$count[b]) }
    }))
    altCount <- unlist(lapply(seq_along(names), function(i) {
        b = as.character(pu[pu$which_label==names[i],]$nucleotide) == as.character(data.frame(gr$ALT)$value[i])
        if(length(b)==0) { return(0) } # neither allele observed
        else if(sum(b)==0) { return(0) } # ref allele observed only
        else { return(pu[pu$which_label==names[i],]$count[b]) }
    }))
    names(refCount) <- names(altCount) <- names
    
    if (verbose) {
        print("Done!")
    }
    return(list(refCount, altCount))
}


#' Get coverage count for positions of interest
#'
#' @param gr GenomicRanges object for positions of interest 
#' @param bamFile bam file
#' @param indexFile bai index file
#' @param verbose Boolean of whether or not to print progress and info
#' @return totCount Total coverage count information for each position of interest
#'
#' @examples
#' \dontrun{
#' # Sites of interest (chr1:4600000, chr2:2000)
#' gr <- GRanges(c('chr1', 'chr2'), IRanges(start=c(4600000, 2000), width=1))
#' # we can get the coverage at these SNP sites from our bams
#' path <- '../data-raw/bams/'
#' files <- list.files(path = path)
#' files <- files[grepl('.bam$', files)]
#' cov <- do.call(cbind, lapply(files, function(f) {
#'     bamFile <- paste0(path, f)
#'     indexFile <- paste0(path, paste0(f, '.bai'))
#'     getCoverage(gr, bamFile, indexFile)
#' }))
#' colnames(cov) <- files
#' }
#'
#' @export
#'
getCoverage <- function (gr, bamFile, indexFile, verbose = FALSE) {
    names <- paste(gr)
    if (verbose) {
        print("Getting coverage for...")
        print(names)
    }

    pp <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = FALSE, max_depth = 1e+07, min_base_quality = 20, min_mapq = 10)
    if (verbose) {
        print("Getting pileup...")
    }
    pu <- pileup(file = bamFile, index = indexFile, scanBamParam = ScanBamParam(which = gr), pileupParam = pp)
    rownames(pu) <- pu$which_label
    if (verbose) {
        print("Getting coverage counts...")
    }
    totCount <- pu[names, ]$count
    totCount[is.na(totCount)] <- 0
    names(totCount) <- names
    if (verbose) {
        print("Done!")
    }
    return(totCount)
}



#' Helper function to get coverage and allele count matrices given a set of putative heterozygous SNP positions
#'
#' @param snps GenomicRanges object for positions of interest
#' @param bamFiles list of bam file
#' @param indexFiles list of bai index file
#' @param n.cores number of cores
#' @param verbose Boolean of whether or not to print progress and info
#' @return
#'   refCount reference allele count matrix for each cell and each position of interest
#'   altCount alternative allele count matrix for each cell and each position of interest
#'   cov total coverage count matrix for each cell and each position of interest
#' 
#' @examples
#' \dontrun{
#' # Get putative hets from ExAC
#' vcfFile <- "../data-raw/ExAC.r0.3.sites.vep.vcf.gz"
#' testRanges <- GRanges(chr, IRanges(start = 1, width=1000))
#' param = ScanVcfParam(which=testRanges)
#' vcf <- readVcf(vcfFile, "hg19", param=param)
#' ## common snps by MAF
#' info <- info(vcf)
#' if(nrow(info)==0) {
#'     if(verbose) {
#'         print("ERROR no row in vcf")
#'     }
#'     return(NA)
#' }
#' maf <- info[, 'AF'] # AF is Integer allele frequency for each Alt allele
#' if(verbose) {
#'     print(paste0("Filtering to snps with maf > ", maft))
#' }
#' vi <- sapply(maf, function(x) any(x > maft))
#' if(verbose) {
#'     print(table(vi))
#' }
#' snps <- rowRanges(vcf)
#' snps <- snps[vi,]
#' ## get rid of non single nucleotide changes
#' vi <- width(snps@elementMetadata$REF) == 1
#' snps <- snps[vi,]
#' ## also gets rid of sites with multiple alt alleles though...hard to know which is in our patient
#' vi <- width(snps@elementMetadata$ALT@partitioning) == 1
#' snps <- snps[vi,]
#' ## Get bams
#' files <- list.files(path = "../data-raw")
#' bamFiles <- files[grepl('.bam$', files)]
#' bamFiles <- paste0(path, bamFiles)
#' indexFiles <- files[grepl('.bai$', files)]
#' indexFiles <- paste0(path, indexFiles)
#' results <- getSnpMats(snps, bamFiles, indexFiles)
#' }
#'
#' @export
#' 
getSnpMats <- function(snps, bamFiles, indexFiles, n.cores=10, verbose=FAlSE) {

    ## loop
    cov <- do.call(cbind, mclapply(seq_along(bamFiles), function(i) {
        bamFile <- bamFiles[i]
        indexFile <- indexFiles[i]
        getCoverage(snps, bamFile, indexFile, verbose)
    }, mc.cores=n.cores))
    colnames(cov) <- bamFiles

    ## any coverage?
    if(verbose) {
        print("Snps with coverage:")
        print(table(rowSums(cov)>0))
    }
    vi <- rowSums(cov)>0; table(vi)
    cov <- cov[vi,]
    snps <- snps[vi,]

    if(verbose) {
        print("Getting allele counts...")
    }
    alleleCount <- mclapply(seq_along(bamFiles), function(i) {
        bamFile <- bamFiles[i]
        indexFile <- indexFiles[i]
        getAlleleCount(snps, bamFile, indexFile, verbose)
    }, mc.cores=n.cores)
    refCount <- do.call(cbind, lapply(alleleCount, function(x) x[[1]]))
    altCount <- do.call(cbind, lapply(alleleCount, function(x) x[[2]]))
    colnames(refCount) <- colnames(altCount) <- bamFiles

    ## check correspondence
    if(verbose) {
        print("altCount + refCount == cov:")
        print(table(altCount + refCount == cov))
        print("altCount + refCount < cov: sequencing errors")
        print(table(altCount + refCount < cov))
        ##vi <- which(altCount + refCount != cov, arr.ind=TRUE)
        ## some sequencing errors evident
        ##altCount[vi]
        ##refCount[vi]
        ##cov[vi]
    }

    results <- list(refCount, altCount, cov)
    return(results)
}


#' Composite Minor Allele Frequency (CLAF) Profile plot
#'
#' @param r Alt allele count in single cells
#' @param n.sc Coverage in single cells
#' @param l Alt allele count in bulk
#' @param n.bulk Coverage in bulk
#' @param region Restrict plotting to select region. Optional. Default: NULL
#' @param filter Boolean of whether to filter for SNPs with coverage. Default: TRUE
#' @param delim Delimiter for names of SNPs as Chromosome[delim]Position. Default: ":" ex. chr1:283838897
#' @param gtf GTF file contents for mapping SNPs to genes. Required if plotGene = TRUE
#' @param plotGene Boolean of whether to plot gene track. Default: FALSE
#'
#'
#' @examples
#' ## Single cell data
#' data(snpsHet_MM16ScSample)
#' ## Bulk exome
#' data(snpsHet_MM16BulkSample)
#' # intersect
#' region <- data.frame('chr'=2, start=0, end=1e9) # deletion region
#' clafProfile(r, cov.sc, l, cov.bulk, region)
#' region <- data.frame('chr'=3, start=0, end=1e9) # neutral region
#' clafProfile(r, cov.sc, l, cov.bulk, region)
#'
#' @export
#'
clafProfile<- function(r, n.sc, l, n.bulk, region=NULL, filter=TRUE, delim=':', gtf=NULL, plotGene=FALSE) {

    if(filter) {
        #####
        # Clean
        #####

        # filter out snps without coverage
        s <- rowSums(n.sc) > 0
        r <- r[s,]
        n.sc <- n.sc[s,]
        l <- l[s]
        n.bulk <- n.bulk[s]
    }

    if(!is.null(region)) {
        ######
        ## Deletion regions
        ######
        print('localing snps to deletion region...')

        snps <- rownames(n.sc)
        vi <- insideCnvs(snp2df(snps, delim=delim), region)
        chr.snps <- snps[vi]

        print('number of snps:')
        print(length(chr.snps))

        # order snps
        chr.snps.order <- snp2df(chr.snps, delim=delim)[,2]
        chr.snps <- chr.snps[order(chr.snps.order)]

        # restrict
        r <- r[chr.snps,]
        n.sc <- n.sc[chr.snps,]
        l <- l[chr.snps]
        n.bulk <- n.bulk[chr.snps]
    }

    if(plotGene) {
        print('map snps to genes')
        snpsName <- snp2df(chr.snps, delim=delim)
        gf <- geneFactors(snpsName, gtf, fill=TRUE) # keep all snvs for now
    }

    ######
    ## Merge bulk and single cell
    ######

    r.tot <- cbind(r, 'Bulk'=l)
    n.tot <- cbind(n.sc, 'Bulk'=n.bulk)

    ## condense
    #no.info <- n.tot==0
    #r.cond <- do.call(cbind, lapply(1:ncol(r.tot), function(i) {
    #    ni <- no.info[,i]
    #    r.cond <- rep(0, nrow(r.tot))
    #    if(sum(!ni)>0) {
    #        r.cond[1:sum(!ni)] <- r.tot[!ni, i]
    #    }
    #    r.cond
    #}))
    #dim(r.cond)
    #n.cond <- do.call(cbind, lapply(1:ncol(n.tot), function(i) {
    #    ni <- no.info[,i]
    #    n.cond <- rep(0, nrow(n.tot))
    #    if(sum(!ni)>0) {
    #        n.cond[1:sum(!ni)] <- n.tot[!ni, i]
    #    }
    #    n.cond
    #    }))
    #dim(n.cond)
    #colnames(r.cond) <- colnames(n.cond) <- colnames(r.tot)

    ######
    ## Convert to frac
    ######

    visualize.mat <- function(r, n.sc, E = l/n.bulk) {
        n <- nrow(r)
        m <- ncol(r)
        mat <- do.call(rbind, lapply(1:n, function(i) {
            do.call(cbind, lapply(1:m, function(j) {
                ri <- r[i,j]
                n.sci <- n.sc[i,j]
                Ei <- E[i]
                if(is.na(Ei)) {
                    mut.frac <- NA
                }
                else if(Ei <= 0.5) {
                    mut.frac <- ri/n.sci
                }
                else if(Ei > 0.5) {
                    mut.frac <- 1-ri/n.sci
                }
                else {
                    mut.frac <- NA
                }

                ## f will be high if inconsistent
                ## f will be low if consistent
                ## f will be NaN if no coverage
                ## use colorRamp from green to red
                f <- mut.frac
                return(f)
            }))
        }))
        rownames(mat) <- rownames(r)
        colnames(mat) <- colnames(r)
        return(mat)
    }
    mat.tot <- visualize.mat(r.tot, n.tot, E = l/n.bulk)
    head(mat.tot)

    ## condence
    #mat.cond <- do.call(cbind, lapply(1:ncol(mat.tot), function(i) {
    #    ni <- no.info[,i]
    #    c <- rep(NA, nrow(mat.tot))
    #    if(sum(!ni)>0) {
    #        c[1:sum(!ni)] <- mat.tot[!ni, i]
    #    }
    #    c
    #}))
    #colnames(mat.cond) <- colnames(mat.tot)
    #rownames(mat.cond) <- rownames(mat.tot)
    #head(mat.cond)
    #dim(mat.cond)

    if(plotGene) {
        # add in gene info
        geneInfo <- as.integer(as.factor(gf[rownames(mat.tot)]))
        # label filled genes
        r <- range(gf[rownames(mat.tot)])
        geneInfo[gf>(r[2]-100)] <- 0
        # hacky method to make alternating gene colors
        is.even <- function(x) x %% 2 == 0
        geneInfo[is.even(geneInfo)] <- geneInfo[is.even(geneInfo)]*10000
        geneInfo[geneInfo==0] <- NA
        geneInfo[geneInfo < 1000] <- 0
        geneInfo[geneInfo > 1000] <- 1
        names(geneInfo) <- rownames(mat.tot)
    }

    # order?
    #order <- names(sort(colSums(mat.tot>0, na.rm=TRUE), decreasing=FALSE))
    #order <- names(sort(colSums(mat.tot, na.rm=TRUE), decreasing=FALSE))
    #mat.tot <- mat.tot[, order]
    #n.tot <- n.tot[, order]

    ######
    ## Plot
    ######
    require(ggplot2)
    require(reshape2)

    m <- melt(t(mat.tot))
    colnames(m) <- c('cell', 'snp', 'alt.frac')
    rownames(m) <- paste(m$cell, m$snp)
    m$alt.frac[is.nan(m$alt.frac)] <- NA
    n <- melt(t(n.tot))
    colnames(n) <- c('cell', 'snp', 'coverage')
    rownames(n) <- paste(n$cell, n$snp)
    n$coverage[n$coverage>30] <- 30  # max for visualization purposes
    #n$coverage <- log10(n$coverage+1)
    n$coverage <- n$coverage^(1/3) # cube root for visualization purposes only
    dat <- cbind(m, coverage=n$coverage)

    # along region
    p <- ggplot(dat, aes(snp, cell)) +
        # geom_tile(alpha=0) +
        geom_point(aes(colour = alt.frac, size = coverage)) +
            scale_size_continuous(range = c(0,3)) +
                # scale_colour_gradientn(colours = rainbow(10)) +
                scale_colour_gradient2(mid="yellow", low = "turquoise", high = "red", midpoint=0.5) +
                    theme(
                        # axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5),lineheight=1),
                        # axis.text.y=element_blank(),
                        axis.title.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        #axis.text.y=element_text(size=rel(0.5))
                        legend.position="bottom"
                        #panel.margin=unit(0 , "lines")
                        )
    #print(p)

    if(plotGene) {
        gdat <- melt(t(geneInfo))
        colnames(gdat) <- c('Var1', 'snp', 'gf')
        g <- ggplot(gdat, aes(snp, Var1)) +
            geom_point(aes(colour = gf, size=2))

        #print(g)

        # plot together
        grid.arrange(p, g, nrow=2, heights=c(5,1))
    } else {
        print(p)
    }

    # stack to more easily visualize sparse regions
    #m <- melt(t(mat.cond))
    #colnames(m) <- c('cell', 'snp', 'alt.frac')
    #rownames(m) <- paste(m$cell, m$snp)
    #m$alt.frac[is.nan(m$alt.frac)] <- NA
    #n <- melt(t(n.cond))
    #colnames(n) <- c('cell', 'snp', 'coverage')
    #rownames(n) <- paste(n$cell, n$snp)
    #n$coverage[n$coverage>30] <- 30
    ##n$coverage <- log10(n$coverage+1)
    #n$coverage <- n$coverage^(1/3)
    #dat2 <- cbind(m, coverage=n$coverage)

    #p2 <- ggplot(dat2, aes(snp, cell)) +
    #    #    geom_tile(alpha=0) +
    #    geom_point(aes(colour = alt.frac, size = coverage, alpha=0.5)) +
    #        scale_size_continuous(range = c(0,3)) +
    #            #                scale_colour_gradientn(colours = rainbow(10)) +
    #            scale_colour_gradient2(mid="yellow", low = "turquoise", high = "red", midpoint=0.5) +
#   #                 theme_bw() +
    #                    theme(
    #                        #                            axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5),lineheight=1),
    #                        axis.text.y=element_blank(),
    #                        axis.title.y=element_blank(),
    #                        axis.ticks.y=element_blank(),
    #                        #axis.text.y=element_text(size=rel(0.5))
    #                        legend.position="bottom",
#   #                         panel.margin=unit(0 , "lines")
    #                        ) + guides(alpha=FALSE)
    #p2

    #grid.arrange(p2, p, ncol=2, widths=c(1,4))

}


#' Run BADGER allele-only model to assess posterior probability of CNVs given allele information only
#'
#' @param r Matrix of alt allele count in single cells
#' @param cov.sc Matrix of coverage in single cells
#' @param l Vector of alt allele count in bulk
#' @param cov.bulk Vector of coverage in bulk
#' @param region Region of interest such as expected CNV boundaries
#' @param gtf GTF file contents for mapping SNPs to genes
#' @param mono Rate of mono-allelic expression. Default: 0.7
#' @param pe Effective error rate to capture error from sequencing, etc. Default: 0.01
#' @param filter Boolean for whether to filter out SNP sites with no coverage. Default: TRUE
#' @param likelihood Boolean for whether to use likelihood based estimate of posterior. Default: FALSE
#' @param n.iter Number of iterations in MCMC. Default: 1000
#' @param quiet Boolean of whether to suppress progress bar. Default: TRUE
#' @param delim Delimiter for names of SNPs as Chromosome[delim]Position. Default: ":" ex. chr1:283838897
#'
#' @return Posterior probability of deletion or LOH
#'
#' @examples
#' ## Single cell data
#' data(snpsHet_MM16ScSample)
#' ## Bulk exome
#' data(snpsHet_MM16BulkSample)
#' \dontrun{
#' gtfFile <- 'data-raw/Homo_sapiens.GRCh37.75.gtf'
#' gtf <- read.table(gtfFile, header=F, stringsAsFactors=F, sep='\t')
#' }
#' region <- data.frame('chr'=2, start=0, end=1e9) # deletion region
#' \dontrun{
#' results <- calcAlleleCnvProb(r, cov.sc, l, cov.bulk, region, gtf)
#' }
#' region <- data.frame('chr'=3, start=0, end=1e9) # neutral region
#' \dontrun{
#' results <- calcAlleleCnvProb(r, cov.sc, l, cov.bulk, region, gtf)
#' }
#'
#' @export
#'
calcAlleleCnvProb <- function(r, cov.sc, l, cov.bulk, region, gtf, mono = 0.7, pe = 0.01, filter=TRUE, likelihood=FALSE, n.iter=1000, quiet=TRUE, delim=':') {

    #####
    # Clean
    #####

    if(filter) {
        # filter out snps without coverage
        s <- rowSums(cov.sc) > 0
        r <- r[s,]
        cov.sc <- cov.sc[s,]
        l <- l[s]
        cov.bulk <- cov.bulk[s]
    }

    ######
    ## Deletion regions
    ######
    if(!is.null(region)) {
        print('localing snps to region...')

        snps <- rownames(cov.sc)
        vi <- insideCnvs(snp2df(snps, delim=delim),region)
        chr.snps <- snps[vi]

        print('number of snps:')
        print(length(chr.snps))
        if(length(chr.snps) <= 1) {
            pm <- rep(NA, ncol(r))
            names(pm) <- colnames(r)
            return(pm)
        }
    }
    else {
        chr.snps <- rownames(cov.sc)
    }

    r <- r[chr.snps,]
    n.sc <- cov.sc[chr.snps,]
    l <- l[chr.snps]
    n.bulk <- cov.bulk[chr.snps]

    ####
    ## Map snps to genes
    ####
    print('mapping snps to genes...')

    # associate each snp with a gene factor
    snpsName <- snp2df(chr.snps, delim=delim)
    snps2genes <- geneFactors(snpsName, gtf, fill=T) # if not found in gtf, assume annotation error; make each unique gene factor
    names(snps2genes) <- chr.snps

    genes.of.interest <- unique(snps2genes)
    print('number of genes:')
    print(length(genes.of.interest))

    # associate each gene factor with a set of snps
    genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
        which(snps2genes %in% genes.of.interest[i])
    })
    names(genes2snps.dict) <- genes.of.interest

    #####
    ## Model
    #####

    print('converting to multi-dimensional arrays...')

    ## Convert to multi-dimensions based on j
    I.j <- unlist(lapply(genes2snps.dict, length))
    numGenes <- length(genes2snps.dict)
    numSnpsPerGene <- max(I.j)
    numCells <- ncol(r)
    ## j, i, k
    r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        snps <- genes2snps.dict[[i]]
        for(s in seq_along(snps)) {
            r.array[i,s,] <- r[snps[s],]
        }
    }
    n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        snps <- genes2snps.dict[[i]]
        for(s in seq_along(snps)) {
            n.sc.array[i,s,] <- n.sc[snps[s],]
        }
    }
    l.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        snps <- genes2snps.dict[[i]]
        for(s in seq_along(snps)) {
            l.array[i,s] <- l[snps[s]]
        }
    }
    n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        snps <- genes2snps.dict[[i]]
        for(s in seq_along(snps)) {
            n.bulk.array[i,s] <- n.bulk[snps[s]]
        }
    }

    print('aggregating data to list...')
    data <- list(
        'l' = l.array,
        'r' = r.array,
        'n.bulk' = n.bulk.array,
        'n.sc' = n.sc.array,
        # 'mu' = mu,
        # 'sigma' = sigma,
        # 'g' = g,
        'J' = length(I.j),  # how many genes
        'K' = ncol(r),  # how many cells
        'I.j' = I.j,
        'pseudo' = pe,
        # 'tau' = 1,
        'mono' = mono)

    modelFile <- system.file("bug", "snpModel.bug", package = "badger")

    print('Initializing model...')
    ## 2 random chains
    #chains <- 2
    #init <- lapply(1:chains, function(i) list('S'=sample(c(0,1), ncol(r), replace=T)))
    ## 2 more chains, 1 with all starting with deletions, 1 with all starting without deletions
    #init <- c(list(list('S'=rep(0, ncol(r))), list('S'=rep(1, ncol(r)))), init)
    #model <- jags.model(modelFile, data=data, inits=init, n.chains=chains+2, n.adapt=0, quiet=quiet)
    model <- jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    # Joe says 4 chains is a standard, so just stick with 4 chains
    update(model, 300, progress.bar=ifelse(quiet,"none","text"))

    print('Running model...')
    print(variable.names(model))
    if(likelihood) {
        parameters <- c('fma', 'h', 'b', 'd')
        samples <- coda.samples(model, parameters, n.iter=n.iter, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine samples across chains

        # likelihood assuming CNV is absent
        pm0 <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
            cnvLik <- sapply(seq_len(numGenes), function(gi) {
                snps <- genes2snps.dict[[gi]]
                geneLik <- sapply(seq_along(snps), function(si) {
                    S <- 0
                    h <- samples[,paste0('h[',gi,',',si,',',ci,']')]
                    b <- samples[,paste0('b[',gi,',',ci,']')]
                    d <- samples[,paste0('d[',gi,',',ci,']')]
                    fma <- samples[,paste0('fma[',gi,',',si,']')]
                    p <- (h*(1-b) + (pe*d + (1-pe)*(1-d))*b)*(1-S) + fma*S
                    snpLik <- dbinom(r.array[gi,si,ci],
                                     n.sc.array[gi,si,ci],
                                     p)
                    mean(snpLik) # take arthmetic mean for now
                })
                #geneLik <- prod(geneLik)
                geneLik <- sum(log(geneLik)) # due to overflow, use sum of log
            })
            #prod(cnvLik)
            exp(sum(cnvLik)) # exponential back
        }))
        pm1 <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
            cnvLik <- sapply(seq_len(numGenes), function(gi) {
                snps <- genes2snps.dict[[gi]]
                geneLik <- sapply(seq_along(snps), function(si) {
                    S <- 1
                    h <- samples[,paste0('h[',gi,',',si,',',ci,']')]
                    b <- samples[,paste0('b[',gi,',',ci,']')]
                    d <- samples[,paste0('d[',gi,',',ci,']')]
                    fma <- samples[,paste0('fma[',gi,',',si,']')]
                    p <- (h*(1-b) + (pe*d + (1-pe)*(1-d))*b)*(1-S) + fma*S
                    snpLik <- dbinom(r.array[gi,si,ci],
                                     n.sc.array[gi,si,ci],
                                     p)
                    mean(snpLik) # take arthmetic mean for now
                })
                #geneLik <- prod(geneLik)
                geneLik <- sum(log(geneLik)) # due to overflow, use sum of log
            })
            #prod(cnvLik)
            exp(sum(cnvLik)) # exponentiate back
        }))
        colnames(pm1) <- colnames(pm0) <- colnames(r)
        pm <- list(pm0, pm1)
    } else {
        parameters <- 'S'
        samples <- coda.samples(model, parameters, n.iter=n.iter, progress.bar=ifelse(quiet,"none","text"))
        samples <- do.call(rbind, samples) # combine samples across chains
        pm <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
            c(mean(samples[,paste("S[",ci,"]",sep="")]))
        }))
        colnames(pm) <- colnames(r)
    }

    print('Complete!')
    return(pm)
}


