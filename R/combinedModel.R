#' Run BADGER combined model to assess posterior probability of CNVs given both allele and expression data
#'
#' @param r Matrix of alt allele count in single cells
#' @param cov.sc Matrix of coverage in single cells
#' @param l Vector of alt allele count in bulk
#' @param cov.bulk Vector of coverage in bulk
#' @param region Region of interest such as expected CNV boundaries
#' @param gtf GTF file contents for mapping SNPs to genes. Required if plotGene = TRUE
#' @param gexp Normalized gene expression matrix
#' @param fits Fit for variance around mean
#' @param gos Gene position table
#' @param m Expected mean deviation due to copy number change
#' @param mono Rate of mono-allelic expression. Default: 0.7
#' @param pe Effective error rate to capture error from sequencing, etc. Default: 0.01
#' @param filter Boolean for whether to filter out SNP sites with no coverage. Default: TRUE
#' @param n.iter Number of iterations in MCMC. Default: 1000
#' @param quiet Boolean of whether to suppress progress bar. Default: TRUE
#' @param delim Delimiter for names of SNPs as Chromosome[delim]Position. Default: ":" ex. chr1:283838897
#' @return List of posterior probabilities for CNV and direction of CNV (deletion vs. amplification)
#'
#' @examples
#' data('MM16.counts')
#' mat <- log2(MM16.counts + 1)
#' data('Normal.counts')
#' mat.ref <- log2(Normal.counts + 1)
#' mats <- normalizedExpression(mat, mat.ref)
#' gexp <- mats[[1]]
#' fits <- mvFit(gexp)
#' region <- data.frame('chr'=1, start=0, end=1e9)
#' data(snpsHet_MM16ScSample)
#' data(snpsHet_MM16BulkSample)
#' \dontrun{
#' gtfFile <- 'data-raw/Homo_sapiens.GRCh37.75.gtf'
#' gtf <- read.table(gtfFile, header=F, stringsAsFactors=F, sep='\t')
#' region <- data.frame('chr'=2, start=0, end=1e9) # deletion region
#' library(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
#' gos <- getBM(values=rownames(mat.tot),attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),filters=c("ensembl_gene_id"),mart=mart.obj)
#' gos$pos <- (gos$start_position + gos$end_position)/2
#' results <- calcCombCnvProb(r, cov.sc, l, cov.bulk, region, gtf, gexp, fits, gos, m=0.15)
#' }
#'
#' @export
calcCombCnvProb <- function(r, cov.sc, l, cov.bulk, region, gtf, gexp, fits, gos, m, filter=TRUE, pe = 0.01, mono = 0.7, n.iter=1000, quiet=TRUE, delim=':') {

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
    ## Region
    ######
    print('localing snps to region...')

    snps <- rownames(cov.sc)
    vi <- insideCnvs(snp2df(snps, delim=delim), region)
    chr.snps <- snps[vi]

    print('number of snps:')
    print(length(chr.snps))
    # temp fix: should be able to accomodate no SNPs if has gene
    if(length(chr.snps) <= 1) {
        pm <- rep(NA, ncol(r)+1)
        names(pm) <- colnames(r)
        return(pm)
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
    print('number of genes with snps:')
    print(length(genes.of.interest))

    # associate each gene factor with a set of snps
    genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
        which(snps2genes %in% genes.of.interest[i])
    })
    names(genes2snps.dict) <- genes.of.interest

    # restrict to genes within region of interest
    restrict <- function(names, region) {
        start <- region$start
        end <- region$end
        gos <- gos[gos$chromosome_name == region$chr,]
        gos$pos <- gos$start_position
        rownames(gos) <- make.unique(gos[,1])
        gos <- gos[names,]
        n <- names[which(gos$pos > start & gos$pos < end)]
        n
    }
    vi <- restrict(rownames(gexp), region)
    print('number of genes expessed:')
    print(length(vi))
    if(length(vi) <= 1) {
        pm <- rep(NA, 2*ncol(r))
        names(pm) <- colnames(r)
        return(pm)
    }
    gexp <- gexp[vi,]
    ng <- nrow(gexp)
    sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))


    ####
    ## Model
    ####

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
        # snp info
        'l' = l.array,
        'r' = r.array,
        'n.bulk' = n.bulk.array,
        'n.sc' = n.sc.array,
        'J' = length(I.j),  # how many genes
        'K' = ncol(r),  # how many cells
        'I.j' = I.j,
        'pseudo' = pe,
        'mono' = mono,
        ## expression
        'gexp' = gexp,
        'JJ' = nrow(gexp),
        'sigma0' = sigma0,
        'mag0' = m
        #'mu0'  = mu0,
        #'sigma0' = sigma0
        #'t' = mu0/2
    )

    modelFile <- system.file("bug", "combinedModel.bug", package = "badger")

    print('Initializing model...')
    # Joe says 4 chains is a standard, so just stick with 4 chains
    model <- jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    update(model, 300, progress.bar=ifelse(quiet,"none","text"))

    parameters <- c('S', 'dd')
    samples <- coda.samples(model, parameters, n.iter=300, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine chains

    samplesSummary <- colMeans(samples)
    return(samplesSummary)
}

