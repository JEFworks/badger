################################################################### Combined
calcCombCnvProb <- function(r, cov.sc, l, cov.bulk, region, gexp, fits, pe = 0.01, mono = 0.7, n.iter=1000, quiet=TRUE, delim=':') {

    #####
    # Clean
    #####

    # filter out snps without coverage
    s <- rowSums(cov.sc) > 0
    r <- r[s,]
    cov.sc <- cov.sc[s,]

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
    snps2genes <- geneFactors(snpsName, gtfFile, fill=T) # if not found in gtf, assume annotation error; make each unique gene factor
    names(snps2genes) <- chr.snps

    genes.of.interest <- unique(snps2genes)
    print('number of genes with snps:')
    print(length(genes.of.interest))

    # associate each gene factor with a set of snps
    genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
        which(snps2genes %in% genes.of.interest[i])
    })
    names(genes2snps.dict) <- genes.of.interest

    # possible that gene does not have snp
    restrict <- function(names, region) {
        start <- region$start
        end <- region$end
        gos <- getBM(values=names,attributes=c("ensembl_transcript_id", "hgnc_symbol", "chromosome_name","start_position","end_position"),filters=c("hgnc_symbol"),mart=mart.obj)
        gos <- gos[gos$chromosome_name == region$chr,]
        gos$pos <- gos$start_position
        rownames(gos) <- make.unique(gos$hgnc_symbol)
        gos <- gos[names,]
        n <- names[which(gos$pos > start & gos$pos < end)]
        n
    }
    # get all genes expressed in region
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

    modelFile <- 'bug/combinedModel.bug'

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

