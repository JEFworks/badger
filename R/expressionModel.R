# Functionalities related to expression model


#' Normalize gene expression to a normal reference for expression-based karyotyping
#'
#' @param mat Test gene expression matrix to be normalized
#' @param mat.ref Normal reference gene expression matrix
#' @param minMeanBoth Filter out genes with mean gene expression below this threshold in both samples. Default: 4.5
#' @param minMeanTest Filter out genes with mean gene expression below this threshold in the test gene expression matrix. Default: 6
#' @param minMeanRef Filter out genes with mean gene expression below this threshold in the reference gene expression matrix. Default: 8
#' @param verbose Boolean for whether to print out verbose warnings. Default: TRUE
#' @return List of two normalized gene expression matrices corresponding to mat and mat.ref
#' @examples
#' data('MM16.counts')
#' mat <- log2(MM16.counts + 1)
#' data('Normal.counts')
#' mat.ref <- log2(Normal.counts + 1)
#' mats <- normalizedExpression(mat, mat.ref)
#'
#' @export
#'
normalizedExpression <- function(mat, mat.ref, minMeanBoth = 4.5, minMeanTest = 6, minMeanRef = 8, verbose=TRUE) {

    ## filter
    vi <- intersect(rownames(mat), rownames(mat.ref))
    if(verbose) {
        print("Number of shared genes:")
        print(length(vi))
        if(length(vi) < 1) {
            print("WARNING: Check if your gene names match")
        }
    }
    mat <- mat[vi,]
    mat.ref <- mat.ref[vi,]

    vi <- (rowMeans(mat) > minMeanBoth & rowMeans(mat.ref) > minMeanBoth) | rowMeans(mat) > minMeanTest | rowMeans(mat.ref) > minMeanRef
    if(verbose) {
        print("Number of genes passed filtering:")
        print(table(vi))
    }
    mat <- mat[vi,]
    mat.ref <- mat.ref[vi,]

    ## normalize
    ## library size
    mat <- scale(mat)
    mat.ref <- scale(mat.ref)
    ## to reference
    ref.means <- rowMeans(mat.ref)
    mat <- mat - ref.means
    mat.ref <- mat.ref - ref.means

    return(list(mat, mat.ref))
}


#' Visualize expression heatmap by chromosome for expression-based karyotyping
#'
#' @param mat.tot Gene expression matrix
#' @param gos Dataframe of positional information for each gene in the gene expression matrix
#' @param zlim Lower and upper bounds for heatmap. Default: c(-2,2)
#' @param window.size Window size for sliding window average. Default: 101
#' @param autosomesOnly Boolean for whether to plot only autosomes. Default: TRUE, if false will also plot X chromosome
#' @param orderCells Booean for whether to automaticaly cluster and order cells. Default: FALSE
#' @return List of sliding window averaged expression by chromosome along with heatmap representation
#' @examples
#' data('MM16.counts')
#' mat <- log2(MM16.counts + 1)
#' data('Normal.counts')
#' mat.ref <- log2(Normal.counts + 1)
#' mats <- normalizedExpression(mat, mat.ref)
#' mat.tot <- cbind(mats[[1]], mats[[2]])
#' library(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
#' gos <- getBM(values=rownames(mat.tot),attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),filters=c("ensembl_gene_id"),mart=mart.obj)
#' gos$pos <- (gos$start_position + gos$end_position)/2
#' tl <- plotExpHeatmap(mat.tot, gos, zlim=c(-0.5,0.5), window.size = 201)
#'
#' @export
#'
plotExpHeatmap <- function(mat.tot, gos, zlim=c(-2,2), window.size = 101, autosomesOnly=TRUE, orderCells=FALSE) {

    ## organize into chromosomes
    tl <- tapply(1:nrow(gos),as.factor(gos$chromosome_name),function(ii) {
        na.omit(mat.tot[gos[,1][ii[order(gos$pos[ii],decreasing=F)]],])
    })
    ## only care about these chromosomes
    tl <- tl[c(as.factor(1:22), "X")]
    if(autosomesOnly) {
        tl <- tl[as.factor(1:22)]
    }

    if(orderCells) {
        avgd <- do.call(rbind, lapply(names(tl),function(nam) {
            d <- tl[[nam]]
            d <- colMeans(d)
            d
        }))
        hc <- hclust(dist(t(avgd)))
        order <- hc$order
    }

    require(RColorBrewer)
    pcol <- rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
    # https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
    chr.sizes <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)
    l <- layout(matrix(seq(1, length(tl)),1,length(tl),byrow=T), widths=chr.sizes/1e7)
    zlim <- zlim
    window.size <- window.size

    lapply(names(tl),function(nam) {
        d <- tl[[nam]]
        if(orderCells) { d <- d[, order] }
        d <- apply(d,2,caTools::runmean,k=window.size, align="center")
        d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
        par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
        image(1:nrow(d),1:ncol(d),d,col=pcol,zlim=zlim,xlab="",ylab="",axes=F,main=nam); box()
    })

    return(tl)
}


#' Estimate expression magnitude variance as a function of number of genes
#'
#' @param gexp Normalized gene expression matrix, preferably restricted to region of known neutral copy number
#' @param num.genes List of number of genes for random sampling. Default: seq(5, 100, by=5)
#' @param rep Number of replicates. Default: 1000
#' @param plot Boolean of whether to plot. Default: FALSE
#' @return Fit for variance around mean normalized gene expression as a function of the number of genes for each cell
#'
#' @examples
#' data('MM16.counts')
#' mat <- log2(MM16.counts + 1)
#' data('Normal.counts')
#' mat.ref <- log2(Normal.counts + 1)
#' mats <- normalizedExpression(mat, mat.ref)
#' gexp <- mats[[1]]
#' fits <- mvFit(gexp)
#'
#' @export
#'
mvFit <- function(gexp, num.genes = seq(5, 100, by=5), rep = 1000, plot=FALSE) {

    mean.var.comp <- lapply(num.genes, function(ng) {
        set.seed(0)
        # 50 replicates
        m <- do.call(rbind, lapply(1:rep, function(i) {
            nrmchr.sub <- gexp[sample(1:nrow(gexp), ng),]
            nm <- apply(nrmchr.sub, 2, mean)
            nm
        }))
        return(m)
    })
    names(mean.var.comp) <- num.genes

    fits <- lapply(1:ncol(gexp), function(k) {
        # for one cell
        mean.comp <- do.call(cbind, lapply(mean.var.comp, function(x) x[,k]))

        if(plot) {
            perf.test <- function(mat) {
                require(ggplot2)
                require(reshape2)
                m <- melt(mat)
                p <- ggplot(m) + geom_boxplot(aes(x = factor(Var2), y = value))
                return(p)
            }
            perf.test(mean.comp)
        }

        if(plot) {
            plot(log10(num.genes),log10(apply(mean.comp, 2, var)), type="l")
        }

        df <- data.frame('x'=num.genes, 'y'=apply(mean.comp, 2, var))
        fit <- lm(log10(y)~log10(x), data=df)

        if(plot) {
            x2 <- log10(num.genes)
            y2 <- predict(fit, x=x2, interval="predict")[, 'fit']
            plot(x2, y2)
            plot(log10(df$x), log10(df$y), type="l")
            points(x2, y2, type="l", col="red")
           plot(df$x, df$y, type="l")
            points(10^x2, 10^y2, type="l", col="red")
        }

        return(fit)
    })
    return(fits)
}

#' Run BADGER expression-only model to assess posterior probability of CNVs given normalized expression data only
#'
#' @param gexp Normalized gene expression matrix
#' @param fits Fit for variance around mean
#' @param m Expected mean deviation due to copy number change
#' @param region Region of interest such as expected CNV boundaries
#' @param gos Gene position table
#' @param quiet Boolean for whether to suppress progress display
#' @return List of posterior probabilities for amplification and deletion
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
#' set.seed(0)
#' library(biomaRt) ## for gene coordinates
#' mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
#' gos <- getBM(values=rownames(gexp),attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),filters=c("ensembl_gene_id"),mart=mart.obj)
#' gos$pos <- (gos$start_position + gos$end_position)/2
#' results <- calcGexpCnvProb(gexp, fits, 0.15, region, gos)
#'
#' @export
#'
calcGexpCnvProb <- function(gexp, fits, m, region, gos, quiet=TRUE) {

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
        pm <- rep(NA, ncol(gexp))
        names(pm) <- colnames(gexp)
        return(list(pm, pm, pm))
    }
    gexp <- gexp[vi,]
    # smooth
    #mat <- apply(gexp, 2, runmean, k=window.size)
    mu0 <- apply(gexp, 2, mean)
    ng <- nrow(gexp)
    sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
    #gexp <- rbind(apply(gexp, 2, mean), apply(gexp, 2, median))

    ####
    ## Model
    ####
    print('aggregating data to list...')
    data <- list(
        'K' = length(mu0),
        'JJ' = nrow(gexp),
        'gexp' = gexp,
        'sigma0' = sigma0,
        'mag0' = m
    )
    modelFile <-  system.file("bug", "expressionModel.bug", package = "badger")

    print('Initializing model...')
    # Joe says 4 chains is a standard, so just stick with 4 chains
    #model <- jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    #update(model, 1000, progress.bar=ifelse(quiet,"none","text"))
    inits <- list(
        list(S = rep(0, ncol(gexp)), dd = 0),
        list(S = rep(1, ncol(gexp)), dd = 0),
        list(S = rep(0, ncol(gexp)), dd = 1),
        list(S = rep(1, ncol(gexp)), dd = 1)
    )
    require(rjags)
    model <- jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
    update(model, 100, progress.bar=ifelse(quiet,"none","text"))

    parameters <- c('S', 'dd', 'mu')
    samples <- coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine chains

    print('...Done!')

    snpLike <- samples
    v <- colnames(snpLike)
    S <- snpLike[,grepl('S', v)]
    dd <- snpLike[,grepl('dd', v)]
    mu <- snpLike[,grepl('mu', v)]
    #plot(mu0, colMeans(mu))
    delcall <- apply(S*(1-dd), 2, mean)
    delcall
    ampcall <- apply(S*dd, 2, mean)
    ampcall
    #plot(mu0, delcall)
    #plot(mu0, ampcall)
    names(ampcall) <- names(delcall) <- colnames(gexp)

    return(list('posterior probability of amplification'=ampcall,
                'posterior probability of deletion'=delcall,
                'estimated mean normalized expression deviation'=mu0))
}


