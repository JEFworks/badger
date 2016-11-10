# HMM boundary detection methods

#' Predict deletion boundaries based on allele patterns across many single cells
#'
#' @param r Alt allele count in single cells
#' @param n.sc Coverage in single cells
#' @param l Alt allele count in bulk
#' @param n.bulk Coverage in bulk
#' @param k Smoothing window. Must be smaller than number of snps ie. nrow(r). Default: 31
#' @param heights List ofraversal heights for the constructed hierarchical dendrogram to identify more pure subpopulations with 1 being the most shallow cut (all cells in one group). Default = 1:10
#' @param plot Boolean of whether to plot. Default: TRUE
#'
#'
#' @examples
#' ## Single cell data
#' data(snpsHet_MM16ScSample)
#' ## Bulk exome
#' data(snpsHet_MM16BulkSample)
#' predBounds(r, cov.sc, l, cov.bulk, k=301)
#'
#' @export
#'
predBounds <- function(r, n.sc, l, n.bulk, k=31, heights=1:10, plot=TRUE) {

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
    ## Minor allele fraction
    mat.tot <- visualize.mat(r, n.sc, E = l/n.bulk)

    ## turn nan to NA
    x <- mat.tot
    x[is.nan(x)] <- NA

    require(caTools)
    mat.smooth <- apply(x, 2, caTools::runmean, k=k)
    maf <- apply(mat.smooth, 1, mean, na.rm=TRUE)

    d <- dist(t(mat.smooth))
    d[is.na(d)] <- 0
    hc <- hclust(d, method='ward.D2')

    ## cut tree at various heights to establish groups
    boundsnps.pred <- lapply(heights, function(h) {

        ct <- cutree(hc, k = h)

        cuts <- unique(ct)

        ## look at each group, if deletion present
        boundsnps.pred <- lapply(cuts, function(group) {
        if(sum(ct==group) > 1) {

            ## minor allele fraction
            mat.tot <- visualize.mat(r[, ct==group], n.sc[, ct==group], E = rowSums(r[, ct==group]>0)/n.bulk)
            ## turn nan to NA
            x <- mat.tot
            x[is.nan(x)] <- NA
            mat.smooth <- apply(x, 2, caTools::runmean, k=k)

            # coverage
            cov.tot <- log10(n.sc+1)
            cov.smooth <- apply(cov.tot, 2, caTools::runmean, k=k)
            cov <- apply(cov.smooth, 1, sum, na.rm=TRUE)
            cov <- na.omit(cov)
            head(cov)

            df <- data.frame('maf'=maf, 'cov'=cov)

            require(depmixS4)
            mod <- depmix(
                response = list(maf ~ 1, cov ~ 1),
                data = df,
                nstates = 2,
                family = list(gaussian(), gaussian())
            )
            mod <- setpars(mod, c(0, 1,
                                  0.90, 0.01, 0.01, 0.99,
                                  0, 0.1,
                                  2.5, 0.75, #stricter coverage requirements for deletion
                                  0.5, 0.25,
                                  4, 3))
            #f <- fit(mod, fixed = c(1, 1, rep(0, 4*3)))
            #f <- fit(mod)
            #summary(f)
            #esttrans <- posterior(f)
            esttrans <- viterbi(mod)

            ## Get boundaries from states
            #plot(1:nrow(r), esttrans$state, type="l")
            boundsnps <- rownames(mat.tot)[esttrans$state == 1]
        }})
    })

    boundsnps <- unique(unlist(boundsnps.pred))

    pred <- rep(1, nrow(mat.tot))
    names(pred) <- rownames(mat.tot)
    if(length(boundsnps>0)) {
        pred[boundsnps] <- 0
    }

    if(plot) {
        par(mfrow=c(3,1), mar=rep(5,4))
        plot(1:nrow(mat.tot), maf, type="l")
        plot(1:nrow(mat.tot), cov, type="l")
        plot(1:nrow(mat.tot), pred, type="l")
    }

    return(boundsnps)
}
