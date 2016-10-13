#devtools::use_package("reshape2")
#devtools::use_package("ggplot2")
#devtools::use_package("gridExtra")

#####
# Plotting and visualization functions
#####

#' Composite Minor Allele Frequency (CLAF) Profile plot
#'
#' @param r alt allele count in single cells
#' @param cov.sc coverage in single cells
#' @param l alt allele count in bulk
#' @param cov coverage in bulk
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
clafProfile<- function(r, n.sc, l, n.bulk, filter=TRUE, region=NULL, delim=':', gtf=NULL, plotGene=FALSE) {

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

        snps <- rownames(cov.sc)
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

