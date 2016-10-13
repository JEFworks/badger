#' Functionalities related to expression model
#'
#' @author: Jean Fan
#'

#load('data/MM16.RData')
#mat <- log2(MM16.counts + 1)
#load('data/Normal.RData')
#mat.ref <- log2(Normal.counts + 1)
#genes.int <- intersect(rownames(mat), rownames(mat.ref))
#mat <- mat[genes.int,]
#mat.ref <- mat.ref[genes.int,]


normalizedExpression <- function(mat, mat.ref) {

    # filter
    vi <- (rowMeans(mat) > 4.5 & rowMeans(mat.ref) > 4.5) | rowMeans(mat) > 6 | rowMeans(mat.ref) > 8
    if(verbose) {
        print("Number of genes passed filtering:")
        print(table(vi))
    }
    mat <- mat[vi,]
    mat.ref <- mat.ref[vi,]

    # normalize
    # library size
    #mat <- t(t(mat) - colMeans(mat))
    mat <- scale(mat)
    #mat.ref <- t(t(mat.ref) - colMeans(mat.ref))
    mat.ref <- scale(mat.ref)
    #hist(mat)
    #hist(mat.ref)
    #hist(colMeans(mat))
    #hist(colMeans(mat.ref))
    # normalize to normal
    #mat.means <- rowMeans(mat)
    #mat <- mat - mat.means
    #mat.ref <- mat.ref - mat.means
    ref.means <- rowMeans(mat.ref)
    mat <- mat - ref.means
    mat.ref <- mat.ref - ref.means
    #hist(colMeans(mat))
    #hist(colMeans(mat.ref))

    return(list(mat, mat.ref))
}

#mat.tot <- cbind(mat, mat.ref)
#library(biomaRt) ## for gene coordinates
#mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
#gos <- getBM(values=rownames(mat.tot),attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),filters=c("ensembl_gene_id"),mart=mart.obj)
#gos$pos <- (gos$start_position + gos$end_position)/2
#head(gos)

plotExpHeatmap <- function(mat.tot, gos, zlim=c(-2,2), window.size = 101, orderCells=FALSE) {

    ## organize into chromosomes
    tl <- tapply(1:nrow(gos),as.factor(gos$chromosome_name),function(ii) {
        na.omit(mat.tot[gos[,1][ii[order(gos$pos[ii],decreasing=F)]],])
    })
    ## only care about these chromosomes
    tl <- tl[c(as.factor(1:22), "X")]
    tl <- tl[as.factor(1:22)]

    if(orderCells) {
        avgd <- do.call(rbind, lapply(names(tl),function(nam) {
            d <- tl[[nam]]
            d <- colMeans(d)
            d
        }))
        hc <- hclust(dist(t(avgd)))
        order <- hc$order
    }

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
        #d <- apply(d,2,rollmean,k=nrow(d)/3, align="center")
        d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
        par(mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
        image(1:nrow(d),1:ncol(d),d,col=pcol,zlim=zlim,xlab="",ylab="",axes=F,main=nam); box()
    })

    return(tl)
}
