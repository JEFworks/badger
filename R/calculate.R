# Calculate various errors and bias metrics


# Helper function Peter likes
sn <- function(x) { names(x) <- x; x }

# Helper to convert from list of "chr:start:end" to data frame
snp2df <- function(pos, delim=':') {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- as.character(posp[,1])
    posp[,2:ncol(posp)] <- as.numeric(as.character(posp[,2:ncol(posp)]))
    return(posp)
}

#' Calculate effective error rate
#'
#' Calculate effective error rate from sequencing, alignment, and other potential error sources
#'
#' @param type Type of SNPs provided in \code{snpsCov}, \code{snpsRef}, and \code{snpsAlt} used to estimate effective error. This should be (an unambiguous abbreviation of) one of 'homo_ref', 'homo_alt', or 'het'.
#' @param snpsCov A matrix of snp positions as rows and cells as columns with each value corresponding to the total read count at the corresponding snp position in the corresponding cell
#' @param snpsRef A matrix of snp positions as rows and cells as columns with each value corresponding to the reference read count at the corresponding snp position in the corresponding cell. Should have the same snp positions and cells as in \code{snpsCov}
#' @param snpsAlt A matrix of snp positions as rows and cells as columns with each value corresponding to the alternative read count at the corresponding snp position in the corresponding cell. Should have the same snp positions and cells as in \code{snpsCov}
#' @param cnvs A vector of copy number variation regions or other regions with entry as a string in format 'chromosome:start:end'. Optional. Default = NULL
#' @param n Minimum fraction of cell in which we must have coverage at a snp position. Useful for eliminating false positives. Default = 0.1
#'
#' @return Effective error rate
#'
#' @examples
#' \dontrun{
#' deletion.region <- data.frame('chr'=2, start=0, end=1e9)
#' calcErrorRate('het', snpsCov, snpsRef, snpsAlt, cnvs=deletion.region) # ignore regions imapcted by deletion
#' calcErrorRate('het', snpsCov, snpsRef, snpsAlt) # use all regions
#' }
calcErrorRate <- function(type, snpsCov, snpsRef, snpsAlt, cnvs=NULL, n=0.1) {

    snpsCovName <- snp2df(rownames(snpsCov))

    calcErrorRateHomoRef <- function(snpsCov, snpsRef, cnvs, n) {
        # Filter for SNPs that are not inside CNV regions and have coverage in more than n*100% of cells
        if (!is.null(cnvs)) {
            vi <- !(insideCnvs(snpsCovName, cnvs)) & rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        } else {
            vi <- rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        }
        # Exclude positions that never showed any ref coverage (mistake? somatic?)
        vi2 <- vi & rowSums(snpsRef) > 0
        # Calculate error as frequency of observing non-reference allele
        er <- sum(rowSums(snpsCov[vi2,]) - rowSums(snpsRef[vi2,]))/ sum(rowSums(snpsCov[vi2,]))
        return(er)
    }

    calcErrorRateHomoAlt <- function(snpsCov, snpsAlt, cnvs, n) {
        # Filter for SNPs that are not inside CNV regions and have coverage in more than n*100% of cells
        if (!is.null(cnvs)) {
            vi <- !(insideCnvs(snpsCovName, cnvs)) & rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        } else {
            vi <- rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        }
        # Exclude positions that never showed any alt coverage (mistake? somatic?)
        vi2 <- vi & rowSums(snpsAlt) > 0
        # Calculate error as frequency of observing non-alternative allele
        er <- sum(rowSums(snpsCov[vi2,]) - rowSums(snpsAlt[vi2,]))/ sum(rowSums(snpsCov[vi2,]))
        return(er)
    }

    calcErrorRateHet <- function(snpsCov, snpsRef, snpsAlt, cnvs, n) {
        # Filter for SNPs that are not inside CNV regions and have coverage in more than n*100% of cells
        if (!is.null(cnvs)) {
            vi <- !(insideCnvs(snpsCovName, cnvs)) & rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        } else {
            vi <- rowSums(snpsCov > 0) > round(ncol(snpsCov)*n)
        }
        # Calculate error as frequency of observing neither known alleles
        er <- sum(rowSums(snpsCov[vi,]) - rowSums(snpsAlt[vi,]) - rowSums(snpsRef[vi,]))/ sum(rowSums(snpsCov[vi,]))
        er*2
    }

    switch(type,
           homo_ref = calcErrorRateHomoRef(snpsCov, snpsRef, cnvs, n),
           homo_alt = calcErrorRateHomoAlt(snpsCov, snpsAlt, cnvs, n),
           het = calcErrorRateHet(snpsCov, snpsRef, snpsAlt, cnvs, n)
    )
}

