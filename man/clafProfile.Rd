% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/snpModel.R
\name{clafProfile}
\alias{clafProfile}
\title{Composite Minor Allele Frequency (CLAF) Profile plot}
\usage{
clafProfile(r, n.sc, l, n.bulk, region = NULL, filter = TRUE, delim = ":",
  gtf = NULL, plotGene = FALSE)
}
\arguments{
\item{r}{Alt allele count in single cells}

\item{n.sc}{Coverage in single cells}

\item{l}{Alt allele count in bulk}

\item{n.bulk}{Coverage in bulk}

\item{region}{Restrict plotting to select region. Optional. Default: NULL}

\item{filter}{Boolean of whether to filter for SNPs with coverage. Default: TRUE}

\item{delim}{Delimiter for names of SNPs as Chromosome[delim]Position. Default: ":" ex. chr1:283838897}

\item{gtf}{GTF file contents for mapping SNPs to genes. Required if plotGene = TRUE}

\item{plotGene}{Boolean of whether to plot gene track. Default: FALSE}
}
\description{
Composite Minor Allele Frequency (CLAF) Profile plot
}
\examples{
## Single cell data
data(snpsHet_MM16ScSample)
## Bulk exome
data(snpsHet_MM16BulkSample)
# intersect
region <- data.frame('chr'=2, start=0, end=1e9) # deletion region
clafProfile(r, cov.sc, l, cov.bulk, region)
region <- data.frame('chr'=3, start=0, end=1e9) # neutral region
clafProfile(r, cov.sc, l, cov.bulk, region)
}

