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
                        min_base_quality=20,
                        min_mapq=10
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

        if(verbose) {
            print("Done!")
        }

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
                        min_base_quality=20,
                        min_mapq=10
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

            if(verbose) {
                print("Done!")
            }
            return(totCount)
    }
