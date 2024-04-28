
## ----------- ancillary.R ------------ ##
#                                        #
#      msa                               #
#      mltree                            #
#      gapless_msa                       #
#                                        #
## ------------------------------------ ##


## ---------------------------------------------------------------- ##
#     msa <- function(sequences, ids, seqtype, sfile, inhouse )      #
## ---------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein, DNA or CDS sequences.
#' @usage msa(sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE, inhouse = FALSE)
#' @param sequences vector containing the sequences as strings.
#' @param seqtype it should be either "prot" of "dna" or "cds" (see details).
#' @param ids character vector containing the sequences' ids.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @param inhouse logical, if TRUE the in-house MUSCLE software is used. It must be installed on your system and in the search path for executables.
#' @details If seqtype is set to "cds" the sequences must not contain stop codons and they will be translated using the standard code. Afterward, the amino acid alignment will be used to lead the codon alignment.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($id) returns the identifiers, the third element ($aln) provides the alignment in fasta format and the fourth element ($ali) gives the alignment in matrix format.
#' @examples \dontrun{msa(sequences = sapply(c("P19446", "P40925", "P40926"), ptm::get.seq),ids = c("wmelon", "cyt", "mit"))}
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @importFrom bio3d write.fasta
#' @export

msa <- function (sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE, inhouse = FALSE){

  if (length(sequences) < 2) {
    stop("At least two sequences are required!")
  } else if (length(sequences) != length(ids)) {
    stop("The number of sequences and sequences' ids doesn't match!")
  }
  if (inhouse) {
    if (seqtype == "cds"){
      dnaSeq <- sequences
      cod <- strsplit(gsub("(.{3})", "\\1 ", dnaSeq), split = " ")
      sequences <- unlist(lapply(sequences, function(x) tr(x)))
    }
    seqs <- lapply(sequences, function(x) strsplit(x, split = "")[[1]])
    sqs <- bio3d::seqbind(seqs[[1]], seqs[[2]], blank = "-")
    c <- 2
    while (c < length(sequences)) {
      c <- c + 1
      sqs <- bio3d::seqbind(sqs, seqs[[c]], blank = "-")
    }
    aln <- bio3d::seqaln(sqs, id = ids, exefile = "muscle")
    aln$seq <- sequences
    if (seqtype == "cds"){
      aln$cod <- aln$ali
      for (i in 1:length(cod)){
        contador <- 1
        for (j in 1:ncol(aln$ali)){
          if (aln$ali[i,j] != "-"){
            aln$cod[i,j] <- cod[[i]][contador]
            contador <- contador + 1
          } else {
            aln$cod[i,j] <- "---"
          }
        }
      }
    }
    if (sfile != FALSE & seqtype == "cds") {
      for (i in 1:length(ids)){
        t <- paste(">", ids[i], sep = "")
        cat(t, "\n", file = sfile, append = TRUE)
        tt <- paste(aln$cod[i,], collapse = "")
        cat(tt, "\n", file = sfile, append = TRUE)
      }
    } else if (sfile != FALSE & seqtype != "cds") {
      bio3d::write.fasta(aln, file = sfile)
    }
    if (file.exists("aln.fa")) {
      system("rm aln.fa")
    }
    return(aln)
  }  else {
    ## --- Using the Biostring and muscle R packages
    # seq <- Biostrings::AAStringSet(sequences)

    if (requireNamespace('Biostrings', quietly = TRUE)){
      if (seqtype == "prot"){
        seq <- Biostrings::AAStringSet(sequences)
      } else if (seqtype == "dna"){
        seq <- Biostrings::DNAStringSet(sequences)
      }
    } else {
      stop("You must install the package Biostrings in order to use this function")
    }

    if (requireNamespace('muscle', quietly = TRUE)){
      aln1 <- muscle::muscle(seq)
    } else {
      stop("You must install the package muscle in order to use this function")
    }
    aln <- list()
    aln$seq <- sequences
    aln$ids <- ids
    aln$aln <- as.character(aln1)
    l <- sapply(aln$aln, function(x) strsplit(x, split = ""))
    aln$ali <- matrix(unlist(l), nrow = length(sequences),
                      byrow = TRUE)
    if (sfile != FALSE) {
      for (i in 1:length(aln$aln)) {
        t <- paste(">", aln$ids[i], sep = "")
        cat(t, file = sfile, append = TRUE)
        if (seqtype == "cds"){
          tt <- paste("\n", aln$cod[i], "\n", sep = "")
        } else {
          tt <- paste("\n", aln$aln[i], "\n", sep = "")
        }
        cat(tt, file = sfile, append = TRUE)
      }
    }
    return(aln)
  }
}

## ---------------------------------------------------------------- ##
#                     mltree <- function()                           #
## ---------------------------------------------------------------- ##
#' Build Up a ML Tree
#' @description Given an alignment builds an ML tree.
#' @usage mltree(msa, df = TRUE, gapl = TRUE, model = "WAG")
#' @param msa input alignment.
#' @param df logical. When TRUE msa should be a dataframe, when FALSE msa should be a string giving the path to a fasta file containing the alignment.
#' @param gal logical, when TRUE a gapless alignment is used.
#' @param model allows to choose an amino acid models (see the function phangorn::as.pml)
#' @details The function makes a NJ tree and then improvove it using an optimization procedure based on ML.
#' @return a ML optimize tree (and parameters)
#' @author Juan Carlos Aledo
#' @seealso gapless_msa
#' @importFrom ape nj
#' @importFrom ape dist.aa
#' @importFrom bio3d read.fasta
#' @importFrom phangorn as.phyDat
#' @importFrom phangorn pml
#' @importFrom phangorn optim.pml
#' @export

mltree <- function(msa, df = TRUE, gapl = TRUE, model = "WAG"){
  if (df == TRUE){
    aln <- msa
  } else {
    aln <- bio3d::read.fasta(msa)$ali
  }
  if (gapl == TRUE){
    data <- gapless_msa(aln)
  } else {
    data <- aln
  }

  tre.ini <- nj(ape::dist.aa(data))
  fit.ini <- pml(tre.ini, as.phyDat(as.matrix(data), type = "AA"), model = model)
  fit <- optim.pml(fit.ini, model = model)
  return(fit)
}

## ---------------------------------------------------------------- ##
#                   gapless_msa <- function()                        #
## ---------------------------------------------------------------- ##
#' Remove Gaps in a MSA
#' @description Removes gaps in a given msa.
#' @usage gapless_msa(msa, seqtype = 'AA', df = TRUE, sfile = FALSE)
#' @param msa input alignment.
#' @param seqtype the nature of the sequences: 'DNA' or 'AA'.
#' @param df logical. When TRUE msa should be a matrix, when FALSE msa should be a string giving the path to a fasta file containing the alignment.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @details It should be noted that this function does not carry out the alignment itself.
#' @return an alignment without gaps in form of matrix or a file containing such an alignment in fasta format.
#' @examples \dontrun{gapless(msa = bla)}
#' @seealso msa
#' @importFrom seqinr read.fasta
#' @export

gapless_msa <- function(msa, seqtype = "AA", df = TRUE, sfile = FALSE){
  if (df == FALSE){
    msa <- read.fasta(msa, seqtype = seqtype)
    M <- matrix(rep(NA, length(msa)*length(msa[[1]])), nrow = length(msa))
    rownames(M) <- attributes(msa)$names
    for (i in 1:length(msa)){
      M[i,] <- as.character(msa[[i]])
    }
  } else {
    M <- as.matrix(msa)
  }
  gapless <- c()
  for (i in 1:ncol(M)){
    if (! "-" %in% M[,i]){
      gapless <- c(gapless, i)
    }
  }
  output <- M[, gapless]

  if (sfile != FALSE) {
    for (i in 1:nrow(output)) {
      t <- paste(">", rownames(output)[i], sep = "")
      cat(t, file = sfile, append = TRUE)
      tt <- paste(output[i,], collapse = "")
      tt <- paste("\n", tt, "\n", sep = "")
      cat(tt, file = sfile, append = TRUE)
    }
  }

  return(as.data.frame(output))
}

