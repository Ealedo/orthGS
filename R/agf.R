#' Angiosperms Gymnosperms Ferns
#'
#' @format
#' A dataframe with 275 rows (GS proteins) and 23 columns:
#' \describe{
#'		\item{n}{Reference number}
#'		\item{phylo_id}{Unique identification label of the protein/gen}
#'		\item{species}{Species}
#'		\item{taxon}{Acrogymnospermae, Angiospermae  or Polypodiopsida}
#'		\item{dna}{CDS sequence}
#'		\item{prot}{Protein sequence}
#'		\item{short}{Unique three letter identification of the species}
#'		\item{gs}{Either GS2, GS1a, GS1b_Gym or GS1b_Ang. Here the ferns proteins have been forced to be either GS1a or GS2}
#'		\item{pI}{Computed isoelectric point}
#'		\item{factor}{Either GS2, GS1a, GS1b_Gym, GS1b_Ang, Ferns}.
#'		\item{size}{Protein size as number of residues}
#'		\item{CSpos}{Position of the signal peptide when inferred}
#'		\item{prediction}{Either mitochondrial transfer peptide, signal peptide, thylakoid luminal transfer peptide, other}
#'		\item{Lk_SP}{Likelihood for signal peptide}
#'		\item{Lk_mTP}{Likelihood for mitochondrial transfer peptide}
#'		\item{Lk_cTP}{Likelihood for chloroplatic transfer peptide}
#'		\item{Lk_Thylak}{Likelihood for thylakoid luminal transfer peptide}
#'		\item{secAa}{Second aminoacid following the initiation methionine}
#'		\item{core}{Protein sequence or GS2 forms after removing the CT and NT not present in other GS lineage}
#'		\item{database}{source of the sequence}
#'		\item{acc}{access identifier}
#'		\item{up_id}{uniprot identifier}
#'		\item{note}{notes}
#' }
#' @source
#' It has been handly curated by the authors
"agf"
