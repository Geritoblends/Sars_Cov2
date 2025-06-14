library(seqinr)

path <- "/home/gero/Documents/projects/tec/biologia/project/alpha.fasta"


extract_surface_glycoproteins <- function(filepath) {
  # Read FASTA file
  fasta_seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE, whole.header = TRUE)
  
  # Get all headers
  headers <- names(fasta_seqs)
  
  # Find surface glycoprotein entries (case sensitive)
  is_surface <- grepl("\\|surface glycoprotein", headers)
  
  if (sum(is_surface) == 0) {
    stop("No surface glycoprotein genes found in the file")
  }
  
  # Extract matching sequences
  surface_seqs <- fasta_seqs[is_surface]
  
  # Return as named character vector
  sequences <- unlist(surface_seqs)
  names(sequences) <- headers[is_surface]
  
  return(sequences)
}

