library(Biostrings)
library(pwalign)

#  Función de alineamiento con limpieza previa
needleman_wunsch <- function(ref_seq, var_seq) {
    dna_mat <- nucleotideSubstitutionMatrix(
                                            match = 1,
                                            mismatch = -1
    )
    alignment <- pairwiseAlignment(
                                   pattern = var_seq,
                                   subject = ref_seq,
                                   substitutionMatrix = dna_mat,
                                   gapOpening = -0.5,
                                   gapExtension = -0.5,
                                   type = "global"
    )
    aln_ref <- as.character(subject(alignment))
    aln_var <- as.character(pattern(alignment))
    print("Secuencia de referencia:")
    print(aln_ref)
    print("Secuencia de variante:")
    print(aln_var)
    return(c(aln_ref, aln_var))
}

# refseq = "ABCIFNN"
# varseq = "ABIFNN"
refseq = "actgg"
varseq = "acgg"

needleman_wunsch(refseq, varseq)[[1]]
