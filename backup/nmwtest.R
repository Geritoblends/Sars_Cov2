library(seqinr)
library(Biostrings)
library(pwalign)

s1 = read.fasta("test.fasta", forceDNAtolower = FALSE, as.string = TRUE)[[1]]
s2 = read.fasta("wuhan.fasta", forceDNAtolower = FALSE, as.string = TRUE)[[3]]


dna_mat <- nucleotideSubstitutionMatrix(
    match = 1,    # Score for matches
    mismatch = -1, # Penalty for mismatches
    baseOnly = TRUE
)

gap_open <- -0.5      # Cost to open a gap
gap_extend <- -0.5   # Cost to extend a gap

# 2. Biological alignment with proper parameters
alignment <- pairwiseAlignment(
    pattern = s1,
    subject = s2,
    substitutionMatrix = dna_mat,
    gapOpening = gap_open,
    gapExtension = gap_extend,
    type = "global"  # Needleman-Wunsch
)

as.character(pattern(alignment))
as.character(subject(alignment))

# writePairwiseAlignments(alignment)  # Shows alignment with mismatches

