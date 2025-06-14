# Required libraries
library(seqinr)
library(dplyr)
library(ggplot2)
library(Biostrings)
library(pwalign)
library(foreach)
library(future)
library(doFuture)
library(progressr)

# Setup parallel backend
plan(multicore)  # More stable than multicore for Biostrings
registerDoFuture()
handlers(global = TRUE)

# File paths
ref_file <- "wuhan.fasta"
variant_files <- c(
  "alpha.fasta",
  "beta.fasta",
  "gamma.fasta",
  "delta.fasta",
  "omicronba1.fasta",
  "omicronba2.fasta",
  "omicronjn1.fasta"
)

# ✅ Clean only problematic characters (not everything!)

clean_seq <- function(seq) {
  gsub("[^ACTGUactgu]", "A", seq)
}


# ✅ Extract surface glycoprotein sequences
extract_surface_glycoproteins <- function(filepath) {
  fasta_seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE, whole.header = TRUE)
  headers <- names(fasta_seqs)
  is_surface <- grepl("\\|surface glycoprotein", headers)
  if (sum(is_surface) == 0) stop("No surface glycoprotein genes found in the file")
  surface_seqs <- fasta_seqs[is_surface]
  sequences <- unlist(surface_seqs)
  names(sequences) <- headers[is_surface]
  return(sequences)
}

# ✅ Needleman-Wunsch with cleaning and DNAString
needleman_wunsch <- function(ref_seq, var_seq) {
  ref_seq <- DNAString(clean_seq(ref_seq))
  var_seq <- DNAString(clean_seq(var_seq))

  dna_mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)

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

# ✅ Compare one sequence vs ref
compare_sequences <- function(ref_seq, var_seq, variant_name) {
  aligned <- needleman_wunsch(ref_seq, var_seq)
  aln_ref <- aligned[1]
  aln_var <- aligned[2]

  len <- min(nchar(aln_ref), nchar(aln_var))
  len <- len - (len %% 3)

  result <- list()
  idx <- 0

  for (i in seq(1, len, by = 3)) {
    ref_codon <- substr(aln_ref, i, i + 2)
    var_codon <- substr(aln_var, i, i + 2)

    if (grepl("-", ref_codon) | grepl("-", var_codon) | ref_codon == var_codon) next

    ref_aa <- seqinr::translate(s2c(ref_codon))
    var_aa <- seqinr::translate(s2c(var_codon))
    mut_type <- ifelse(ref_aa == var_aa, "synonymous", "non-synonymous")

    idx <- idx + 1
    result[[idx]] <- data.frame(
      codon_position = i / 3,
      ref_codon = ref_codon,
      var_codon = var_codon,
      ref_aa = ref_aa,
      var_aa = var_aa,
      aa_change = paste0(ref_aa, i / 3, var_aa),
      type = mut_type,
      variant = variant_name,
      stringsAsFactors = FALSE
    )
  }

  if (idx == 0) return(NULL)
  bind_rows(result)
}

# ✅ Load reference
ref_seq <- extract_surface_glycoproteins(ref_file)[[1]]

# ✅ Run in parallel with progress bar
mutations_list <- list()

with_progress({
  p <- progressor(steps = length(variant_files))

  for (file in variant_files) {
    variant_name <- gsub(".fasta$", "", basename(file))
    variant_data <- extract_surface_glycoproteins(file)

    res <- foreach(i = seq_along(variant_data), .combine = rbind) %dofuture% {
      var_seq <- variant_data[[i]]
      compare_sequences(ref_seq, var_seq, variant_name)
    }

    mutations_list[[variant_name]] <- res
    p(sprintf("Done: %s", variant_name))
  }
})

# ✅ Combine all mutation data
mutations <- bind_rows(mutations_list)

# ✅ Plot top mutations
top_mutations <- mutations %>%
  group_by(variant, aa_change) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(variant) %>%
  slice_max(order_by = count, n = 10) %>%
  ungroup()

ggplot(top_mutations, aes(x = reorder(aa_change, count), y = count, fill = variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap(~variant, scales = "free_y") +
  labs(
    title = "Top 10 mutaciones más frecuentes por variante",
    x = "Cambio aminoacídico",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
