library(seqinr)
library(dplyr)
library(ggplot2)
library(foreach)
library(future)
library(doFuture)
library(progressr)
library(Biostrings)
library(pwalign)

# Enable progressr
handlers(global = TRUE)
plan(multicore)  # Use multicore or multisession depending on OS
registerDoFuture()

# File paths
ref_file <- "wuhan.fasta"
variant_files <- c(
                   "alpha.fasta", "beta.fasta", "gamma.fasta", 
                   "delta.fasta", "omicronba1.fasta", 
                   "omicronba2.fasta", "omicronjn1.fasta"
)

needleman_wunsch <- function(ref_seq, var_seq) {
    dna_mat <- nucleotideSubstitutionMatrix(
        match = 1,    # Score for matches
        mismatch = -1, # Penalty for mismatches
        baseOnly = TRUE
    )

    gap_open <- -0.5      # Cost to open a gap
    gap_extend <- -0.5   # Cost to extend a gap

    # 2. Biological alignment with proper parameters
    alignment <- pairwiseAlignment(
       pattern = var_seq,
       subject = ref_seq,
       substitutionMatrix = dna_mat,
       gapOpening = gap_open,
       gapExtension = gap_extend,
       type = "global"  # Needleman-Wunsch
    )

    aln_var <- as.character(pattern(alignment))
    aln_ref <- as.character(subject(alignment))

    return(list(aligned_seqs = matrix(c(aln_ref, aln_var), byrow = TRUE, nrow = 2)))
}

# Initialize results dataframe
mutations <- data.frame(
                        codon_position = integer(),
                        ref_codon = character(),
                        var_codon = character(),
                        ref_aa = character(),
                        var_aa = character(),
                        aa_change = character(),
                        type = character(),
                        variant = character(),
                        stringsAsFactors = FALSE
)

compare_sequences <- function(ref_seq, var_seq, variant_name) {
    aln <- needleman_wunsch(ref_seq, var_seq)
    
    # Work directly with character vectors - no string conversion needed
    ref_aln <- aln$aligned_seqs[1,]
    var_aln <- aln$aligned_seqs[2,]
    
    # Calculate number of complete codons
    min_len <- min(length(ref_aln), length(var_aln))
    n_codons <- min_len %/% 3
    
    if (n_codons == 0) return(NULL)
    
    # Pre-allocate vectors for results
    mutations <- vector("list", n_codons)
    n_mutations <- 0
    
    # Process codons in chunks of 3
    for (codon_pos in seq_len(n_codons)) {
        start_idx <- (codon_pos - 1) * 3 + 1
        end_idx <- start_idx + 2
        
        # Extract codons directly from vectors
        ref_codon <- ref_aln[start_idx:end_idx]
        var_codon <- var_aln[start_idx:end_idx]
        
        # Skip if gaps present or incomplete codons
        if (any(ref_codon == "-") || any(var_codon == "-")) next
        if (length(ref_codon) < 3 || length(var_codon) < 3) next
        
        # Skip if codons are identical
        if (identical(ref_codon, var_codon)) next
        
        # Convert to amino acids
        ref_aa <- translate(ref_codon)  # Assuming translate() can handle character vectors
        var_aa <- translate(var_codon)
        
        # Store mutation data
        n_mutations <- n_mutations + 1
        mutations[[n_mutations]] <- list(
            codon_position = codon_pos,
            ref_codon = paste(ref_codon, collapse = ""),
            var_codon = paste(var_codon, collapse = ""),
            ref_aa = ref_aa,
            var_aa = var_aa,
            aa_change = paste0(ref_aa, codon_pos, var_aa),
            type = ifelse(ref_aa == var_aa, "synonymous", "non-synonymous"),
            variant = variant_name
        )
    }
    
    # Return early if no mutations found
    if (n_mutations == 0) return(NULL)
    
    # Create single data frame efficiently
    mutations <- mutations[seq_len(n_mutations)]  # Remove unused slots
    
    # Extract vectors for data frame construction
    data.frame(
        codon_position = sapply(mutations, `[[`, "codon_position"),
        ref_codon = sapply(mutations, `[[`, "ref_codon"),
        var_codon = sapply(mutations, `[[`, "var_codon"),
        ref_aa = sapply(mutations, `[[`, "ref_aa"),
        var_aa = sapply(mutations, `[[`, "var_aa"),
        aa_change = sapply(mutations, `[[`, "aa_change"),
        type = sapply(mutations, `[[`, "type"),
        variant = sapply(mutations, `[[`, "variant"),
        stringsAsFactors = FALSE
    )
}

compare_sequences <- function(ref_seq, var_seq, variant_name) {
    aln <- needleman_wunsch(ref_seq, var_seq)
    print("Aligned two sequences successfully") 
    print("Reference sequence:")
    print(ref_seq)
    print("Variant sequence:")
    print(var_seq)

    # Work directly with character vectors
    ref_aln <- aln$aligned_seqs[1,]
    var_aln <- aln$aligned_seqs[2,]
    
    min_len <- min(length(ref_aln), length(var_aln))
    n_codons <- min_len %/% 3
    
    if (n_codons == 0) return(NULL)
    
    # Pre-allocate vectors directly - no intermediate lists!
    codon_positions <- integer(n_codons)
    ref_codons <- character(n_codons) 
    var_codons <- character(n_codons)
    ref_aas <- character(n_codons)
    var_aas <- character(n_codons)
    aa_changes <- character(n_codons)
    types <- character(n_codons)
    
    n_mutations <- 0
    
    for (codon_pos in seq_len(n_codons)) {
        start_idx <- (codon_pos - 1) * 3 + 1
        end_idx <- start_idx + 2
        
        ref_codon <- ref_aln[start_idx:end_idx]
        var_codon <- var_aln[start_idx:end_idx]
        
        # Skip invalid or identical codons
        if (any(ref_codon == "-") || any(var_codon == "-") || 
            identical(ref_codon, var_codon)) next
        
        # Translate codons
        ref_aa <- translate(ref_codon)
        var_aa <- translate(var_codon)
        
        # Store directly in pre-allocated vectors
        n_mutations <- n_mutations + 1
        codon_positions[n_mutations] <- codon_pos
        ref_codons[n_mutations] <- paste(ref_codon, collapse = "")
        var_codons[n_mutations] <- paste(var_codon, collapse = "")
        ref_aas[n_mutations] <- ref_aa
        var_aas[n_mutations] <- var_aa
        aa_changes[n_mutations] <- paste0(ref_aa, codon_pos, var_aa)
        types[n_mutations] <- ifelse(ref_aa == var_aa, "synonymous", "non-synonymous")
    }
    
    if (n_mutations == 0) return(NULL)
    
    # Single data frame creation with pre-sized vectors
    data.frame(
        codon_position = codon_positions[seq_len(n_mutations)],
        ref_codon = ref_codons[seq_len(n_mutations)],
        var_codon = var_codons[seq_len(n_mutations)],
        ref_aa = ref_aas[seq_len(n_mutations)],
        var_aa = var_aas[seq_len(n_mutations)],
        aa_change = aa_changes[seq_len(n_mutations)],
        type = types[seq_len(n_mutations)],
        variant = variant_name,
        stringsAsFactors = FALSE
    )
}

# Alternative: Skip data frame creation entirely if you just need the data
compare_sequences_list <- function(ref_seq, var_seq, variant_name) {
    aln <- needleman_wunsch(ref_seq, var_seq)
    print("needleman wunsched once")
    ref_aln <- aln$aligned_seqs[1,]
    var_aln <- aln$aligned_seqs[2,]
    
    min_len <- min(length(ref_aln), length(var_aln))
    n_codons <- min_len %/% 3
    
    if (n_codons == 0) return(list())
    
    mutations <- vector("list", n_codons)
    n_mutations <- 0
    
    for (codon_pos in seq_len(n_codons)) {
        start_idx <- (codon_pos - 1) * 3 + 1
        ref_codon <- ref_aln[start_idx:(start_idx + 2)]
        var_codon <- var_aln[start_idx:(start_idx + 2)]
        
        if (any(ref_codon == "-") || any(var_codon == "-") || 
            identical(ref_codon, var_codon)) next
        
        ref_aa <- translate(ref_codon)
        var_aa <- translate(var_codon)
        
        n_mutations <- n_mutations + 1
        # Store as named vector instead of data frame
        mutations[[n_mutations]] <- c(
            codon_position = codon_pos,
            ref_codon = paste(ref_codon, collapse = ""),
            var_codon = paste(var_codon, collapse = ""),
            ref_aa = ref_aa,
            var_aa = var_aa,
            aa_change = paste0(ref_aa, codon_pos, var_aa),
            type = ifelse(ref_aa == var_aa, "synonymous", "non-synonymous"),
            variant = variant_name
        )
    }
    
    mutations[seq_len(n_mutations)]
}


# Read S gene only from each file (every 12th gene starting from 3)
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



# Main logic
ref_data <- extract_surface_glycoproteins(ref_file)
ref_seq <- ref_data[[1]]

mutations <- list()

with_progress({
    p <- progressor(steps = length(variant_files))

    for (file in variant_files) {
        variant_data <- extract_surface_glycoproteins(file)
        var_name <- tools::file_path_sans_ext(basename(file))

        message(sprintf("Processing S gene variant: %s (%d sequences)", var_name, length(variant_data)))

        res <- foreach(j = seq_along(variant_data), .combine = rbind) %dofuture% {
            compare_sequences(ref_seq, variant_data[[j]], var_name)
            # p(sprintf("Processed %s", variant_data$names[j]))
        }

        mutations[[var_name]] <- res
        p(sprintf("Completed %s", var_name))
    }
})

# Combine and plot
mutations_df <- bind_rows(mutations)

top_mutations <- mutations_df %>%
    group_by(variant, aa_change) %>%
    tally(sort = TRUE) %>%
    group_by(variant) %>%
    top_n(10, n)

ggplot(top_mutations, aes(x = reorder(aa_change, n), y = n, fill = variant)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~variant, scales = "free") +
    coord_flip() +
    labs(title = "Top 10 S Gene Amino Acid Changes per Variant", x = "AA Change", y = "Count")
