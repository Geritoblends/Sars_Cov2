library(seqinr)
library(dplyr)
library(ggplot2)
library(foreach)
library(future)
library(doFuture)
library(progressr)

# Enable progressr
handlers(global = TRUE)
plan(multicore)
registerDoFuture()

# Your existing needleman_wunsch function
needleman_wunsch<- function(seq1, seq2, gap = -1, mismatch = -1, match = 1) {
  len1 <- length(seq1)
  len2 <- length(seq2)
  
  # Score matrix
  M <- matrix(0, nrow = len1 + 1, ncol = len2 + 1)
  M[1, ] <- gap * 0:len2
  M[, 1] <- gap * 0:len1
  
  # Fill matrix
  for (i in 2:(len1 + 1)) {
    for (j in 2:(len2 + 1)) {
      score_diag <- M[i - 1, j - 1] + if (seq1[i - 1] == seq2[j - 1]) match else mismatch
      score_up   <- M[i - 1, j] + gap
      score_left <- M[i, j - 1] + gap
      M[i, j] <- max(score_diag, score_up, score_left)
    }
  }
  
  # Backtrace
  i <- len1 + 1
  j <- len2 + 1
  align1 <- character(len1 + len2)
  align2 <- character(len1 + len2)
  k <- len1 + len2
  
  while (i > 1 || j > 1) {
    current <- M[i, j]
    if (i > 1 && j > 1 && current == M[i - 1, j - 1] + if (seq1[i - 1] == seq2[j - 1]) match else mismatch) {
      align1[k] <- seq1[i - 1]
      align2[k] <- seq2[j - 1]
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && current == M[i - 1, j] + gap) {
      align1[k] <- seq1[i - 1]
      align2[k] <- "-"
      i <- i - 1
    } else {
      align1[k] <- "-"
      align2[k] <- seq2[j - 1]
      j <- j - 1
    }
    k <- k - 1
  }
  
  # Drop unused prefix of vectors
  alignment <- rbind(align1[(k + 1):(len1 + len2)], align2[(k + 1):(len1 + len2)])
  
  return(list(score = M[len1 + 1, len2 + 1], aligned_seqs = alignment))
}

# FIXED compare_sequences function with proper debugging
compare_sequences <- function(ref_seq_str, var_seq_str, variant_name) {
    print(paste("Processing variant:", variant_name))
    
    # CRITICAL FIX: Convert strings to character vectors
    ref_seq <- unlist(strsplit(ref_seq_str, ""))
    var_seq <- unlist(strsplit(var_seq_str, ""))
    
    print(paste("Ref seq length:", length(ref_seq), "Var seq length:", length(var_seq)))
    
    aln <- needleman_wunsch(ref_seq, var_seq)
    print(paste("Alignment score:", aln$score))
    
    ref_aln <- aln$aligned_seqs[1,]
    var_aln <- aln$aligned_seqs[2,]
    
    min_len <- min(length(ref_aln), length(var_aln))
    n_codons <- min_len %/% 3
    
    print(paste("Alignment length:", min_len, "Complete codons:", n_codons))
    
    if (n_codons == 0) {
        print("WARNING: No complete codons found!")
        return(NULL)
    }
    
    # Pre-allocate vectors
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
        
        # CRITICAL FIX: Convert vectors to strings for translate function
        tryCatch({
            ref_codon_str <- paste(ref_codon, collapse = "")
            var_codon_str <- paste(var_codon, collapse = "")
            
            ref_aa <- translate(s2c(ref_codon_str))  # seqinr translate needs s2c
            var_aa <- translate(s2c(var_codon_str))
            
            # Store directly in pre-allocated vectors
            n_mutations <- n_mutations + 1
            codon_positions[n_mutations] <- codon_pos
            ref_codons[n_mutations] <- ref_codon_str
            var_codons[n_mutations] <- var_codon_str  
            ref_aas[n_mutations] <- ref_aa
            var_aas[n_mutations] <- var_aa
            aa_changes[n_mutations] <- paste0(ref_aa, codon_pos, var_aa)
            types[n_mutations] <- ifelse(ref_aa == var_aa, "synonymous", "non-synonymous")
            
        }, error = function(e) {
            print(paste("Translation error at codon", codon_pos, ":", e$message))
        })
    }
    
    print(paste("Mutations found:", n_mutations))
    
    if (n_mutations == 0) {
        print("WARNING: No mutations detected!")
        return(NULL)
    }
    
    # Return data frame with only the mutations found
    result <- data.frame(
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
    
    print(paste("Returning data frame with", nrow(result), "mutations"))
    return(result)
}

# Your existing read_s_genes function
read_s_genes <- function(filepath) {
    seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE)
    s_gene_indices <- seq(1, length(seqs), by = 1)
    list(sequences = unlist(seqs[s_gene_indices]), names = names(seqs)[s_gene_indices])
}

# Main logic with debugging
variant_files <- c("test.fasta")
ref_file <- "wuhan.fasta"

print("=== READING REFERENCE FILE ===")
ref_data <- read_s_genes(ref_file)
ref_seq <- ref_data$sequences[1]
print(paste("Reference sequence length:", nchar(ref_seq)))

mutations <- list()

print("=== PROCESSING VARIANTS ===")
with_progress({
    p <- progressor(steps = length(variant_files))

    for (file in variant_files) {
        print(paste("Processing file:", file))
        variant_data <- read_s_genes(file)
        var_name <- tools::file_path_sans_ext(basename(file))
        
        print(paste("Variant sequences found:", length(variant_data$sequences)))

        message(sprintf("Processing S gene variant: %s (%d sequences)", var_name, length(variant_data$sequences)))

        res <- foreach(j = seq_along(variant_data$sequences), .combine = rbind) %dofuture% {
            compare_sequences(ref_seq, variant_data$sequences[j], var_name)
        }

        mutations[[var_name]] <- res
        print(paste("Results for", var_name, ":", ifelse(is.null(res), "NULL", paste(nrow(res), "mutations"))))
        p(sprintf("Completed %s", var_name))
    }
})

print("=== COMBINING RESULTS ===")
print("Mutations list contents:")
print(sapply(mutations, function(x) {
  if(is.null(x)) return("NULL")
  if(is.data.frame(x)) return(paste("data.frame:", nrow(x), "rows"))
  return(class(x))
}))

# Combine results
mutations_df <- bind_rows(mutations)

print(paste("Final mutations_df dimensions:", nrow(mutations_df), "x", ncol(mutations_df)))

if(nrow(mutations_df) > 0) {
    print("Column names:")
    print(colnames(mutations_df))
    print("First few mutations:")
    print(head(mutations_df))
    
    # Create plot
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
        
} else {
    print("ERROR: No mutations found!")
    print("Check your FASTA files - they might be identical or have formatting issues")
}
