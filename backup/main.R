# Set working directory
setwd("/home/gero/Documents/projects/tec/biologia/project")

# Translation table (codon to amino acid)
trad <- c(UUU="F", UUC="F", UUA="L", UUG="L",
          UCU="S", UCC="S", UCA="S", UCG="S",
          UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
          UGU="C", UGC="C", UGA="STOP", UGG="W",
          CUU="L", CUC="L", CUA="L", CUG="L",
          CCU="P", CCC="P", CCA="P", CCG="P",
          CAU="H", CAC="H", CAA="Q", CAG="Q",
          CGU="R", CGC="R", CGA="R", CGG="R",
          AUU="I", AUC="I", AUA="I", AUG="M",
          ACU="T", ACC="T", ACA="T", ACG="T",
          AAU="N", AAC="N", AAA="K", AAG="K",
          AGU="S", AGC="S", AGA="R", AGG="R",
          GUU="V", GUC="V", GUA="V", GUG="V",
          GCU="A", GCC="A", GCA="A", GCG="A",
          GAU="D", GAC="D", GAA="E", GAG="E",
          GGU="G", GGC="G", GGA="G", GGG="G")

# Function to create progress bar
create_progress_bar <- function(total, prefix = "") {
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  return(list(
    update = function(i) setTxtProgressBar(pb, i),
    close = function() close(pb)
  ))
}

# Function to read FASTA files and extract S gene sequences
read_s_gene_sequences <- function(file) {
  lines <- readLines(file)
  headers <- grep("^>", lines, value = TRUE)
  sequences <- character(length(headers))
  
  # Extract sequences
  seq_start <- which(grepl("^>", lines))
  for (i in seq_along(seq_start)) {
    start <- seq_start[i] + 1
    end <- ifelse(i < length(seq_start), seq_start[i+1] - 1, length(lines))
    sequences[i] <- paste(toupper(lines[start:end]), collapse = "")
  }
  
  # Select only S gene sequences (every 12th sequence starting at 3)
  s_indices <- seq(3, length(headers), by = 12)
  s_headers <- headers[s_indices]
  s_sequences <- sequences[s_indices]
  
  return(list(headers = s_headers, sequences = s_sequences))
}

# Function to compare sequences without full alignment
compare_sequences <- function(ref_seq, var_seq) {
  # Convert to RNA (T to U)
  ref_seq <- gsub("T", "U", ref_seq)
  var_seq <- gsub("T", "U", var_seq)
  
  # Split into codons
  ref_codons <- substring(ref_seq, seq(1, nchar(ref_seq), 3), seq(3, nchar(ref_seq), 3))
  var_codons <- substring(var_seq, seq(1, nchar(var_seq), 3), seq(3, nchar(var_seq), 3))
  
  # Determine the minimum number of codons to compare
  min_codons <- min(length(ref_codons), length(var_codons))
  
  # Initialize results
  mutations <- data.frame(
    codon_position = integer(),
    ref_codon = character(),
    var_codon = character(),
    ref_aa = character(),
    var_aa = character(),
    aa_change = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Compare codon by codon
  for (i in 1:min_codons) {
    ref_codon <- ref_codons[i]
    var_codon <- var_codons[i]
    
    # Skip incomplete codons
    if (nchar(ref_codon) < 3 || nchar(var_codon) < 3) next
    
    # Translate codons
    ref_aa <- ifelse(ref_codon %in% names(trad), trad[ref_codon], "?")
    var_aa <- ifelse(var_codon %in% names(trad), trad[var_codon], "?")
    
    if (ref_codon != var_codon) {
      type <- ifelse(ref_aa == var_aa, "synonymous", "nonsynonymous")
      aa_change <- ifelse(type == "nonsynonymous", 
                         paste0(ref_aa, i, var_aa),
                         paste0(ref_aa, i))
      
      mutations <- rbind(mutations, data.frame(
        codon_position = i,
        ref_codon = ref_codon,
        var_codon = var_codon,
        ref_aa = ref_aa,
        var_aa = var_aa,
        aa_change = aa_change,
        type = type,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(mutations)
}

# Main analysis function
analyze_variants <- function(ref_file, variant_files) {
  # Read reference sequence (assuming it's the S gene)
  ref_data <- readLines(ref_file)
  ref_seq <- paste(toupper(ref_data[-1]), collapse = "")
  
  # Initialize results
  all_mutations <- data.frame(
    variant = character(),
    sequence_id = character(),
    codon_position = integer(),
    ref_codon = character(),
    var_codon = character(),
    ref_aa = character(),
    var_aa = character(),
    aa_change = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each variant file
  for (file in variant_files) {
    cat("\nProcessing file:", file, "\n")
    variant_name <- sub(".fasta$", "", basename(file))
    
    # Read and extract S gene sequences
    s_data <- read_s_gene_sequences(file)
    total_sequences <- length(s_data$sequences)
    
    # Setup progress bar
    pb <- create_progress_bar(total_sequences, prefix = paste0(variant_name, ": "))
    
    # Process each S gene sequence
    for (i in seq_along(s_data$sequences)) {
      # Compare with reference
      mutations <- compare_sequences(ref_seq, s_data$sequences[i])
      
      if (nrow(mutations) > 0) {
        mutations$variant <- variant_name
        mutations$sequence_id <- sub("^>", "", s_data$headers[i])
        all_mutations <- rbind(all_mutations, mutations)
      }
      
      # Update progress
      pb$update(i)
    }
    
    pb$close()
  }
  
  return(all_mutations)
}

# Run the analysis
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

results <- analyze_variants(ref_file, variant_files)

# Analyze and plot most prevalent mutations by variant
if (nrow(results) > 0) {
  # Calculate mutation frequencies
  mutation_freq <- results %>%
    group_by(variant, codon_position, aa_change, type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(variant) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup()
  
  # Get top 10 most frequent mutations per variant
  top_mutations <- mutation_freq %>%
    group_by(variant) %>%
    arrange(desc(count)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Plot with base R graphics
  par(mar = c(10, 4, 4, 2))
  layout(matrix(1:length(unique(top_mutations$variant)), ncol = 1))
  
  for (v in unique(top_mutations$variant)) {
    variant_data <- top_mutations %>% filter(variant == v)
    barplot(variant_data$count, names.arg = variant_data$aa_change, 
            las = 2, main = paste("Top Mutations -", v),
            ylab = "Count", col = "steelblue")
  }
  
  # Print summary statistics
  cat("\n=== MUTATION SUMMARY ===\n")
  cat("Total mutations detected:", nrow(results), "\n")
  cat("\nBreakdown by variant:\n")
  print(table(results$variant))
  cat("\nBreakdown by mutation type:\n")
  print(table(results$type))
  
  # Show most common mutations overall
  cat("\nMost frequent amino acid changes overall:\n")
  overall_counts <- results %>% 
    group_by(aa_change) %>% 
    summarise(count = n()) %>% 
    arrange(desc(count))
  print(head(overall_counts, 20))
} else {
  cat("No mutations found in any sequences.\n")
}
