# Set working directory
setwd("/home/gero/Documents/projects/tec/biologia/project")

# Load required packages for parallel processing (included in base R)
library(parallel)
library(utils)
library(doFuture)
library(foreach)
library(future)
library(progressr)


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

# Worker function for parallel processing
process_variant_sequence <- function(seq_data, ref_seq, variant_name) {
  mutations <- compare_sequences(ref_seq, seq_data$sequence)
  if (nrow(mutations) > 0) {
    mutations$variant <- variant_name
    mutations$sequence_id <- seq_data$header
  }
  return(mutations)
}

# Base R alternative to dplyr's group_by and summarise
summarize_mutations <- function(results) {
  # Create unique identifier for each mutation type
  results$mutation_id <- paste(results$variant, results$codon_position, 
                               results$aa_change, results$type, sep = "|")
  
  # Count occurrences
  counts <- table(results$mutation_id)
  
  # Convert to data frame
  freq_df <- data.frame(
    mutation_id = names(counts),
    count = as.numeric(counts),
    stringsAsFactors = FALSE
  )
  
  # Split back into columns
  split_cols <- strsplit(freq_df$mutation_id, "\\|")
  freq_df$variant <- sapply(split_cols, `[`, 1)
  freq_df$codon_position <- as.integer(sapply(split_cols, `[`, 2))
  freq_df$aa_change <- sapply(split_cols, `[`, 3)
  freq_df$type <- sapply(split_cols, `[`, 4)
  
  # Calculate frequencies per variant
  variant_totals <- table(results$variant)
  freq_df$freq <- freq_df$count / as.numeric(variant_totals[freq_df$variant])
  
  # Remove temporary column
  freq_df$mutation_id <- NULL
  
  return(freq_df)
}

# Get top N mutations per variant
get_top_mutations <- function(freq_df, n = 10) {
  variants <- unique(freq_df$variant)
  top_muts <- data.frame()
  
  for (v in variants) {
    variant_data <- freq_df[freq_df$variant == v, ]
    variant_data <- variant_data[order(-variant_data$count), ]
    top_muts <- rbind(top_muts, head(variant_data, n))
  }
  
  return(top_muts)
}

# Main analysis function with parallel processing
analyze_variants_parallel <- function(ref_file, variant_files, workers = parallelly::availableCores() - 1) {
  # Setup future plan
  future::plan(multisession, workers = workers)
  doFuture::registerDoFuture()
  
  ref_data <- readLines(ref_file)
  ref_seq <- paste(toupper(ref_data[-1]), collapse = "")
  
  job_list <- list()
  for (file in variant_files) {
    variant_name <- sub("\\.fasta$", "", basename(file))
    s_data <- read_s_gene_sequences(file)
    for (i in seq_along(s_data$sequences)) {
      job_list[[length(job_list) + 1]] <- list(
        sequence = s_data$sequences[[i]],
        header = sub("^>", "", s_data$headers[[i]]),
        variant = variant_name
      )
    }
  }
  
  total_jobs <- length(job_list)
  cat("Prepared", total_jobs, "sequence comparisons.\n")
  
  # Optional: Enable progress bar
  progressr::handlers(global = TRUE)
  results <- list()
  with_progress({
    p <- progressor(along = job_list)
    results <- foreach(job = job_list, .combine = rbind, .packages = c()) %dofuture% {
      mut <- compare_sequences(ref_seq, job$sequence)
      if (nrow(mut) > 0) {
        mut$variant <- job$variant
        mut$sequence_id <- job$header
      }
      p()  # Progress step
      mut
    }
  })
  
  return(results)
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

# Use parallel processing (adjust cores as needed)
cat("Starting SARS-CoV-2 variant analysis...\n")
cat("Using", 4, "cores for parallel processing\n")

results <- analyze_variants_parallel(ref_file, variant_files, workers = 7)

# Analyze and plot most prevalent mutations by variant
if (nrow(results) > 0) {
  cat("\nAnalyzing mutation patterns...\n")
  
  # Calculate mutation frequencies using base R
  mutation_freq <- summarize_mutations(results)
  
  # Get top 10 most frequent mutations per variant
  top_mutations <- get_top_mutations(mutation_freq, 10)
  
  # Create output directory for plots
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  
  # Plot with base R graphics
  variants <- unique(top_mutations$variant)
  
  # Create individual plots for each variant
  for (v in variants) {
    variant_data <- top_mutations[top_mutations$variant == v, ]
    
    # Create filename
    plot_file <- paste0("plots/", v, "_top_mutations.png")
    png(plot_file, width = 1200, height = 800, res = 120)
    
    par(mar = c(10, 4, 4, 2))
    barplot(variant_data$count, 
            names.arg = variant_data$aa_change, 
            las = 2, 
            main = paste("Top 10 Mutations -", toupper(v)),
            ylab = "Count", 
            col = "steelblue",
            cex.names = 0.8)
    
    dev.off()
    cat("Saved plot:", plot_file, "\n")
  }
  
  # Create combined plot
  png("plots/all_variants_comparison.png", width = 1600, height = 1200, res = 120)
  par(mfrow = c(ceiling(length(variants)/2), 2), mar = c(8, 4, 3, 1))
  
  for (v in variants) {
    variant_data <- top_mutations[top_mutations$variant == v, ]
    barplot(variant_data$count, 
            names.arg = variant_data$aa_change, 
            las = 2, 
            main = paste("Top Mutations -", toupper(v)),
            ylab = "Count", 
            col = "steelblue",
            cex.names = 0.7)
  }
  dev.off()
  cat("Saved combined plot: plots/all_variants_comparison.png\n")
  
  # Print summary statistics
  cat("\n=== MUTATION SUMMARY ===\n")
  cat("Total mutations detected:", nrow(results), "\n")
  cat("\nBreakdown by variant:\n")
  print(table(results$variant))
  cat("\nBreakdown by mutation type:\n")
  print(table(results$type))
  
  # Show most common mutations overall
  cat("\nMost frequent amino acid changes overall:\n")
  aa_counts <- table(results$aa_change)
  print(head(sort(aa_counts, decreasing = TRUE), 20))
  
  # Save detailed results
  write.csv(results, "detailed_mutations.csv", row.names = FALSE)
  write.csv(mutation_freq, "mutation_frequencies.csv", row.names = FALSE)
  write.csv(top_mutations, "top_mutations_per_variant.csv", row.names = FALSE)
  
  cat("\nResults saved to CSV files:\n")
  cat("- detailed_mutations.csv\n")
  cat("- mutation_frequencies.csv\n") 
  cat("- top_mutations_per_variant.csv\n")
  
} else {
  cat("No mutations found in any sequences.\n")
}
