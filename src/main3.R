# Working directory
setwd("/home/gero/Documents/projects/tec/biologia/project")

# Cargar librerías
library(seqinr)
library(dplyr)
library(ggplot2)
library(Biostrings)
library(pwalign)
library(foreach)
library(doFuture)
library(future)

plan(multisession)  
registerDoFuture()
# handlers(global = TRUE)

# Archivos
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

#  Función para extraer glicoproteína de superficie
# extract_surface_glycoproteins <- function(filepath) {
#     fasta_seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)
#     headers <- names(fasta_seqs)
#     is_surface <- grepl("\\|surface glycoprotein", headers)
#     if (sum(is_surface) == 0) stop("No se ha encontrado ningún gen S en el archivo fasta.")
#     surface_seqs <- fasta_seqs[is_surface]
#     sequences <- unlist(surface_seqs)
#     names(sequences) <- headers[is_surface]
#     return(sequences)
# }

extract_surface_glycoproteins <- function(filepath, max_seqs = 10) {
    # Leer archivo FASTA
    fasta_seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE, whole.header = TRUE, forceDNAtolower = FALSE)
    headers <- names(fasta_seqs)

    # Encontrar secuencias que contengan "surface glycoprotein" en el encabezado
    is_surface <- grepl("\\|surface glycoprotein", headers)

    if (sum(is_surface) == 0) {
        stop("No surface glycoprotein genes found in the file")
    }

    # Extraer secuencias de glicoproteína de superficie
    surface_seqs <- fasta_seqs[is_surface]

    # Limitar a las primeras max_seqs secuencias
    if (length(surface_seqs) > max_seqs) {
        surface_seqs <- surface_seqs[1:max_seqs]
        cat(sprintf("Limited to first %d sequences from %s\n", max_seqs, basename(filepath)))
    }

    sequences <- unlist(surface_seqs)
    names(sequences) <- names(surface_seqs)
    return(sequences)
}


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
                                   gapOpening = -2,
                                   gapExtension = -2,
                                   type = "global"
    )
    aln_ref <- as.character(subject(alignment))
    aln_var <- as.character(pattern(alignment))
    return(c(aln_ref, aln_var))
}

#  Inicializar tabla de mutaciones
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

#  Comparar una secuencia de variante contra la referencia
compare_sequences <- function(ref_seq, var_seq, variant_name) {
    mutations_local <- data.frame(
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

    if (nchar(ref_seq) != nchar(var_seq)) {
        aligned <- needleman_wunsch(ref_seq, var_seq)
        aln_ref <- aligned[1]
        aln_var <- aligned[2]
    } else {
        aln_ref <- ref_seq
        aln_var <- var_seq
    }

    len <- min(nchar(aln_ref), nchar(aln_var))
    len <- len - (len %% 3)

    for (i in seq(1, len, by = 3)) {
        ref_codon <- substr(aln_ref, i, i+2)
        var_codon <- substr(aln_var, i, i+2)
        if (grepl("N", var_codon) || grepl("-", ref_codon) || grepl("-", var_codon)) next

        ref_aa <- seqinr::translate(s2c(ref_codon))
        var_aa <- seqinr::translate(s2c(var_codon))
        aa_change <- paste0(ref_aa, ceiling(i/3), var_aa)

        if (ref_codon != var_codon) {
            mut_type <- ifelse(ref_aa == var_aa, "synonymous", "non-synonymous")
            mutations_local <- rbind(mutations_local, data.frame(
                codon_position = i/3,
                ref_codon = ref_codon,
                var_codon = var_codon,
                ref_aa = ref_aa,
                var_aa = var_aa,
                aa_change = aa_change,
                type = mut_type,
                variant = variant_name,
                stringsAsFactors = FALSE
            ))
        }
    }

    return(mutations_local)
}

#  Cargar secuencia de referencia
ref_s_seq <- extract_surface_glycoproteins(ref_file)
ref_seq <- ref_s_seq[[1]]

#  Procesar variantes

mutation_results <- foreach(file = variant_files, .combine = rbind) %do% {
    variant_name <- gsub(".fasta", "", basename(file))
    var_seqs <- extract_surface_glycoproteins(file)

    foreach(i = seq_along(var_seqs), .combine = rbind) %dofuture% {
        var_seq <- var_seqs[[i]]
        compare_sequences(ref_seq, var_seq, variant_name)
    }
}
mutations <- mutation_results


#  Guardar resultados
write.csv(mutations, "mutations_summary.csv", row.names = FALSE)

#  Graficar las 10 mutaciones más frecuentes por variante
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
