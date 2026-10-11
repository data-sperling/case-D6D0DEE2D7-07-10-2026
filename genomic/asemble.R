# Sample Nr. 03-26/10/05 genome visualisation


cran_packages <- c("ggplot2", "dplyr", "readr")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(ggplot2)
library(dplyr)
library(readr)
library(Biostrings)


# should be in root
fasta_files <- list.files(
  pattern = "\\.(fna|fa|fasta)(\\.gz)?$",
  ignore.case = TRUE
)

gff_files <- list.files(
  pattern = "\\.(gff|gff3)(\\.gz)?$",
  ignore.case = TRUE
)

if (length(fasta_files) == 0) {
  stop("No FASTA file found. Check your working directory.")
}

if (length(gff_files) == 0) {
  stop("No GFF/GFF3 file found. Check your working directory.")
}

# explicitly select teh correct file
fasta_file <- fasta_files[1]
gff_file <- gff_files[1]

cat("Using FASTA:", fasta_file, "\n")
cat("Using GFF:", gff_file, "\n")


# read FAST file
genome <- readDNAStringSet(fasta_file, format = "fasta")

# add FASTA sequence identifiers
seq_ids <- names(genome)
seq_lengths <- width(genome)

genome_summary <- data.frame(
  sequence_id = seq_ids,
  length_bp = seq_lengths
) |>
  arrange(desc(length_bp))

write_csv(genome_summary, "beetroot_sequence_lengths.csv")


# plot chromosome vs sequence lengths
p1 <- ggplot(
  genome_summary,
  aes(x = reorder(sequence_id, length_bp), y = length_bp / 1e6)
) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(
    title = "Red beetroot genome sequence lengths",
    x = "Chromosome / sequence",
    y = "Length (Mbp)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "beetroot_genome_lengths.png",
  p1, width = 9, height = 6, dpi = 300
)

print(p1)


# read GFF/GFF3 annotation file as:
# seqid, source, type, start, end, score, strand, phase, attributes
gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = c(
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
  ),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) |>
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  ) |>
  filter(!is.na(start), !is.na(end))


# id records
genes <- gff |>
  filter(tolower(type) == "gene")

if (nrow(genes) == 0) {
  warning(
    "No feature type 'gene' found. Plotting mRNA records instead."
  )
  genes <- gff |>
    filter(tolower(type) %in% c("mrna", "transcript"))
}

if (nrow(genes) == 0) {
  stop("No gene or mRNA/transcript records found in the GFF file.")
}


# plot gene distribution along chromosome
gene_counts <- genes |>
  count(seqid, name = "gene_count") |>
  arrange(desc(gene_count))

write_csv(gene_counts, "beetroot_gene_counts.csv")

p2 <- ggplot(
  gene_counts,
  aes(x = reorder(seqid, gene_count), y = gene_count)
) +
  geom_col(fill = "darkseagreen4") +
  coord_flip() +
  labs(
    title = "Annotated genes per genomic sequence",
    x = "Chromosome / sequence",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "beetroot_gene_counts.png",
  p2, width = 9, height = 6, dpi = 300
)

print(p2)


# plot GC content in windows along each sequence
# def. window size: 1 million bp.
# change window_size to 100k for 100 kb windows.
window_size <- 1000000

gc_list <- lapply(seq_along(genome), function(i) {
  dna <- genome[[i]]
  dna_length <- length(dna)

  starts <- seq.int(1, dna_length, by = window_size)

  do.call(rbind, lapply(starts, function(s) {
    e <- min(s + window_size - 1, dna_length)
    fragment <- subseq(dna, start = s, end = e)

    gc <- letterFrequency(
      fragment,
      letters = c("G", "C"),
      as.prob = FALSE
    )

    valid_bases <- sum(letterFrequency(
      fragment,
      letters = c("A", "C", "G", "T"),
      as.prob = FALSE
    ))

    gc_total <- sum(gc)

    data.frame(
      sequence_id = names(genome)[i],
      start = s,
      end = e,
      midpoint = (s + e) / 2,
      gc_percent = if (valid_bases > 0) {
        100 * gc_total / valid_bases
      } else {
        NA_real_
      }
    )
  }))
})

gc_data <- bind_rows(gc_list)

write_csv(gc_data, "beetroot_gc_content.csv")

p3 <- ggplot(
  gc_data,
  aes(x = midpoint / 1e6, y = gc_percent)
) +
  geom_line(color = "steelblue", linewidth = 0.35) +
  facet_wrap(~ sequence_id, scales = "free_x") +
  labs(
    title = "GC content across the red beetroot genome",
    subtitle = paste("Window size:", format(window_size, big.mark = ","), "bp"),
    x = "Position along sequence (Mbp)",
    y = "GC content (%)"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  "beetroot_gc_content.png",
  p3, width = 12, height = 8, dpi = 300
)

print(p3)

# read from CLI
cat("\nAnalysis complete!\n")
cat("Generated plots:\n")
cat("- beetroot_genome_lengths.png\n")
cat("- beetroot_gene_counts.png\n")
cat("- beetroot_gc_content.png\n")
cat("\nGenerated CSV files:\n")
cat("- beetroot_sequence_lengths.csv\n")
cat("- beetroot_gene_counts.csv\n")
cat("- beetroot_gc_content.csv\n")
