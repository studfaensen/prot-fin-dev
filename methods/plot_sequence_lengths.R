args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    message("Usage: Rscript plot_sequence_lengths.R <fasta-file> <out-file.png>")
    q()
}

# Load necessary libraries
library(seqinr)
library(ggplot2)
library(ggdist)

# Read the FASTA file
fasta_file <- args[1]
sequences <- read.fasta(file = fasta_file, seqtype = "AA", as.string = TRUE)

# Extract sequence lengths
sequence_lengths <- sapply(sequences, nchar)

plt <- ggplot(data.frame(length = sequence_lengths), aes(x = "", y = length)) +
  stat_halfeye(fill = "darkgray", width = 2) +
  geom_boxplot(width = 0, outlier.colour = "black", outlier.size = 0.1) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank()) +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "Frequencies", y = "Lengths") +
  ggtitle("Raincloud Plot of Sequence Lengths")

ggsave(args[2], plot = plt, width = 8, height = 3, dpi = 300)
