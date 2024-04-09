args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    message("Usage: Rscript plot_sequence_lengths.R <fasta-file> <out-file.png>")
    q()
}

# Load necessary libraries
library(seqinr)
library(ggplot2)
library(ggdist)
library(tibble)

# Read the FASTA file
fasta_files <- args[-length(args)]
sequence_lengths <- lapply(fasta_files, function(f) {
  sequences <- read.fasta(file = f, seqtype = "AA", as.string = TRUE, seqonly=TRUE)
  sapply(sequences, nchar)
})

groups <- lapply(seq_along(sequence_lengths), function(i) {tibble(value = sequence_lengths[[i]], group = fasta_files[i])})

data <- dplyr::bind_rows(groups)

plt <- ggplot(data, aes(x = group, y = value)) +
  geom_jitter(position = position_jitter(width = 0.1, seed = 1), size = 0.1, colour = "blue") +
  stat_halfeye(fill = "darkgray", width = 2, position = position_nudge(x = 0.11, y=0)) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank()) +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "Frequencies", y = "Lengths") +
  ggtitle("Raincloud Plot of Sequence Lengths")

ggsave(args[length(args)], plot = plt, width = 8, height = 3, dpi = 300)
