# Load necessary libraries
library(dplyr)
library(ggseqlogo)

# Filter out any rows with NA in motif_sequence_RNA
candidate_sequences_clean <- candidate_sequences_new_corrected %>%
  filter(!is.na(motif_sequence_RNA))

# Split sequences by nuclease class
sequences_by_class <- candidate_sequences_clean %>%
  group_by(nuclease_class) %>%
  summarise(seqs = list(motif_sequence_RNA), .groups = "drop")

# Function to generate cleavage-centered motif logo
plot_cleavage_logo <- function(class_name, output_file) {
  
  # Extract sequences for this class
  seqs <- sequences_by_class %>%
    filter(nuclease_class == class_name) %>%
    pull(seqs) %>%
    unlist()
  
  # Determine motif length and cleavage site index (middle of motif)
  motif_len <- nchar(seqs[1])          # assume all motifs same length
  cleavage_index <- ceiling(motif_len / 2)  # cleavage site is center
  rel_positions <- 1:motif_len - cleavage_index  # -n ... 0 ... +n
  
  # Generate logo
  p <- ggseqlogo(seqs, method = "bits") +
    scale_x_continuous(
      breaks = seq(1, motif_len, by = 5),
      labels = rel_positions[seq(1, motif_len, by = 5)]
    ) +
    ggtitle(paste("Motif Logo:", class_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10)
    )
  
  # Save as PNG
  ggsave(output_file, plot = p, width = 8, height = 4)
  
  return(p)
}

# Generate logos for each nuclease class and save as PNGs
logo_RNase_Y_high <- plot_cleavage_logo("RNase_Y_high_confidence", "logo_RNase_Y_high_confidence.png")
logo_RNase_Y_cand <- plot_cleavage_logo("RNase_Y_candidate", "logo_RNase_Y_candidate.png")
logo_other_nuc <- plot_cleavage_logo("other_endonuclease", "logo_other_endonuclease.png")

