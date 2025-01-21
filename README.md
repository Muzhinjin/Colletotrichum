# Colletotrichum

CAZ enzymes
# Load necessary libraries
library(ggplot2)
library(extrafont)
library(reshape2)

# Create the data
data <- data.frame(
  Isolates = c("6Ca", "34Ca", "40Ca", "28Cg", "57Cg", "58Cg", "84Cg"),
  Species = c("C. nymphaeae", "C. nymphaeae", "C. nymphaeae", 
              "C. siamense", "C. siamense", "C. siamense", "C. siamense"),
  AA = c(196, 198, 166, 178, 187, 185, 202),
  CBM = c(9, 9, 9, 13, 12, 12, 13),
  CE = c(70, 73, 66, 60, 57, 63, 64 ),
  GH = c(379, 378, 368, 404, 410, 413, 413),
  GT = c(103, 101, 102, 102, 102, 105, 110),
  PL = c(43, 42, 42, 46, 47, 47, 47)
)

# Convert data to long format
data_long <- pivot_longer(data, cols = c(AA, CBM, CE, GH, GT, PL), 
                          names_to = "Compartment", values_to = "Value")

# Customize colors
compartment_colors <- c("AA" = "orange", 
                        "CBM" = "brown", 
                        "CE" = "cyan",
                        "GH" = "green", 
                        "GT" = "lightcoral",
                        "PL" = "gray")

# Create the stacked bar plot
plot3 <- ggplot(data_long, aes(x = Isolates, y = Value, fill = Compartment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = compartment_colors) +
  facet_wrap(~Species, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 12, colour = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 12, colour = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 12, colour = "black"),
    legend.text = element_text(family = "Times New Roman", size = 12, colour = "black"),
    legend.title = element_text(family = "Times New Roman", size = 12, colour = "black"),
    strip.text = element_text(face = "italic", family = "Times New Roman", size = 12, colour = "black"),
    plot.title = element_text(family = "Times New Roman", size = 12, colour = "black", hjust = 0.5)
  ) +
  labs(x = "Isolates", y = "CAZymes number", fill = "CAZyme type",
       title = "A")
plot3

library(pheatmap)
library(extrafont)  # Font support
library(grid)       # For grid text adjustments

# Load fonts (ensure Times New Roman is installed)
loadfonts(device = "win")  # Use for Windows, adjust for macOS/Linux

# Create the data frame
data <- data.frame(
  CAZyme_type = c("AA1:CW(Lignin)", "AA2:CW(Lignin)", "AA3:CW(Cellulose-Lignin)", "AA7", 
                  "AA9: CW(Cellulose-Hemicellulose-Lignin)", "CBM1-CW(Cellulose)", "CBM18", 
                  "CBM91", "CE3:Hemicellulose", "CE4-Hemicellulose", "CE5: Hemicellulose", 
                  "GH13:FCW+ESR(α-glucans)", "GH16:FCW(β-glycans)", "GH18-FCW(Chitin)", 
                  "GH2:CW( Cellulose)", "GH28-PCW(Pectin)", "GH3: CW (β-glycans)", 
                  "GH31:PG+ESR+PCW(hemicellulose)", "GH43: PCW(Hemicellulose + pectin)", 
                  "GH47:PG(N-/O-glycans)", "GH5-PCW(Cellulose+ Hemicellulose)", "GH71: FCW(β-1,3-glucan)", 
                  "GH76-PCW(Chitin)", "GH78-PCW(Pectin)", "GT1", "GT2", "PL1:PCW(Pectin)", "PL3:PCW(Pectin)"),
  `6Ca` = c(22, 11, 43, 7, 25, 11, 16, 17, 1, 6, 15, 13, 20, 16, 7, 16, 30, 11, 17, 11, 19, 4, 11, 15, 16, 20, 11, 11),
  `34Ca` = c(22, 7, 43, 7, 21, 1, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 20, 4, 11, 15, 3, 20, 11, 11),
  `40Ca` = c(21, 7, 42, 7, 21, 7, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 18, 4, 11, 15, 17, 20, 11, 11),
  `57Cg` = c(25, 8, 52, 7, 18, 15, 25, 16, 0, 4, 20, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 4, 10, 17, 16, 17, 13, 11),
  `58Cg` = c(25, 8, 52, 7, 20, 15, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 5, 10, 17, 16, 18, 13, 11),
  `28Cg` = c(25, 8, 52, 7, 17, 2, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 10, 21, 10, 18, 5, 10, 18, 23, 18, 13, 11),
  `84Cg` = c(27, 12, 49, 79, 28, 0, 0, 0, 11, 14, 18, 13, 26, 25, 10, 21, 31, 10, 43, 9, 22, 12, 11, 10, 25, 20, 17, 14)
)

# Prepare the data for pheatmap
rownames(data) <- data$CAZyme_type  # Use CAZyme types as row names
data_clean <- data[, -1]  # Exclude the first column

# Define species annotations (group isolates)
species_annotation <- data.frame(
  Species = factor(c("C. nymphaeae", "C. nymphaeae", "C. nymphaeae", 
                     "C. siamense", "C. siamense", "C. siamense", "C. siamense")),
  row.names = colnames(data_clean)  # Explicitly set row names
)

# Define annotation colors
annotation_colors <- list(
  Species = c("C. nymphaeae" = "brown", "C. siamense" = "lightblue")
)
# Remove the "X" prefix from isolate column names (if it exists)
colnames(data_clean) <- gsub("^X", "", colnames(data_clean))

# Define the pheatmap plot
Plot4 <- pheatmap(
  as.matrix(data_clean),                # Heatmap data
  annotation_col = species_annotation, # Column annotations
  annotation_colors = annotation_colors, # Colors for annotations
  fontsize = 12,                        # Font size
  color = colorRampPalette(c("white", "orange"))(50), # Heatmap colors
  cluster_cols = FALSE,                 # Disable column clustering
  cluster_rows = FALSE,                 # Disable row clustering
  display_numbers = formatted_data,     # Display formatted numbers
  annotation_legend = TRUE,             # Show legend for annotations
  number_color = "black",               # Color for numbers in cells
  fontsize_row = 12,                    # Row font size
  fontsize_col = 12,                    # Column font size
  angle_col = "45",                     # Rotate column names
  fontfamily = "Times New Roman"        # Use Times New Roman font
)

# Save the plot (optional)
ggsave("heatmap_with_annotations.png", plot = Plot4, width = 10, height = 8)


# Add species names above the heatmap panels
grid.text(expression(italic("C. nymphaeae")), 
          x = unit(0.1, "npc"), y = unit(1, "npc"), 
          gp = gpar(fontsize = 12, fontfamily = "Times New Roman", col = "black"))
grid.text(expression(italic("C. siamense")), 
          x = unit(0.35, "npc"), y = unit(1, "npc"), 
          gp = gpar(fontsize = 12, fontfamily = "Times New Roman", col = "black"))


library(pheatmap)
library(extrafont)  # Font support
library(grid)       # For grid text adjustments

# Load fonts (ensure Times New Roman is installed)
loadfonts(device = "win")  # Use for Windows, adjust for macOS/Linux

# Create the data frame
data <- data.frame(
  CAZyme_type = c("AA1:CW(Lignin)", "AA2:CW(Lignin)", "AA3:CW(Cellulose-Lignin)", "AA7", 
                  "AA9: CW(Cellulose-Hemicellulose-Lignin)", "CBM1-CW(Cellulose)", "CBM18", 
                  "CBM91", "CE3:Hemicellulose", "CE4-Hemicellulose", "CE5: Hemicellulose", 
                  "GH13:FCW+ESR(α-glucans)", "GH16:FCW(β-glycans)", "GH18-FCW(Chitin)", 
                  "GH2:CW( Cellulose)", "GH28-PCW(Pectin)", "GH3: CW (β-glycans)", 
                  "GH31:PG+ESR+PCW(hemicellulose)", "GH43: PCW(Hemicellulose + pectin)", 
                  "GH47:PG(N-/O-glycans)", "GH5-PCW(Cellulose+ Hemicellulose)", "GH71: FCW(β-1,3-glucan)", 
                  "GH76-PCW(Chitin)", "GH78-PCW(Pectin)", "GT1", "GT2", "PL1:PCW(Pectin)", "PL3:PCW(Pectin)"),
  `6Ca` = c(22, 11, 43, 7, 25, 11, 16, 17, 1, 6, 15, 13, 20, 16, 7, 16, 30, 11, 17, 11, 19, 4, 11, 15, 16, 20, 11, 11),
  `34Ca` = c(22, 7, 43, 7, 21, 1, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 20, 4, 11, 15, 3, 20, 11, 11),
  `40Ca` = c(21, 7, 42, 7, 21, 7, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 18, 4, 11, 15, 17, 20, 11, 11),
  `57Cg` = c(25, 8, 52, 7, 18, 15, 25, 16, 0, 4, 20, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 4, 10, 17, 16, 17, 13, 11),
  `58Cg` = c(25, 8, 52, 7, 20, 15, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 5, 10, 17, 16, 18, 13, 11),
  `28Cg` = c(25, 8, 52, 7, 17, 2, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 10, 21, 10, 18, 5, 10, 18, 23, 18, 13, 11),
  `84Cg` = c(27, 12, 49, 79, 28, 0, 0, 0, 11, 14, 18, 13, 26, 25, 10, 21, 31, 10, 43, 9, 22, 12, 11, 10, 25, 20, 17, 14)
)

data_clean <- data[, -1]  # Exclude 'Species' column for the heatmap
species_annotation <- data.frame(Species = data$Species, row.names = colnames(data_clean))

# Define annotation colors
annotation_colors <- list(Species = c("C. nymphaeae" = "brown", "C. siamense" = "lightblue"))

# Create the heatmap
pheatmap(
  as.matrix(data_clean),
  annotation_col = species_annotation,
  annotation_colors = annotation_colors,
  fontsize = 12,
  color = colorRampPalette(c("white", "red"))(50),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = "45",
  fontfamily = "Times New Roman"
)

#Effectors graph

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(ggtext) # Ensure this is loaded for element_markdown
library(gridExtra)

# Create the data
data <- data.frame(
  Isolates = c("6Ca", "34Ca", "40Ca", "28Cg", "57Cg", "58Cg", "84Cg"),
  Species = c("C. nymphaeae", "C. nymphaeae", "C. nymphaeae", 
              "C. siamense", "C. siamense", "C. siamense", "C. siamense"),
  Apoplastic = c(250, 250, 245, 315, 328, 329, 321),
  Cytoplasmic = c(45, 48, 42, 72, 69, 69, 75),
  Apoplastic_cytoplasmic = c(52, 54, 52, 84, 81, 82, 83)
)

# Convert data to long format
data_long <- pivot_longer(data, cols = c(Apoplastic, Cytoplasmic, Apoplastic_cytoplasmic), 
                          names_to = "Compartment", values_to = "Value")

# Customize colors
compartment_colors <- c("Apoplastic" = "orange", 
                        "Cytoplasmic" = "brown", 
                        "Apoplastic_cytoplasmic" = "cyan")

# Create the stacked bar plot
plot1 <- ggplot(data_long, aes(x = Isolates, y = Value, fill = Compartment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = compartment_colors) +
  facet_wrap(~Species, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.text = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.title = element_text(family = "Times New Roman", size = 14, colour = "black"),
    strip.text = element_text(face = "italic", family = "Times New Roman", size = 14, colour = "black"),
    plot.title = element_text(family = "Times New Roman", size = 14, colour = "black", hjust = 0.5)
  ) +
  labs(x = "Isolates", y = "Number of effectors", fill = "Type of effectors",
       title = "A")
plot1

#Secondary Metabolites

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Create the data frame with isolate names and species
data <- data.frame(
  Isolates = c("6Ca", "34Ca", "40Ca", "28Cg", "57Cg", "58Cg", "84Cg"),
  Species = c("C. nymphaeae", "C. nymphaeae", "C. nymphaeae", 
              "C. siamense", "C. siamense", "C. siamense", "C. siamense"),
  Fungal_RiPP_like = c(13, 12, 12, 16, 13, 14, 15),
  Fungal_RiPP_like_NRPS_like = c(0, 1, 1, 2, 2, 1, 1),
  Fungal_RiPP_like_T1PKS = c(0, 4, 4, 3, 4, 4, 4),
  Fungal_RiPP_like_T3PKS  = c(0, 1, 1, 0, 0, 0, 0),
  T3PKS_T1PKS = c(0, 0, 0, 1, 1, 1, 1),
  Betalactone = c(1, 1, 1, 0, 0, 0, 0),
  T3PKS = c(0, 0, 0, 1, 1, 1, 1),
  Indole = c(6, 6, 6, 5, 5, 6, 5),
  Isocyanide = c(1, 1, 1, 1, 1, 1, 1),
  NRPS = c(7, 7, 7, 7, 7, 7, 7),
  NRPS_betalactone = c(1, 1, 1, 1, 1, 0, 0),
  NRPS_T1PKS = c(3, 3, 2, 4, 7, 7, 7),
  NRPS_terpene = c(2, 1, 1, 1, 1, 1, 1),
  NRPS_like = c(7, 7, 7, 9, 9, 9, 9),
  T1PKS = c(15, 16, 16, 25, 25, 24, 26),
  Terpene = c(11, 11, 11, 14, 10, 12, 13)
)

# Reshape the data to long format for ggplot
data_long <- data %>%
  gather(key = "Biosynthetic_class", value = "Count", -Isolates, -Species)

# Remove isolates with zero counts across all biosynthetic classes
data_long <- data_long %>%
  filter(Count > 0)

# Manually define colors for each biosynthetic cluster
colors <- c(
  Fungal_RiPP_like = "blue",
  Fungal_RiPP_like_NRPS_like = "red",
  Fungal_RiPP_like_T1PKS = "green",
  Fungal_RiPP_like_T3PKS = "purple",
  T3PKS_T1PKS = "orange",
  Betalactone = "darkgreen",
  T3PKS = "yellow",
  Indole = "brown",
  Isocyanide = "gray",
  NRPS = "cyan",
  NRPS_betalactone = "mistyrose",
  NRPS_T1PKS = "magenta",
  NRPS_terpene = "ivory",
  NRPS_like = "salmon3",
  T1PKS = "tan4",
  Terpene = "wheat"
)


# Create the stacked bar plot with facet wrap by species
plot2 <-ggplot(data_long, aes(x = Isolates, y = Count, fill = Biosynthetic_class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Species, scales = "free_x") +  # Facet by species and free the x-axis scale
  scale_fill_manual(values = colors) +  # Use the defined colors for biosynthetic classes
  theme_minimal(base_family = "Times New Roman", base_size = 14) +  # Set font and size
  labs(title = "", 
       x = "Isolates", 
       y = "Number of SM backbone genes") +
  
  theme(
    axis.text.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.text.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.text = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.title = element_text(family = "Times New Roman", size = 14, colour = "black"),
    strip.text = element_text(face = "italic", family = "Times New Roman", size = 14, colour = "black"),
    plot.title = element_text(family = "Times New Roman", size = 14, colour = "black", hjust = 0.5)
  )

plot2

ggsave("plot2.pdf", plot = plot2, width = 12, height = 8, dpi = 300)

#Arrange the plots in a 2-column, 1-row layout
grid.arrange(plot1, plot2, ncol = 2)



colors()
# Load required libraries
library(pheatmap)
library(extrafont)  # Font support
library(grid)       # For grid text adjustments

# Load fonts (ensure Times New Roman is installed)
loadfonts(device = "win")  # Use for Windows, adjust for macOS/Linux

# Create the data frame
data <- data.frame(
  CAZyme_type = c("AA1:CW(Lignin)", "AA2:CW(Lignin)", "AA3:CW(Cellulose-Lignin)", "AA7", 
                  "AA9: CW(Cellulose-Hemicellulose-Lignin)", "CBM1-CW(Cellulose)", "CBM18", 
                  "CBM91", "CE3:Hemicellulose", "CE4-Hemicellulose", "CE5: Hemicellulose", 
                  "GH13:FCW+ESR(α-glucans)", "GH16:FCW(β-glycans)", "GH18-FCW(Chitin)", 
                  "GH2:CW( Cellulose)", "GH28-PCW(Pectin)", "GH3: CW (β-glycans)", 
                  "GH31:PG+ESR+PCW(hemicellulose)", "GH43: PCW(Hemicellulose + pectin)", 
                  "GH47:PG(N-/O-glycans)", "GH5-PCW(Cellulose+ Hemicellulose)", "GH71: FCW(β-1,3-glucan)", 
                  "GH76-PCW(Chitin)", "GH78-PCW(Pectin)", "GT1", "GT2", "PL1:PCW(Pectin)", "PL3:PCW(Pectin)"),
  `6Ca` = c(22, 11, 43, 7, 25, 11, 16, 17, 1, 6, 15, 13, 20, 16, 7, 16, 30, 11, 17, 11, 19, 4, 11, 15, 16, 20, 11, 11),
  `34Ca` = c(22, 7, 43, 7, 21, 1, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 20, 4, 11, 15, 3, 20, 11, 11),
  `40Ca` = c(21, 7, 42, 7, 21, 7, 1, 17, 1, 6, 15, 13, 20, 11, 7, 16, 30, 11, 17, 11, 18, 4, 11, 15, 17, 20, 11, 11),
  `57Cg` = c(25, 8, 52, 7, 18, 15, 25, 16, 0, 4, 20, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 4, 10, 17, 16, 17, 13, 11),
  `58Cg` = c(25, 8, 52, 7, 20, 15, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 11, 22, 12, 18, 5, 10, 17, 16, 18, 13, 11),
  `28Cg` = c(25, 8, 52, 7, 17, 2, 21, 16, 0, 6, 7, 12, 21, 10, 10, 19, 30, 10, 21, 10, 18, 5, 10, 18, 23, 18, 13, 11),
  `84Cg` = c(27, 12, 49, 79, 28, 0, 0, 0, 11, 14, 18, 13, 26, 25, 10, 21, 31, 10, 43, 9, 22, 12, 11, 10, 25, 20, 17, 14)
)

# Prepare the data for pheatmap
rownames(data) <- data$CAZyme_type  # Use CAZyme types as row names
data_clean <- data[, -1]  # Exclude the first column

# Define species annotations (group isolates)
species_annotation <- data.frame(
  Species = factor(c("C. nymphaeae", "C. nymphaeae", "C. nymphaeae", 
                     "C. siamense", "C. siamense", "C. siamense", "C. siamense"))
)

rownames(species_annotation) <- colnames(data_clean)

# Italicize species names
species_legend <- list(Species = c("C. nymphaeae" = "brown", "C. siamense" = "lightblue"))
names(species_legend$Species) <- expression(italic("C. nymphaeae"), italic("C. siamense"))

# Remove the "X" from isolate column names
colnames(data)[-1] <- gsub("^X", "", colnames(data)[-1])

# Format numbers without .00
formatted_data <- apply(data_clean, c(1, 2), function(x) sprintf("%.0f", x))


# Generate the heatmap
pheatmap(as.matrix(data_clean),
         annotation_col = species_annotation, # Add species grouping
         annotation_colors = list(Species = c("C. nymphaeae" = "brown", "C. siamense" = "lightblue")),
         fontsize = 12,
         color = colorRampPalette(c("white", "blue", "red"))(50),
         cluster_cols = F,
         cluster_rows = FALSE,
         display_numbers = formatted_data, # Display formatted numbers
         number_color = "black",
         fontsize_row = 12,
         fontsize_col = 12,
         angle_col = "45",
         main = "CAZyme Distribution Heatmap",
         fontfamily = "Times New Roman")









# Load necessary libraries
library(ggplot2)
library(tidyr)

# Create the data frame
data <- data.frame(
  CAZyme = c("AA", "CBM", "CE", "GH", "GT", "PL"),
  `6Ca` = c(196, 9, 70, 379, 103, 43),
  `34Ca` = c(198, 9, 73, 378, 101, 42),
  `40Ca` = c(166, 9, 66, 368, 102, 42),
  `28Cg` = c(178, 13, 60, 404, 102, 46),
  `57Cg` = c(187, 12, 57, 410, 102, 47),
  `58Cg` = c(185, 12, 63, 413, 105, 47),
  `84Cg` = c(202, 15, 64, 413, 110, 47)
)

# Reshape the data to long format
data_long <- data %>%
  gather(key = "Isolate", value = "Count", -CAZyme)

# Create a new column to classify isolates into species groups
data_long$Species <- ifelse(data_long$Isolate %in% c("6Ca", "34Ca", "40Ca"), "C. nymphaeae", "C. siamense")

# Plot the stacked bar plot with isolates clustered by species grouping
ggplot(data_long, aes(x = Isolate, y = Count, fill = CAZyme)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Plot of CAZymes by Species and Isolate",
       x = "Isolate", y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("AA" = "#1f77b4", "CBM" = "#ff7f0e", "CE" = "#2ca02c",
                               "GH" = "#d62728", "GT" = "#9467bd", "PL" = "#8c564b")) +
  scale_x_discrete(labels = c("6Ca" = "C. nymphaeae", "34Ca" = "C. nymphaeae", "40Ca" = "C. nymphaeae", 
                              "28Cg" = "C. siamense", "57Cg" = "C. siamense", "58cg" = "C. siamense", "84Cg" = "C. siamense")) +
  facet_wrap(~Species, scales = "free_x", ncol = 2)  # Facet by species grouping (C. nymphaeae vs. C. siamense)




# Sample gene count data based on the provided figure (full dataset)
input <- c(
  "84Cg&58Cg&57Cg&28Cg&40Ca&34Ca&6Ca" = 10538,               # 84Cg only
  "84Cg&58Cg&57Cg&28Cg" = 5463,           # 84Cg & 58Cg
  "40Ca&34Ca&6Ca" = 3210,      # 84Cg, 57Cg & 58Cg
  "57Cg&58Cg" = 1085, # 84Cg, 28Cg, 57Cg & 58Cg
  "84Cg&28Cg" = 803,       # 84Cg, 28Cg & 40Ca
  "34Ca&40Ca&57Cg&58Cg" = 101,       # 84Cg, 58Cg & 40Ca
  "84Cg&28Cg&34Ca&48Ca" = 75,             # 57Cg & 28Cg
  "34Ca&6Ca" = 46,                  # 40Ca only
  "28Cg&57Cg&84Cg" = 39,                  # 34Ca only
  "57Cg&28Cg&84Cg" = 33,                   # 6Ca only
  "58Cg&57Cg&28Cg" = 32,        # 84Cg, 58Cg & 28Cg
  "84Cg&58Cg&34Ca&6Ca&40Ca" = 28,        # 84Cg, 57Cg & 34Ca
  "57Cg&58Cg&84Cg" = 23, # 84Cg, 58Cg, 57Cg, 28Cg & 40Ca
  "28Cg&57Cg" = 21,         # 84Cg, 57Cg & 6Ca
  "34Ca" = 19,        # 84Cg, 40Ca & 34Ca
  "6Ca" = 18,# 58Cg, 28Cg & 6Ca
  "40Ca" = 18,
  "57Cg" = 18,
  "58Cg" = 18,
  "84Cg" = 18
)


# Load the UpSetR package
library(UpSetR)

# Generate the UpSet plot
upset(fromExpression(input), 
      nintersects = 20,               # Show all 16 intersections
      nsets = 7,                      # Total number of sets
      order.by = "freq",              # Order intersections by frequency
      decreasing = TRUE,              # Show in descending order
      mb.ratio = c(0.6, 0.4),         # Adjust ratio of main bar to sets
      number.angles = 0,              # Angle of number labels
      text.scale = 2,               # Scale the text size
      point.size = 3.8,               # Adjust size of dots in intersections
      line.size = 1                  # Adjust line size connecting dots
)


# Input list
input <- c(
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 6364,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 1077,
  "A.alternata SRC1lrK2f&&A.arborescens BMP0308" = 791,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 784,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 432,
  "&A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 352,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25" = 235,
  "A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 179,
  "A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 162,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43" = 144,
  "A.dauci A2016&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 133,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 110,
  "A.linariae 25&A.solani NL03003" = 106,
  "A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 105,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.tomatophila BMP2032" = 97,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 83,
  "A.dauci A2016&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 75,
  "A.brassicicola Abra43&A.tomatophila BMP2032" = 69,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308" = 64,
  "A.solani NL03003" = 58,
)


# Convert list to fromExpression format
upset_data <- UpSetR::fromExpression(input)

# Plotting with UpSetR
upset(
  upset_data, 
  nintersects = 20,               # Show all intersections
  nsets = 7,                       # Total number of sets
  order.by = "freq",               # Order intersections by frequency
  decreasing = TRUE,               # Show in descending order
  mb.ratio = c(0.6, 0.4),          # Adjust ratio of main bar to sets
  number.angles = 0,               # Angle of number labels
  text.scale = 2,                  # Scale the text size
  point.size = 3.8,                # Adjust size of dots in intersections
  line.size = 1                    # Adjust line size connecting dots
)



# Load necessary library
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

# Input list with cleaned duplicates
input <- c(
  "A.dauci&A.alternata&A.brassicicola&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 6753,  
  "A.dauci&A.alternata&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 1055,
  "A.alternata&A.arborescens BMP0308&A.burnsii CBS10738" = 484,
  "A.alternata&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 454,
  "A.dauci&A.alternata&A.brassicicola&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003" = 360,  
  "A.alternata&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 333,  
  "A.alternata&A.arborescens BMP0308" = 314,
  "A.linariae AL-1&A.solani NL03003" = 245,             
  "A.dauci&A.alternata&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003" = 239,
  "A.alternata&A.burnsii CBS10738" = 202,
  "A.alternata&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003" = 127,  
  "A.dauci&A.linariae AL-1&A.solani NL03003" = 124,
  "A.arborescens BMP0308&A.alternata&A.burnsii CBS10738" = 117,
  "A.alternata&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 89,
  "A.dauci&A.alternata&A.brassicicola&A.arborescens BMP0308&A.burnsii CBS10738" = 89,
  "A.alternata&A.brassicicola&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003" = 80,
  "A.dauci&A.alternata&A.brassicicola&A.arborescens BMP0308&A.linariae AL-1&A.burnsii CBS10738" = 74,
  "A.arborescens BMP0308&A.burnsii CBS10738" = 70,  
  "A.dauci&A.alternata&A.brassicicola&A.arborescens BMP0308" = 70,  
  "A.dauci&A.brassicicola&A.arborescens BMP0308&A.linariae AL-1&A.solani NL03003&A.burnsii CBS10738" = 67
)

# Convert list to fromExpression format
upset_data <- UpSetR::fromExpression(input)

# Plotting with UpSetR
upset(
  upset_data, 
  nintersects = 20,               # Show all intersections
  nsets = 7,                       # Total number of sets
  order.by = "freq",               # Order intersections by frequency
  decreasing = TRUE,               # Show in descending order
  mb.ratio = c(0.6, 0.4),          # Adjust ratio of main bar to sets
  number.angles = 0,               # Angle of number labels
  text.scale = 2,                  # Scale the text size
  point.size = 3.8,                # Adjust size of dots in intersections
  line.size = 1                    # Adjust line size connecting dots
)

# Save the plot as a PNG
png("upset_plot.png", width = 800, height = 600)
upset(
  upset_data, 
  nintersects = 20,
  nsets = 7,
  order.by = "freq",
  decreasing = TRUE,
  mb.ratio = c(0.6, 0.4),
  number.angles = 0,
  text.scale = 2,
  point.size = 3.8,
  line.size = 1
)
dev.off()  # Close the device


# Install UpSetR package if not already installed
if (!require("UpSetR")) install.packages("UpSetR")

# Load required library
library(UpSetR)

# Define the input data as a named vector
input <- c(
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 6364,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 1077,
  "A.alternata SRC1lrK2f&A.arborescens BMP0308" = 791,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 784,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 432,
  "A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 352,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 235,
  "A.linariae 25&A.tomatophila BMP2032" = 179,
  "A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 162,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43" = 144,
  "A.dauci A2016&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 133,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 110,
  "A.linariae 25&A.solani NL03003" = 106,
  "A.alternata SRC1lrK2f&A.arborescens BMP0308&A.linariae 25&A.solani NL03003" = 105,
  "A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308&A.tomatophila BMP2032" = 97,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 83,
  "A.dauci A2016&A.brassicicola Abra43&A.arborescens BMP0308&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 75,
  "A.brassicicola Abra43&A.tomatophila BMP2032" = 69,
  "A.dauci A2016&A.alternata SRC1lrK2f&A.brassicicola Abra43&A.arborescens BMP0308" = 64,
  "A.solani NL03003" = 58
)

# Convert list to fromExpression format
upset_data <- UpSetR::fromExpression(input)

# Plotting with UpSetR
upset(
  upset_data, 
  nintersects = 20,               # Show all intersections
  nsets = 7,                       # Total number of sets
  matrix.color = "#800080", main.bar.color = "#800000",
  sets.bar.color = "#008080", sets.x.label = "Cluster count",
  mainbar.y.label = "Cluster count", mainbar.y.max = NULL,
  order.by = "freq",               # Order intersections by frequency
  decreasing = TRUE,               # Show in descending order
  mb.ratio = c(0.6, 0.4),          # Adjust ratio of main bar to sets
  number.angles = 0,               # Angle of number labels
  text.scale = 2,                  # Scale the text size
  point.size = 3.8,                # Adjust size of dots in intersections
  line.size = 1                    # Adjust line size connecting dots
)



# Define the input data as a named vector
input <- c(
  "A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 83,
  "A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 162,
  "A.alternata SRC1lrK2f&A.linariae 25" = 27,
  "A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003" = 48,
  "A.alternata SRC1lrK2f&A.solani NL03003" = 18,
  "A.alternata SRC1lrK2f&A.solani NL03003&A.tomatophila BMP2032" = 8,
  "A.linariae 25&A.solani NL03003" = 106,
  "A.alternata SRC1lrK2f&A.tomatophila BMP2032" = 45,
  "A.alternata SRC1lrK2f&A.linariae 25&A.tomatophila BMP2032" = 25,
  "A.linariae 25&A.tomatophila BMP2032" = 179,
  "A.solani NL03003&A.tomatophila BMP2032" = 29,
  "A.solani NL03003" = 0,
  "A.linariae 25" = 20,
  "A.alternata SRC1lrK2f" = 16,
  "A.tomatophila BMP2032" = 2
)









# Load necessary libraries
library(ggplot2)
library(reshape2)

# Load necessary libraries
library(ggplot2)
library(reshape2)
# Load necessary libraries
library(ggplot2)
library(reshape2)

# Create the data frame
data <- data.frame(
  Species = c("A. linariae 25", "A. tomatophila BMP2032", "A. solani NL03003", "A. alternata SRC1lrK2f", 
              "A. arborescens FERA675", "A. brassicicola Abra43", "A. dauci A2016"),
  AA = c(114, 82, 105, 114, 112, 102, 117),
  CBM = c(40, 34, 38, 49, 50, 39, 42),
  CE = c(42, 23, 41, 44, 45, 39, 39),
  GH = c(241, 213, 235, 236, 235, 236, 246),
  GT = c(117, 92, 112, 103, 109, 101, 101),
  PL = c(19, 16, 22, 25, 23, 22, 25)
)

# Reshape data for ggplot
data_long <- melt(data, id.vars = "Species", variable.name = "Treatment", value.name = "Count")

# Plot with customized colors, black axis labels, and italicized "Alternaria" in x-axis title
ggplot(data_long, aes(x = Species, y = Count, fill = Treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("AA" = "#1f77b4", "CBM" = "#ff7f0e", "CE" = "#2ca02c", 
                               "GH" = "lightblue", "GT" = "#9467bd", "PL" = "#8c564b")) +
  labs(
    x = expression(italic("Alternaria") ~ "species"),
    y = "Number of CAZymes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5)
  )



# Load necessary libraries
library(ggplot2)
library(reshape2)

# Create the data frame with manually formatted expressions for Species
data <- data.frame(
  Species = factor(c(
    expression(italic("A. linariea") ~ plain("25")),
    expression(italic("A. tomatophila") ~ plain("BMP2032")),
    expression(italic("A. solani") ~ plain("NL03003")),
    expression(italic("A. alternata") ~ plain("SRC1lrK2f")),
    expression(italic("A. arborescens") ~ plain("FERA675")),
    expression(italic("A. brassicicola") ~ plain("Abra43")),
    expression(italic("A. dauci") ~ plain("A2016"))
  )),
  AA = c(114, 82, 105, 114, 112, 102, 117),
  CBM = c(40, 34, 38, 49, 50, 39, 42),
  CE = c(42, 23, 41, 44, 45, 39, 39),
  GH = c(241, 213, 235, 236, 235, 236, 246),
  GT = c(117, 92, 112, 103, 109, 101, 101),
  PL = c(19, 16, 22, 25, 23, 22, 25)
)

# Reshape data for ggplot
data_long <- melt(data, id.vars = "Species", variable.name = "Treatment", value.name = "Count")

# Plot with customized colors, black axis labels, and formatted Species labels
ggplot(data_long, aes(x = Species, y = Count, fill = Treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("AA" = "#1f77b4", "CBM" = "#ff7f0e", "CE" = "#2ca02c", 
                               "GH" = "lightblue", "GT" = "#9467bd", "PL" = "#8c564b")) +
  labs(
    x = expression(italic("Alternaria") ~ "species"),
    y = "Number of CAZymes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 12),
    plot.title = element_text(hjust = 0.5)
  )







# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggtext) # For element_markdown

# Create the data frame
data <- data.frame(
  Species = c("A. linariae 25", "A. tomatophila BMP2032", "A. solani NL03003", "A. alternata SRC1lrK2f", 
              "A. arborescens FERA675", "A. brassicicola Abra43", "A. dauci A2016"),
  AA = c(114, 82, 105, 114, 112, 102, 117),
  CBM = c(40, 34, 38, 49, 50, 39, 42), 
  CE = c(42, 23, 41, 44, 45, 39, 39),
  GH = c(241, 213, 235, 236, 235, 236, 246),
  GT = c(117, 92, 112, 103, 109, 101, 101),
  PL = c(19, 16, 22, 25, 23, 22, 25)
)

# Reshape data for ggplot
data_long <- melt(data, id.vars = "Species", variable.name = "CAZyme", value.name = "Count")

# Calculate total for each species
data_totals <- data_long %>%
  group_by(Species) %>%
  summarise(Total = sum(Count))

# Convert species names to italicize the species part after "A."
data_long$Species <- gsub("^(A\\.) (\\w+)", "*\\1* *\\2*", data_long$Species)
data_totals$Species <- gsub("^(A\\.) (\\w+)", "*\\1* *\\2*", data_totals$Species)

# Plot with totals on the bars and italicized species names
ggplot(data_long, aes(y = Species, x = Count, fill = CAZyme)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = data_totals, aes(x = Total + 20, y = Species, label = Total), # Adjust '20' for spacing
            color = "black", size = 5, family = "Times New Roman", inherit.aes = FALSE) +
  scale_fill_manual(values = c("AA" = "#1f77b4", "CBM" = "#ff7f0e", "CE" = "#2ca02c", 
                               "GH" = "lightblue", "GT" = "#9467bd", "PL" = "#8c564b")) +
  labs(
    y = expression(italic("Alternaria") ~ "species"),
    x = "Number of CAZymes"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(family = "Times New Roman", size = 14, colour = "black"), # Italicize species names
    axis.text.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.text = element_text(family = "Times New Roman", size = 14, colour = "black"), # Legend text styling
    plot.title = element_text(hjust = 0.5)
  )
library(ggtext)

# Install and load pheatmap package if you haven't already
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# Load your data into R
data <- read.csv("path/to/your_file.csv", row.names = 1)  # Make sure the path is correct
# If your data is separated by tabs instead of commas, use sep = "\t" in read.csv()

# Create the heatmap
pheatmap(data,
         cluster_rows = TRUE,        # Cluster the rows
         cluster_cols = TRUE,        # Cluster the columns
         scale = "none",             # Scale data; options are "row", "column", "none"
         color = colorRampPalette(c("white", "blue", "red"))(100),  # Customize colors
         display_numbers = TRUE,     # Display values in the cells
         fontsize_number = 8,        # Adjust font size for numbers
         fontsize_row = 8,           # Font size for row labels
         fontsize_col = 8,           # Font size for column labels
         main = "Gene Cluster Heatmap")  # Title for the heatmap


# Sample data (provided)
data <- data.frame(
  Gene_cluster = c("T1PKS", "T1PKS", "T1PKS", "T1PKS", "T1PKS", "T1PKS/NRPS", "T1PKS/NRPS", 
                   "T1PKS/NRPS", "T1PKS/NRPS", "NRPS-like/indole/T1PKS", "T1PKS", "T1PKS", 
                   "Terpene", "T1PKS", "Terpene", "NRPS-like", "T1PKS/NRPS", "NRPS", 
                   "T1PKS", "T1PKS/NRPS", "T1PKS", "NRPS-like", "T1PKS", "Fungal-RiPP-like",
                   "NRPS", "Terpene", "Fungal-RiPP-like", "T1PKS", "T1PKS"),
  Putative_product = c("Alternariol", "1,3,6,8-tetrahydroxynaphthalene", "Alternapyrone", 
                       "Ankaflavin/monascin", "Betaenone", "FR901512", "Asperthecin", 
                       "Equisetin", "Metachelin C/dimerumic acid", "Ustethylin A", 
                       "Solanapyrone D", "Abscisic acid", "Squalestatin1", 
                       "Dehydrocurvularin", "Giberellic acid", "Choline", 
                       "Desmethlbassianin", "Pseurotin/azaspiren", "Gregatin", 
                       "ACT Toxin", "Zearalenone", "Cichorine", "Nivalenol/vomitoxin", 
                       "Oryzine/A/oryzine B", "Cyclo-(D-Phe-L-Phe-D-Val-L-Val)", "PR toxin", 
                       "Polymyxin", "Scytalone/T3HN", "Depudecin"),
  A_linariea_25 = c(100, 100, 80, 8, 0, 100, 28, 0, 0, 20, 0, 25, 40, 0, 0, 100, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_tomatophila_BMP2032 = c(100, 100, 80, 0, 0, 0, 0, 0, 37, 0, 0, 25, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_solani_NL03003 = c(100, 100, 100, 50, 0, 0, 0, 0, 50, 20, 100, 25, 40, 50, 28, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_alternata_SRC1lrK2f = c(100, 100, 80, 0, 85, 0, 0, 36, 50, 0, 0, 25, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_arborescens_FERA675 = c(100, 100, 80, 0, 0, 0, 0, 36, 50, 0, 0, 37, 40, 0, 0, 100, 0, 40, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 66),
  A_brassicicola_Abra43 = c(0, 100, 0, 0, 85, 0, 0, 0, 37, 0, 0, 0, 40, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 12, 100, 50, 40, 40, 0),
  A_dauci_A2016 = c(100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 40, 0, 0, 100, 0, 0, 22, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0)
)
# Remove 'Gene_Cluster' and 'Putative_Product' columns to create a numeric matrix
heatmap_data <- as.matrix(data[, -c(1, 2)])
rownames(heatmap_data) <- paste(data$Gene_Cluster, data$Putative_Product, sep = ": ")

# Create annotations for gene clusters as a facet
gene_cluster_annotation <- data.frame(Gene_Cluster = data$Gene_Cluster)
rownames(gene_cluster_annotation) <- rownames(heatmap_data)

# Define custom colors for each gene cluster
cluster_colors <- c("T1PKS" = "orange", 
                    "T1PKS/NRPS" = "cyan", 
                    "NRPS-like/indole/T1PKS" = "pink", "Terpene" = "yellow", 
                    "NRPS-like" = "brown", "NRPS" = "red", "Fungal-RiPP-like" = "darkgrey")

# Create the annotation and assign the colors for gene clusters
gene_cluster_annotation <- data.frame(Gene_Cluster = data$Gene_Cluster)
rownames(gene_cluster_annotation) <- rownames(heatmap_data)
annotation_colors <- list(Gene_Cluster = cluster_colors)
color_palette <- colorRampPalette(c("white", "blue", "red"))(100)
# Generate the heatmap
pheatmap(heatmap_data,
         annotation_row = gene_cluster_annotation,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = color_palette,
         fontsize_row = 14,
         fontsize_col = 14,
         labels_col = italic_colnames,
         fontsize = 14,
         fontfamily = "Times", # Times New Roman font
         display_numbers = T) # No numbers in heatmap cells



# Load necessary library
library(pheatmap)

# Sample data (provided)
data <- data.frame(
  Gene_cluster = c("T1PKS", "T1PKS", "T1PKS", "T1PKS", "T1PKS", "T1PKS/NRPS", "T1PKS/NRPS", 
                   "T1PKS/NRPS", "T1PKS/NRPS", "NRPS-like/indole/T1PKS", "T1PKS", "T1PKS", 
                   "Terpene", "T1PKS", "Terpene", "NRPS-like", "T1PKS/NRPS", "NRPS", 
                   "T1PKS", "T1PKS/NRPS", "T1PKS", "NRPS-like", "T1PKS", "Fungal-RiPP-like",
                   "NRPS", "Terpene", "Fungal-RiPP-like", "T1PKS", "T1PKS"),
  Putative_product = c("Alternariol", "1,3,6,8-tetrahydroxynaphthalene", "Alternapyrone", 
                       "Ankaflavin/monascin", "Betaenone", "FR901512", "Asperthecin", 
                       "Equisetin", "Metachelin C/dimerumic acid", "Ustethylin A", 
                       "Solanapyrone D", "Abscisic acid", "Squalestatin1", 
                       "Dehydrocurvularin", "Giberellic acid", "Choline", 
                       "Desmethlbassianin", "Pseurotin/azaspiren", "Gregatin", 
                       "ACT Toxin", "Zearalenone", "Cichorine", "Nivalenol/vomitoxin", 
                       "Oryzine/A/oryzine B", "Cyclo-(D-Phe-L-Phe-D-Val-L-Val)", "PR toxin", 
                       "Polymyxin", "Scytalone/T3HN", "Depudecin"),
  A_linariea_25 = c(100, 100, 80, 8, 0, 100, 28, 0, 0, 20, 0, 25, 40, 0, 0, 100, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_tomatophila_BMP2032 = c(100, 100, 80, 0, 0, 0, 0, 0, 37, 0, 0, 25, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_solani_NL03003 = c(100, 100, 100, 50, 0, 0, 0, 0, 50, 20, 100, 25, 40, 50, 28, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_alternata_SRC1lrK2f = c(100, 100, 80, 0, 85, 0, 0, 36, 50, 0, 0, 25, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  A_arborescens_FERA675 = c(100, 100, 80, 0, 0, 0, 0, 36, 50, 0, 0, 37, 40, 0, 0, 100, 0, 40, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 66),
  A_brassicicola_Abra43 = c(0, 100, 0, 0, 85, 0, 0, 0, 37, 0, 0, 0, 40, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 12, 100, 50, 40, 40, 0),
  A_dauci_A2016 = c(100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 40, 0, 0, 100, 0, 0, 22, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0)
)

# Remove 'Gene_Cluster' and 'Putative_Product' columns to create a numeric matrix
heatmap_data <- as.matrix(data[, -c(1, 2)])
rownames(heatmap_data) <- paste(data$Gene_cluster, data$Putative_product, sep = ": ")

# Create annotations for gene clusters as a facet
gene_cluster_annotation <- data.frame(Gene_cluster = data$Gene_cluster)
rownames(gene_cluster_annotation) <- rownames(heatmap_data)

# Define custom colors for each gene cluster
cluster_colors <- c("T1PKS" = "orange", 
                    "T1PKS/NRPS" = "cyan", 
                    "NRPS-like/indole/T1PKS" = "pink", 
                    "Terpene" = "yellow", 
                    "NRPS-like" = "brown", 
                    "NRPS" = "red", 
                    "Fungal-RiPP-like" = "darkgrey")

# Create the annotation and assign the colors for gene clusters
annotation_colors <- list(Gene_cluster = cluster_colors)

# Define color palette for the heatmap
color_palette <- colorRampPalette(c("white", "blue", "red"))(100)

italic_colnames <- c(
  expression(italic("A. linariea") ~ 25),
  expression(italic("A. tomatophila") ~ BMP2032),
  expression(italic("A. solani") ~ NL03003),
  expression(italic("A. alternata") ~ SRC1lrK2f),
  expression(italic("A. arborescens") ~ FERA675),
  expression(italic("A. brassicicola") ~ Abra43),
  expression(italic("A. dauci") ~ A2016)
)

# Generate the heatmap
pheatmap(heatmap_data,
         annotation_row = gene_cluster_annotation,
         annotation_colors = annotation_colors,
         cluster_rows = T,
         cluster_cols = T,
         scale = "none",
         color = color_palette,
         fontsize_row = 14,
         fontsize_col = 14,
         labels_col = italic_colnames,
         fontsize = 14,
         cellwidth = 12,
         angle_col = 45,
         cellheight = 6,
         fontfamily = "Times", # Times New Roman font
         display_numbers = FALSE) # Set to FALSE to hide numbers in heatmap cells







# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Sample data frame (replace with actual data loading if needed)
data <- data.frame(
  Gene_cluster = c("NRPS", "TIPKS", "NRP-metallophore", "fungal-RiPP-like", "NRPS-T1PKS", "NRPS-like,indole,T1PKS", 
                   "Terpene", "isocyanide-nrp", "T3PKS", "NRPS-like", "NAPAA"),
  A_solani_NL03003 = c(2, 10, 1, 10, 2, 1, 4, 1, 1, 6, 0),
  A_linariea_25 = c(5, 13, 0, 6, 4, 0, 4, 1, 0, 4, 0),
  A_alternata_SRC1lrK2f = c(3, 8, 0, 4, 2, 0, 3, 1, 1, 4, 1),
  A_arborescens_FERA675 = c(5, 7, 0, 4, 4, 0, 3, 1, 1, 5, 0),
  A_tomatophila_BMP2032 = c(8, 12, 0, 4, 0, 1, 3, 1, 0, 3, 0),
  A_brassicicola_Abra43 = c(4, 3, 0, 8, 0, 0, 6, 1, 1, 4, 0),
  A_dauci_A2016 = c(4, 8, 0, 4, 0, 1, 3, 1, 0, 5, 0)
)

# Reshape the data to long format for ggplot
data_long <- data %>%
  pivot_longer(cols = -Gene_cluster, names_to = "Species", values_to = "Count")

# Format species names in italics for plotting
data_long$Species <- factor(data_long$Species, levels = unique(data_long$Species))

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Sample data (replace this with your actual data or load from a CSV)
data <- data.frame(
  Cluster = c("NRPS", "TIPKS", "NRP_metallophore", "fungal_RiPP_like", "NRPS_T1PKS", 
              "NRPS_like_indole_T1PKS", "Terpene", "isocyanide_nrp", "T3PKS", 
              "NRPS_like", "NAPAA"),
  A_solani_NL03003 = c(2, 10, 1, 10, 2, 1, 4, 1, 1, 6, 0),
  A_linariea_25 = c(5, 13, 0, 6, 4, 0, 4, 1, 0, 4, 0),
  A_alternata_SRC1lrK2f = c(3, 8, 0, 4, 2, 0, 3, 1, 1, 4, 1),
  A_arborescens_FERA675 = c(5, 7, 0, 4, 4, 0, 3, 1, 1, 5, 0),
  A_tomatophila_BMP2032 = c(8, 12, 0, 4, 0, 1, 3, 1, 0, 3, 0),
  A_brassicicola_Abra43 = c(4, 3, 0, 8, 1, 0, 6, 1, 1, 4, 0),
  A_dauci_A2016 = c(4, 8, 0, 4, 0, 1, 3, 1, 0, 5, 0)
)

# Melt data to long format
data_long <- melt(data, id.vars = "Cluster", variable.name = "Species", value.name = "Count")

# Custom colors for each gene cluster
custom_colors <- c("NRPS" = "#1f77b4", "TIPKS" = "#ff7f0e", "NRP_metallophore" = "#2ca02c",
                   "fungal_RiPP_like" = "#d62728", "NRPS_T1PKS" = "#9467bd", "NRPS_like_indole_T1PKS" = "#8c564b",
                   "Terpene" = "#e377c2", "isocyanide_nrp" = "#7f7f7f", "T3PKS" = "#bcbd22",
                   "NRPS_like" = "#17becf", "NAPAA" = "#aec7e8")

# Create stacked bar plot
ggplot(data_long, aes(x = Species, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
labs(
  x = expression(italic("Alternaria") ~ "species"),
  y = "Number of CAZymes"
) +
  scale_fill_manual(values = custom_colors) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 14),
    plot.title = element_text(hjust = 0.5)
  )

labs(title = "Secondary Metabolites Gene Clusters Across Species",
     x = "Species",
     y = "Count") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 12, color = "black"),  # Set font style and size
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggtext) # For element_markdown

# Create the data frame
data <- data.frame(
  Species = c("A. linariae 25", "A. tomatophila BMP2032", "A. solani NL03003", "A. alternata SRC1lrK2f", 
              "A. arborescens FERA675", "A. brassicicola Abra43", "A. dauci A2016"),
  AA = c(114, 82, 105, 114, 112, 102, 117),
  CBM = c(40, 34, 38, 49, 50, 39, 42), 
  CE = c(42, 23, 41, 44, 45, 39, 39),
  GH = c(241, 213, 235, 236, 235, 236, 246),
  GT = c(117, 92, 112, 103, 109, 101, 101),
  PL = c(19, 16, 22, 25, 23, 22, 25)
)



library(ggplot2)
library(dplyr)

# Create the data frame
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create the data frame
data <- data.frame(
  Species = c("A. solani NL03003", "A. linariea 25", "A. alternata SRC1lrK2f", 
              "A. arborescens FERA675", "A. tomatophila BMP2032", 
              "A. brassicicola Abra43", "A. dauci A2016"),
  NRPS = c(2, 5, 3, 5, 8, 4, 4),
  TIPKS = c(10, 13, 8, 7, 12, 3, 8),
  NRP_metallophore = c(1, 0, 0, 0, 0, 0, 0),
  Fungal_RiPP_like = c(10, 6, 4, 4, 4, 8, 4),
  NRPS_T1PKS = c(2, 0, 2, 4, 0, 0, 0),
  NRPS_like_indole_T1PKS = c(1, 4, 0, 0, 1, 0, 1),
  Terpene = c(4, 0, 3, 3, 3, 6, 3),
  Isocyanide_nrp = c(1, 4, 1, 1, 1, 1, 0),
  T3PKS = c(1, 0, 0, 0, 0, 4, 5),
  NRPS_like = c(6, 4, 4, 5, 3, 0, 0),
  NAPAA = c(0, 0, 1, 0, 0, 0, 0)
)
# Reshape data for ggplot
data_long <- melt(data, id.vars = "Species", variable.name = "SM", value.name = "Count")


# Calculate total for each species
data_totals <- data_long %>%
  group_by(Species) %>%
  summarise(Total = sum(Count))

# Convert species names to italicize the species part after "A."
data_long$Species <- gsub("^(A\\.) (\\w+)", "*\\1* *\\2*", data_long$Species)
data_totals$Species <- gsub("^(A\\.) (\\w+)", "*\\1* *\\2*", data_totals$Species)

# Plot with totals on the bars and italicized species names
ggplot(data_long, aes(y = Species, x = Count, fill = SM)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = data_totals, aes(x = Total + 20, y = Species, label = Total), # Adjust '20' for spacing
            color = "black", size = 5, family = "Times New Roman", inherit.aes = FALSE) +
  scale_fill_manual(values = c("NRPS" = "#1f77b4", "TIPKS" = "#ff7f0e", "NRP_metallophore" = "#2ca02c",
                               "Fungal_RiPP_like" = "#d62728", "NRPS_T1PKS" = "#9467bd", "NRPS_like_indole_T1PKS" = "#8c564b",
                               "Terpene" = "#e377c2", "Isocyanide_nrp" = "#7f7f7f", "T3PKS" = "#bcbd22",
                               "NRPS_like" = "#17becf", "NAPAA" = "#aec7e8")) +
  labs(
    y = expression(italic("Alternaria") ~ "species"),
    x = "Secondary  metabolites (SM) gene clusters"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(family = "Times New Roman", size = 14, colour = "black"), # Italicize species names
    axis.text.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black"),
    axis.title.x = element_text(family = "Times New Roman", size = 14, colour = "black"),
    legend.text = element_text(family = "Times New Roman", size = 14, colour = "black"), # Legend text styling
    plot.title = element_text(hjust = 0.5)
  )

# Load the VennDiagram package
library(VennDiagram)

# Define the input data as a named vector
input <- c(
  "A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 83,
  "A.linariae 25&A.solani NL03003&A.tomatophila BMP2032" = 162,
  "A.alternata SRC1lrK2f&A.linariae 25" = 27,
  "A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003" = 48,
  "A.alternata SRC1lrK2f&A.solani NL03003" = 18,
  "A.alternata SRC1lrK2f&A.solani NL03003&A.tomatophila BMP2032" = 8,
  "A.linariae 25&A.solani NL03003" = 106,
  "A.alternata SRC1lrK2f&A.tomatophila BMP2032" = 45,
  "A.alternata SRC1lrK2f&A.linariae 25&A.tomatophila BMP2032" = 25,
  "A.linariae 25&A.tomatophila BMP2032" = 179,
  "A.solani NL03003&A.tomatophila BMP2032" = 29,
  "A.solani NL03003" = 0,
  "A.linariae 25" = 20,
  "A.alternata SRC1lrK2f" = 16,
  "A.tomatophila BMP2032" = 2
)

# Create the Venn diagram using the specified counts for each set and intersection
venn.plot <- venn.diagram(
  x = list(
    "A.alternata SRC1lrK2f" = input["A.alternata SRC1lrK2f"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25"] + 
      input["A.alternata SRC1lrK2f&A.solani NL03003"] + 
      input["A.alternata SRC1lrK2f&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"],
    
    "A.linariae 25" = input["A.linariae 25"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25"] + 
      input["A.linariae 25&A.solani NL03003"] + 
      input["A.linariae 25&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.tomatophila BMP2032"] + 
      input["A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"],
    
    "A.solani NL03003" = input["A.solani NL03003"] + 
      input["A.alternata SRC1lrK2f&A.solani NL03003"] + 
      input["A.linariae 25&A.solani NL03003"] + 
      input["A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003"] + 
      input["A.alternata SRC1lrK2f&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"],
    
    "A.tomatophila BMP2032" = input["A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.tomatophila BMP2032"] + 
      input["A.linariae 25&A.tomatophila BMP2032"] + 
      input["A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"] + 
      input["A.alternata SRC1lrK2f&A.linariae 25&A.solani NL03003&A.tomatophila BMP2032"]
  ),
  
  category.names = c("A.alternata SRC1lrK2f", "A.linariae 25", "A.solani NL03003", "A.tomatophila BMP2032"),
  filename = NULL,
  output = TRUE,
  
  # Set text size for numbers
  cex = 1.5,
  fill = c("skyblue", "pink1", "mediumorchid", "orange"),
  alpha = 0.5,
  lty = "blank"
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Create the data frame
data <- data.frame(
  Sample = c("AA0", "AA6", "AA8", "CBM42", "CBM43", "GH45", "GH53", "GH55", "GT25", 
             "GT50", "GT59", "GT76", "PL11", "PL26", "PL42", "CBM36", "CBM37", "CBM38", 
             "CBM39", "CBM63", "GH33"),
  A_linariea = c(5, 0, 1, 0, 0, 2, 1, 1, 4, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Tomatophila = c(0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0),
  A_silani = c(3, 1, 1, 1, 2, 2, 1, 4, 4, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
  A_alternata = c(2, 1, 1, 1, 2, 2, 1, 4, 4, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1),
  A_bauci = c(7, 0, 1, 1, 2, 2, 1, 4, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
  A_brasicola = c(5, 1, 1, 1, 2, 2, 1, 4, 3, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0),
  Aborecsebe = c(5, 1, 1, 1, 2, 2, 1, 4, 3, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0)
)

# Calculate the differences between each species and A. linariea
difference_data <- data
difference_data[-1] <- sweep(data[-1], 1, data$A_linariea, "-")

# Melt the data for ggplot2
difference_melted <- melt(difference_data, id.vars = "Sample")

# Create a heatmap to show the differences
ggplot(difference_melted, aes(x = variable, y = Sample, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Difference in Abundance Compared to A. linariea",
       x = "Species",
       y = "Sample") +
  theme_minimal()



# Display the differences
print(difference_data)
