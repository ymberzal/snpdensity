# snpdensity
# Load required libraries
library(vcfR)
library(ggplot2)
library(dplyr)
library(gtools)  # For mixedsort

# Step 1: Load the VCF file
vcf <- read.vcfR("/media/yami/Elements/DDRAD_SAHIWAL/arsucd/3vcf/bcftools/vcf_files_sahiwal/sahiwal20snponlyfilternamed.vcf")

# Step 2: Extract chromosome and position data
chrom <- vcf@fix[, 1]  # Chromosome column
pos <- as.numeric(vcf@fix[, 2])  # Position column

# Step 3: Create a data frame with chromosome and position
snp_data <- data.frame(chrom = chrom, pos = pos)

# Step 4: Create bins (e.g., 1 Mb bins) for each chromosome
bin_size <- 1e6  # 1 Mb bins

# Create bins for each chromosome
snp_data <- snp_data %>%
  group_by(chrom) %>%
  mutate(bin = floor(pos / bin_size) * bin_size)

# Step 5: Ensure chromosomes are ordered naturally
# First, make sure chromosomes are treated as factors and ordered correctly
snp_data$chrom <- factor(snp_data$chrom, levels = mixedsort(unique(snp_data$chrom)))

# Step 6: Count SNPs per chromosome bin
snp_density <- snp_data %>%
  group_by(chrom, bin) %>%
  summarize(snp_count = n())

# Step 7: Plot SNP density as a heatmap with chromosomes in natural order
ggplot(snp_density, aes(x = bin, y = chrom, fill = snp_count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Color scale for SNP density
  labs(x = "Genomic Position (bins)", y = "Chromosome", fill = "SNP Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Load the required libraries
library(vcfR)
library(ggplot2)
library(dplyr)
library(gtools)

# Step 1: Load the VCF file
# Replace with your actual VCF file path
vcf <- read.vcfR("/media/yami/Elements/DDRAD_SAHIWAL/arsucd/3vcf/bcftools/vcf_files_sahiwal/sahiwal20snponlyfilternamed.vcf")

# Step 2: Extract chromosome and position data from the VCF
chrom <- vcf@fix[, 1]  # Chromosome column
pos <- as.numeric(vcf@fix[, 2])  # Position column

# Step 3: Create a data frame with chromosomes
snp_data <- data.frame(chrom = chrom, pos = pos)

# Step 4: Convert chromosome names to a factor and use mixedsort for natural ordering
snp_data$chrom <- factor(snp_data$chrom, levels = mixedsort(unique(snp_data$chrom)))

# Step 5: Count SNPs per chromosome
snp_density <- snp_data %>%
  group_by(chrom) %>%
  summarize(snp_count = n())

# Step 6: Plot SNP density by chromosome
ggplot(snp_density, aes(x = chrom, y = snp_count, fill = chrom)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "SNP Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
