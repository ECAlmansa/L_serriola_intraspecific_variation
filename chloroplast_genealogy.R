# Generate chloroplast genealogy data
# Elvira Castillo Almansa (castilloalmansae@gmail.com)

# Load packages
library(geneHapR)
library(vcfR)
library(ggplot2)
library(scatterpie)
library(rnaturalearth)

## HAPLOTYPE IDENTIFICATION ##

# Import genotype data and annotation file
data_vcf <- read.vcfR("/genotype_data.vcf.gz", convertNA = TRUE)

# Identify haplotypes
haplo <- vcf2hap(data_vcf, na_drop=TRUE)
hapSummary <- hap_summary(haplo)

# Generate haplotype summary
plotHapTable(hapSummary)

## HAPLOTYPE DISTRIBUTION PLOT ##

# Load data
df <- read.csv("/data.csv", header = FALSE)

# Identify haplotype columns
haplos <- colnames(df)[3:ncol(df)]

# Sum haplotype counts per row
df$total <- rowSums(df[haplos])

# Load base map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Color range
okabe_ito <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#999999"
)
color_range <- colorRampPalette(okabe_ito)(38) # for 38 haplotypes

# Plot
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") + # base map
  geom_scatterpie(
    data = df,
    aes(x = Lon, y = Lat, r = sqrt(total) * 0.7), # pie size and location
    cols = haplos # pie sclices
  ) +
  scale_fill_manual(values = color_range, name = "Legend") + # pie slices colors
  coord_sf(xlim = c(-15, 80), ylim = c(30, 60), expand = FALSE) + # extent
  theme_minimal() +
  theme(
    axis.title = element_blank(), # remove title
    axis.text = element_blank(), #remove text
    axis.ticks = element_blank() # remove ticks
  )
