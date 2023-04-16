# Load the installed Package
library(readr)
library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(viridis)
library(ggsci)
library(RColorBrewer)
library(factoextra)


# Reading the contents of TSV file using read_tsv() method
otu_table <-readr::read_tsv("C:/Users/jens_/OneDrive/Dokumente/Master/Master Thesis/16S sequencing/Analysis/SCHE23-1.asvs.tsv")

#Deleting, reordering and renaming of columns
otu_clean = subset(otu_table, select = -c(otu, asv, seq, uparse_info))

colnames(otu_clean) <- c('tax', 'DMG1629', 'DMG1630', 'DMG1631',
                         'CTRL1638', 'CTRL1639', 'CTRL1640', 'CTRL1641',
                         'DMG1642', 'DMG1643', 'DMG1644')

col_order <-c('tax', 'DMG1629', 'DMG1630', 'DMG1631','DMG1642', 'DMG1643', 'DMG1644',
              'CTRL1638', 'CTRL1639', 'CTRL1640', 'CTRL1641')
otu_clean <- otu_clean[, col_order]

# Extract taxonomic information from the tax column (change last expression for desired group)
otu_clean$tax <- str_extract(otu_clean$tax, "(?<=\\|family;)[^|]+")
otu_clean$tax <- str_replace(otu_clean$tax, "unclassified_", "")
# Replace NA values with "unclassified" and remove scrap
otu_clean$tax[is.na(otu_clean$tax)] <- "unclassified"
otu_clean$tax <- gsub(';.*', '', otu_clean$tax)
otu_clean$tax <- str_extract(otu_clean$tax, "\\w+")

#Calculate relative abundances
otu_clean$DMG1629 <- (otu_clean$DMG1629/sum(otu_clean$DMG1629))
otu_clean$DMG1630 <- (otu_clean$DMG1630/sum(otu_clean$DMG1630))
otu_clean$DMG1631 <- (otu_clean$DMG1631/sum(otu_clean$DMG1631))
otu_clean$CTRL1638 <- (otu_clean$CTRL1638/sum(otu_clean$CTRL1638))
otu_clean$CTRL1639 <- (otu_clean$CTRL1639/sum(otu_clean$CTRL1639))
otu_clean$CTRL1640 <- (otu_clean$CTRL1640/sum(otu_clean$CTRL1640))
otu_clean$CTRL1641 <- (otu_clean$CTRL1641/sum(otu_clean$CTRL1641))
otu_clean$DMG1642 <- (otu_clean$DMG1642/sum(otu_clean$DMG1642))
otu_clean$DMG1643 <- (otu_clean$DMG1643/sum(otu_clean$DMG1643))
otu_clean$DMG1644 <- (otu_clean$DMG1644/sum(otu_clean$DMG1644))

#Sum up different OTUs with the same name
otu_summed <- otu_clean %>%
  group_by(tax) %>%
  summarize(across(everything(), sum))


# Convert the OTU table into a long format for ggplot2
otu_table_long <- melt(otu_summed)
otu_table_long <- otu_table_long[otu_table_long$value >= 0.001, ]

#Calculate number of groups displayed in the plot
num_nonredundant <- length(unique(otu_table_long$tax))
print(num_nonredundant)

# Create a bar plot
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(data = otu_table_long, aes(x = variable, y = value, fill = tax)) +
        geom_bar(stat = "identity") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
        scale_fill_manual(values = col_vector) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x ='', y = "Relative Abundance", fill = 'Family')+
        theme(legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.key.size = unit(0.25, "cm"))



library(ggfortify)

otu <- data.frame(t(otu_summed[-1]))$
print(c('tax'))
colnames(otu[,2:118]) <- otu_summed[, 1]

pca_res <- prcomp(otu_summed[,2:11], scale. = TRUE)

p <- autoplot(pca_res, data= otu_summed, colour = 'tax')

p +   theme(legend.position = "none")

