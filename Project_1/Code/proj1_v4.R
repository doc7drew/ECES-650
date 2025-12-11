# rm(list = ls(all = T)) # clears the env in R
#install.packages("vegan")
#install.packages("patchwork")

# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("remotes")
# BiocManager::install("biobakery/maaslin3")

# load required packages
library(vegan)
library(maaslin3)
library(tidyverse)
library(readxl)
library(patchwork)

###################
### file import
###################

# directory variables
file_location <- "C:\\Users\\a-rei\\OneDrive\\Grad school\\Drexel\\4 - FA 2025\\ECES 650 - Statistical Analysis of Genomics\\Project 1\\"
nfc_qiime2 <- "project1PCROutput\\project1PCROutput\\qiime2\\"
nfc_abund_feat_file <- "abundance_tables\\feature-table2.tsv"
nfc_rel_abund_file <- "rel_abundance_tables\\rel-table-ASV_with-DADA2-tax.tsv"
qiime2_abund_file <- "qiime2_final\\feature-table2.tsv"
qiime2_tax_file <- "qiime2_final\\data_taxonomy.tsv"
clinical <- "diet_data\\diet_data\\ABdatabase_correctedJune25_updated.xlsx"

# bring in files
nfc_rel_abund <- read_tsv(paste0(file_location, nfc_qiime2, nfc_rel_abund_file), col_names = TRUE)
nfc_abund <- read_tsv(paste0(file_location, nfc_qiime2, nfc_abund_feat_file), col_names = TRUE) |> select(-PCR89)
qiime2_abund <- read_tsv(paste0(file_location, qiime2_abund_file), col_names = TRUE) |> select(-last_col())
qiime2_tax <- read_tsv(paste0(file_location, qiime2_tax_file), col_names = TRUE)
clin_data <- read_excel(paste0(file_location, clinical), sheet = 2) |> slice(1:n() - 1)

# create datasets for taxonomy & relative abundance from the feature file
str(nfc_rel_abund)
str(nfc_abund)
str(qiime2_abund) # change 

colnames(nfc_rel_abund) <- gsub("a$", "", colnames(nfc_rel_abund)) #remove a at the end of string
colnames(nfc_abund) <- gsub("a$", "", colnames(nfc_abund))

nfc_tax_df <- nfc_rel_abund |>
  tibble::column_to_rownames(var = names(nfc_rel_abund)[1]) |>
  select(-starts_with("PCR"))

nfc_rel_abund_df <- nfc_rel_abund |>
  tibble::column_to_rownames(var = names(nfc_rel_abund)[1]) |>
  select(-c(Kingdom:sequence))

qiime2_tax_df <- qiime2_tax |>
  tibble::column_to_rownames(var = names(qiime2_tax)[1]) |>
  separate(
    col = Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";") |>
  mutate(across(
    .cols = Kingdom:Species,
    .fns = ~ {
      .x |>
        trimws() |> # Remove any leading/trailing spaces
        str_remove(pattern = "^[a-z]__")
    }
  ))


# move OTU column to actual row name
nfc_abund_df <- nfc_abund |>
  tibble::column_to_rownames(var = names(nfc_abund)[1]) |>
  rename(PCR01 = PCR1,
         PCR02 = PCR2,
         PCR03 = PCR3,
         PCR04 = PCR4,
         PCR05 = PCR5,
         PCR06 = PCR6,
         PCR07 = PCR7,
         PCR08 = PCR8,
         PCR09 = PCR9)

rename_map <- clin_data |>
  select(
    new_name = Name,
    current_name = Sample_name) |>
  deframe()

qiime2_abund_df <-qiime2_abund |>
  tibble::column_to_rownames(var = names(qiime2_abund)[1]) |>
  rename(!!!rename_map)

colnames(qiime2_abund_df) <- gsub("_", "", colnames(qiime2_abund_df))

qiime2_abund_df <- qiime2_abund_df |>
  rename(PCR01 = PCR1,
         PCR02 = PCR2,
         PCR03 = PCR3,
         PCR04 = PCR4,
         PCR05 = PCR5,
         PCR06 = PCR6,
         PCR07 = PCR7,
         PCR08 = PCR8,
         PCR09 = PCR9)
  

# transpose so PCR is obs
nfc_rel_abund_df_t <- t(nfc_rel_abund_df)
nfc_abund_df_t <- t(nfc_abund_df)
qiime2_abund_df_t <- t(qiime2_abund_df) 

#qiime2_abund_df_t <- qiime2_abund_df_t[-nrow(qiime2_abund_df_t),]

# inspect and clean datasets
str(nfc_rel_abund_df_t)
str(nfc_abund_df_t)
str(qiime2_abund_df_t)

str(clin_data)

variable.names(clin_data)

clin_data <- clin_data |>
  filter(Group != "discarted") |>
  mutate(across(11:31, as.numeric),
         across(c("Time", "TimeAB", "Diet", "DietAB", "Group", "Group2", "Gender"), as.factor),
         Name = gsub("_", "", Name))

str(clin_data)
# levels(clin_data$DietAB)
# levels(clin_data$Group)
# levels(clin_data$Group2)

# diversity indices
nfc_q2_simpson <- diversity(nfc_abund_df_t, "simpson")
nfc_q2_simpson_df <- data.frame(value = nfc_q2_simpson)

nfc_q2_shannon <- diversity(nfc_abund_df_t, "shannon")
nfc_q2_shannon_df <- data.frame(value = nfc_q2_shannon)

q2_simpson <- diversity(qiime2_abund_df_t, "simpson")
q2_simpson_df <- data.frame(value = q2_simpson)

q2_shannon <- diversity(qiime2_abund_df_t, "shannon")
q2_shannon_df <- data.frame(value = q2_shannon)

nfc_q2_simp_p1 <- ggplot(nfc_q2_simpson_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Simpson Diversity Index", 
       subtitle = "NF-Core QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

nfc_q2_shan_p2 <- ggplot(nfc_q2_shannon_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Shannon Diversity Index - NF-Core QIIME2",
       subtitle = "NF-Core QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

q2_simp_p1 <- ggplot(q2_simpson_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Simpson Diversity Index - QIIME2", 
       subtitle = "QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

q2_shan_p2 <- ggplot(q2_shannon_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Shannon Diversity Index - QIIME2", 
       subtitle = "QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

nfc_q2_simp_p1 + q2_simp_p1
nfc_q2_shan_p2 + q2_shan_p2

# new dataset for indices
clin_data2 <- clin_data |>
  # The Name column is already cleaned (no underscore) from previous steps
  mutate(Sample_ID = Name)
levels(clin_data2$DietAB) <- c("Premium", "Conventional", "None")


# Create PCR ID column for merging
nfc_q2_simpson_df_merged <- nfc_q2_simpson_df |>
  tibble::rownames_to_column(var = "Sample_ID") # Move the Sample IDs (PCR#) from the row names into a column for merging
nfc_q2_shannon_df_merged <- nfc_q2_shannon_df |>
  tibble::rownames_to_column(var = "Sample_ID")
q2_simpson_df_merged <- q2_simpson_df |>
  tibble::rownames_to_column(var = "Sample_ID") # Move the Sample IDs (PCR#) from the row names into a column for merging
q2_shannon_df_merged <- q2_shannon_df |>
  tibble::rownames_to_column(var = "Sample_ID")


# Combining All T1 and T2 samples
nfc_q2_simp_Data <- clin_data2 |>
  left_join(nfc_q2_simpson_df_merged, by = "Sample_ID") |> # Join the clinical data and the diversity scores by the common Sample_ID (PCR#)
  na.omit() # Remove any rows where the diversity score was not found (e.g., sample dropped out)
nfc_q2_shan_Data <- clin_data2 |>
  left_join(nfc_q2_shannon_df_merged, by = "Sample_ID") |>
  na.omit()
q2_simp_Data <- clin_data2 |>
  left_join(q2_simpson_df_merged, by = "Sample_ID") |> # Join the clinical data and the diversity scores by the common Sample_ID (PCR#)
  na.omit() # Remove any rows where the diversity score was not found (e.g., sample dropped out)
q2_shan_Data <- clin_data2 |>
  left_join(q2_shannon_df_merged, by = "Sample_ID") |>
  na.omit()


# boxplot for all 3 categories
nfc_q2_simp_boxplot <- ggplot(nfc_q2_simp_Data, aes(x = DietAB, y = value, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Simpson Diversity Index by Diet Category",
    subtitle = "NF-Core QIIME2",
    x = "Diet Category",
    y = "Simpson Index (Higher = More Diverse)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Tilt labels to avoid overlap
  )

nfc_q2_shan_boxplot <- ggplot(nfc_q2_shan_Data, aes(x = DietAB, y = value, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Shannon Diversity Index by Diet Category",
    subtitle = "NF-Core QIIME2",
    x = "Diet Category",
    y = "Shannon Index (Higher = More Diverse)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

q2_simp_boxplot <- ggplot(q2_simp_Data, aes(x = DietAB, y = value, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Simpson Diversity Index by Diet Category",
    subtitle = "QIIME2",
    x = "Diet Category",
    y = "Simpson Index (Higher = More Diverse)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Tilt labels to avoid overlap
  )

q2_shan_boxplot <- ggplot(q2_shan_Data, aes(x = DietAB, y = value, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Shannon Diversity Index by Diet Category",
    subtitle = "QIIME2",
    x = "Diet Category",
    y = "Shannon Index (Higher = More Diverse)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
nfc_q2_simp_boxplot + q2_simp_boxplot
nfc_q2_shan_boxplot + q2_shan_boxplot


# bray & gower dissimilarity indices
nfc_q2_bray <- vegdist(nfc_abund_df_t, "bray")
nfc_q2_bray_df <- data.frame(value = nfc_q2_bray)

nfc_q2_gower <- vegdist(nfc_abund_df_t, "gower")
nfc_q2_gower_df <- data.frame(value = nfc_q2_gower)

q2_bray <- vegdist(qiime2_abund_df_t, "bray")
q2_bray_df <- data.frame(value = nfc_q2_bray)

q2_gower <- vegdist(qiime2_abund_df_t, "gower")
q2_gower_df <- data.frame(value = q2_gower)

nfc_q2_b_p3 <- ggplot(nfc_q2_bray_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Bray Dissimilarity Index", 
       subtitle = "NF-Core QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

nfc_q2_g_p4 <- ggplot(nfc_q2_gower_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Gower Dissimilarity Index", 
       subtitle = "NF-Core QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

q2_b_p3 <- ggplot(q2_bray_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Bray Dissimilarity Index", 
       subtitle = "QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

q2_g_p4 <- ggplot(q2_gower_df, aes(x = value)) +
  geom_histogram() +
  labs(title = "Gower Dissimilarity Index", 
       subtitle = "QIIME2", x = "Index Value", y = "Count") +
  theme_minimal()

nfc_q2_b_p3 + q2_b_p3
nfc_q2_g_p4 + q2_g_p4


# Avg Dissimilarity per sample (bray)
nfc_q2_bray_matrix <- as.matrix(nfc_q2_bray)

# Calculate the mean dissimilarity of each sample to all others
# (rowMeans is a vectorized way to get the average distance for each sample)
nfc_q2_bray_Avg_Dissim <- data.frame(
  Sample_ID = rownames(nfc_q2_bray_matrix),
  bray_dissimilarity = rowMeans(nfc_q2_bray_matrix)
)

# Avg Dissimilarity per sample (gower)
nfc_q2_gower_matrix <- as.matrix(nfc_q2_gower)
nfc_q2_gower_Avg_Dissim <- data.frame(
  Sample_ID = rownames(nfc_q2_gower_matrix),
  gower_dissimilarity = rowMeans(nfc_q2_gower_matrix)
)

q2_bray_matrix <- as.matrix(q2_bray)

# Calculate the mean dissimilarity of each sample to all others
# (rowMeans is a vectorized way to get the average distance for each sample)
q2_bray_Avg_Dissim <- data.frame(
  Sample_ID = rownames(q2_bray_matrix),
  bray_dissimilarity = rowMeans(q2_bray_matrix)
)

# Avg Dissimilarity per sample (gower)
q2_gower_matrix <- as.matrix(q2_gower)
q2_gower_Avg_Dissim <- data.frame(
  Sample_ID = rownames(q2_gower_matrix),
  gower_dissimilarity = rowMeans(q2_gower_matrix)
)


# merge (bray)
nfc_q2_bray_data <- clin_data2 |>
  left_join(nfc_q2_bray_Avg_Dissim, by = "Sample_ID") |>
  na.omit()

# merge (gower)
nfc_q2_gower_data <- clin_data2 |>
  left_join(nfc_q2_gower_Avg_Dissim, by = "Sample_ID") |>
  na.omit()

q2_bray_data <- clin_data2 |>
  left_join(q2_bray_Avg_Dissim, by = "Sample_ID") |>
  na.omit()

q2_gower_data <- clin_data2 |>
  left_join(q2_gower_Avg_Dissim, by = "Sample_ID") |>
  na.omit()


# boxplot (bray & gower)
nfc_q2_bray_boxplot <- ggplot(nfc_q2_bray_data, aes(x = DietAB, y = bray_dissimilarity, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Bray Dissimilarity by Diet Category",
    subtitle = "NF-core QIIME2",
    x = "Diet Category",
    y = "Bray Dissimilarity (Higher = Less Similar)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

nfc_q2_gower_boxplot <- ggplot(nfc_q2_gower_data, aes(x = DietAB, y = gower_dissimilarity, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Gower Dissimilarity by Diet Category",
    subtitle = "NF-core QIIME2",
    x = "Diet Category",
    y = "Gower Dissimilarity (Higher = Less Similar)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

q2_bray_boxplot <- ggplot(q2_bray_data, aes(x = DietAB, y = bray_dissimilarity, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Bray Dissimilarity by Diet Category",
    subtitle = "QIIME2",
    x = "Diet Category",
    y = "Bray Dissimilarity (Higher = Less Similar)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

q2_gower_boxplot <- ggplot(q2_gower_data, aes(x = DietAB, y = gower_dissimilarity, fill = DietAB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(
    title = "Gower Dissimilarity by Diet Category",
    subtitle = "QIIME2",
    x = "Diet Category",
    y = "Gower Dissimilarity (Higher = Less Similar)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Display the Beta Diversity plots side-by-side
nfc_q2_bray_boxplot + q2_bray_boxplot
nfc_q2_gower_boxplot + q2_gower_boxplot

# rarefaction
nfc_q2_spAbund <- rowSums(nfc_abund_df_t) # gives the number of individuals found in each plot
#nfc_q2_spAbund

q2_spAbund <- rowSums(qiime2_abund_df_t) # gives the number of individuals found in each plot
#q2_spAbund

nfc_q2_raremin <- min(rowSums(nfc_abund_df_t)) # rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
nfc_q2_raremin

q2_raremin <- min(rowSums(qiime2_abund_df_t)) # rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
q2_raremin

nfc_q2_sRare <- rarefy(nfc_abund_df_t, nfc_q2_raremin) # expected rarefied # of species
q2_sRare <- rarefy(qiime2_abund_df_t, q2_raremin) # expected rarefied # of species

nfc_q2_rcurve <- rarecurve(nfc_abund_df_t, step = 5000, col = "blue", tidy = TRUE)
q2_rcurve <- rarecurve(qiime2_abund_df_t, step = 5000, col = "blue", tidy = TRUE)

nfc_q2_rcurve_label <- nfc_q2_rcurve |>
  group_by(Site) |>
  filter(Sample == max(Sample))

q2_rcurve_label <- q2_rcurve |>
  group_by(Site) |>
  filter(Sample == max(Sample))

# https://www.youtube.com/watch?v=_OEdFjc1D9I&t=1s
ggplot(nfc_q2_rcurve, aes(x = Sample, y = Species, group = Site)) +
  geom_line(aes(color = Site),
            linewidth = 0.5) +
  geom_text(data = nfc_q2_rcurve_label, aes(label = Site, color = Site),
            hjust = 0, nudge_x = 500, size = 3) +
  labs(x = "Subsample Size (Individuals/Reads)", 
       y = "Rarefied Species Richness",
       title = "Rarefaction Curves by Sample",
       subtitle = "NF-Core QIIME2") +
  theme_minimal()

ggplot(q2_rcurve, aes(x = Sample, y = Species, group = Site)) +
  geom_line(aes(color = Site),
            linewidth = 0.5) +
  geom_text(data = q2_rcurve_label, aes(label = Site, color = Site),
            hjust = 0, nudge_x = 500, size = 3) +
  labs(x = "Subsample Size (Individuals/Reads)", 
       y = "Rarefied Species Richness",
       title = "Rarefaction Curves by Sample",
       subtitle = "QIIME2") +
  theme_minimal()
  


########
# stats
########
#install.packages("corrplot")
library(corrplot)

str(clin_data)

clin_data_sub <- clin_data |>
  select(Age:MGPI)

summary(clin_data_sub)

cor_mat <- cor(clin_data_sub, method = "pearson", use = "pairwise.complete.obs")
corrplot(cor_mat, type = 'lower', tl.col = 'black', tl.srt = 45)
title(main = 'Correlation matrix of clinical values')

# values greater than 0.9
high_cor <- abs(cor_mat) >= 0.9 # T/F matrix
high_val <- cor_mat[high_cor] # extraction
high_val

# show only high correlations
cor2 <- cor_mat
cor2[abs(cor2) < 0.9] <- NA
round(cor2, 4)

#drop BMC, vLDL, cLDL, Insuline, MGBI, MGT, MGPD, MGPI
table(clin_data$Group)


clin_align <- clin_data |>
  #rename(SampleID = "Name") |> 
  mutate(DietAB_Char = trimws(as.character(DietAB))) |>
  
  mutate(Diet_Group = factor(case_when(
    # relabel
    DietAB_Char == "a" ~ "ProbioticYogurt", 
    DietAB_Char == "b" ~ "ControlYogurt",
    TRUE ~ "NonYogurt"
  ))) |>
  filter(Time == "T2") |>
  # This makes the PCR IDs (from the 'Name' column) the row names.
  tibble::column_to_rownames(var = "Name")

colnames(clin_align)


clin_pred <- clin_align |>
  select(
    Sample_name, CP, Time, TimeAB, Diet, DietAB, Diet_Group, Gender, Age, BMC,
    HbA1c, HOMA_IR, TG, cHDL, cLDL, PCR_hs, MG, MM, MGBD
  )




# align datasets

# ID common samples
nfc_common <- intersect(rownames(nfc_abund_df_t), rownames(clin_pred))
q2_common <- intersect(rownames(qiime2_abund_df_t), rownames(clin_pred))

# subset and reorder

# features - NF Core & QIIME2: subset using the common sample
final_ncf <- nfc_abund_df_t[nfc_common, ]
final_q2 <- qiime2_abund_df_t[q2_common, ]

# clinical data - subset using the common sample
final_clin_nfc <- clin_pred[nfc_common, ]
final_clin_q2 <- clin_pred[q2_common, ]

summary(final_clin_nfc)
summary(final_clin_q2)



# CCA model
nfc_cca <- cca(final_ncf ~ Diet_Group + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD,
               data = final_clin_nfc)
q2_cca <- cca(final_q2 ~ Diet_Group + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD,
               data = final_clin_q2)

summary(nfc_cca)
summary(q2_cca)


plot(nfc_cca, 
     type = "t", 
     display = c("bp", "cn"), 
     main = "CCA Biplot: Microbial Community Constrained by Clinical Factors - NF-Core",
     scaling = 1)

plot(q2_cca, 
     type = "t", 
     display = c("bp", "cn"), 
     main = "CCA Biplot: Microbial Community Constrained by Clinical Factors - QIIME2",
     scaling = 1)


anova.cca(nfc_cca, permutations = 999)
anova.cca(nfc_cca, by = "term", permutations = 999)
anova.cca(q2_cca, permutations = 999)
anova.cca(q2_cca, by = "term", permutations = 999)


maas_mod_nfc <- maaslin3(
  input_data = final_ncf,
  input_metadata = final_clin_nfc,
  output = "maas_ncf",
  formula = "~Diet_Group + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD +
    Age + Gender",
  normalization = 'TSS',
  transform = 'LOG',
)
summary(maas_mod_nfc)

maas_mod_q2 <- maaslin3(
  input_data = final_q2,
  input_metadata = final_clin_q2,
  output = "maas_q2",
  formula = "~Diet_Group + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD +
    Age + Gender",
  normalization = 'TSS',
  transform = 'LOG',
)
summary(maas_mod_nfc)


maas_mod_q2

maas_mod_nfc_group <- maaslin3(
  input_data = final_ncf,
  input_metadata = final_clin_nfc,
  output = "maas_ncf_group",
  formula = "~group(Diet_Group) + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD +
    Age + Gender",
  normalization = 'TSS',
  transform = 'LOG',
)

maas_mod_q2_group <- maaslin3(
  input_data = final_q2,
  input_metadata = final_clin_q2,
  output = "maas_q2_group",
  formula = "~group(Diet_Group) + BMC + HbA1c + HOMA_IR + TG + cHDL + cLDL + PCR_hs + MG + MM + MGBD +
    Age + Gender",
  normalization = 'TSS',
  transform = 'LOG',
)
summary(maas_mod_nfc2)
maas_mod2$fit_data_prevalence
maas_mod2$fit_data_abundance



# different options to view the ""variables" from the model
# its the same as going to the output folder if you dont want to leave R
# maas_mod_q2$data
# maas_mod_q2$normalized_data
# maas_mod_q2$filtered_data
# maas_mod_q2$transformed_data
# maas_mod_q2$metadata
# maas_mod_q2$standardized_metadata
# maas_mod_q2$formula
# maas_mod_q2$fit_data_abundance
# maas_mod_q2$fit_data_prevalence








# diversity per group, intra group differences @ intervention
#random effect

# this would've made life so much easier
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# 
# library(phyloseq)