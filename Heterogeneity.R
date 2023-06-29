library(tidyverse)
library(dplyr)
library(Biobase)
library(ggplot2)
library(ggpubr)

TPM_data = read.csv("./WGCNA/rsubread_GENE_tpm.csv", row.names = 1)


# 06 Ctrl Male only
c06_M_C_tpm.df = TPM_data[,1:3]
#c06_M_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c06_M_C_tpm.df = c06_M_C_tpm.df[rowSums(c06_M_C_tpm.df[])>0,]

# New column with:
c06_M_C_tpm.df$MeanTPM = apply(c06_M_C_tpm.df, 1, mean) # rowmean
c06_M_C_tpm.df$SdTPM = apply(c06_M_C_tpm.df, 1, sd) # SD
c06_M_C_tpm.df$CvTPM = c06_M_C_tpm.df$SdTPM / c06_M_C_tpm.df$MeanTPM # CV as CV = SD/mean

c06_M_C_tpm.ggScat = ggplot(c06_M_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) + 
  geom_point(size = 0.2) + theme_bw()


# 22 Ctrl Male only
c22_M_C_tpm.df = TPM_data[,c(13,15,17)]
#c22_M_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c22_M_C_tpm.df = c22_M_C_tpm.df[rowSums(c22_M_C_tpm.df[])>0,]

# New column with:
c22_M_C_tpm.df$MeanTPM = apply(c22_M_C_tpm.df, 1, mean) # rowmean
c22_M_C_tpm.df$SdTPM = apply(c22_M_C_tpm.df, 1, sd) # SD
c22_M_C_tpm.df$CvTPM = c22_M_C_tpm.df$SdTPM / c22_M_C_tpm.df$MeanTPM # CV as CV = SD/mean

c22_M_C_tpm.ggScat = ggplot(c22_M_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) +
  geom_point(size = 0.2) + theme_bw()


# 28 Ctrl Male only
c28_M_C_tpm.df = TPM_data[,c(25,27,29)]
#c28_M_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c28_M_C_tpm.df = c28_M_C_tpm.df[rowSums(c28_M_C_tpm.df[])>0,]

# New column with:
c28_M_C_tpm.df$MeanTPM = apply(c28_M_C_tpm.df, 1, mean) # rowmean
c28_M_C_tpm.df$SdTPM = apply(c28_M_C_tpm.df, 1, sd) # SD
c28_M_C_tpm.df$CvTPM = c28_M_C_tpm.df$SdTPM / c28_M_C_tpm.df$MeanTPM # CV as CV = SD/mean

c28_M_C_tpm.ggScat = ggplot(c28_M_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) +
  geom_point(size = 0.2) + theme_bw()

################################

# 06 Ctrl Female only
c06_F_C_tpm.df = TPM_data[,7:9]
#c06_F_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c06_F_C_tpm.df = c06_F_C_tpm.df[rowSums(c06_F_C_tpm.df[])>0,]

# New column with:
c06_F_C_tpm.df$MeanTPM = apply(c06_F_C_tpm.df, 1, mean) # rowmean
c06_F_C_tpm.df$SdTPM = apply(c06_F_C_tpm.df, 1, sd) # SD
c06_F_C_tpm.df$CvTPM = c06_F_C_tpm.df$SdTPM / c06_F_C_tpm.df$MeanTPM # CV as CV = SD/mean

c06_F_C_tpm.ggScat = ggplot(c06_F_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) +
  geom_point(size = 0.2) + theme_bw()


# 22 Ctrl Female only
c22_F_C_tpm.df = TPM_data[,c(19,21,23)]
#c22_F_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c22_F_C_tpm.df = c22_F_C_tpm.df[rowSums(c22_F_C_tpm.df[])>0,]

# New column with:
c22_F_C_tpm.df$MeanTPM = apply(c22_F_C_tpm.df, 1, mean) # rowmean
c22_F_C_tpm.df$SdTPM = apply(c22_F_C_tpm.df, 1, sd) # SD
c22_F_C_tpm.df$CvTPM = c22_F_C_tpm.df$SdTPM / c22_F_C_tpm.df$MeanTPM # CV as CV = SD/mean

c22_F_C_tpm.ggScat = ggplot(c22_F_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) +
  geom_point(size = 0.2) + theme_bw()


# 28 Ctrl Female only
c28_F_C_tpm.df = TPM_data[,c(32,34,35)]
#c28_F_C_tpm.mt = as.matrix(TPM_data[,1:3])

# Remove rows with only zeros
c28_F_C_tpm.df = c28_F_C_tpm.df[rowSums(c28_F_C_tpm.df[])>0,]

# New column with:
c28_F_C_tpm.df$MeanTPM = apply(c28_F_C_tpm.df, 1, mean) # rowmean
c28_F_C_tpm.df$SdTPM = apply(c28_F_C_tpm.df, 1, sd) # SD
c28_F_C_tpm.df$CvTPM = c28_F_C_tpm.df$SdTPM / c28_F_C_tpm.df$MeanTPM # CV as CV = SD/mean

c28_F_C_tpm.ggScat = ggplot(c28_F_C_tpm.df, aes(x = log2(MeanTPM), y = CvTPM)) +
  geom_point(size = 0.2) + theme_bw()

# Combine figures
Cv_vs_log2TPM_ScatPlots = ggarrange(c06_M_C_tpm.ggScat, c06_F_C_tpm.ggScat,
                                     c22_M_C_tpm.ggScat, c22_F_C_tpm.ggScat,
                                     c28_M_C_tpm.ggScat, c28_F_C_tpm.ggScat,
                                     nrow = 3, ncol = 2)
Cv_vs_log2TPM_ScatPlots



max(c06_M_C_tpm.df$CvTPM)
max(c22_M_C_tpm.df$CvTPM)
max(c28_M_C_tpm.df$CvTPM)
max(c06_F_C_tpm.df$CvTPM)
max(c22_F_C_tpm.df$CvTPM)
max(c28_F_C_tpm.df$CvTPM)

bin_size = 10

c06_M_C_tpm.df %>%
  mutate(bin_CvTPM = factor(CvTPM%/%bin_size*10)) %>%
  ggplot(aes(x = bin_CvTPM, y = log2(MeanTPM))) +
  geom_boxplot()

# Boxplots (3 bins)
c06_M_C_tpm.df$bin = cut(c06_M_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c06_M_C_tpm.ggBox = ggplot(c06_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c22_M_C_tpm.df$bin = cut(c22_M_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c22_M_C_tpm.ggBox = ggplot(c22_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c28_M_C_tpm.df$bin = cut(c28_M_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c28_M_C_tpm.ggBox = ggplot(c28_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c06_F_C_tpm.df$bin = cut(c06_F_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c06_F_C_tpm.ggBox = ggplot(c06_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c22_F_C_tpm.df$bin = cut(c22_F_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c22_F_C_tpm.ggBox = ggplot(c22_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c28_F_C_tpm.df$bin = cut(c28_F_C_tpm.df$CvTPM, c(0, .5, 1, 1.5))
c28_F_C_tpm.ggBox = ggplot(c28_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

Cv_vs_log2TPM_BoxPlots = ggarrange(c06_M_C_tpm.ggBox, c06_F_C_tpm.ggBox,
                                    c22_M_C_tpm.ggBox, c22_F_C_tpm.ggBox,
                                    c28_M_C_tpm.ggBox, c28_F_C_tpm.ggBox,
                                    nrow = 3, ncol = 2)
Cv_vs_log2TPM_BoxPlots

# Boxplots (6 bins)
c06_M_C_tpm.df$bin = cut(c06_M_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c06_M_C_tpm.ggBox = ggplot(c06_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c22_M_C_tpm.df$bin = cut(c22_M_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c22_M_C_tpm.ggBox = ggplot(c22_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c28_M_C_tpm.df$bin = cut(c28_M_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c28_M_C_tpm.ggBox = ggplot(c28_M_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c06_F_C_tpm.df$bin = cut(c06_F_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c06_F_C_tpm.ggBox = ggplot(c06_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c22_F_C_tpm.df$bin = cut(c22_F_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c22_F_C_tpm.ggBox = ggplot(c22_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

c28_F_C_tpm.df$bin = cut(c28_F_C_tpm.df$CvTPM, c(0, .25, .5, .75, 1, 1.25, 1.5))
c28_F_C_tpm.ggBox = ggplot(c28_F_C_tpm.df) + geom_boxplot(aes(bin, log2(MeanTPM))) + ylim(-15, 17)

Cv_vs_log2TPM_BoxPlots = ggarrange(c06_M_C_tpm.ggBox, c06_F_C_tpm.ggBox,
                                   c22_M_C_tpm.ggBox, c22_F_C_tpm.ggBox,
                                   c28_M_C_tpm.ggBox, c28_F_C_tpm.ggBox,
                                   nrow = 3, ncol = 2)
Cv_vs_log2TPM_BoxPlots



cars$bin <- cut(cars$dist, c(1, 10, 30, 50, 200))
ggplot(cars) + geom_boxplot(aes(bin, speed))

table(c28_M_C_tpm.df$bin)
cars %>% 
  mutate(bin_dist = factor(dist%/%bin_size*10)) %>% 
  ggplot(aes(x = bin_dist, y = speed)) +
  geom_boxplot()


# CV Analysis starting here
setName <- function(df, CVname) {
  curated_df <- data.frame(GENE = rownames(df), CvTPM = df[,"CvTPM"])
  names(curated_df)[names(curated_df) == "CvTPM"] <- CVname
  return(curated_df)
}

curated_06FC <- setName(c06_F_C_tpm.df, "CV_06")
curated_22FC <- setName(c22_F_C_tpm.df, "CV_22")
curated_28FC <- setName(c28_F_C_tpm.df, "CV_28")

curated_06MC <- setName(c06_M_C_tpm.df, "CV_06")
curated_22MC <- setName(c22_M_C_tpm.df, "CV_22")
curated_28MC <- setName(c28_M_C_tpm.df, "CV_28")

# Create a df which contains only the CV
CV_F_C <- curated_06FC %>%
  left_join(curated_22FC, by = "GENE") %>%
  left_join(curated_28FC, by = "GENE") %>%
  na.omit()

CV_M_C <- curated_06MC %>%
  left_join(curated_22MC, by = "GENE") %>%
  left_join(curated_28MC, by = "GENE") %>%
  na.omit()

# Compute the absolute change of CV along lifespan
CV_Change <- function(df) {
  change <- df %>%
    mutate(CV_ABS = abs(CV_22 - CV_06) + abs(CV_28 - CV_22)) %>%
    mutate(Trend = case_when(
      (CV_22 - CV_06)>=0 & (CV_28 - CV_22)>=0 ~ "/",
      (CV_22 - CV_06)<0 & (CV_28 - CV_22)<0 ~ "\\",
      (CV_22 - CV_06)>=0 & (CV_28 - CV_22)<0 ~ "A",
      (CV_22 - CV_06)<0 & (CV_28 - CV_22)>=0  ~ "V"
    )) %>% 
    arrange(desc(CV_ABS))
  return(change)
}

CV_FC_ABS <- CV_Change(CV_F_C)
CV_MC_ABS <- CV_Change(CV_M_C)