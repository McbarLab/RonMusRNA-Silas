# This analysis is coded by Di "Silas" Kuang
# Cited from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0

workingdir <- "C:/Users/skuang/Documents/RonMusRNA-Silas/WGCNA"

if (getwd() != workingdir) {
  # Set the working directory to "./WGCNA/"
  setwd(workingdir)
}

currentwd <- getwd()
options(stringsAsFactors = FALSE)

#install BiocManager dependent packages
library(BiocManager)
if (!require("WGCNA", quietly = TRUE))
  BiocManager::install("WGCNA")
if (!require("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")
if (!require("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
if (!require("impute", quietly = TRUE))
  BiocManager::install("impute")
if(!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if(!require("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")
if(!require("EnhancedVolcano", quietly = TRUE))
  BiocManager::install('EnhancedVolcano')
if(!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")


# Load packages
library(tidyverse)
library(WGCNA)
library(matrixStats)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(survival)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)