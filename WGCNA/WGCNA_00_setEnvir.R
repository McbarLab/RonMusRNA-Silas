# This analysis is coded by Di "Silas" Kuang
# Cited from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Email: dkuang5@wisc.edu
# RStudio version: 2023.03.0 Build 386
# R version: 4.3.0

currentwd <- "."

if (getwd() != currentwd) {
  # Set the working directory to "./WGCNA/"
  setwd("./WGCNA/")
}

currentwd <- getwd()

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