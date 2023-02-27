# Load the readr and tools packages
library(readr)
library(tools)
library(dplyr)

# set number of digits
options(digits = 22)

# Get a list of all csv files in the directory
csv_files <- list.files("./DGE", pattern = "*.csv", full.names = TRUE)

# Loop over each csv file in the list and apply the code
for (csv_file in csv_files) {
  # Read in the csv file
  data <- read_csv(csv_file)

    
  # Mutate logFC and direction
    colnames(data) <- c("ID","Ensembl", "Symbol", "logFC", "logCPM", "LR", "PValue", "FDR", "Status", "Description")
    data$logFC = -1*data$logFC 
    data$Status = case_when(data$Status == "UP" ~ "DOWN",
                            data$Status == "DOWN" ~ "UP")
    
  # Extract the file name
  file_name <- basename(csv_file)
  file_name_no_ext <- file_path_sans_ext(file_name)
  
  # Split the file name into five parts
  name_parts <- strsplit(file_name_no_ext, "_|-")
  
  numerator <- paste(name_parts[[1]][3], name_parts[[1]][4], name_parts[[1]][5], sep = "_")
  denominator <- paste(name_parts[[1]][6], name_parts[[1]][7], name_parts[[1]][8], sep = "_")
  trunk_name <- paste(denominator, numerator, sep = "-")
  
  # Swap the two strings and create the new file name
  new_name <- paste(name_parts[[1]][1], "_", name_parts[[1]][2], "_", trunk_name, ".csv", sep = "")
  
  col4 <- paste(trunk_name, "logFC", sep = "_")
  col5 <- paste(trunk_name, "logCPM", sep = "_")
  col6 <- paste(trunk_name, "LR", sep = "_")
  col7 <- paste(trunk_name, "PValue", sep = "_")
  col8 <- paste(trunk_name, "FDR", sep = "_")
  col9 <- paste(trunk_name, "Status", sep = "_")
    
  # Change the header names
  colnames(data) <- c("","Ensembl", "Symbol", col4, col5, col6, col7, col8, col9, "Description")

  # Save the data as a csv file with the new name
  write_csv(data, paste("Corrected_DGE/", new_name, sep = ""))
}
