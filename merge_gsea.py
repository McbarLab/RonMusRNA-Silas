import pandas as pd
import glob
import os

# Path to the directory containing CSV files
csv_files_path = './GSEA_Pathway_csv/'

# List all CSV files in the directory
csv_files = glob.glob(csv_files_path + '*.csv')

# Check if there are any CSV files in the directory
if not csv_files:
    print("No CSV files found in the specified directory.")
else:
    # Create a dictionary to store DataFrames, with modified file names as keys
    data_frames = {}

    # Loop through each CSV file and read it into a DataFrame
    for csv_file in csv_files:
        # Extract the file name without extension
        file_name = os.path.basename(csv_file).split('.')[0]

        # Remove the " gseaKEGG_all" from the file name
        sheet_name = file_name.replace(" gseaKEGG_all", "")

        # Read CSV file into a DataFrame
        df = pd.read_csv(csv_file)

        # Add the DataFrame to the dictionary with the modified file name as the key
        data_frames[sheet_name] = df

    # Save the dictionary of DataFrames to an Excel file with each DataFrame as a separate tab
    with pd.ExcelWriter('master_GSEA_20240306.xlsx', engine='xlsxwriter') as writer:
        for sheet_name, df in data_frames.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
