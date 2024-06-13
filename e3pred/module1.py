import os
import pandas as pd
import subprocess
import joblib

# Function to prompt for the input FASTA filename
def get_input_file():
    fasta_file = input("Enter the input FASTA filename (e.g., 'XYZ.fasta'): ")
    base_filename = os.path.splitext(fasta_file)[0]
    return fasta_file, base_filename

# Function to run pfeature_comp.py command
def run_pfeature_comp_command(fasta_file, base_filename, feature):
    command = f"python3 -m e3pred.pfeature_comp -i {fasta_file} -o {base_filename}_{feature}.csv -j {feature}"
    subprocess.run(command, shell=True)

# Columns to extract after combining
columns_to_extract = [
    'BTC_T', 'BTC_S', 'BTC_H', 'BTC_D', 'SOC1_G1', 
    'AAI_NADH010107', 'AAI_JUNJ780101', 'AAI_BIOV880102', 
    'AAI_NADH010102', 'AAI_BIOV880101', 'AAI_RADA880106', 
    'APAAC1_D', 'PAAC1_D', 'AAC_D', 'AAI_NADH010103', 
    'AAI_HUTJ700103', 'AAI_KHAG800101', 'AAI_NADH010101', 
    'AAI_NADH010104', 'AAI_FASG760103', 'AAI_CHOC760101', 
    'AAI_ZIMJ680103', 'AAI_DAYM780201', 'PAAC1_K', 'AAC_K', 
    'APAAC1_K', 'AAI_LEVM760103', 'AAI_GRAR740103', 'APAAC1_S', 'AAC_S'
]

# Function to combine CSV files horizontally
def combine_csv_files_horizontally(input_files):
    combined_df = pd.read_csv(input_files[0])
    for file in input_files[1:]:
        df = pd.read_csv(file)
        combined_df = pd.concat([combined_df, df], axis=1)
    print("All files combined")
    return combined_df

# Function to extract specified columns from a DataFrame
def extract_columns(dataframe, columns, output_file):
    extracted_data = dataframe[columns]
    extracted_data.to_csv(output_file, index=False)
    print(f"Extracted columns saved to '{output_file}'")

# Function to load a specific model
def load_model(model_filename):
    model_path = os.path.join(os.path.dirname(__file__), 'model', model_filename)
    return joblib.load(model_path)

# Function to prompt the user to select a model
def select_model(model_filenames):
    print("Select a model to run:")
    for i, model_filename in enumerate(model_filenames, 1):
        print(f"{i}. {model_filename}")
    while True:
        choice = input("Enter the number of the model you want to use: ")
        try:
            choice_index = int(choice) - 1
            if 0 <= choice_index < len(model_filenames):
                return model_filenames[choice_index]
            else:
                print("Invalid choice. Please enter a valid number.")
        except ValueError:
            print("Invalid input. Please enter a number.")

# Function to load models and make predictions on the test data
def run_models_on_test_data(test_data_path, selected_model, scaler_filename):
    test_data = pd.read_csv(test_data_path)
    X_test = test_data
    scaler_path = os.path.join(os.path.dirname(__file__), 'model', scaler_filename)
    scaler = joblib.load(scaler_path)
    X_test_scaled = scaler.transform(X_test)
    model = load_model(selected_model)
    y_pred = model.predict(X_test_scaled)
    model_name = selected_model.replace('_model.pkl', '').replace('_', ' ').title()
    return y_pred, model_name

# Function to parse the FASTA file and extract sequence information
def parse_fasta(fasta_file):
    sequences = []
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    # Initialize variables to hold the current header and sequence
    header = None
    sequence = []

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Process the previous header and sequence if they exist
            if header and sequence:
                sequences.append(parse_header_and_sequence(header, ''.join(sequence)))

            # Start a new header and reset the sequence
            header = line
            sequence = []
        else:
            # Accumulate the sequence lines
            sequence.append(line)

    # Process the last header and sequence
    if header and sequence:
        sequences.append(parse_header_and_sequence(header, ''.join(sequence)))

    return sequences

def parse_header_and_sequence(header, sequence):
    if header.startswith('>sp|') or header.startswith('>tr|'):
        parts = header.split("|")
        sp = parts[0][1:]  # Remove '>' character
        uniprot_id = parts[1]
        entry_name = parts[2].split()[0]
        # Capture description up to 'OS='
        if 'OS=' in parts[2]:
            description = parts[2].split('OS=')[0].strip()
        else:
            description = parts[2].strip()
        os = "Unknown"
        ox = "Unknown"
        pe = "Unknown"
        sv = "Unknown"
        gene_name = "Unknown"
        if "OS=" in header:
            os = header.split("OS=")[1].split(" OX=")[0].strip()
        if "OX=" in header:
            ox = header.split("OX=")[1].split(" ")[0]
        if "PE=" in header:
            pe = header.split("PE=")[1].split(" ")[0]
        if "SV=" in header:
            sv = header.split("SV=")[1].split(" ")[0]
        if "GN=" in header:
            gene_name = header.split("GN=")[1].split(" ")[0]
    elif header.startswith('>gi|'):
        parts = header.split("|")
        sp = parts[0][1:]  # Remove '>' character
        uniprot_id = parts[1]
        entry_name = parts[2].split()[0]
        # Capture description up to '['
        if '[' in parts[2]:
            description = parts[2].split('[')[0].strip()
        else:
            description = parts[2].strip()
        os = "Unknown"
        ox = "Unknown"
        pe = "Unknown"
        sv = "Unknown"
        gene_name = "Unknown"
        if "[" in header:
            os = header.split("[")[1].split("]")[0]
        if "OX=" in header:
            ox = header.split("OX=")[1].split(" ")[0]
    else:
        sp = "Unknown"
        uniprot_id = "Unknown"
        entry_name = "Unknown"
        description = header[1:]  # Remove '>' character
        os = "Unknown"
        ox = "Unknown"
        pe = "Unknown"
        sv = "Unknown"
        gene_name = "Unknown"

    sequence_length = len(sequence)
    return (sp, uniprot_id, entry_name, description, os, ox, pe, sv, gene_name, sequence_length, sequence)

def create_output_csv(sequences, labels, output_file):
    output_df = pd.DataFrame({
        'UniProt ID': [seq[1] for seq in sequences],
        'Entry Name': [seq[2] for seq in sequences],
        'Description': [seq[3] for seq in sequences],
        'Species': [seq[4] for seq in sequences],
        'NCBI Taxonomic ID': [seq[5] for seq in sequences],
        'Protein Existence Evidence': [seq[6] for seq in sequences],
        'Sequence Version No.': [seq[7] for seq in sequences],
        'Gene Name': [seq[8] for seq in sequences],
        'Sequence Length': [seq[9] for seq in sequences],
        'Predicted Label': labels,
        'UniProt Link': [f"https://www.uniprot.org/uniprot/{seq[1]}" for seq in sequences]
    })
    output_df.to_csv(output_file, index=False)
    print(f"Output saved to '{output_file}'")

# Main function to execute the whole process
def main():
    fasta_file, base_filename = get_input_file()
    sequences = parse_fasta(fasta_file)
    for feature in ["aac", "btc", "paac", "apaac", "soc", "aai"]:
        run_pfeature_comp_command(fasta_file, base_filename, feature)
    input_files = [f"{base_filename}_{feature}.csv" for feature in ["aac", "aai", "apaac", "btc", "paac", "soc"]]
    combined_df = combine_csv_files_horizontally(input_files)
    output_file = f"{base_filename}_extracted_features.csv"
    extract_columns(combined_df, columns_to_extract, output_file)
    for file in input_files:
        os.remove(file)
    print("Intermediate files have been deleted.")
    model_filenames = [
        "logistic_regression_model.pkl",
        "random_forest_model.pkl",
        "gradient_boosting_model.pkl",
        "adaboost_model.pkl",
        "bagging_model.pkl",
        "extra_trees_model.pkl",
        "hist_gradient_boosting_model.pkl",
        "support_vector_machine_model.pkl",
        "mlp_model.pkl",
        "xgboost_model.pkl"
    ]
    selected_model = select_model(model_filenames)
    labels, model_name = run_models_on_test_data(output_file, selected_model, 'scaler.pkl')
    output_csv_file = f"{base_filename}_predictions.csv"
    create_output_csv(sequences, labels, output_csv_file)

if __name__ == "__main__":
    main()

