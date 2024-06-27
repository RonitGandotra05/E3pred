import os
import pandas as pd
import subprocess
import joblib
import requests
import xml.etree.ElementTree as ET
import csv

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
        'label': labels,
        'UniProt Link': [f"https://www.uniprot.org/uniprot/{seq[1]}" for seq in sequences]
    })
    output_df.to_csv(output_file, index=False)
    print(f"Output saved to '{output_file}'")

# Function to call UbiBrowser API
def call_ubibrowser_api(url):
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.content
        else:
            print(f"Failed to fetch data from URL: {url}")
            print(f"HTTP Status Code: {response.status_code}")
            print(f"Response Text: {response.text}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None

def parse_xml_response(xml_content):
    try:
        root = ET.fromstring(xml_content)
        interactions = root.findall('.//Interaction')
        
        data = []
        for interaction in interactions:
            enzyme_uniprot_AC = interaction.find('enzyme_uniprot_AC').text if interaction.find('enzyme_uniprot_AC') is not None else ''
            enzyme_gene_name = interaction.find('enzyme_gene_name').text if interaction.find('enzyme_gene_name') is not None else ''
            substrate_uniprot_AC = interaction.find('substrate_uniprot_AC').text if interaction.find('substrate_uniprot_AC') is not None else ''
            substrate_gene_name = interaction.find('substrate_gene_name').text if interaction.find('substrate_gene_name') is not None else ''
            species = interaction.find('species').text if interaction.find('species') is not None else ''
            confidence_score = interaction.find('confidence_score').text if interaction.find('confidence_score') is not None else ''
            likelihood_ratio = interaction.find('likelihood_ratio').text if interaction.find('likelihood_ratio') is not None else ''
            rank = interaction.find('rank').text if interaction.find('rank') is not None else ''
            p_value = interaction.find('P-value').text if interaction.find('P-value') is not None else ''
            level = interaction.find('Level').text if interaction.find('Level') is not None else ''
            enzyme_type = interaction.find('enzyme_type').text if interaction.find('enzyme_type') is not None else ''

            domain_lr = interaction.find('.//DOMAIN_LR').text if interaction.find('.//DOMAIN_LR') is not None else ''
            domain_enzyme = interaction.find('.//domain_enzyme').text if interaction.find('.//domain_enzyme') is not None else ''
            domain_substrate = interaction.find('.//domain_substrate').text if interaction.find('.//domain_substrate') is not None else ''

            go_lr = interaction.find('.//GO_LR').text if interaction.find('.//GO_LR') is not None else ''
            go_enzyme = interaction.find('.//GO_enzyme').text if interaction.find('.//GO_enzyme') is not None else ''
            go_substrate = interaction.find('.//GO_substrate').text if interaction.find('.//GO_substrate') is not None else ''

            net_lr = interaction.find('.//NET_LR').text if interaction.find('.//NET_LR') is not None else ''
            net_number = interaction.find('.//net_number').text if interaction.find('.//net_number') is not None else ''
            net_score = interaction.find('.//net_score').text if interaction.find('.//net_score') is not None else ''

            motif_lr = interaction.find('.//MOTIF_LR').text if interaction.find('.//MOTIF_LR') is not None else ''
            motif = interaction.find('.//motif').text if interaction.find('.//motif') is not None else ''
            motif_others_elements = interaction.findall('.//motif_others')
            motif_others = '; '.join([elem.text if elem is not None else '' for elem in motif_others_elements])

            data.append({
                'enzyme_uniprot_AC': enzyme_uniprot_AC,
                'enzyme_gene_name': enzyme_gene_name,
                'substrate_uniprot_AC': substrate_uniprot_AC,
                'substrate_gene_name': substrate_gene_name,
                'species': species,
                'confidence_score': confidence_score,
                'likelihood_ratio': likelihood_ratio,
                'rank': rank,
                'P-value': p_value,
                'Level': level,
                'enzyme_type': enzyme_type,
                'domain_lr': domain_lr,
                'domain_enzyme': domain_enzyme,
                'domain_substrate': domain_substrate,
                'go_lr': go_lr,
                'go_enzyme': go_enzyme,
                'go_substrate': go_substrate,
                'net_lr': net_lr,
                'net_number': net_number,
                'net_score': net_score,
                'motif_lr': motif_lr,
                'motif': motif,
                'motif_others': motif_others
            })
        
        return data
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return []

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
    # Assuming sequences and y_pred are defined earlier in your script
    
    # Separate URL templates for ubiquitination and deubiquitination
    ubiquitination_url_template = "http://ubibrowser.bio-it.cn/v2/home/API?process=ubiquitination&type=enzyme&term={}"
    deubiquitination_url_template = "http://ubibrowser.bio-it.cn/v2/home/API?process=deubiquitination&type=enzyme&term={}"

    # Initialize lists to store parsed data
    ubiquitination_data = []
    deubiquitination_data = []

    # Initialize a flag to check if any entries with label == 1 were processed
    processed_label_1 = False

    for sequence, label in zip(sequences, labels):
        entry_name = sequence[2]  # Assuming sequence[2] contains the entry name

        if label == 1:
            processed_label_1 = True  # Flag that at least one label == 1 entry was processed

            # Call UbiBrowser API for ubiquitination data
            ubiquitination_url = ubiquitination_url_template.format(entry_name)
            print(f"Calling UbiBrowser API for ubiquitination data: {entry_name}")
            ubiquitination_xml_content = call_ubibrowser_api(ubiquitination_url)
            if ubiquitination_xml_content:
                ubiquitination_parsed_data = parse_xml_response(ubiquitination_xml_content)
                if ubiquitination_parsed_data:
                    ubiquitination_data.extend(ubiquitination_parsed_data)
                else:
                    print(f"No ubiquitination data parsed for entry '{entry_name}'")
            else:
                print(f"Failed to fetch ubiquitination data for entry '{entry_name}'")

            # Call UbiBrowser API for deubiquitination data
           # deubiquitination_url = deubiquitination_url_template.format(entry_name)
           # print(f"Calling UbiBrowser API for deubiquitination data: {entry_name}")
           # deubiquitination_xml_content = call_ubibrowser_api(deubiquitination_url)
           # if deubiquitination_xml_content:
            #    deubiquitination_parsed_data = parse_xml_response(deubiquitination_xml_content)
             #   if deubiquitination_parsed_data:
              #      deubiquitination_data.extend(deubiquitination_parsed_data)
               # else:
                #    print(f"No deubiquitination data parsed for entry '{entry_name}'")
           # else:
            #    print(f"Failed to fetch deubiquitination data for entry '{entry_name}'")

        elif label == 0:
            # Skip processing for label == 0 entries
            continue

        else:
            # Handle unexpected labels if needed
            print(f"Unexpected label value: {label}")

    # After processing all entries, print messages
    if processed_label_1:
        print("Pipeline execution completed.")
    else:
        print("No entries with label == 1 found or all were skipped.")
        
   # if not ubiquitination_data and not deubiquitination_data:
        #print("No data found for ubiquitination and deubiquitination.")
    if not ubiquitination_data:
        print("No data found for ubiquitination.")
   # elif not deubiquitination_data:
        #print("No data found for deubiquitination.")

    # Write ubiquitination data to CSV
    ubiquitination_csv_file = f"{base_filename}_ubiquitination_data.csv"
    with open(ubiquitination_csv_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=ubiquitination_data[0].keys())
        writer.writeheader()
        writer.writerows(ubiquitination_data)

    print(f"Ubiquitination data written to {ubiquitination_csv_file}")

    # Write deubiquitination data to CSV
   # deubiquitination_csv_file = f"{base_filename}_deubiquitination_data.csv"
  #  with open(deubiquitination_csv_file, 'w', newline='', encoding='utf-8') as f:
    #    writer = csv.DictWriter(f, fieldnames=deubiquitination_data[0].keys())
     #   writer.writeheader()
      #  writer.writerows(deubiquitination_data)

 #   print(f"Deubiquitination data written to {deubiquitination_csv_file}")

    print("Pipeline execution completed.")

if __name__ == "__main__":
    main()
