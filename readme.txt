README.txt
==========

e3pred - A Protein Sequence Analysis Pipeline
=============================================

Overview
--------
This Python package provides a comprehensive pipeline for analyzing protein sequences from FASTA files, extracting features, making predictions using pre-trained machine learning models, and fetching protein interaction data from the UbiBrowser API. The pipeline is designed to automate various bioinformatics tasks and streamline the process of protein sequence analysis.

Features
--------
- Parses input FASTA files and extracts sequence information.
- Extracts various features using the `pfeature_comp.py` script from the `e3pred` module.
- Combines extracted features into a single CSV file.
- Selects relevant columns from the combined feature set.
- Allows user to choose from multiple pre-trained machine learning models for predictions.
- Makes predictions on the extracted features and assigns labels to sequences.
- Interacts with the UbiBrowser API to fetch protein interaction data for sequences predicted to be ubiquitinated.
- Generates output CSV files with sequence information, predicted labels, and interaction data.

Requirements
------------
- Python 3.x
- pandas
- joblib
- requests
- xml.etree.ElementTree (standard library)
- subprocess (standard library)
- os (standard library)
- csv (standard library)

Installation
------------
1. Ensure you have Python 3.x installed on your system.
2. Install required Python packages using pip:
3. Ensure the `pfeature_comp.py` script from the `e3pred` module is available and correctly configured.

Usage
-----
1. Place your input FASTA file in the same directory as the script.
2. Run the script: use the command "e3pred" to run the script after installing.
3. Follow the on-screen prompts:
- Enter the input FASTA filename (e.g., 'XYZ.fasta').
- Select a machine learning model from the list.
4. The script will perform the following steps automatically:
- Parse the input FASTA file and extract sequence information.
- Run feature extraction commands for different feature sets.
- Combine the extracted feature CSV files into one DataFrame.
- Extract specified columns from the combined DataFrame and save them to a CSV file.
- Load the selected machine learning model, make predictions, and assign labels to the sequences.
- Create an output CSV file with the sequences and their predicted labels.
- Call the UbiBrowser API for entries labeled as ubiquitinated to retrieve interaction data.
- Parse the XML responses and save the data to CSV files.

Output
------
The script generates several output files:
- `<base_filename>_extracted_features.csv`: Combined and extracted feature set.
- `<base_filename>_predictions.csv`: Sequences with predicted labels and additional information.
- `<base_filename>_ubiquitination_data.csv`: Protein interaction data for ubiquitinated sequences.

Files and Functions
-------------------
- `script.py`: The main script that orchestrates the entire pipeline.
- `pfeature_comp.py`: External script for feature extraction (part of the `e3pred` module).

Functions in `script.py`:
- `get_input_file()`: Prompts for the input FASTA filename.
- `run_pfeature_comp_command(fasta_file, base_filename, feature)`: Runs the feature extraction command.
- `combine_csv_files_horizontally(input_files)`: Combines CSV files horizontally.
- `extract_columns(dataframe, columns, output_file)`: Extracts specified columns from a DataFrame.
- `load_model(model_filename)`: Loads a specific model.
- `select_model(model_filenames)`: Prompts the user to select a model.
- `run_models_on_test_data(test_data_path, selected_model, scaler_filename)`: Runs models on test data.
- `parse_fasta(fasta_file)`: Parses the FASTA file and extracts sequence information.
- `parse_header_and_sequence(header, sequence)`: Parses individual FASTA headers and sequences.
- `create_output_csv(sequences, labels, output_file)`: Creates an output CSV with sequences and labels.
- `call_ubibrowser_api(url)`: Calls the UbiBrowser API.
- `parse_xml_response(xml_content)`: Parses the XML response from the UbiBrowser API.

Notes
-----
- Ensure all required dependencies are installed and accessible.
- Ensure the `pfeature_comp.py` script is available and correctly configured.
- The UbiBrowser API interaction assumes internet connectivity.

Contact
-------
For any issues or questions, please contact Ronit Gandotra at ronitgandotra@gmail.com.


