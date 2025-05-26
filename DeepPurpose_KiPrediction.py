##Homework PYTHON SCRIPT - DeepPurpose Ki prediction
#Explanations line by line made with Gemini

import pandas as pd
from DeepPurpose import DTI
from DeepPurpose import utils
from Bio import SeqIO

# This Python script is designed to predict Drug-Target Interactions (DTI)
# using a pre-trained model from the DeepPurpose library.
# It reads drug data (in SMILES format) from an Excel file
# and a protein sequence (target) from a FASTA file,
# then uses the model to predict interaction scores.

# === 1. Read Analogues from Excel ===
# This section handles extracting drug (chemical compound) data
# from your Excel file.

df = pd.read_excel("ANALOGI COLOCVIU.xlsx")
# `pd.read_excel()` is a function from the `pandas` library (aliased as `pd`).
# It reads the content of the Excel file named "ANALOGI COLOCVIU.xlsx"
# and loads it into a pandas DataFrame object.
# A DataFrame is a tabular data structure, similar to a spreadsheet.

# Find the SMILES column (assuming it's named "SMILES" or similar)
smiles_col = [col for col in df.columns if 'smiles' in col.lower()]
# This line iterates through all column names in the `df` DataFrame.
# For each column name (`col`), it checks if the string 'smiles' (in lowercase)
# is present within the column name (also converted to lowercase for case-insensitive search).
# The result is a list (`smiles_col`) containing the names of columns that meet this condition.
# Example: if you have a column named "SMILES", "SMILES_STRING", or "drug_smiles", it will be detected.

if not smiles_col:
    raise ValueError("SMILES column not found in the Excel file.")
# This `if` statement checks if the `smiles_col` list is empty.
# If it's empty (meaning no column containing 'smiles' in its name was found),
# it raises a `ValueError` with the specified message.
# This ensures the script stops if it cannot find the necessary data.

smiles_list = df[smiles_col[0]].astype(str).tolist()
# `smiles_col[0]` takes the first column name from the `smiles_col` list (assuming you have one relevant column).
# `df[smiles_col[0]]` selects that column from the DataFrame.
# `.astype(str)` converts all values in that column to string data type.
# This is important because SMILES (Simplified Molecular-Input Line-Entry System) are textual representations of chemical structures.
# `.tolist()` converts the pandas Series (the selected column) into a standard Python list.
# Thus, `smiles_list` will now contain a list of SMILES strings from your Excel file.

# === 2. Read Protein FASTA Sequence ===
# This section handles extracting the target protein sequence from a FASTA file.

record = next(SeqIO.parse("tinta.fasta", "fasta"))
# `SeqIO.parse()` is a function from the `BioPython` library (aliased as `Bio.SeqIO`).
# It's used to parse biological sequence files.
# "tinta.fasta" is the name of your FASTA file containing the protein sequence.
# "fasta" specifies the file format.
# `SeqIO.parse()` returns an iterator, and `next()` retrieves the first "record" from that file.
# A FASTA file can contain multiple sequences, but here we assume you have a single target sequence.

protein_sequence = str(record.seq)
# `record.seq` accesses the sequence object within the FASTA record.
# `str()` converts the sequence object into a standard Python string.
# `protein_sequence` will now store the amino acid sequence of your target protein as a string.

# === 3. Load Pre-trained DeepPurpose Model ===
# This section loads a pre-trained deep learning model from the DeepPurpose library.

print("Loading pre-trained model...")
# Displays an informative message to the console.

model = DTI.model_pretrained(model='MPNN_CNN_DAVIS')
# This is a crucial line.
# `DTI.model_pretrained()` is a function from the `DTI` module of the DeepPurpose library.
# It is used to download and load a pre-trained Drug-Target Interaction model.
# `model='MPNN_CNN_DAVIS'` specifies the exact name of the pre-trained model you want to use.
# 'MPNN_CNN_DAVIS' is a common DTI model in DeepPurpose, trained on the DAVIS dataset.
# `model` will now store the DeepPurpose model instance, ready for use.

# === 4. Prepare Inputs for Prediction ===
# This section organizes the data into a format that the model can understand for prediction.

X_drugs = smiles_list
# Creates a variable `X_drugs` which is simply your list of SMILES strings.

X_targets = [protein_sequence] * len(smiles_list)
# Creates an `X_targets` variable.
# Since we want to predict the interaction of each drug in `smiles_list` with the *same* target protein,
# this line creates a list where the single protein sequence is replicated `len(smiles_list)` times.
# This ensures you have a list of protein sequences of the same length as your drug list,
# allowing for a one-to-one correspondence for prediction.

dummy_labels = [0] * len(X_drugs)
# DeepPurpose expects a "Label" column in the input DataFrame, even for prediction
# (for cases where you don't have actual labels).
# This line creates a list of zeros, of the same length as `X_drugs`.
# These will serve as "dummy" labels for the prediction process.

# Define encoding types based on the pre-trained model (MPNN_CNN_DAVIS)
drug_encoding_type = 'MPNN'
target_encoding_type = 'CNN'
# These lines define the encoding types that the 'MPNN_CNN_DAVIS' pre-trained model expects.
# "MPNN" (Message Passing Neural Network) is an encoding method for drugs (molecular graphs).
# "CNN" (Convolutional Neural Network) is an encoding method for protein sequences.
# These are essential for `utils.data_process` to know how to handle your data.

# Process the data using utils.data_process for prediction
# Use split_method='no_split' as we are only predicting, not splitting data for training/validation/test
print("Processing data for prediction...")
# Displays an informative message.

data_for_prediction = utils.data_process(X_drugs, X_targets, dummy_labels,
                                               drug_encoding=drug_encoding_type,
                                               target_encoding=target_encoding_type,
                                               split_method='no_split')
# This is another crucial line.
# `utils.data_process()` is a function from the `utils` module of DeepPurpose.
# It takes your raw data (SMILES, protein sequences, and labels)
# and processes it into the internal format required by the DeepPurpose model.
# - `X_drugs`: the list of SMILES.
# - `X_targets`: the list of protein sequences.
# - `dummy_labels`: the list of dummy labels.
# - `drug_encoding=drug_encoding_type`: specifies the encoding type for drugs.
# - `target_encoding=target_encoding_type`: specifies the encoding type for proteins.
# - `split_method='no_split'`: indicates that we do not want to split the data into training/validation/test sets.
#   We just want to process all data into a single set for prediction.
# The result, `data_for_prediction`, will be an internal DeepPurpose object (typically a pandas DataFrame with columns like 'drug_encoding', 'target_encoding', 'Label')
# that is ready to be fed to the model.

# === 5. Run Prediction ===
# This section uses the loaded model to make predictions on the prepared data.

print("Making predictions...")
# Displays an informative message.

predictions = model.predict(data_for_prediction)
# `model.predict()` is the method that performs the actual predictions.
# We pass it the `data_for_prediction` object (which contains the processed drugs and targets).
# The model processes this data and returns an array (or list) of numerical values,
# where each value represents the predicted interaction score for a drug-target pair.
# `predictions` will store these scores.

# === 6. Save Results to CSV ===
# This section saves the obtained predictions to a CSV file.

results_df = pd.DataFrame({
    'SMILES': smiles_list,
    'Predicted_Score': predictions
})
# Creates a new pandas DataFrame named `results_df`.
# It will have two columns:
# - 'SMILES': contains your original list of SMILES for the drugs.
# - 'Predicted_Score': contains the prediction scores obtained from the model.

results_df.to_csv("predictii_Ki_DeepPurpose.csv", index=False)
# `results_df.to_csv()` saves the `results_df` DataFrame to a CSV file.
# "predictii_Ki_DeepPurpose.csv" is the name of your output file.
# `index=False` prevents pandas from writing the DataFrame index as a separate column in the CSV file.

print("✅ Predictions saved to 'predictii_Ki_DeepPurpose.csv'")
# Displays a confirmation message after saving the results.

---

### How to Reproduce and Use It:

To run this script, you will need the following:

1.  **Python Environment:** Ensure you have a Python environment (preferably Python 3.7+).
2.  **Installed Libraries:**
    * `pandas`: `pip install pandas`
    * `DeepPurpose`: `pip install DeepPurpose`
    * `biopython`: `pip install biopython`
    * **Note:** DeepPurpose depends on PyTorch. Installing DeepPurpose should install most dependencies, but if you encounter PyTorch-related errors, install it separately following instructions from the PyTorch website (choose the correct version for your system and GPU, if you have one). For CPU-only: `pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu`
3.  **Input Files:**
    * An Excel file named `ANALOGI COLOCVIU.xlsx` in the same directory as your script. This file must contain a column with SMILES. The column name can be "SMILES", "smiles", "SMILES_STRING", etc., as long as it contains "smiles" (case-insensitive).
    * A FASTA file named `tinta.fasta` in the same directory as your script. This file should contain your target protein sequence.

**Steps to Run:**

1.  **Save the Code:** Copy the code above into a text file and save it as `script.py` (or any name you prefer, e.g., `dti_prediction.py`).
2.  **Place Input Files:** Ensure `ANALOGI COLOCVIU.xlsx` and `tinta.fasta` are in the same directory as your `script.py`.
3.  **Open a Terminal:** Navigate to the directory where you saved your files.
4.  **Run the Script:** Execute the command: `python script.py`

After running, a file named `predictii_Ki_DeepPurpose.csv` will be created in the same directory, containing your original drug SMILES and their predicted interaction scores.

**To Reproduce with Other Data:**

* **Drugs:** Edit the `ANALOGI COLOCVIU.xlsx` file with your new SMILES. Ensure the format of the SMILES column remains consistent.
* **Target Protein:** Edit the `tinta.fasta` file with your new target protein sequence.
