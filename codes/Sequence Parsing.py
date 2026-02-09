from urllib.request import urlopen  # To fetch UniProt records from the web
from Bio import SwissProt           # To parse UniProt/Swiss-Prot flat files
from io import StringIO             # To treat downloaded text as a file object

# ---- FETCH FROM UNIPROT ----
# UniProt accession for given sequence
# (corresponding UniProt ID from BLAST against nr database)
uniprot_id = "P60568"
url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"

# Fetch UniProt data by API call and parse it using Biopython
with urlopen(url) as handle:
    record = SwissProt.read(StringIO(handle.read().decode()))

# Extract the full protein sequence
sequence = record.sequence

# ---- FIND THE SIGNAL PEPTIDE ----
# Logic: iterate through annotated features to find 'SIGNAL'
signal_peptide_end = None
for feature in record.features:
    if feature.type == "SIGNAL":
        signal_peptide_end = int(feature.location.end)
        break

# ---- FIND THE MATURE PROTEIN (CHAIN) ----
mature_start = None
mature_end = None

# Logic: iterate through annotated features to find 'CHAIN'
for feature in record.features:
    if feature.type == "CHAIN":
        # UniProt CHAIN feature defines the mature protein
        mature_start = int(feature.location.start) + 1
        mature_end = int(feature.location.end)
        break

# ---- EXTRACT AND PRINT RESULTS ----
if signal_peptide_end:
    print(f"Signal peptide: residues 1–{signal_peptide_end}")    
else:
    print("Signal peptide is missing")

if mature_start and mature_end and signal_peptide_end:
    # Extract mature protein sequence
    mature_sequence = sequence[mature_start - 1 : mature_end]
    print(f"Mature protein region: residues {mature_start}–{mature_end}")
    print("Mature protein sequence:")
    print(mature_sequence)
else:
    print("Mature protein is missing")
