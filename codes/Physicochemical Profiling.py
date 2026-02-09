from Bio.SeqUtils.ProtParam import ProteinAnalysis

# This code allows the user to calculate the physicochemical properties of a protein sequence.
# Input - protein structure
# Output - the code prints the MW, pI, GRAVY, and Instability Index 


# function to calc the values of the physiochemial properties
def physicochemical_profile(sequence: str):
   
    sequence = sequence.replace(" ", "").replace("\n", "").upper()   # remove any whitespace and convert everything to upper

    analysis = ProteinAnalysis(sequence)   
    #ProteinAnalysis function is based on ProtParam tools, refrence : https://biopython.org/wiki/ProtParam

    results = {
        "Molecular_Weight (Da)": round(analysis.molecular_weight(), 2),
        "Isoelectric_Point (pI)": round(analysis.isoelectric_point(), 2),
        "GRAVY": round(analysis.gravy(), 3),
        "Instability_Index": round(analysis.instability_index(), 2),
        "Stability": "Stable" if analysis.instability_index() < 40 else "Unstable"
    }

    return results


# executing the function created to calc the values
protein_sequence = """
MYRMQLLSCIALSLALVTNSAPTSSSTKKTQLQLEHLLLDLQMVILNGINNYKNPKLTRMLTFKFYMPKKATELKHLQCLEEELKPLEEVLNLAQSKNFHLRPRDLISNINVIVLELKGSETTFMCEYADEKTATIVEFLNRWITFCQSIISTLT
"""

profile = physicochemical_profile(protein_sequence)

print("Physicochemical Results")
print("------------------------------------------------------------------------")
for key, value in profile.items():
    print(f"{key}: {value}")
print("------------------------------------------------------------------------")
