import pandas as pd

def detect_aggregation_propensity(sequence, netsurfp_file):
    """
    Step 1: Identify hydrophobic regions based on sequence only
            (>4 consecutive hydrophobic residues using sliding window)
    Step 2: From those, identify surface-exposed regions using NetSurfP RSA
    """

    hydrophobic = "AVILMFWY"     # hydrophobic amino acids
    window_size = 4             # as given in the question
    rsa_cutoff = 0.25            # surface-exposed threshold

    # -----------------------------
    # LOAD NETSURFP OUTPUT
    # -----------------------------
    df = pd.read_csv(netsurfp_file)
    df.columns = df.columns.str.strip().str.lower()

    # Identify residue position column robustly
    pos_col = next(
        c for c in df.columns
        if c in ["n", "pos", "residue", "residue_number", "#"]
    )

    # Build RSA lookup dictionary
    rsa_dict = dict(zip(df[pos_col], df["rsa"]))

    found_any = False
    found_surface = False

    print("Hydrophobic regions (sequence-only):")

    # -----------------------------
    # SLIDING WINDOW SCAN
    # -----------------------------
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]

        # Check if all residues are hydrophobic
        is_hydrophobic = True
        for aa in window:
            if aa not in hydrophobic:
                is_hydrophobic = False
                break

        if is_hydrophobic:
            start = i + 1
            end = i + window_size
            print(f"Residues {start}–{end}: {window}")
            found_any = True

            # -----------------------------
            # NETSURFP SURFACE CHECK
            # -----------------------------
            if all(rsa_dict.get(pos, 0) >= rsa_cutoff for pos in range(start, end + 1)):
                print(f"  → Surface-exposed hydrophobic patch")
                found_surface = True

    if not found_any:
        print("No hydrophobic regions found.")

    if not found_surface:
        print("\nNo surface-exposed hydrophobic patches detected.")




sequence = 'MYRMQLLSCIALSLALVTNSAPTSSSTKKTQLQLEHLLLDLQMVILNGINNYKNPKLTRMLTFKFYMPKKATELKHLQCLEEELKPLEEVLNLAQSKNFHLRPRDLISNINVIVLELKGSETTFMCEYADEKTATIVEFLNRWITFCQSIISTLT'

detect_aggregation_propensity(sequence, "netsurfp_output.csv")