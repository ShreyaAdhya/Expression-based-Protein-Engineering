from Bio.Seq import Seq
from Bio.Data import CodonTable

# -----------------------------
# INPUTS
# -----------------------------
optimized_dna = seq   
original_protein = prot_seq 

bamHI = "GGATCC"
xhoI  = "CTCGAG"
stop_codons = ["TAA", "TAG", "TGA"]
polyA_signals = ["AATAAA"]



# -----------------------------
# Check restriction sites
# -----------------------------
def check_restriction_sites(seq):
    if bamHI in seq:
        return False
    if xhoI in seq:
        return False
    return True

# -----------------------------
# Check start and stop codon
# -----------------------------
def check_start_stop(seq):
    start_ok = False
    stop_ok = False

    if seq.startswith("ATG"):
        start_ok = True

    last_codon = seq[-3:]
    if last_codon in stop_codons:
        stop_ok = True

    if start_ok and stop_ok:
        return True
    else:
        return False

# -----------------------------
# Check reading frame
# -----------------------------
def check_frame(seq):
    length = len(seq)
    if length % 3 == 0:
        return True
    else:
        return False

# -----------------------------
# Check premature stop codons
# -----------------------------
def check_premature_stops(seq):
    internal_seq = seq[3:-3]

    for i in range(0, len(internal_seq), 3):
        codon = internal_seq[i:i+3]
        if codon in stop_codons:
            return False

    return True

# -----------------------------
# Calculate GC content
# -----------------------------
def gc_content(seq):
    g_count = seq.count("G")
    c_count = seq.count("C")
    gc = g_count + c_count
    gc_percent = (gc / len(seq)) * 100
    return gc_percent

# -----------------------------
# Check GC range
# -----------------------------
def check_gc_range(seq, low=40, high=50):
    gc = gc_content(seq)

    if gc >= low and gc <= high:
        return True
    else:
        return False


# -----------------------------
# Check polyA signals
# -----------------------------
def check_polyA(seq):
    for signal in polyA_signals:
        if signal in seq:
            return False
    return True

# -----------------------------
# Translate DNA to protein
# -----------------------------
def translate(seq):
    dna_seq = Seq(seq)
    protein = dna_seq.translate(to_stop=False)
    return str(protein)


# -----------------------------
# RUN ALL CHECKS
# -----------------------------
results = {
    "No internal BamHI/XhoI": check_restriction_sites(optimized_dna),
    "Correct start/stop codon": check_start_stop(optimized_dna),
    "Correct reading frame": check_frame(optimized_dna),
    "No premature stop codons": check_premature_stops(optimized_dna),
    "GC content (40–50%) ": check_gc_range(optimized_dna),
    "No polyA signals": check_polyA(optimized_dna),
    "Protein sequence match": translate(optimized_dna) == original_protein
}

# -----------------------------
# REPORT
# -----------------------------
print("VALIDATION\n" + "-"*30)
for k, v in results.items():
    print(f"{k}: {'PASS' if v else 'FAIL'}")

print("\nGC Content:", round(gc_content(optimized_dna), 2), "%")
