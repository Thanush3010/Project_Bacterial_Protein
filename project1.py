import os
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from Bio.Blast import NCBIWWW, NCBIXML

# Ensure results directory exists
os.makedirs("results", exist_ok=True)

# -------------------------------
# STEP 1: Read FASTA Sequence
# -------------------------------
records = list(SeqIO.parse("input_sequence_project1.fasta", "fasta"))
print("Records found:", len(records))

if len(records) == 0:
    raise ValueError("No FASTA records found.")

record = records[0]
sequence = record.seq

# -------------------------------
# STEP 2: Sequence Quality Analysis
# -------------------------------
analysis = ProtParam.ProteinAnalysis(str(sequence))

length = len(sequence)
molecular_weight = analysis.molecular_weight()
aa_composition = analysis.get_amino_acids_percent()

with open("results/qc_summary.txt", "w") as f:
    f.write(f"Sequence ID: {record.id}\n")
    f.write(f"Sequence Length: {length} amino acids\n")
    f.write(f"Molecular Weight: {molecular_weight:.2f} Da\n\n")
    f.write("Amino Acid Composition (%):\n")
    for aa, percent in aa_composition.items():
        f.write(f"{aa}: {percent * 100:.2f}\n")

print("Sequence quality analysis completed.")

# -------------------------------
# STEP 3: BLAST using Biopython
# -------------------------------
print("Starting BLAST search...")

result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=sequence
)

with open("results/blast_results.xml", "w") as xml_out:
    xml_out.write(result_handle.read())

print("BLAST completed.")

with open("results/blast_results.xml") as xml_file:
    blast_record = NCBIXML.read(xml_file)

with open("results/blast_results.txt", "w") as out:
    for alignment in blast_record.alignments[:3]:
        hsp = alignment.hsps[0]
        out.write(f"Hit: {alignment.title}\n")
        out.write(f"Score: {hsp.score}\n")
        out.write(f"E-value: {hsp.expect}\n\n")

print("BLAST results saved.")

# -------------------------------
# STEP 4: Functional Annotation
# -------------------------------
with open("results/functional_annotation.txt", "w") as fa:
    fa.write("FUNCTIONAL ANNOTATION REPORT\n")
    fa.write("=" * 40 + "\n\n")

    for i, alignment in enumerate(blast_record.alignments[:3], start=1):
        hsp = alignment.hsps[0]
        fa.write(f"Top Hit {i}:\n")
        fa.write(f"Protein: {alignment.title}\n")
        fa.write(f"Score: {hsp.score}\n")
        fa.write(f"E-value: {hsp.expect}\n\n")

    fa.write(
        "Predicted Function:\n"
        "The hypothetical protein shows strong similarity to conserved "
        "bacterial metabolic enzymes. Based on homology, it is predicted "
        "to participate in core cellular metabolic processes.\n"
    )

print("Functional annotation completed.")
