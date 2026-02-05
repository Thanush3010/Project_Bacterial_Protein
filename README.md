Project 1: In-silico Functional Characterization of a Hypothetical Bacterial Protein
📌 Project Overview

This project demonstrates the use of Biopython for computational analysis of a hypothetical bacterial protein sequence. A complete bioinformatics pipeline was implemented to perform sequence validation, homology search using BLAST, and functional annotation based on sequence similarity.

The project was developed as part of a Biopython course assignment, emphasizing programmatic biological data analysis rather than web-based tools.

🎯 Objectives

To parse and validate a protein sequence in FASTA format

To perform sequence quality analysis using Biopython

To execute BLAST programmatically using Biopython

To identify homologous sequences and infer protein function

To generate reproducible and well-documented analysis outputs

🧪 Methodology

Sequence Input
A hypothetical bacterial protein sequence was provided in FASTA format.

Sequence Quality Analysis
Protein properties such as sequence length, molecular weight, and amino acid composition were calculated using Bio.SeqUtils.ProtParam.

Homology Search (BLAST)
BLASTP was performed programmatically using Biopython’s NCBIWWW.qblast() function against the NCBI non-redundant protein database.

BLAST Result Parsing
The BLAST output was saved in XML format and parsed to extract the top homologous sequences based on score and E-value.

Functional Annotation
The protein’s biological function was inferred using homology-based annotation from the top BLAST hits.

🛠️ Tools & Technologies

Python 3.x

Biopython

NCBI BLAST (programmatic access)

Visual Studio Code

Windows PowerShell
