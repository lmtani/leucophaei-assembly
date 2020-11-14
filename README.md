# *S. leucophaei* Genome Assembly

This repository contains a workflow that automates the following tasks:

1. Download Illumina and Nanopore reads from SRA.
1. Assemble Nanopore reads into contigs with Canu.
1. Map Illumina reads into canu contigs.
1. Run Pilon to improve the assembly.
1. Evaluate the assembly with BUSCO.

![diagram](https://user-images.githubusercontent.com/12699242/99152175-e08cf100-267e-11eb-9c3e-a1ff8f4a5c3b.png)

## Results

