# *S. leucophaei* Genome Assembly

This repository contains a workflow that automates the following tasks:

1. Download Illumina and Nanopore reads from SRA.
1. Assemble Nanopore reads into contigs with Canu.
1. Map Illumina reads into canu contigs.
1. Run Pilon to improve the assembly.
1. Evaluate the assembly completeness with BUSCO.

![diagram](https://user-images.githubusercontent.com/12699242/99152175-e08cf100-267e-11eb-9c3e-a1ff8f4a5c3b.png)

## Inputs

- SRA code for illumina reads
- SRA code for nanopore reads
- http url for the submitted genome (.gz)
- Canu parameters and docker image - fast, genome_size, corrected error rate, correction min coverage.

## Outputs

- Contigs, assembled using Canu and Nanopore reads and corrected with Pilon and Illumina reads.
- BUSCO statistics.
- Dotplot (Last) comparing the assembly with the NCBI submitted genome.


## Benchmark

Here some metrics obtained using different Canu parameters.

| inputs*  | contigs | total length | L50 | pilon assembly | BUSCO                                          |
|----------|---------|--------------|-----|----------------|------------------------------------------------|
| 1 (.144) | 57      | 19177454     | 12  | 19300039       | `C:97.0%[S:96.9%,D:0.1%],F:0.7%,M:2.3%,n:1764` |
| 2 (.105) | 49      | 19075034     | 12  | 19197306       | `C:97.0%[S:96.9%,D:0.1%],F:0.7%,M:2.3%,n:1764` |

> The full list of inputs are defined in the respective `inputs_N.json` file, also in this repository.


- [inputs_1.json](https://1drv.ms/u/s!AiqTy8_f1TgKuV_NMVES9wlEeoSW?e=oaTff6) and [outputs](https://1drv.ms/u/s!AiqTy8_f1TgKuVy5oowCsH8M5lEK?e=7HVsro)
- [inputs_2.json](https://1drv.ms/u/s!AiqTy8_f1TgKuV64knyq9dXUqV2W?e=tajhuL) and [outputs](https://1drv.ms/u/s!AiqTy8_f1TgKuV3Y1uhIzD0qBi1o?e=u8v5y9)
