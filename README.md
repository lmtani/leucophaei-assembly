# *S. leucophaei* Genome Assembly
[![DOI:10.1094/MPMI-08-20-0218-A](https://zenodo.org/badge/DOI/10.1094/MPMI-08-20-0218-A.svg)](https://doi.org/10.1094/MPMI-08-20-0218-A)


This repository contains a workflow that automates the following tasks:

1. Download Illumina and Nanopore reads from SRA.
1. Assemble Nanopore reads into contigs with Canu.
1. Map Illumina reads into canu contigs.
1. Run Pilon to improve the assembly.
1. Evaluate the assembly completeness with BUSCO.
1. Compare contigs with the submitted assembly using LAST

![diagram](https://user-images.githubusercontent.com/12699242/99152175-e08cf100-267e-11eb-9c3e-a1ff8f4a5c3b.png)

## About WDL format

More information about this [workflow language can be found here](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md).

To run it you'll need to use a tool like [cromwell](https://cromwell.readthedocs.io/en/stable/)

This workflow has [Docker](https://www.docker.com/products/docker-desktop) as a dependency. It allows us to use images with pre-installed bioinformatics programs.

Example:

```bash
java -jar <path-to-cromwell.jar> run -i CanuAssembly.inputs.json CanuAssembly.wdl
```

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

\* The full list of inputs are defined in the respective `inputs_N.json` file, available in the following links.


- [inputs_1.json](https://1drv.ms/u/s!AiqTy8_f1TgKuV_NMVES9wlEeoSW?e=oaTff6) and [outputs](https://1drv.ms/u/s!AiqTy8_f1TgKuVy5oowCsH8M5lEK?e=7HVsro)
- [inputs_2.json](https://1drv.ms/u/s!AiqTy8_f1TgKuV64knyq9dXUqV2W?e=tajhuL) and [outputs](https://1drv.ms/u/s!AiqTy8_f1TgKuV3Y1uhIzD0qBi1o?e=u8v5y9)


## Cite us

This repository contains the assembly workflow used in:

>Crestana GS, Taniguti LM, Dos Santos CP, Benevenuto J, Ceresini P, Carvalho G, Kitajima JP, Monteiro-Vitorello CBB. **Complete Chromosome-Scale Genome Sequence Resource for Sporisorium panici-leucophaei, the Causal Agent of Sourgrass Smut Disease. Mol Plant Microbe Interact.** 2020 Dec 28. doi: 10.1094/MPMI-08-20-0218-A. PMID: 33369501.

