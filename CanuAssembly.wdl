version 1.0


struct CanuResult {
    String name
    File? dotplot
    File? fasta
}

struct IlluminaReads {
    File r1
    File r2
}

task Unzip {
    input {
        File gzipped_file
        String outfile
    }

    command <<<
        gunzip -c ~{gzipped_file} > ~{outfile}

    >>>

    runtime {
        docker: "actions/gzip:1.6-5"
    }

    output {
        File out = "~{outfile}"
    }
}

task GetSraIllumina {
    input {
        String sra_illumina_paired
    }

    command <<<
        fastq-dump --gzip --split-files ~{sra_illumina_paired}
    >>>

    runtime {
        docker: "cyverseuk/fastq-dump:latest"
    }

    output {
        File r1 = "~{sra_illumina_paired}_1.fastq.gz"
        File r2 = "~{sra_illumina_paired}_2.fastq.gz"
    }
}

task GetSraNanopore {
    input {
        String sra_reads
    }

    command <<<
        fastq-dump -I --gzip ~{sra_reads}
    >>>

    runtime {
        docker: "cyverseuk/fastq-dump:latest"
    }

    output {
        File fastq = "~{sra_reads}.fastq.gz"
    }
}

task LastdbIndex {
    input {
        String index_name
        File reference_fasta
    }

    command <<<
        lastdb ~{index_name} ~{reference_fasta}
    >>>

    runtime {
        docker: "biocontainers/last-align:v963-2-deb_cv1"
    }

    output {
        Array[File] files = glob("~{index_name}.*")
    }
}

task LastAlign {

    input {
        String index_name
        String test_name
        Array[File] index_files
        File new_assembly
    }

    command <<<
        set -e
        for i in ~{sep=" " index_files}; do ln -s $i .; done
        lastal ~{index_name} ~{new_assembly} > ~{test_name}.paf
        last-dotplot --rot2=h -x 1500 -y 1500 ~{test_name}.paf ~{test_name}.png
    >>>

    runtime {
        docker: "biocontainers/last-align:v963-2-deb_cv1"
    }
    output {
        File paf = "~{test_name}.paf"
        File dotplot = "~{test_name}.png"
    }
}

task CanuAssemble {
    input {
        File nanopore_reads
        CanuParameters p
    }

    String corrected_error_rate = if (defined(p.corrected_error_rate)) then "correctedErrorRate=~{p.corrected_error_rate}"  else ""
    String correction_min_coverage = if (defined(p.correction_min_coverage)) then "corMinCoverage=~{p.correction_min_coverage}"  else ""

    command <<<
        canu -p ~{p.prefix} -d result \
            genomeSize=~{p.genome_size} \
            ~{true="-fast" false="" p.fast} \
            ~{corrected_error_rate} \
            ~{correction_min_coverage} \
            -nanopore-raw ~{nanopore_reads}
    >>>

    runtime {
        docker: "~{p.container_image}"
        memory: "16 GB"
        cpu: "8"
        disks: "local-disk " + 20 + " HDD"
        preemptible: "3"
    }

    output {
        File contigs_fasta = "result/${p.prefix}.contigs.fasta"
        File report = "result/~{p.prefix}.report"
    }
}


task AlignToAssemblyBwa {
    input {
        File r1
        File r2
        File reference
    }

    command <<<
        set -e
        bwa index ~{reference}
        bwa mem -t 7 ~{reference} ~{r1} ~{r2} | samtools sort -@ 4 -o align.bam -
        samtools index align.bam
    >>>

    runtime {
        docker: "taniguti/bwa-samtools"
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + 200 + " HDD"
        preemptible: "3"
    }

    output {
        File bam = "align.bam"
        File bai = "align.bam.bai"
    }
}

task Pilon {
    input {
        String project_name
        File alignments
        File align_idxs
        File reference
    }

    command <<<
        mkdir resultados
        java -Xmx5G -jar /workspace/pilon-1.23.jar --genome ~{reference} \
             --frags ~{alignments} --outdir resultados/ --output ~{project_name} \
             --changes --fix bases,gaps --tracks --threads 5
    >>>

    runtime {
        docker: "taniguti/pilon"
        memory: "16 GB"
        cpu: "4"
        disks: "local-disk " + 200 + " HDD"
        preemptible: "3"
    }

    output {
        File changes = "resultados/~{project_name}Changes.wig"
        File copy_number = "resultados/~{project_name}CopyNumber.wig"
        File coverage = "resultados/~{project_name}Coverage.wig"
        File fasta = "resultados/~{project_name}.fasta"
    }
}


task Busco4 {
    input {
        File assembly
        String lineage

    }

    command <<<
        busco --in ~{assembly} \
              --lineage_dataset ~{lineage} \
              --mode geno \
              -o results \
              -c 7 \
              --augustus_species ustilago_maydis
    >>>

    runtime {
        docker: "ezlabgva/busco:v4.1.2_cv1"
    }

    output {
        File short_table = "results/short_summary.specific.basidiomycota_odb10.results.txt"
    }
}


struct CanuParameters {
    String genome_size
    String prefix
    Boolean fast
    String container_image
    String? corrected_error_rate
    String? correction_min_coverage
}

workflow CanuAssemblyNanopore {
    input {
        String project_name
        # Assembly parameters
        CanuParameters canu_parameters
        # Inputs from NCBI
        File reference_gz
        String sra_illumina
        String sra_nanopore
    }

    String index_name = "my-index"

    call GetSraNanopore {
        input:
            sra_reads=sra_nanopore
    }

    call CanuAssemble {
        input:
            nanopore_reads=GetSraNanopore.fastq,
            p=canu_parameters
    }

    call Unzip as UnzipFasta {
        input:
            gzipped_file=reference_gz,
            outfile="reference.fa"
    }

    call LastdbIndex {
        input:
            index_name=index_name,
            reference_fasta=UnzipFasta.out
    }

    call LastAlign {
        input:
            index_name=index_name,
            index_files=LastdbIndex.files,
            new_assembly=CanuAssemble.contigs_fasta,
            test_name=project_name
    }

    call GetSraIllumina {
        input:
            sra_illumina_paired=sra_illumina
    }

    call AlignToAssemblyBwa {
        input:
            r1=GetSraIllumina.r1,
            r2=GetSraIllumina.r2,
            reference=CanuAssemble.contigs_fasta
    }

    call Pilon {
        input:
            project_name=project_name,
            reference=CanuAssemble.contigs_fasta,
            alignments=AlignToAssemblyBwa.bam,
            align_idxs=AlignToAssemblyBwa.bai
    }

    call Busco4 {
        input:
            assembly=Pilon.fasta,
            lineage="basidiomycota_odb10"
    }

    output {
        CanuResult results =  {
            "name": project_name,
            "fasta": CanuAssemble.contigs_fasta,
            "dotplot": LastAlign.dotplot
        }
        File pilon_fasta = Pilon.fasta
        File busco_short_table = Busco4.short_table
    }
}
