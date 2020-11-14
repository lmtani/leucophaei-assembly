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

task DownloadSraPairedReads {
    input {
        String sra_illumina_paired
    }

    command <<<
        fastq-dump -I --gzip --split-files ~{sra_illumina_paired}
    >>>

    runtime {
        docker: "cyverseuk/fastq-dump:latest"
    }

    output {
        File r1 = "~{sra_illumina_paired}_1.fastq.gz"
        File r2 = "~{sra_illumina_paired}_2.fastq.gz"
    }
}

task DownloadSraReads {
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

task CreateLastdbIndex {
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

task RunCanuAssemble {
    input {
        String prefix
        String genome_size
        String corrected_error_rate
        File nanopore_reads
        Boolean fast
    }

    command <<<
    canu -p ~{prefix} -d result \
        genomeSize=~{genome_size} \
        ~{true="-fast" false="" fast} \
        correctedErrorRate=~{corrected_error_rate} \
        -nanopore-raw ~{nanopore_reads}
    >>>

    runtime {
        docker: "taniguti/canu:1.8"
        memory: "16 GB"
        cpu: "8"
        disks: "local-disk " + 20 + " HDD"
        preemptible: "3"
    }

    output {
        File contigs_fasta = "result/${prefix}.contigs.fasta"
        File report = "result/~{prefix}.report"
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
        bwa mem -t 7 ~{reference} ~{r1} ~{r2} > align.sam
        samtools view -u align.sam | samtools sort -@ 4 - -o align.bam
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

task RunPilon {
    input {
        File alignments
        File align_idxs
        File reference
    }

    command <<<
        mkdir resultados
        java -Xmx5G -jar /workspace/pilon-1.23.jar --genome ~{reference} \
             --frags ~{alignments} --outdir resultados/ \
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
        File changes = "resultados/pilonChanges.wig"
        File copy_number = "resultados/pilonCopyNumber.wig"
        File coverage = "resultados/pilonCoverage.wig"
        File fasta = "resultados/pilon.fasta"
    }
}


task RunBusco4 {
    input {
        File assembly
        String lineage

    }

    command <<<
        busco --in ~{assembly} \
              --lineage_dataset ~{lineage} \
              --mode prot \
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


workflow Canu2AssemblyNanopore {
    input {
        String genome_size
        String project_name
        String corrected_error_rate
        File reference_gz
        String sra_illumina
        String sra_nanopore
    }

    String index_name = "my-index"

    call DownloadSraReads {
        input:
            sra_reads=sra_nanopore
    }

    String test_name = project_name + "_" + corrected_error_rate
    call RunCanuAssemble {
        input:
            prefix=test_name,
            genome_size=genome_size,
            nanopore_reads=DownloadSraReads.fastq,
            corrected_error_rate=corrected_error_rate,
            fast=false
    }

    call Unzip as UnzipFasta {
        input:
            gzipped_file=reference_gz,
            outfile="reference.fa"
    }

    call CreateLastdbIndex {
        input:
            index_name=index_name,
            reference_fasta=UnzipFasta.out
    }

    call LastAlign {
        input:
            index_name=index_name,
            index_files=CreateLastdbIndex.files,
            new_assembly=RunCanuAssemble.contigs_fasta,
            test_name=test_name
    }

    call DownloadSraPairedReads {
        input:
            sra_illumina_paired=sra_illumina
    }

    call AlignToAssemblyBwa {
        input:
            r1=DownloadSraPairedReads.r1,
            r2=DownloadSraPairedReads.r2,
            reference=UnzipFasta.out
    }

    call RunPilon {
        input:
            reference=UnzipFasta.out,
            alignments=AlignToAssemblyBwa.bam,
            align_idxs=AlignToAssemblyBwa.bai
    }

    call RunBusco4 {
        input:
            assembly=RunPilon.fasta,
            lineage="basidiomycota_odb10"
    }

    output {
        CanuResult results =  {
            "name": test_name,
            "fasta": RunCanuAssemble.contigs_fasta,
            "dotplot": LastAlign.dotplot
        }
        File pilon_fasta = RunPilon.fasta
        File busco_short_table = RunBusco4.short_table
    }
}
