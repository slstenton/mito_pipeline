version 1.0

workflow mitopipeline {

    input {
  	    File cram
	    File refFasta
        File refFastaIndex
        File MToolBox_config
        File MitoSAlt_config
        String name
        Int runtime_cpu
        Int runtime_preemptible
    }
        
    call cram_to_fastq {
        input:
            cram = cram,
            refFasta = refFasta,
            refFastaIndex = refFastaIndex,
            name = name,
            runtime_cpu = runtime_cpu,
            collate_threads = runtime_cpu * 2,
			runtime_memory = floor(runtime_cpu * 6.5 * 0.25),
            runtime_disk_gb = 10 + ceil(size(cram, 'GB') * 12),
            runtime_preemptible = runtime_preemptible
    }
    
    call mtDNA_SNV_indel {
        input:
            MToolBox_config = MToolBox_config,
            fastq_R1 = cram_to_fastq.fastq_R1,
            fastq_R2 = cram_to_fastq.fastq_R2,
            name = name,
            runtime_cpu = runtime_cpu,
            runtime_memory = floor(runtime_cpu * 6.5),
            runtime_disk_gb = 10 + ceil(size(cram, 'GB') * 24),
            runtime_preemptible = runtime_preemptible
    }

    call mtDNA_deletion {
        input:
        	MitoSAlt_config = MitoSAlt_config,
			fastq_mt1 = mtDNA_SNV_indel.fastq_mt1,
            fastq_mt2 = mtDNA_SNV_indel.fastq_mt2,
            name = name,
            runtime_cpu = runtime_cpu,
            runtime_memory = floor(runtime_cpu * 6.5 * 0.25),
            runtime_disk_gb = 5 + ceil(size(cram, 'GB') * 24),
            runtime_preemptible = runtime_preemptible
    }

    output {
  	    File MToolBox_vcf = mtDNA_SNV_indel.MToolBox_vcf
	    File MToolBox_calls = mtDNA_SNV_indel.MToolBox_calls
	    File MToolBox_anno = mtDNA_SNV_indel.MToolBox_anno
	    File MitoSAlt_bam = mtDNA_deletion.MitoSAlt_bam
        File MitoSAlt_bw = mtDNA_deletion.MitoSAlt_bw
        File MitoSAlt_breakpoint = mtDNA_deletion.MitoSAlt_breakpoint
        File MitoSAlt_cluster = mtDNA_deletion.MitoSAlt_cluster
        File MitoSAlt_summary = mtDNA_deletion.MitoSAlt_summary
    }
}

task cram_to_fastq {
    input {
        File cram
        File refFasta
        File refFastaIndex
        String name
        Int runtime_cpu
        Int collate_threads
        Int runtime_memory
        Int runtime_disk_gb
        Int runtime_preemptible
    }
    command <<<
      
        set -euxo pipefail
      
        # convert cram to bam
        samtools view -b -T ~{refFasta} -o ~{name}.bam ~{cram}
      
        # convert bam to fastq
        samtools flagstat -@ ~{runtime_cpu} ~{name}.bam > ~{name}_flagstat.txt
        samtools collate -O -u -n 128 --threads ~{collate_threads} ~{name}.bam temp \
        | samtools fastq -F 0x900 -n -t --threads 2 \
              -1 ~{name}.R1.fastq.gz \
              -2 ~{name}.R2.fastq.gz \
              -0 >(tee >(expr `wc -l` / 4 > count_0.txt) | gzip > ~{name}_0.fq.gz) \
              -s >(tee >(expr `wc -l` / 4 > count_single.txt) | gzip > ~{name}_single.fq.gz) \
              -
    >>>
    output {
        File fastq_R1 = "~{name}.R1.fastq.gz"
        File fastq_R2 = "~{name}.R2.fastq.gz"
    }
    runtime {
        disks: 'local-disk ~{runtime_disk_gb} HDD'
        cpu: '~{runtime_cpu}'
        memory: '~{runtime_memory} GB'
        docker: 'quay.io/biocontainers/samtools:1.9--h8571acd_11'
        preemptible: '~{runtime_preemptible}'
    }
}

task mtDNA_SNV_indel {
    input {
        File MToolBox_config
        File fastq_R1
        File fastq_R2
        String name
        Int runtime_cpu
        Int runtime_memory
        Int runtime_disk_gb
        Int runtime_preemptible
    }
    command <<<
        
        set -euxo pipefail
  
        # run MToolBox
        mv ~{fastq_R1} .
        mv ~{fastq_R2} .
        ../MToolBox/MToolBox/MToolBox.sh -i ~{MToolBox_config}
        mv OUT_~{name}/* .
    >>>
    output {
        File MToolBox_vcf = 'sample.vcf'
        File MToolBox_calls = '~{name}-table.txt'
        File MToolBox_anno = glob("*.annotation.csv")[0]
        File fastq_mt1 = 'outmt1.fastq'
        File fastq_mt2 = 'outmt2.fastq'
    }
    runtime {
        disks: 'local-disk ~{runtime_disk_gb} HDD'
        cpu: '~{runtime_cpu}'
        memory: '~{runtime_memory} GB'
        docker: 'sstenton/mtoolbox'
        preemptible: '~{runtime_preemptible}'
  }
}

task mtDNA_deletion {
    input {
        File MitoSAlt_config
        File fastq_mt1
        File fastq_mt2
        String name
        Int runtime_cpu
        Int runtime_memory
        Int runtime_disk_gb
        Int runtime_preemptible
    }
    command <<<
    
        set -euxo pipefail
    
        cd /
        mv MitoSAlt_1.1.1 cromwell_root
        cd cromwell_root/MitoSAlt_1.1.1
        perl MitoSAlt1.1.1.pl ~{MitoSAlt_config} ~{fastq_mt1} ~{fastq_mt2} ~{name}
        cd /
        MitoSAlt_indel_PATH=$(find . -name 'indel' -d)
        MitoSAlt_bam_PATH=$(find . -name 'bam' -d | egrep -v internal/bam)
        MitoSAlt_bw_PATH=$(find . -name 'bw' -d)
        mv "$MitoSAlt_indel_PATH" ./cromwell_root
        mv "$MitoSAlt_bam_PATH" ./cromwell_root
        mv "$MitoSAlt_bw_PATH" ./cromwell_root
        cd cromwell_root/
        touch ~{name}.tsv
        mv indel/* .
        mv bam/* .
        mv bw/* .
    >>>
    output {
        File MitoSAlt_bam = '~{name}.bam'
        File MitoSAlt_bw = '~{name}.bw'
        File MitoSAlt_breakpoint = '~{name}.breakpoint'
        File MitoSAlt_cluster = '~{name}.cluster'
        File MitoSAlt_summary = '~{name}.tsv'
    }
    runtime {
        disks: 'local-disk ~{runtime_disk_gb} HDD'
        cpu: '~{runtime_cpu}'
        memory: '~{runtime_memory} GB'
        docker: 'sstenton/mitosalt'
        preemptible: '~{runtime_preemptible}'
    }
}