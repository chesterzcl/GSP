params.genome_dir="/home/zl436/palmer_scratch/ref"
params.genome="ROS.fna"
params.cohort_name="23K9"
params.data_dir="/home/zl436/palmer_scratch/fastq"
params.reads="${params.data_dir}/*_{1,2}.fastq.gz"
params.group_label_file="${params.data_dir}/hoh_test_breed.txt"

process prepare_genome_samtools{

        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        input:
        path genome

	output:
	path "${genome}.fai"

        script:
        """
        samtools faidx ${genome}
        """
}

process prepare_genome_bwa{

        container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

        input:
        path genome

	output:
	path "${genome.baseName}.*"

        script:
        """
        bwa index ${genome} -p ${genome.baseName}
        """
}

process prepare_genome_picard {

	container 'quay.io/biocontainers/picard:1.141--hdfd78af_6'
	
	input:
	path genome

	output:
	path "${genome.baseName}.dict"
	
	script:
	"""
	picard CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
	"""
}

process extract_chromosome_info{
	input:
	path genome

	output:
	path "chromosome_list.txt"

	script:
	"""
	grep '^>NC' ${genome}|cut -d " " -f 1|cut -c 2- >"chromosome_list.txt"
	"""
}

process reads_quality_control{

	container 'quay.io/biocontainers/fastqc:0.11.9--0'
	
	cpus 5

	input:
	tuple val(replicateID),path(reads)
	
	output:
	tuple val(replicateID),path("*_1_fastqc.zip"),path("*_2_fastqc.zip")

	script:
	"""
	fastqc -t ${task.cpus} ${reads}
	"""
}

process get_trimming_parameters{
	
	container 'chesterli214/gsexp_amd64'

	input:
	tuple val(replicateID),path(QC_results_1),path(QC_results_2)

	output:
	tuple val(replicateID),env(trimming_params)

	script:
	"""
	unzip ${QC_results_1}
	unzip ${QC_results_2}
	qual1=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
	qual2=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
	trimming_params="LEADING:8 TRAILING:10 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:10:\${qual1} MINLEN:50"
	"""
}

process quality_trimming{

	container 'quay.io/biocontainers/trimmomatic:0.36--3'
	
	cpus 5

	input:
	tuple val(replicateID),path(reads)
	tuple val(replicateID),val(trimming_params)

	output:
	tuple val(replicateID),path("*_{1,2}_trimmed.fastq.gz")	

	script:
	"""
	trimmomatic PE -threads ${task.cpus} ${reads} ${replicateID}_1_trimmed.fastq.gz ${replicateID}_1_unpaired.fastq.gz ${replicateID}_2_trimmed.fastq.gz ${replicateID}_2_unpaired.fastq.gz ${trimming_params}
	"""
}

process reads_mapping{
	container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

	cpus 5

	input:
	path genome_index
	tuple val(replicateID),path(trimmed_reads)

	output:
	tuple val(replicateID),path("*.sam")

	script:
	"""
	bwa mem -t ${task.cpus} -R '@RG\\tID:${replicateID}\\tSM:${replicateID}.SM\\tLB:${replicateID}.LB' -o ${replicateID}.sam ${genome_index[0].baseName} ${trimmed_reads}
	rm ${trimmed_reads}
	"""
}

process post_alignment_processing{
	container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

	cpus 5

	input:
	tuple val(replicateID),path(samfile)

	output:
	tuple val(replicateID),path("${replicateID}_st_dr.bam"),path("${replicateID}_st_dr.bam.bai")

	script:
	"""
	samtools fixmate -@ ${task.cpus} -O BAM ${replicateID}.sam ${replicateID}.bam
	rm ${replicateID}.sam
	samtools sort -@ ${task.cpus} -O BAM -o ${replicateID}_sorted.bam ${replicateID}.bam
	rm ${replicateID}.bam
	samtools rmdup ${replicateID}_sorted.bam ${replicateID}_st_dr.bam 
	rm ${replicateID}_sorted.bam
	samtools index ${replicateID}_st_dr.bam 
	"""
}

process gatk_variant_calling{
	container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

	cpus 2

	input:
	tuple val(replicateID),path(bamfile),path(bamindex)
        path ref_fasta
        path ref_index
        path ref_dict

	output:
	tuple val(replicateID),path("${replicateID}.g.vcf")

	script:
	"""
	gatk HaplotypeCaller -R ${ref_fasta} -I ${bamfile} -O ${replicateID}.g.vcf -ERC GVCF
	"""

}

process gatk_joint_genotyping{
	container 'quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1'

	cpus 2
	memory '12 GB'	

	input:
	val cohort_name
        path sample_map
        path ref_fasta
        path ref_index
        path ref_dict
	val interval_name

	output:
        tuple val(interval_name),path("${cohort_name}_${interval_name}.vcf"),path("${cohort_name}_${interval_name}.vcf.idx")

	script:
	"""
	gatk GenomicsDBImport --sample-name-map ${sample_map} --genomicsdb-workspace-path ${cohort_name}_${interval_name}_db -L ${interval_name}

    	gatk GenotypeGVCFs -R ${ref_fasta} -V gendb://${cohort_name}_${interval_name}_db -O ${cohort_name}_${interval_name}.vcf -L ${interval_name}

	"""
}


process signature_profiling{
	container 'chesterli214/gsexp_amd64'
	
	cpus 4

	input:
        tuple val(interval_name),path(vcf_file),path(vcf_index)
	path group_label

	output:
	tuple val(interval_name),path("${params.cohort_name}_${interval_name}_signature.txt")

	script:
	"""
	GSEXP SigFreq -t ${task.cpus} -i ${vcf_file} -o ${params.cohort_name}_${interval_name}_signature.txt -p ${group_label} --var-type biSNP --group-num 1 --tar-lower 0.9 --ref-upper 0.1 --min-sample 1
	"""
}


workflow{
        
	//***************************************************************
	//Generating genome indices

        prepare_genome_samtools("${params.genome_dir}/${params.genome}")

	prepare_genome_picard("${params.genome_dir}/${params.genome}")

	prepare_genome_bwa("${params.genome_dir}/${params.genome}")

	extract_chromosome_info("${params.genome_dir}/${params.genome}")

	//***************************************************************
	//Prepare input data

	reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
	reads_ch.view()
	intvl_ch=extract_chromosome_info.out.splitText().map{it->it.trim()}
	intvl_ch.view()

	//***************************************************************
	//Preprocessing

	reads_quality_control(reads_ch)
	
	get_trimming_parameters(reads_quality_control.out)
	//get_trimming_parameters.out.view()

	reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)

	quality_trimming(reads_ch,get_trimming_parameters.out)
	//quality_trimming.out.view()

	//***************************************************************
	//Reads alignment

	reads_mapping(prepare_genome_bwa.out,quality_trimming.out)

	//***************************************************************
	//Post processing
	
	post_alignment_processing(reads_mapping.out)

	//***************************************************************
	//Variant calling
	
	gatk_variant_calling(post_alignment_processing.out,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out)

	//Generate sample map

	sample_map = gatk_variant_calling.out.collectFile(){id,gvcf->["${params.cohort_name}_map.tsv","${id}\t${gvcf}\n"]}.first()

	//***************************************************************
	//Joint variant calling	
	
	gatk_joint_genotyping(params.cohort_name,sample_map,"${params.genome_dir}/${params.genome}",prepare_genome_samtools.out,prepare_genome_picard.out,intvl_ch)
	gatk_joint_genotyping.out.view()

	//***************************************************************
	//Genomic signature profiling 

	signature_profiling(gatk_joint_genotyping.out,params.group_label_file)
	signature_profiling.out.view()

	
}


