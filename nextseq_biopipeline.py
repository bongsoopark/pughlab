# Title: NextSeq Mapping and Downstream Pipeline in Pugh Lab
# Center for Eukaryotic Gene Regulation
# Department of Biochemistry and Molecular Biology
# The Genome Institutes, Penn State University, USA
# The version update will be archived in PSU GitHub (git.psu.edu)
# All right reserved 2016, PSU Bioinformaitcs & Genomics

# NextSeq automated pipeline was developed by Bongsoo Park, iGenomics Initiative
# Research Infrastructure Team (CEGR)
# Shaun Mahony, Bongsoo Park, William Lai, and Gretta Armstrong
# Import Python Libraries

import os, sys
import os.path
import time
import MySQLdb
import csv
import argparse

# Main biopipeline script
# This script will generata 48-96 pbs scripts for submitting jobs into LionX clusters
# The default parameters will be used for the automated pipeline.

def __main__():
	parser = argparse.ArgumentParser()
	
	# The first parameters were designed in Verion 0.1 (2014)
	parser.add_argument('--configure', dest='configure', action="store")
	parser.add_argument('--nextseq', dest='nextseqname', action="store")
	parser.add_argument('--source', dest='source', action="store")
	parser.add_argument('--readtype', dest='readtype', action="store")
	parser.add_argument('--fastqc', dest='fastqc', action="store")
	parser.add_argument('--mapping', dest='mapping', action="store")
	parser.add_argument('--mappingstat', dest='mappingstat', action="store")
	parser.add_argument('--analysis', dest='analysis', action="store")
	parser.add_argument('--errorcheck', dest='errorcheck', action="store")
	parser.add_argument('--samplesheet', dest='samplesheet', action="store")
	
	# The below parameters were added in Verion 0.3 (2015)
	parser.add_argument('--rmdup', dest='rmdup', action="store")
	parser.add_argument('--reanalysis', dest='reanalysis', action="store")
	parser.add_argument('--reffeature', dest='reffeature', action="store")
	parser.add_argument('--dataset', dest='dataset', action="store")
	parser.add_argument('--nonuniq', dest='nonuniq', action="store")
	args = parser.parse_args()

	# Core Pipeline v0.1-v0.3 : Default parameters
	configure = args.configure
	nextseqname = args.nextseqname
	source_command = args.source
	readtype = args.readtype
	fastqc = args.fastqc
	mapping = args.mapping
	mappingstat = args.mappingstat
	analysis = args.analysis
	errorcheck = args.errorcheck
	samplesheet = args.samplesheet
	
	# Core Pipeline v0.2 : Variable parameters
	f = open(configure,"r")
	variable_parameters = {}
	for line in f:
		line = line.strip()
		data = line.split("\t")
		variable_parameters.update({data[0]:data[1]})	
	f.close()

	# Assignment to Variable Parameters!
	for the_key in variable_parameters:
		exec(the_key+" = variable_parameters[\""+the_key+"\"]")

	# Core Pipeline v0.3 : Remove duplicate reads and ReAnalysis option	
	rmdup = args.rmdup
	reanalysis = args.reanalysis
	reffeature = args.reffeature
	dataset = args.dataset
	nonuniq = args.nonuniq

	# The default source files are FastQ files from NextSeq
	if source != source_command:
		source = source_command

	fn = open(samplesheet,"r")

	# The introduction of the scripts below
	print "============================================================================="
	print " NextSeq Mapping and Downstream Pipeline, Pugh Lab"
	print " Center for Eukaryotic Gene Regulation"
	print " Department of Biochemistry and Molecular Biology"
	print " The Pennsylvania State University, 2014-2016"
	print " Motivated and Inspired by Yeast, Mouse, Human ENCODE Projects"
	print "----------------------------------------------------------------------------"
	print " Developed by Bongsoo Park, iGenomics Initiative"
	print " Designed by Research Infrastructure Team (Shaun, Bongsoo, Will, and Gretta)"
	print " Version "+version
	print " Last update : "+last_update
	print "============================================================================"

	# The sample names are stored in the run_name_array
	# Initialize the variables
	run_name_array = []
	sample_cnt = 0
	current_directory = os.getcwd()

	# Archive Sample Counts into the sample archive
	# The running time on the cluster will be modified if the number of tags is huge!
	# In NextSeq 500 sequencer, the maximum tags would be closed to 400M tags
	try:
		f_cnt = open("SampleCount.csv","r")
		sample_cnt_list = {}
		for ele in f_cnt:
			ele = ele.strip()
			tmp = ele.split("\t") 
			sample_cnt_list.update({tmp[0]:int(tmp[1])})
		f_cnt.close()
	except:
		# If there no SampleCount.csv information, 10M tags will be assigned.
		f_cnt = open(samplesheet,"r")
		sample_cnt_list = {}
		for ele in f_cnt:
			ele = ele.strip()
			tmp = ele.split(",")
			sample_cnt_list.update({tmp[0]:10000000})
		f_cnt.close()	

	# Retrieve mysql key : pughlab db in Francline VM Host
	# MySQL connection will be assigned to variables.
	# The summary report will be updated after completing the downstream analysis
	f_sql = open(script_path+"/mysql_key.txt","r")
	mysql_key = []
	for ele in f_sql:
		ele = ele.strip()
		mysql_key.append(ele)
	f_sql.close()

	# Followed the lines of SampleSheet.csv file 0.1-0.2
	# In Version0.3, SampleMappingSheet.csv will be loaded for secondary mappings.
	# The Ref1, Ref2 were assigned in the previous pipeline (fastq_ready_ref2.py).
	for line in fn:
		line = line.strip()
		# Three columns (ID, SampleID, Index sequences v0.1-0.2|Ref1,or 2 v0.3-0.4)
		data = line.split(",")
	
		# The first and second lines will be ignored. 
		# The header line will be excluded.
		if data[0] != "[Data]" and data[0] != "SampleID":
			ele = data[0]
			run_name = data[1]
			uniq_id = data[1][:5]
			# The first five digits is PughLab uniq ID -SeqExperiment ID
			ref_genome = data[1][5:]
			ref_genome_id = data[2]
			
			# Assay ID has been assigned. XO - Deduplication applied
			# Other assays are not applied
			try:
				assay_id = data[3]
				if assay_id == "MN" or assay_id == "IN":
					rmdup_flag = "No"
				else:
					rmdup_flag = "Yes"
			except:
				assay_id = "NA"
				rmdup_flag = rmdup
			run_name_array.append(run_name)
			sample_cnt = sample_cnt + 1
						
			try:
				# source:fastq
				run_name2 = data[3]
			except:
				run_name2 = "No READ2"

			print "INFO :: "+run_name+".pbs script has been generated."
			
			# walltime and memory will be assigned. (pbs script)
			# francline db will be served as HPC information archive for each genome
			db = MySQLdb.connect(host=mysql_key[0], user=mysql_key[1], passwd=mysql_key[2], db=mysql_key[3])
			cursor = db.cursor()
			try:
				cursor.execute("SELECT HPC_INFO, GENOME_SIZE from PughLabRefGenomeInfo where NAME = '"+ref_genome+"';")
				fetch_data = cursor.fetchone()
				hpc_info = fetch_data[0].split(",")
				gsize = fetch_data[1]
				walltime = hpc_info[0]
				pmem = hpc_info[1]
				ppn = hpc_info[2]
			# If there is no hpc wall time information, the default number will be assigned. phiX mapping is assigned with no Index.
			except:
				if ref_genome == "phiX":
					walltime = "1"
					pmem = "1gb"
					ppn = "1"
					gsize = str(100000)
				else:
					walltime = "24"
					pmem = "8gb"
					ppn = "4"
					gsize = str(1000000)
			db.close()

			# If the number of tags is more than 60 millon, the wall time will be 36 hours
			try:
				if sample_cnt_list[run_name] > 60000000:
					walltime = "36"
					pmem = "10gb"
					ppn = "4"
			except:
				tmp = 0
			# ppn=4 is the optimized for the mapping, but it will increase the queue time
			# ppn=4 is only for mouse and human genomes
			# ppn=2 was assigned. 
			pbs_script = open(run_name+".pbs","w")
			if reanalysis=="yes":
				# Francline will support the ReAnalysis function
				pbs_script.write("#PBSFREE: ReAnalysis\n");
			pbs_script.write("#PBS -l walltime="+walltime+":00:00\n")
			pbs_script.write("#PBS -l nodes=1:ppn="+ppn+"\n")
			pbs_script.write("#PBS -l pmem="+pmem+"\n")
			pbs_script.write("#PBS -j oe\n")
			pbs_script.write("#PBS -q "+server+"-pughhpc\n")
			pbs_script.write("\n")
			pbs_script.write("# Pugh Lab, Mapping and Downstream Pipeline : "+version+"\n\n")
			pbs_script.write("cd "+current_directory+"\n")
			pbs_script.write("\n")
			
			# The available modules for the pipeline
			pbs_script.write("# Step 0 : Module Loading\n")
			pbs_script.write("module load bedtools/2.22.1\n")
			pbs_script.write("module load samtools/1.1\n")
			pbs_script.write("module load java/1.7.0_21\n")
			pbs_script.write("echo 'Step 0 : Module load completed...'\n")
			pbs_script.write("\n")
			
			# Fastq.gz to FastQC report
			if fastqc == "yes":
				pbs_script.write("# Step 1 : FastQC Report\n")
				if ref_genome_id == "Ref1":
					pbs_script.write(python_path+" "+script_path+"/fastq_merge.py "+ele+" "+uniq_id+"\n")
					pbs_script.write(python_path+" "+script_path+"/fastq_report.py "+ele+" "+uniq_id+"\n")
					pbs_script.write("echo 'Step 1 : FastQ reports completed...'\n")
					pbs_script.write("\n# FastQ Reports have been generated\n\n")
				else:
					pbs_script.write("\n# FastQ Reports process will be skipped in the second reference genome\n\n")
			if mapping == "yes":
				pbs_script.write("# Step 2 : Mapping - BWA-MEM 0.7.12 as a default option\n")
				mapper_name = mapper.split("-")
				
				# Currently, BWA and BOWTIE2 are available for mappers
				# In version 0.3, the spike-in control samples will be added.
				# Only BWA will support Spike-in Control samples
				# Picard MarkDuplicate will be applied.
				if readtype == "SR":
					if mapper_name[0] == "bwa":
						# Spike-In Control by Ecoli samples
						pbs_script.write(share_path+"/"+mapper+"/bwa mem -M -t "+ppn+" "+genome_path+"/recommend_"+mapper+"/ec2/ec2_all.fa "+uniq_id+".fq > "+uniq_id+"ec2.sam\n")
						pbs_script.write("samtools view -Sbh "+uniq_id+"ec2.sam > "+uniq_id+"ec2.bam\n")
						pbs_script.write("samtools sort "+uniq_id+"ec2.bam "+uniq_id+"ec2.sorted\n")
						# Real genome mappings
						pbs_script.write(share_path+"/"+mapper+"/bwa mem -M -t "+ppn+" "+genome_path+"/recommend_"+mapper+"/"+ref_genome+"/"+ref_genome+"_all.fa "+uniq_id+".fq > "+run_name+".sam\n")
					elif mapper_name[0] == "bowtie2":
						pbs_script.write(share_path+"/"+mapper+"/bowtie2 -p "+ppn+" "+genome_path+"/recommend_"+mapper+"/"+ref_genome+"/"+ref_genome+" "+uniq_id+".fq -S "+run_name+".sam\n")
					else:
						print "Error"
						exit(0)	
				elif readtype == "PE":
					if mapper_name[0] == "bwa":
						# Spike-In Control by Ecoli samples
						pbs_script.write(share_path+"/"+mapper+"/bwa mem -M -t "+ppn+" "+genome_path+"/recommend_"+mapper+"/ec2/ec2_all.fa "+uniq_id+"_R1.fq "+uniq_id+"_R2.fq > "+uniq_id+"ec2.sam\n")
						pbs_script.write("samtools view -Sbh "+uniq_id+"ec2.sam > "+uniq_id+"ec2.bam\n")
						pbs_script.write("samtools sort "+uniq_id+"ec2.bam "+uniq_id+"ec2.sorted\n")
						# Real genome mappings
						pbs_script.write(share_path+"/"+mapper+"/bwa mem -M -t "+ppn+" "+genome_path+"/recommend_"+mapper+"/"+ref_genome+"/"+ref_genome+"_all.fa "+uniq_id+"_R1.fq "+uniq_id+"_R2.fq > "+run_name+".sam\n")
					elif mapper_name[0] == "bowtie2":
						pbs_script.write(share_path+"/"+mapper+"/bowtie2 -p "+ppn+" "+genome_path+"/recommend_"+mapper+"/"+ref_genome+"/"+ref_genome+" -1 "+run_name+" -2 "+run_name2+" -S "+run_name+".sam\n")

				else:
					print "Error"
				pbs_script.write("samtools view -Sbh "+run_name+".sam > "+run_name+".origin.bam\n")

				pbs_script.write("\n#BWA MEM Mapping has been completed.\n")
				#Uniquely mapped reads
				pbs_script.write("# Step 3 : Sorting and Extract the uniquely mapped region\n")
				pbs_script.write("mkdir -p TotBam\n")
				
				#Picard input file generation by using samtools 1.2
				pbs_script.write("samtools sort "+run_name+".origin.bam "+run_name+".sorted\n")
				
				#Uniquely mapped reads will be retrieved from the sorted bam file
				pbs_script.write("samtools view -bh -F 4 -q 1 "+run_name+".sorted.bam > "+run_name+".uniq.bam\n")
				pbs_script.write("samtools view -bh -F 4 -q 5 "+run_name+".sorted.bam > "+run_name+".bam\n")
				pbs_script.write("echo 'Step 3: Sorting bam file and extract uniquely mapped bam completed...'\n")
				pbs_script.write("\n\n")
				
				# Extract BED file, however there was an issue on the memory size in this process
				pbs_script.write("# Step 4 : Deduplication process\n")
				# Check list #1
				# -b : bam format
				# -h : header
				# -f 0x2 : properly paired read only
				# -f 0x40 : READ1 only
				if readtype == "PE":
					# Mark Duplicate using PiCard (0.3 pipeline)
					# 0x400 (PCR Duplicates Mark)
					# 0x2 (Properly paired reads)
					# MAQ>20 filter for HiSeq Single Read
					# MAQ>5 filter for NextSeq Pair-End Reads
					pbs_script.write("# 4-1 : Deduplicated bam, bed files by PiCard\n")
					pbs_script.write("java -Xmx3g -jar "+picard+"picard.jar MarkDuplicates I="+run_name+".sorted.bam"+" O="+run_name+".picard.bam M="+run_name+".picard.txt\n")
					pbs_script.write("java -Xmx3g -jar "+picard+"picard.jar MarkDuplicates I="+uniq_id+"ec2.sorted.bam"+" O="+uniq_id+"ec2.dedup.bam M="+uniq_id+"ec2.picard.txt REMOVE_DUPLICATES=true\n")
					pbs_script.write("samtools view -bh -F 0x400 -f 0x2 "+run_name+".picard.bam > "+run_name+".dedup.bam\n")
					pbs_script.write("samtools view -bh -F 4 -f 0x40 -q 5 "+run_name+".dedup.bam | bedtools bamtobed -i stdin > "+run_name+".dedup.bed\n")
					pbs_script.write("\n")
					pbs_script.write("# 4-2 : Deduplicated bam, bed files by Samtools\n")
					pbs_script.write(samtools+" rmdup "+run_name+".uniq.bam "+run_name+".rmdup.bam\n")
					pbs_script.write("samtools view -bh -F 4 -f 0x42 -q 5 "+run_name+".rmdup.bam | bedtools bamtobed -i stdin > "+run_name+".rmdup.bed\n")
					
					pbs_script.write("samtools view -bh -F 4 -f 0x40 -q 5 "+run_name+".bam | bedtools bamtobed -i stdin > "+run_name+".temp.bed\n")
					pbs_script.write("\n")

				elif readtype == "SR":
					# MAQ>20 filter
					pbs_script.write("samtools view -bh -F 4 -q 20 "+run_name+".bam | bedtools bamtobed -i stdin > "+run_name+".temp.bed\n")
				else:
					print "Error : Check your read type (PE, or SR)"
					exit(0)
				# Step 4 
				pbs_script.write("echo 'Step 4 : BAM to BED completed...'\n")
				pbs_script.write("\n\n")
				
				# Non-uniquely mapped reads - multi.bam
				pbs_script.write("# Step 5 : Non-uniquely mapped reads check\n")
				pbs_script.write("samtools view -F 4 -h "+run_name+".sorted.bam | awk '{if($1~/^@/ || $5==0){print $0}}' > "+run_name+".multi.sam\n")
				pbs_script.write("samtools view -Sbh "+run_name+".multi.sam > "+run_name+".multi.bam\n")
				pbs_script.write("samtools view -bh "+run_name+".multi.bam | bedtools bamtobed -i stdin > "+run_name+".multi.bed\n")
				# Step 5 
				pbs_script.write("echo 'Step 5 : Non-uniquely mapped reads check completed...'\n")
				pbs_script.write("\n\n")

				# BAM to TAB converting
				pbs_script.write("# Step 6 : Index bam files and BAM to SCIDX converting using ScriptManager\n")
				pbs_script.write("samtools index "+run_name+".bam\n")
				pbs_script.write("samtools idxstats "+run_name+".bam > "+run_name+".bam.idxstats\n")
				pbs_script.write("java -Xmx2048M -Xms512M -jar "+script_path+"/9_ScriptManager/PEHistogram.jar -B "+run_name+".bam\n")
				pbs_script.write("java -Xmx2048M -Xms512M -jar "+script_path+"/9_ScriptManager/BAMtoIDX.jar -b "+run_name+".bam -i "+run_name+".bam.bai -o "+run_name+".idx\n")
				
				# Step 6 
				pbs_script.write("echo 'Step 6 : BAM to SCIDX converting has been completed...'\n")
				pbs_script.write("\n\n")
				
				# Downstream black list filter - dm3, mm9, mm10, hg19
				pbs_script.write("# Step 7: Black list filter process\n")
				if ref_genome == "hg19" or ref_genome == "mm10" or ref_genome == "mm9" or ref_genome == "dm3" or ref_genome == "sacCer3":
					pbs_script.write("# Detected genome : "+ref_genome+"\n")
					if rmdup == "no" or rmdup_flag == "no":
						pbs_script.write("bedtools intersect -v -a "+run_name+".temp.bed -b "+script_path+"/x_Blacklist/"+ref_genome+"-blacklist.bed > "+run_name+".bed\n")
					else:
						pbs_script.write("bedtools intersect -v -a "+run_name+".dedup.bed -b "+script_path+"/x_Blacklist/"+ref_genome+"-blacklist.bed > "+run_name+".bed\n")
				else:
					if rmdup == "no" or rmdup_flag == "no":
						pbs_script.write("mv "+run_name+".temp.bed "+run_name+".bed\n")
					else:
						pbs_script.write("mv "+run_name+".dedup.bed "+run_name+".bed\n")
				
				# Repeated Masked region reads will be tested.
				if nonuniq == "yes":
					pbs_script.write("rm -rf "+run_name+".bed\n")			
					pbs_script.write("mv "+run_name+".multi.bed "+run_name+".bed\n")			
				
				# Reference feature (promoter) will be tested.
				if reffeature == "promoter":
					if ref_genome == "hg19" or ref_genome == "mm10" or ref_genome == "sacCer3":
						pbs_script.write("mv "+run_name+".bed "+run_name+".uniq.bed\n")			
						pbs_script.write("bedtools intersect -a "+run_name+".uniq.bed -b "+script_path+"/w_GenomicFeatures/"+ref_genome+"_TSS_promoter.bed > "+run_name+".bed\n")
				pbs_script.write("echo 'Step 7 : Black list filter completed...'\n")
				pbs_script.write("\n\n")
				
				# BED(0 based) to TAB(1 based) converting
				pbs_script.write("# Step 8 : BED to TAB (READ1 only, properly paired reads)\n")
				pbs_script.write(python_path+" "+share_path+"/tab2genetrack/tabs2genetrack.py -s 1 -i "+run_name+".bed -f bed -o "+run_name+".tab\n")
				if errorcheck == "yes":
					pbs_script.write("echo 'Step 8 : Format converting completed...'\n")
					pbs_script.write("\n# Error Reporting : Format converting\n")				
					pbs_script.write(python_path+" "+script_path+"/error_reporting.py --step 2 --sampleid "+run_name+"\n")
				pbs_script.write("\n")
				
				#Mapping Statistics
				#[0] Index count
				#[1] Mapped count
				#[2] Uniquely mapped count
				#[3] Mapped coordinate
				#[4] Rmdup by Samtools count
				#[5] Non-uniquely mapped count
				#[6] Spike-in Tag count
				#[7] Adapter dimer tag count 
				#[8] Deduplicated tags of READ1 and READ2
				#[9] Uniquely Mapped, but non-properly paired
				#[10] Uniquely Mapped, but duplicated
				#[11] Uniquely Mapped, deduplicated
				pbs_script.write("# Step 9 : Generate the basic mapping statistics (READ1)\n")
				if readtype == "PE":
					pairend_option = " -f 0x40"
				else:
					pairend_option = ""
				pbs_script.write("samtools view"+pairend_option+" -c "+run_name+".sorted.bam > "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view"+pairend_option+" -F 4 -c "+run_name+".sorted.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view"+pairend_option+" -F 4 -q 1 -c "+run_name+".sorted.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write("wc -l "+run_name+".tab | awk '{print $1}' >> "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view -c "+pairend_option+" -c "+run_name+".rmdup.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view -F 4 "+pairend_option+" -c "+run_name+".multi.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view -F 4 -q 1 "+pairend_option+" -c "+uniq_id+"ec2.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write(python_path+" "+script_path+"/fastq_parse.py ./"+uniq_id+"_R1/"+uniq_id+"_R1.fq_fastqc.html >> "+run_name+"_stat.txt\n")
				pbs_script.write("samtools view "+pairend_option+" -c "+run_name+".dedup.bam >> "+run_name+"_stat.txt\n")
				# CorePipeline v0.4
				if readtype == "PE":
					pbs_script.write("######\n")
					#Unpaired read check script
					#Uniquely Mapped, but non-properly paired
					pbs_script.write("samtools view -F 0x6 -f 0x40 -q 1 -c "+run_name+".picard.bam >> "+run_name+"_stat.txt\n")
					#Uniquely Mapped, but duplicated
					pbs_script.write("samtools view -F 4 -f 0x442 -q 1 -c "+run_name+".picard.bam >> "+run_name+"_stat.txt\n")
					#Uniquely Mapped, deduplicated
					pbs_script.write("samtools view -F 0x404 -f 0x42 -q 1 -c "+run_name+".picard.bam >> "+run_name+"_stat.txt\n")
				pbs_script.write("echo 'Step 7 : The basic mapping stat completed...'\n")
				pbs_script.write("\n\n")
				
				#Index bam files and generate the historgram
				pbs_script.write("# Step 10: Mapping results (BAM,TAB) into the group folder\n")
				if source == "nextseq":
					# Bam files, uniquely mapped reads MAQ>5
					pbs_script.write("cp "+run_name+".bam "+mapping_result+"/"+nextseqname+"/.\n")
					# Tab files, uniquely mapped read1 MAQ>5
					pbs_script.write("cp "+run_name+".idx "+mapping_result+"/"+nextseqname+"/"+run_name+".tab\n")
				pbs_script.write("echo 'Step 10 : Mapping results moving has been completed...'\n")
				pbs_script.write("\n\n")
			
			# The below was a part of the downstream pipeline in the previous version.
			if analysis == "yes":
				# Shift tag option check
				if int(shift_tag_xo) > 0:
					#Peak Calling using GeneTrack
					pbs_script.write("# Shift tag option\n\n")
					pbs_script.write(python_path+" "+script_path+"/shift_tag.py "+run_name+".tab "+str(shift_tag_xo)+" > "+run_name+".tmp\n")
					pbs_script.write("mv "+run_name+".tab "+run_name+".tab2\n")
					pbs_script.write("mv "+run_name+".tmp "+run_name+".tab\n\n")
				
				# Peak Calling using GeneTrack
				if genetrack_param1_F == "0":
					genetrack_param1 = "s"+genetrack_param1_s+"e"+genetrack_param1_e
				else:
					genetrack_param1 = "s"+genetrack_param1_s+"e"+genetrack_param1_e+"F"+genetrack_param1_F
				pbs_script.write("# Step 11: Peak Calling by GeneTrack\n")
				pbs_script.write("cd "+current_directory+"\n")
				pbs_script.write("# Fine grain - default:s5e10F0\n")
				pbs_script.write(python_path+" "+script_path+"/pugh_scripts/genetrack.py -s "+genetrack_param1_s+" -e "+genetrack_param1_e+" -F "+genetrack_param1_F+" "+run_name+".tab\n")
				pbs_script.write("# Course grain - default:s20e40F0\n")
				pbs_script.write(python_path+" "+script_path+"/pugh_scripts/genetrack.py -s "+genetrack_param2_s+" -e "+genetrack_param2_e+" -F "+genetrack_param2_F+" "+run_name+".tab\n")
				pbs_script.write("mkdir -p "+run_name+"_downstream\n")
				pbs_script.write("cp "+current_directory+"/genetrack_s"+genetrack_param1_s+"e"+genetrack_param1_e+"/"+run_name+"*.gff "+current_directory+"/"+run_name+"_downstream/.\n")
				pbs_script.write("cp "+current_directory+"/genetrack_s"+genetrack_param2_s+"e"+genetrack_param2_e+"/"+run_name+"*.gff "+current_directory+"/"+run_name+"_downstream/.\n")
				pbs_script.write("echo 'Step 11 : Peak calling completed...'\n")
				pbs_script.write("\n\n")

				# Singleton filter - all genomes
				pbs_script.write("# Step 12: Singleton peak filter\n")
				pbs_script.write("cd "+run_name+"_downstream\n")
				pbs_script.write(python_path+" "+script_path+"/singleton_filter.py\n")
				pbs_script.write("echo 'Step 12 : Singleton filter completed...'\n")
				pbs_script.write("\n\n")

				# Peak-paring step	
				pbs_script.write("# Step 13 : Peak-pairing process\n")
				pbs_script.write(python_path+" "+script_path+"/pugh_scripts/cwpair2.py -q -f "+cwpair_param_f+" -u "+cwpair_param_u+" -d "+cwpair_param_d+" -b "+cwpair_param_b+" -c asc -s desc .\n")
				pbs_script.write("mkdir -p peak_signals\n")
				pbs_script.write("mkdir -p peak_pair_signals\n")
				pbs_script.write("mv *.gff ./peak_signals/.\n")
				pbs_script.write("cp ./cwpair_output_mode_f"+cwpair_param_f+"u"+cwpair_param_u+"d"+cwpair_param_d+"b"+cwpair_param_b+"/S_"+run_name+"_s*_??S.gff ./peak_pair_signals/.\n")
				if errorcheck == "yes":
					pbs_script.write("\n# Error Reporting : Peak Calling - Pairing\n")				
					pbs_script.write(python_path+" "+script_path+"/error_reporting.py --step 3 --sampleid "+run_name+"\n")
				pbs_script.write("echo 'Step 13 : Peak-paring step completed...'\n")
				pbs_script.write("\n")
				
				# Genomic regions extraction step
				pbs_script.write("# Step 14 : Genomic region extraction\n")
				if ref_genome == "hg19" or ref_genome == "mm10" or ref_genome == "sacCer3":
					pbs_script.write("cd peak_pair_signals\n")
					pbs_script.write("bedtools intersect -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/w_GenomicFeatures/"+ref_genome+"_TSS_proximal.bed > S_"+run_name+"_TSS_proximal.gff\n")
					pbs_script.write("bedtools intersect -v -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/w_GenomicFeatures/"+ref_genome+"_TSS_proximal.bed > S_"+run_name+"_TSS_distal.gff\n")
					pbs_script.write("bedtools intersect -v -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/x_Blacklist/"+ref_genome+"-RepeatMasker.bed > S_"+run_name+"_Repeat_masked.gff\n")
					# Yeast ENCODE dataset for Telomere masking
					if ref_genome == "sacCer3" or ref_genome == "hg19":
						pbs_script.write("bedtools intersect -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/w_GenomicFeatures/"+ref_genome+"_ORF.gff > S_"+run_name+"_ORF_tag.gff\n")
						pbs_script.write("bedtools intersect -v -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/w_GenomicFeatures/"+ref_genome+"_ORF.gff > S_"+run_name+"_Intergenic_tag.gff\n")
						pbs_script.write("bedtools intersect -v -a S_"+run_name+"_s5e10_noS.gff -b "+script_path+"/x_Blacklist/"+ref_genome+"-TelomereMasker.bed > S_"+run_name+"_Telomere_masked.gff\n")
					pbs_script.write("cd ..\n")
				pbs_script.write("echo 'Step 14 : Genomic region extraction completed...'\n")
				pbs_script.write("\n")
				
				# Fasta extraction 
				pbs_script.write("# Step 15 : Fasta Extraction Step\n")
				pbs_script.write(python_path+" "+script_path+"/pugh_scripts/fastaextract.py -u "+fastaextract_param_u+" -d "+fastaextract_param_d+" -q -g "+genome_path+"/recommend_"+mapper+"/"+ref_genome+"/"+ref_genome+"_all.fa ./peak_pair_signals\n")
				pbs_script.write("echo 'Step 15 : Fasta Extract completed...'\n")
				pbs_script.write("\n\n")
				
				# MEME & TOMTOM Analysis - noS & wiS
				pbs_script.write("# Step 16 : MEME & TOMTOM Analysis\n")
				pbs_script.write("cd ./peak_pair_signals/fastaextract-u"+fastaextract_param_u+"d"+fastaextract_param_d+"\n")
				pbs_script.write("rm -rf *.gff\n")
				pbs_script.write(python_path+" "+script_path+"/split_top500_fa.py --configure "+configure+" --option 2 --meme_param_top "+meme_param_top+" --meme_param_min "+meme_param_min+" --meme_param_max "+meme_param_max+" --meme_param_site "+meme_param_site+" --tomtom_db "+tomtom_db+"\n")
				pbs_script.write("sh run.sh\n")
				pbs_script.write("rm -rf run.sh\n")
				if errorcheck == "yes":
					pbs_script.write("\n# Error Reporting : MEME Analysis\n")				
					pbs_script.write(python_path+" "+script_path+"/error_reporting.py --step 4 --sampleid "+run_name+"\n")
				pbs_script.write("echo 'Step 16 : MEME Suite completed...'\n")
				pbs_script.write("\n")
				pbs_script.write("\n")

				# Composite plot & HeatMap Generation (Forward : Blue, Reverse : Red)
				pbs_script.write("# Step 17 : Composite Plot & HeatMap Generation\n")
				pbs_script.write("cd "+current_directory+"/"+run_name+"_downstream\n")
				pbs_script.write("mkdir 0_Bam\n")
				pbs_script.write("mkdir 1_Bed\n")
				pbs_script.write("mkdir 2_Tags\n")
				pbs_script.write("mkdir 3_MEME\n")
				pbs_script.write("cd 3_MEME\n")
				pbs_script.write("mkdir TSS_proximal\n")
				pbs_script.write("mkdir TSS_distal\n")
				pbs_script.write("mkdir Repeat_masked\n")
				pbs_script.write("mkdir Telomere_masked\n")
				pbs_script.write("mkdir ORF_tag\n")
				pbs_script.write("mkdir Intergenic_tag\n")
				pbs_script.write("cd ../2_Tags\n")
				pbs_script.write("cp ../../"+run_name+".tab .\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/TSS_proximal/.\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/TSS_distal/.\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/Repeat_masked/.\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/Telomere_masked/.\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/ORF_tag/.\n")
				pbs_script.write("ln "+run_name+".tab ../3_MEME/Intergenic_tag/.\n")
				
				# MEME motif TagPileUp (Core Pipeline v0.6)
				# Genome Partition has been tested.
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../meme_refs -u 52 -d 52\n")	
				pbs_script.write("cd ../3_MEME/TSS_proximal\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs2 -u 52 -d 52\n")	
				pbs_script.write("cd ../TSS_distal\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs3 -u 52 -d 52\n")	
				pbs_script.write("cd ../Repeat_masked\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs4 -u 52 -d 52\n")	
				pbs_script.write("cd ../Telomere_masked\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs5 -u 52 -d 52\n")	
				pbs_script.write("cd ../ORF_tag\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs6 -u 52 -d 52\n")	
				pbs_script.write("cd ../Intergenic_tag\n")	
				pbs_script.write("perl "+script_path+"/distance_from_ref.pl -i . -r ../../meme_refs7 -u 52 -d 52\n")	
				pbs_script.write("cd ../../2_Tags\n")
				
				if reffeature == "no" or reffeature == "promoter":
					pbs_script.write("sh "+script_path+"/distance_from_ref.sh "+ref_genome+" "+tsspileup_param_u+" "+tsspileup_param_d+"\n")
					pbs_script.write(python_path+" "+script_path+"/tag_pileup.py "+ref_genome+" "+tsspileup_param_u+" "+tsspileup_param_d+" "+tsspileup_param_s+" "+configure+" "+run_name+" > run.sh\n")
				else:
					pbs_script.write("sh "+script_path+"/distance_from_ref.sh "+reffeature+" "+tsspileup_param_u+" "+tsspileup_param_d+"\n")
					pbs_script.write(python_path+" "+script_path+"/tag_pileup.py "+reffeature+" "+tsspileup_param_u+" "+tsspileup_param_d+" "+tsspileup_param_s+" "+configure+" "+run_name+" > run.sh\n")

				pbs_script.write("\n")	
				pbs_script.write("sh run.sh\n")
				pbs_script.write("mv meme_*.png ../peak_pair_signals/fastaextract-u"+fastaextract_param_u+"d"+fastaextract_param_d+"/S_"+run_name+"_"+genetrack_param1+"_noS-u"+fastaextract_param_u+"d"+fastaextract_param_d+"_Top"+meme_param_top+"/.\n")	
				pbs_script.write("echo 'Step 17 : TSS PileUp completed...'\n")
				pbs_script.write("\n")	

				# FIMO Analysis for three candidate motifs
				pbs_script.write("# Step 18 : FIMO analysis for three motifs\n")
				pbs_script.write("echo 'Step 18: FIMO analysis for three motifs'\n")
				pbs_script.write("cd "+current_directory+"/"+run_name+"_downstream/peak_pair_signals/fastaextract-u"+fastaextract_param_u+"d"+fastaextract_param_d+"\n")
				pbs_script.write("cp ./S_"+run_name+"_s"+genetrack_param1_s+"e"+genetrack_param1_e+"_noS-"+"u"+fastaextract_param_u+"d"+fastaextract_param_d+"_Top500/meme.txt .\n")
				pbs_script.write("rm -rf fimo_out_1\n")
				pbs_script.write("rm -rf fimo_out_2\n")
				pbs_script.write("rm -rf fimo_out_3\n")
				pbs_script.write("fimo --motif 1 -o fimo_out_1 meme.txt S_"+run_name+"_s"+genetrack_param1_s+"e"+genetrack_param1_e+"_noS-"+"u"+fastaextract_param_u+"d"+fastaextract_param_d+".fa\n")
				pbs_script.write("fimo --motif 2 -o fimo_out_2 meme.txt S_"+run_name+"_s"+genetrack_param1_s+"e"+genetrack_param1_e+"_noS-"+"u"+fastaextract_param_u+"d"+fastaextract_param_d+".fa\n")
				pbs_script.write("fimo --motif 3 -o fimo_out_3 meme.txt S_"+run_name+"_s"+genetrack_param1_s+"e"+genetrack_param1_e+"_noS-"+"u"+fastaextract_param_u+"d"+fastaextract_param_d+".fa\n")
				pbs_script.write("echo 'Step 18 : FIMO analysis completed...'\n\n")


				# Four Color Plot generation and Heatmap generation
				pbs_script.write("# Step 19 : Four Color Plot generation and Heatmap generation\n")
				pbs_script.write("echo 'Step 19: Four Color Plot generation and Heatmap generation'\n")
				for i in range(1,4):
					pbs_script.write("cd "+current_directory+"/"+run_name+"_downstream/peak_pair_signals/fastaextract-u"+fastaextract_param_u+"d"+fastaextract_param_d+"/fimo_out_"+str(i)+"\n")
					pbs_script.write(python_path+" "+script_path+"/generate_additional_diagrams.py "+run_name+" "+ref_genome+" > run.sh\n")
					pbs_script.write("sh run.sh\n")
				pbs_script.write("echo 'Step 19 : Four Color Plot generation and Heatmap generation...'\n\n")
			
				# Follow-up: DNA Shape and Annotation Pipeline: ToDo List

				# MultiGPS 0.5 Integrative Analysis
				# MACS2 Broad peak callings
				pbs_script.write("# Step 20: MultiGPS & MACS2 analyusis, CorePipeilne v0.6\n")
				pbs_script.write("echo 'Step 20: MultiGPS and MACS2 analyusis, CorePipeilne v0.6'\n")
				# MACS2 and MultiGPS will be applied to only yeast, mouse, human
				if ref_genome == "hg19" or ref_genome == "mm10" or ref_genome == "sacCer3":
					pbs_script.write("cd "+current_directory+"/"+run_name+"_downstream/2_Tags\n")
					pbs_script.write("echo '"+run_name+".tab\tSignal\tIDX\t"+run_name+"\t1' >> target.design\n")
					# MultiGPS 0.5 applied
					pbs_script.write("java -Xmx3G -jar ~/software/multigps/multigps_v0.5.jar --geninfo ~/software/multigps/"+ref_genome+".info --threads 4 --design target.design --verbose --gaussmodelsmoothing --gausssmoothparam 3 --out multigps_out --memepath ~/bin --mememinw "+meme_param_min+" --mememaxw "+meme_param_max+" --seq "+genome_path+"/recommend_"+mapper+"/"+ref_genome+" >multiGPS.out 2>&1\n")

					# MACS2 applied
					if ref_genome == "hg19":
						ref_genome_size = "2.7e9"
					elif ref_genome == "mm10":
						ref_genome_size = "1.8e9"
					elif ref_genome == "sacCer3":
						ref_genome_size = "1.2e7"
					else:
						ref_genome_size = "N/A"
					pbs_script.write("cd "+current_directory+"/"+run_name+"_downstream/1_Bed\n")
					pbs_script.write("ln -s ../../"+run_name+".dedup.bed .\n")
					# MACS2 applied without control sample
					pbs_script.write("/gpfs/group/pughhpc/share/bin/macs2 callpeak -t "+run_name+".dedup.bed -f BED --gsize "+ref_genome_size+" -n "+run_name+"_MACS2 --nomodel\n")

				if errorcheck == "yes":
					pbs_script.write("\n# Error Reporting : Tag PileUp\n\n")				
					pbs_script.write(python_path+" "+script_path+"/error_reporting.py --step 5 --sampleid "+run_name+"\n")
				pbs_script.write("echo 'Step 20 : MultiGPS & MACS2 completed...'\n")
				pbs_script.write("\n")
				pbs_script.write("\n")

	fn.close()

def addslashes(s):
    l = ["\\", '"', "'", "\0", ]
    for i in l:
        if i in s:
            s = s.replace(i, '\\'+i)
    return s

# Main def will be executed.
if __name__ == "__main__" : __main__()
