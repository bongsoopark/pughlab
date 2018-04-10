# NextSeq FastQ file check script
# Written by Bongsoo Park
# Pugh Lab, Center for Eukaryotic Gene Regulation
# Penn State University, 2016

import sys
import os

# The parameter is the RunID (PEGR)
run_id = sys.argv[1]

# Prep dir - New Galaxy/PEGR pipeline
# Old dir - Old Shell Script pipeline

prep_dir = "/gpfs/cyberstar/pughhpc/galaxy-cegr/files/prep/prep_dir"
old_dir = "/storage/home/bxp12/storage/illumina/illuminaNextSeq/NSQData_PughLab"

# Retrieve the folder name
f = open("nextseq_run_list.csv","r") 
folder_name = {}

for line in f:
    line = line.strip()
    data = line.split(",")
    print data
    folder_name.update({data[0]:data[1]})
f.close()

f = open("check_fastq_files.txt","r")
new_fastq = []
old_fastq = []
for line in f:
    line = line.strip()
    data = line.split(",")
    # Run164-216 were curated.
    if data[0] != "run_id" and data[0] == run_id and int(run_id) >= 164:
        # Copy new FastQ files 
        os.system("cp "+prep_dir+"/"+folder_name[data[0]]+"/"+data[0]+"-"+data[2]+"_S"+data[1]+"_R1_001.fastq.gz .\n")
        os.system("gunzip "+data[0]+"-"+data[2]+"_S"+data[1]+"_R1_001.fastq.gz\n")
        os.system("cp "+prep_dir+"/"+folder_name[data[0]]+"/"+data[0]+"-"+data[2]+"_S"+data[1]+"_R2_001.fastq.gz .\n")
        os.system("gunzip "+data[0]+"-"+data[2]+"_S"+data[1]+"_R2_001.fastq.gz\n")
        # Copy old fastQ files
        for i in range(1,5):
            os.system("cp "+old_dir+"/"+folder_name[data[0]]+"/FastQ/"+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R1_001.fastq.gz .\n") 
            os.system("gunzip "+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R1_001.fastq.gz\n") 
            os.system("cp "+old_dir+"/"+folder_name[data[0]]+"/FastQ/"+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R2_001.fastq.gz .\n") 
            os.system("gunzip "+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R2_001.fastq.gz\n") 
        # Remove temporary R1, R2 fastq files
        os.system("rm -rf tmp*.fq\n")

        # Concatenate four separate files into one fastq file (READ1, READ2)
        for i in range(1,5):
            os.system("cat "+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R1_001.fastq >> tmp_R1.fq\n") 
            os.system("cat "+data[3]+"_S"+str(int(data[4]))+"_L00"+str(i)+"_R2_001.fastq >> tmp_R2.fq\n") 
        # Create output file
        os.system("echo 'SampleID: "+data[2]+","+data[3]+"' >> check_fastq_files.output\n")
        # Count line numbers
        try:
            with open(data[0]+"-"+data[2]+"_S"+data[1]+"_R1_001.fastq") as myfile:
                new_read1 = sum(1 for line in myfile if line.rstrip('\n'))
            with open("tmp_R1.fq") as myfile:
                old_read1 = sum(1 for line in myfile if line.rstrip('\n'))
            if new_read1 == old_read1:
                os.system("echo 'Success: READ1 is matched' >> check_fastq_files.output\n")
            else:
                os.system("echo 'Failed: READ1 is not matched ("+str(new_read1)+" vs "+str(old_read1)+")' >> check_fastq_files.output\n")
        except:
            os.system("echo 'Failed: READ1 file missing' >> check_fastq_files.output\n")

        try:
            with open(data[0]+"-"+data[2]+"_S"+data[1]+"_R2_001.fastq") as myfile:
                new_read2 = sum(1 for line in myfile if line.rstrip('\n'))
            with open("tmp_R2.fq") as myfile:
                old_read2 = sum(1 for line in myfile if line.rstrip('\n'))
            if new_read2 == old_read2:
                os.system("echo 'Success: READ2 is matched' >> check_fastq_files.output\n")
            else:
                os.system("echo 'Failed: READ2 is not matched ("+str(new_read1)+" vs "+str(old_read1)+")' >> check_fastq_files.output\n")
        except:
            os.system("echo 'Failed: READ2 file missing' >> check_fastq_files.output\n")

f.close()
