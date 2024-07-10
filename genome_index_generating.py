import subprocess
import os
import sys


class GenomeIndexGenerator:
    """
    In this class we generate the human and mouse index files. For RNA-seq data, we are using STAR for indexing.
    There are two input files for generating human index file: hg38.fa.gz from
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"
    and
    gencode.v36.annotation.gtf.gz downloaded from Gencode (http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz).
    We chose the file from Gencode because of the issues with counting reads with HTSeq-count with annotation files from UCSC. However, we are using Version 36 of Gencode which is the one that
    UCSC uses. See below for more information about the annotations in UCSC:
    http://genome.ucsc.edu/cgi-bin/hgTables with this set up:
    Clade : Mammal
    Genome : Human
    Assembly : GRCh38/hg38
    Group : Genes and Gene Predictions
    Track : Gencode v36
    Table : knownGene
    region : Check genome
    Output format : GTF - gene transfer format

    There are two input files for generating mouse index file: mm39.fa.gz from
    https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/  and
    gencode.vM26.annotation.gtf.gz from Gencode (http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz)
    We are using Version 26 of Gencode which is the one that
    UCSC uses. See below for more information about the annotations in UCSC:
    http://genome.ucsc.edu/cgi-bin/hgTables with this set up:
    Clade : Mammal
    Genome : Mouse
    Assembly : GRCm39/mm39
    Group : Genes and Gene Predictions
    Track : ALL Gencode vm26
    Table : Basic
    region : Check genome
    Output format : GTF - gene transfer format

    Decompress the gzip files using the command: gunzip "filename"

    For other data such as ATAC_seq and MINT_ChIP, we use BWA for generating index files.
    """
    def __init__(self, threads, input_dir, genome_directory_human,
                 genome_directory_mouse):
        self.threads = threads
        self.input_dir = input_dir
        self.human_genome_dir = genome_directory_human
        self.mouse_genome_dir = genome_directory_mouse

    def star_indexing_human(self):
        """
        Generating human index file from human genome reference hg38.fa and human annotation file
        gencode.v36.annotation.gtf.
        :return: human genome directory
        """
        human_genome_fasta_path = os.path.join(self.input_dir, "hg38.fa")
        if not os.path.exists(human_genome_fasta_path):
            sys.exit("The human genome reference (hg38.fa) does not exist in input directory. Follow the description"
                     "and put the file in input directory.")
        #human_gtf_path = os.path.join(self.input_dir, "hg_ucsc.gtf")
        human_gtf_path = os.path.join(self.input_dir, "gencode.v36.annotation.gtf")  #Using annotation file from Gencode
        if not os.path.exists(human_gtf_path):
            sys.exit("The human annotation file does not exist in input directory. Follow the description"
                     "and put the file in input directory.")
        if not os.path.exists(self.human_genome_dir):
            os.makedirs(self.human_genome_dir)
        command = "STAR --runThreadN {} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 100" \
            .format(self.threads, self.human_genome_dir, human_genome_fasta_path, human_gtf_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            print("Error in running STAR!")
        return

    def star_indexing_mouse(self):
        """
        Generating mouse index file from mouse genome reference mm39.fa and mouse annotation file
        gencode.vM26.annotation.gtf.

        """
        mouse_genome_fasta_path = os.path.join(self.input_dir, "mm39.fa")
        if not os.path.exists(mouse_genome_fasta_path):
            sys.exit("The mouse genome reference (mm39.fa) does not exist in input directory. Follow the description"
                     "and put the file in input directory.")

        #mouse_gtf_path = os.path.join(self.input_dir, "mm_ucsc.gtf")
        mouse_gtf_path = os.path.join(self.input_dir, "gencode.vM26.annotation.gtf")

        if not os.path.exists(mouse_gtf_path):
            sys.exit("The mouse annotation file does not exist in input directory. Follow the description"
                     "and put the file in input directory.")
        if not os.path.exists(self.mouse_genome_dir):
            os.makedirs(self.mouse_genome_dir)
        command = "STAR --runThreadN {} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 100" \
            .format(self.threads, self.mouse_genome_dir, mouse_genome_fasta_path, mouse_gtf_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            print("Error in running STAR!")
        return

    def bwa_indexing_human(self):
        human_genome_fasta_path = os.path.join(self.input_dir, "Homo_sapiens_assembly38.fasta")
        if not os.path.exists(human_genome_fasta_path):
            sys.exit("The human genome reference (homo_sapiens_assembly38.fasta) does not exist in input directory. Follow the instruction"
                     "and put the file in input directory.")

        command = "bwa index {} -p human_genome".format(human_genome_fasta_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            print("Error in running BWA!")

    def bwa_indexing_mouse(self):
        mouse_genome_fasta_path = os.path.join(self.input_dir, "GRCm39.primary_assembly.genome.fa")
        if not os.path.exists(mouse_genome_fasta_path):
            sys.exit("The mouse genome reference (GRCm39.primary_assembly.genome.fa) does not exist in input directory. Follow the instruction"
                     "and put the file in input directory.")

        command = "bwa index {} -p mouse_genome".format(mouse_genome_fasta_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            print("Error in running BWA!")