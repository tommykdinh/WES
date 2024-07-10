import subprocess
import os
import sys
import shutil

import pandas as pd


class WESPreprocess:
    """
    In this class, we perform preprocessing steps on PDX WES data.
    These steps include quality control with FASTQC, trimming with trimmomatic, alignment(mapping)
    with BWA-mem.
    The rest of the pipeline:
    1. Removing duplicate reads
    2. Get the read group information from fast files
    3. Create target intervals
    4. Indel Realignment
    5. Base quality recalibration
    6. Converting bed files of interval files to interval files (This is done just once. After you get interval files, you can comment out this step in run_preprocessing.py)
    7. Get the metrics of the bam files
    8. Create a mpileup file from a bam file
    9. Calling variants (SNP and InDel variants) with Varscan
    10. I performed one more step to remove variants that did not have the acceptable format by GATK.
    The problem happened when VarScan returned "Y" as a ref.
    11. Update sequence dictionary of the variant files.
    12. Merging the SNP and InDel variants
    13. Variants are annotated by annovar
    14. Filtering repeats

    """

    def __init__(self, read_r1, read_r2, sample, output, input_directory, sequence, wes_type, model):
        self.read_r1 = read_r1
        self.read_r2 = read_r2
        self.output = output
        self.sample = sample
        self.input = input_directory
        self.sequence = sequence
        self.type = wes_type
        self.model = model

    def fastqc(self, fastq_path, output_path):
        """
        Run the fastqc program on a specified fastq file.
        Parameters
        ----------
        fastq_path : str
            Path to the fastq file to analyze.
        output_path : str
            The parent directory where output will be written.
        """
        # print(output_path)
        output_fastqc_path = os.path.join(output_path, self.sample, "output_fastqc")
        if not os.path.exists(output_fastqc_path):
            os.makedirs(output_fastqc_path)
        print(output_fastqc_path)
        command = "fastqc -o {} {}".format(output_fastqc_path, fastq_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running fastqc!")

        dir = os.listdir(output_fastqc_path)

        for item in dir:
            if item.endswith(".zip"):
                os.remove(os.path.join(output_fastqc_path, item))

        return

    def trimmomatic(self, num_threads):
        """
        Run the trimmomatic program on  fastq files 1 and 2.
        Parameters
        ----------
        num_threads : int
        Number of threads.

        """
        output_trimmomatic_path = os.path.join(self.output, self.sample, "output_trimmomatic")
        if not os.path.exists(output_trimmomatic_path):
            os.makedirs(output_trimmomatic_path)

        output_forward_paired_path = os.path.join(output_trimmomatic_path, "paired_{}_1.fastq.gz".format(self.sample))
        output_forward_unpaired_path = os.path.join(output_trimmomatic_path,
                                                    "unpaired_{}_1.fastq.gz".format(self.sample))
        output_reverse_paired_path = os.path.join(output_trimmomatic_path, "paired_{}_2.fastq.gz".format(self.sample))
        output_reverse_unpaired_path = os.path.join(output_trimmomatic_path,
                                                    "unpaired_{}_2.fastq.gz".format(self.sample))
        illuminaclip_path = os.path.join(self.input, "TruSeq3-PE-2.fa")

        command = "java -jar /risapps/rhel7/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads {} {} {} {} {} {} {} ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
            .format(num_threads, self.read_r1, self.read_r2, output_forward_paired_path, output_forward_unpaired_path,
                    output_reverse_paired_path, output_reverse_unpaired_path, illuminaclip_path)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running trimmomatic!")

        os.remove(output_forward_unpaired_path)
        os.remove(output_reverse_unpaired_path)

        # Performing quality control after trimming
        self.fastqc(output_forward_paired_path, output_trimmomatic_path)
        self.fastqc(output_reverse_paired_path, output_trimmomatic_path)

        return

    def create_read_group_info(self):
        trimmed_reads_path = os.path.join(self.output, self.sample, "output_trimmomatic")
        if not os.path.exists(trimmed_reads_path):
            print("Error! the path to trimmed files does not exist.")

        forward_paired_path = os.path.join(trimmed_reads_path, "paired_{}_1.fastq.gz".format(self.sample))
        id_file = os.path.join(trimmed_reads_path, "id_file.txt")
        command = "zcat {} |  head -n 1 | cut -f 1-4 -d\":\" | sed 's/@//' | sed 's/:/_/g' > {}".format(
            forward_paired_path, id_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error to find id!")

        p_file = os.path.join(trimmed_reads_path, "p_file.txt")
        command_2 = "zcat {} | head -n 1 | cut -f 10 -d\":\" > {}".format(forward_paired_path, p_file)

        try:
            subprocess.check_call(command_2, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error to find sm!")

    def bwa_mem_mapping_to_human_mouse(self, threads):
        if self.model == "PDX":
            self.bwa_mem_mapping(threads, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64")
            self.bwa_mem_mapping(threads, "mouse_genome")
        elif self.model == "human":
            self.bwa_mem_mapping(threads, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64")

    def bwa_mem_mapping(self, num_threads, genome):
        trimmed_reads_path = os.path.join(self.output, self.sample, "output_trimmomatic")
        if not os.path.exists(trimmed_reads_path):
            print("Error! the path to trimmed files does not exist.")

        forward_paired_path = os.path.join(trimmed_reads_path, "paired_{}_1.fastq.gz".format(self.sample))
        reverse_paired_path = os.path.join(trimmed_reads_path, "paired_{}_2.fastq.gz".format(self.sample))

        output_dir = os.path.join(self.output, self.sample)
        output_aligned_reads_path = os.path.join(output_dir, "output_alignment_{}".format(genome))

        if not os.path.exists(output_aligned_reads_path):
            os.makedirs(output_aligned_reads_path)

        command = "bwa mem -t {}  {} {} {}  | samtools sort -o {}_{}_aligned.bam".format(num_threads, genome,
                                                                                         forward_paired_path,
                                                                                         reverse_paired_path,
                                                                                         self.sample, genome)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running bwa-mem!")
        shutil.move("{}_{}_aligned.bam".format(self.sample, genome), output_aligned_reads_path)

    def xenofilter(self, num_threads):
        """
        Run XenofilteR to filter mouse reads from human reads
        """
        if self.model == "PDX":
            input_aligned_bam_human_path = os.path.join(self.output, self.sample,
                                                        "output_alignment_resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64",
                                                        "{}_resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64_aligned.bam".format(
                                                            self.sample))
            if not os.path.exists(input_aligned_bam_human_path):
                print("Error! The aligned human file does noe exist!")
            input_aligned_bam_mouse_path = os.path.join(self.output, self.sample, "output_alignment_mouse_genome",
                                                        "{}_mouse_genome_aligned.bam".format(self.sample))
            if not os.path.exists(input_aligned_bam_mouse_path):
                print("Error! The aligned mouse file does not exist!")

            output_xenofilter_path = os.path.join(self.output, self.sample, "output_xenofilter")
            if not os.path.exists(output_xenofilter_path):
                os.makedirs(output_xenofilter_path)
            xenofilter_rscript = os.path.join(self.output, self.sample, "output_xenofilter", "Xenofilter.R")
            with open(xenofilter_rscript, 'wt') as r1_fp:
                r1_fp.write('library(XenofilteR)\n')
                r1_fp.write(
                    'bp.param <- SnowParam(workers={}, type="SOCK")\n'.format(num_threads))
                r1_fp.write('sample.list <- data.frame(\n')
                r1_fp.write('graft.bam = \"{}\", '.format(input_aligned_bam_human_path))
                r1_fp.write('host.bam = \"{}\") \n'.format(input_aligned_bam_mouse_path))
                r1_fp.write('XenofilteR(sample.list, ')
                r1_fp.write('destination.folder = \"{}\", '.format(output_xenofilter_path))
                r1_fp.write('bp.param=bp.param)')
            try:
                subprocess.check_call("Rscript {}".format(xenofilter_rscript), shell=True)
            except subprocess.CalledProcessError as error:
                sys.exit("Error in running XenofilteR!")

        return

    def sort_bam(self):
        """
        Sort filtered file by coordinate.
        """
        if self.model == "PDX":
            input_file = os.path.join(self.output, self.sample, "output_xenofilter", "Filtered_bams",
                                      "{}_resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64_aligned_Filtered.bam".format(
                                          self.sample))
            if not os.path.exists(input_file):
                sys.exit("Error! The path to the filtered bam file does not exist!")

            output_path = os.path.join(self.output, self.sample, "sorted_filtered_bam")
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            output_file = os.path.join(output_path, "sorted_{}_filtered_human.bam".format(self.sample))
            command = "samtools sort {} -o {}".format(input_file, output_file)  # sort by coordinate
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as error:
                sys.exit("Error in running samtools to sort a bam file!")
            os.remove(input_file)
        return

    def remove_duplicates(self):
        """
        Duplicate reads are removed from the data with picard.
        """
        if self.model == "PDX":
            input_file = os.path.join(self.output, self.sample, "sorted_filtered_bam",
                                  "sorted_{}_filtered_human.bam".format(self.sample))
        else:
            input_file = os.path.join(self.output, self.sample,
                                      "output_alignment_resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64",
                                      "{}_resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64_aligned.bam".format(
                                          self.sample))
        output_path = os.path.join(self.output, self.sample, "remove_duplicates")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_removed_duplicate_file = os.path.join(output_path, "{}_removed_duplicates.bam".format(self.sample))
        output_metrics_file = os.path.join(output_path, "{}_marked_dup_metrics.txt".format(self.sample))
        command = "java -jar $PICARDHOME/picard.jar MarkDuplicates I={} O={} M={} REMOVE_DUPLICATES=true".format(
            input_file, output_removed_duplicate_file, output_metrics_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running picard for removing duplicate reads!")
        os.remove(input_file)

    def add_replace_groups(self):
        input_bam_file = os.path.join(self.output, self.sample, "remove_duplicates",
                                      "{}_removed_duplicates.bam".format(self.sample))
        output_bam_file = os.path.join(self.output, self.sample, "remove_duplicates",
                                       "{}_removed_duplicates_RG.bam".format(self.sample))

        id_file = os.path.join(self.output, self.sample, "output_trimmomatic", "id_file.txt")
        p_file = os.path.join(self.output, self.sample, "output_trimmomatic", "p_file.txt")
        with open(id_file, "r") as f_read_id:
            for readline in f_read_id:
                id = readline.strip()
            # id = str(f_read_id.readline())
            #print(id)

        with open(p_file, "r") as f_read_p:
            for readline in f_read_p:
                p = readline.strip()
            #print(p)
        pu = id + "_" + p
        #print(pu)

        command = "java -jar $PICARDHOME/picard.jar AddOrReplaceReadGroups I={} O={} RGID={} RGLB={} RGPL=Illumina RGPU={} RGSM={}".format(
            input_bam_file, output_bam_file, id, pu, pu, self.sample)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to add groups!")
        os.remove(input_bam_file)
        return

    def index_bam_file(self):
        input_bam_file = os.path.join(self.output, self.sample, "remove_duplicates",
                                      "{}_removed_duplicates_RG.bam".format(self.sample))
        command = "samtools index {}".format(input_bam_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running samtools to index the bam file!")
        return

    def target_create(self):
        """
        Create target intervals
        gold standard indels file was downloaded from:
        https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
        """
        input_bam_file = os.path.join(self.output, self.sample, "remove_duplicates",
                                      "{}_removed_duplicates_RG.bam".format(self.sample))
        output_path = os.path.join(self.output, self.sample, "target_interval")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "realigner.intervals")
        human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
        known_sites = os.path.join(self.input,
                                   "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
        command = "java -jar /risapps/rhel7/gatk/3.7/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {} -I {} -known {} -o {}".format(
            human_fasta_file, input_bam_file, known_sites, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to run realigner target creator!")
        return

    def indel_realigner(self):
        """
        gold standard indels file was downloaded from:
        https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
        """
        input_bam_file = os.path.join(self.output, self.sample, "remove_duplicates",
                                      "{}_removed_duplicates_RG.bam".format(self.sample))
        target_intervals_file = os.path.join(self.output, self.sample, "target_interval", "realigner.intervals")
        if not os.path.exists(target_intervals_file):
            sys.exit("Error! The target interval file does not exist!")

        human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
        known_sites = os.path.join(self.input,
                                   "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
        if not os.path.exists(known_sites):
            sys.exit("Error! The known site file does not exist!")
        output_path = os.path.join(self.output, self.sample, "indel_realign")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "{}_indel_realigned.bam".format(self.sample))
        command = "java -jar /risapps/rhel7/gatk/3.7/GenomeAnalysisTK.jar -T IndelRealigner -R {} -I {} " \
                  "-known {} -targetIntervals {} -o {}".format(human_fasta_file, input_bam_file, known_sites,
                                                               target_intervals_file,
                                                               output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to run indel realigner!")
        os.remove(input_bam_file)
        return

    def base_quality_score_recalibration(self):
        """
        Base quality recalibration
        gold standard Mills indels file was downloaded from:
        https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
        """
        input_bam_file = os.path.join(self.output, self.sample, "indel_realign",
                                      "{}_indel_realigned.bam".format(self.sample))
        human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
        if not os.path.exists(human_fasta_file):
            sys.exit("Error! Fasta file does not exist")
        output_path = os.path.join(self.output, self.sample, "base_recalibation")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "recal_data.table")
        known_sites_1 = os.path.join(self.input,
                                     "resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
        if not os.path.exists(known_sites_1):
            sys.exit("Error! The known site file {} does not exist!".format(known_sites_1))

        command = "singularity run /rsrch3/home/lym_myl_rsch/vravanmehr/WES_gatk3/gatk_4.1.9.0.sif gatk BaseRecalibrator -I {} -R {} --known-sites {} -O {}".format(
            input_bam_file, human_fasta_file, known_sites_1, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to perform base quality score recalibrator!")
        return

    def apply_bqsr(self):
        input_bam_file = os.path.join(self.output, self.sample, "indel_realign",
                                      "{}_indel_realigned.bam".format(self.sample))
        human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")

        recal_data = os.path.join(self.output, self.sample, "base_recalibation", "recal_data.table")
        if not os.path.exists(recal_data):
            sys.exit("Error! The known site file {} does not exist!".format(recal_data))

        output_path = os.path.join(self.output, self.sample, "after_base_recalibation")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "{}_applied_bqsr.bam".format(self.sample))
        command = "singularity run /rsrch3/home/lym_myl_rsch/vravanmehr/WES_gatk3/gatk_4.1.9.0.sif gatk ApplyBQSR  -I {} -R {} --bqsr-recal-file {} -O {}".format(
            input_bam_file, human_fasta_file, recal_data, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to apply bqsr!")
        os.remove(input_bam_file)
        return

    def bedtointerval(self):
        interval_list = os.path.join(self.input, "hglft_genome_9ea6_27a6c0_primary_target_hg38.bed")
        output_file = os.path.join(self.input, "hglft_genome_9ea6_27a6c0_primary_target_hg38.interval")
        dictionary = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")
        command = "java -jar $PICARDHOME/picard.jar BedToIntervalList I={} O={} SD={}".format(interval_list,
                                                                                              output_file, dictionary)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running piccard bedtointerval!")
        return

    def bedtointerval_2(self):
        interval_list = os.path.join(self.input, "hglft_genome_3d975_26b980_capture_target_hg38.bed")
        output_file = os.path.join(self.input, "hglft_genome_3d975_26b980_capture_target_hg38.interval")
        dictionary = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")
        command = "java -jar $PICARDHOME/picard.jar BedToIntervalList I={} O={} SD={}".format(interval_list,
                                                                                              output_file, dictionary)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running piccard bedtointerval!")
        return

    def bedtointerval_3(self):
        interval_list = os.path.join(self.input, "xgen-exome-hyb-panel-v2-probes-hg38.bed")
        output_file = os.path.join(self.input, "xgen-exome-hyb-panel-v2-probes-hg38.interval")
        dictionary = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")
        command = "java -jar $PICARDHOME/picard.jar BedToIntervalList I={} O={} SD={}".format(interval_list,
                                                                                              output_file, dictionary)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running picard bedtointerval!")
        return

    def bedtointerval_4(self):
        interval_list = os.path.join(self.input, "S31285117_hs_hg38", "S31285117_Covered.bed")
        output_file = os.path.join(self.input, "S31285117_hs_hg38", "S31285117_Covered.interval")
        dictionary = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")
        command = "java -jar $PICARDHOME/picard.jar BedToIntervalList I={} O={} SD={}".format(interval_list,
                                                                                              output_file, dictionary)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running picard bedtointerval!")
        return

    def bedtointerval_5(self):
        interval_list = os.path.join(self.input, "S31285117_hs_hg38", "S31285117_Regions.bed")
        output_file = os.path.join(self.input, "S31285117_hs_hg38", "S31285117_Regions.interval")
        dictionary = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.dict")
        command = "java -jar $PICARDHOME/picard.jar BedToIntervalList I={} O={} SD={}".format(interval_list,
                                                                                              output_file, dictionary)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running picard bedtointerval!")
        return

    def metrics(self):
        """
        info: https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-
        Bait_interval: An interval list file that contains the locations of the baits used.
        Target_interval: An interval list file that contains the locations of the targets.
        """
        if self.type == "Nimblegen":
            input_bam_file = os.path.join(self.output, self.sample, "after_base_recalibation", "{}_applied_bqsr.bam".format(self.sample))
            human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
            output_path = os.path.join(self.output, self.sample, "metrics")
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            bait_interval_list = os.path.join(self.input, "hglft_genome_9ea6_27a6c0_primary_target_hg38")
            target_interval_list = os.path.join(self.input, "hglft_genome_3d975_26b980_capture_target_hg38")
            output_file = os.path.join(output_path, "{}_output_hs_matrix.txt".format(self.sample))
            command = "java -jar $PICARDHOME/picard.jar  CollectHsMetrics I={} O={} R={} BAIT_INTERVALS={} TARGET_INTERVALS={}".format(
                input_bam_file, output_file, human_fasta_file, bait_interval_list, target_interval_list)
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as error:
                sys.exit("Error in running piccard!")

        elif self.type == "IDT":
            input_bam_file = os.path.join(self.output, self.sample, "after_base_recalibation", "{}_applied_bqsr.bam".format(self.sample))
            human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
            output_path = os.path.join(self.output, self.sample, "metrics")
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            bait_interval_list = os.path.join(self.input, "xgen-exome-hyb-panel-v2-probes-hg38")
            target_interval_list = os.path.join(self.input, "xgen-exome-hyb-panel-v2-targets-hg38")
            output_file = os.path.join(output_path, "{}_output_hs_matrix.txt".format(self.sample))
            command = "java -jar $PICARDHOME/picard.jar  CollectHsMetrics I={} O={} R={} BAIT_INTERVALS={} TARGET_INTERVALS={}".format(
                input_bam_file, output_file, human_fasta_file, bait_interval_list, target_interval_list)
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as error:
                sys.exit("Error in running picard!")
        elif self.type == "Agilent": # For DLBCL WES.
            input_bam_file = os.path.join(self.output, self.sample, "after_base_recalibation", "{}_applied_bqsr.bam".format(self.sample))
            human_fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
            output_path = os.path.join(self.output, "metrics")
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            bait_interval_list = os.path.join(self.input, "S31285117_hs_hg38", "S31285117_Covered.interval")
            target_interval_list = os.path.join(self.input, "S31285117_hs_hg38","S31285117_Regions.interval")
            output_file = os.path.join(output_path, "{}_hs_matrix.txt".format(self.sample))
            command = "java -jar $PICARDHOME/picard.jar  CollectHsMetrics I={} O={} R={} BAIT_INTERVALS={} TARGET_INTERVALS={}".format(
                input_bam_file, output_file, human_fasta_file, bait_interval_list, target_interval_list)
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as error:
                sys.exit("Error in running picard!")
        else:
            return "Error!"

        return

    def create_mpileup(self):
        input_bam_file = os.path.join(self.output, self.sample, "after_base_recalibation", "{}_applied_bqsr.bam".format(self.sample))

        output_path = os.path.join(self.output, self.sample, "mpileup")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        fasta_file = os.path.join(self.input, "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")
        output_file = os.path.join(output_path, "{}.mpileup".format(self.sample))
        command = "samtools mpileup -B -f {} -q 10 -o {} {}".format(fasta_file, output_file, input_bam_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running samtools to create mpileup of the BAM file!")
        return

    def snp_variant_calling_varscan(self):
        """
        call SNP variants with varscan
        """
        input_file = os.path.join(self.output, self.sample, "mpileup", "{}.mpileup".format(self.sample))

        output_path = os.path.join(self.output, self.sample, "variant_calling_varscan")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "{}_snp_variants.vcf".format(self.sample))
        command = "java -jar /risapps/noarch/varscan/2.4.2/VarScan.v2.4.2.jar mpileup2snp {} --p-value 0.05 --output-vcf 1 --variants > {}".format(
            input_file, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running varscan to call snp variants!")
        return

    def indel_variant_calling_varscan(self):
        """
        call InDel variants with varscan
        """
        input_file = os.path.join(self.output, self.sample, "mpileup", "{}.mpileup".format(self.sample))

        output_path = os.path.join(self.output, self.sample, "variant_calling_varscan")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "{}_indel_variants.vcf".format(self.sample))

        command = "java -jar /risapps/noarch/varscan/2.4.2/VarScan.v2.4.2.jar mpileup2indel {} --p-value 0.05 --output-vcf 1 --variants > {}".format(
            input_file, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running varscan to call indel variants!")
        return

    def snp_varscan_update(self):
        snp_varscan_input_file = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                              "{}_snp_variants.vcf".format(self.sample))
        corrected_snp_varscan = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                             "{}_corr_snp_variants.vcf".format(self.sample))
        command_1 = "awk '($4 == \"A\") || ($4 == \"C\") || ($4 ==\"T\") || ($4 == \"G\")' {} > {}".format(
            snp_varscan_input_file, corrected_snp_varscan)
        try:
            subprocess.check_call(command_1, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in correcting snp variants!")

        corrected_with_header_snp = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                                 "{}_corrected_snp_variants.vcf".format(self.sample))

        command_2 = "cat {} | head -n 24 > {}".format(snp_varscan_input_file, corrected_with_header_snp)
        try:
            subprocess.check_call(command_2, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in getting the header of snp variants!")

        command_3 = "cat {} >> {}".format(corrected_snp_varscan, corrected_with_header_snp)
        try:
            subprocess.check_call(command_3, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in concatenating header and corrected vcf file!")
        os.remove(corrected_snp_varscan)
        # os.remove(snp_varscan_input_file)

    def indel_varscan_update(self):
        indel_varscan_input_file = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                              "{}_indel_variants.vcf".format(self.sample))
        corrected_with_header_indel = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                                   "{}_corrected_indel_variants.vcf".format(self.sample))
        command_1 = "awk '($4 != \"Y\")' {} > {}".format(
            indel_varscan_input_file, corrected_with_header_indel)
        try:
            subprocess.check_call(command_1, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in correcting snp variants!")

        # os.remove(indel_varscan_input_file)

    def update_seq_dict_vcf(self):
        variant_input_file1 = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                           "{}_corrected_indel_variants.vcf".format(self.sample))
        variant_input_file2 = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                           "{}_corrected_snp_variants.vcf".format(self.sample))
        input_bam_file = os.path.join(self.output, self.sample, "after_base_recalibation", "{}_applied_bqsr.bam".format(self.sample))
        output_path = os.path.join(self.output, self.sample, "variant_calling_varscan")

        output_file_1 = os.path.join(output_path, "indel_variants_replacedcontiglines.vcf")
        output_file_2 = os.path.join(output_path, "snp_variants_replacedcontiglines.vcf")
        command_1 = "singularity run /rsrch3/home/lym_myl_rsch/vravanmehr/WES_gatk3/gatk_4.1.9.0.sif gatk " \
                    "UpdateVCFSequenceDictionary -V {} --source-dictionary {} --output {} --replace true".format(
            variant_input_file1,
            input_bam_file, output_file_1)
        try:
            subprocess.check_call(command_1, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to update the indel vcf file!")

        command_2 = "singularity run /rsrch3/home/lym_myl_rsch/vravanmehr/WES_gatk3/gatk_4.1.9.0.sif gatk " \
                    "UpdateVCFSequenceDictionary -V {} --source-dictionary {} --output {} --replace true".format(
            variant_input_file2, input_bam_file, output_file_2)
        try:
            subprocess.check_call(command_2, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running gatk to update the snp vcf file!")
        return

    def merge_variants(self):
        variant_input_file1 = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                           "indel_variants_replacedcontiglines.vcf")
        variant_input_file2 = os.path.join(self.output, self.sample, "variant_calling_varscan",
                                           "snp_variants_replacedcontiglines.vcf")

        output_path = os.path.join(self.output, self.sample, "combined_variants_varscan")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file = os.path.join(output_path, "{}_combined_varscan_variants.vcf".format(self.sample))
        command = "java -jar $PICARDHOME/picard.jar MergeVcfs I={} I={} O={}".format(variant_input_file1,
                                                                                     variant_input_file2, output_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running picard to merge vcf files!")
        return


    def annotate_variants_annovar(self):
        """
        Annotate variants with Annovar:
        """
        input_vcf_file = os.path.join(self.output, self.sample, "combined_variants_varscan",
                                      "{}_combined_varscan_variants.vcf".format(self.sample))
        avinput_file = os.path.join(self.output, self.sample, "combined_variants_varscan",
                                       "{}_combined_varscan_variants.avinput".format(self.sample))

        command_1 = "convert2annovar.pl -format vcf4 {} -includeinfo > {}".format(input_vcf_file, avinput_file)
        try:
            subprocess.check_call(command_1, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running annovar to convert to annovar input format!")

        output_annotated_path = os.path.join(self.output, "annotation_annovar")
        if not os.path.exists(output_annotated_path):
            os.makedirs(output_annotated_path)

        output_annotated_file = os.path.join(output_annotated_path, "{}_annotations".format(self.sample))

        command_2 = "table_annovar.pl {} /risapps/rhel7/annovar/2019.10.24/humandb_hg38/ -buildver hg38 \
        -out {} -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . \
        -csvout -polish".format(avinput_file, output_annotated_file)
        try:
            subprocess.check_call(command_2, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running annovar!")

        output_gene_annotated_path = os.path.join(self.output, "gene_annotation_annovar")
        if not os.path.exists(output_gene_annotated_path):
            os.makedirs(output_gene_annotated_path)

        output_annotated_gene_file = os.path.join(output_gene_annotated_path, "{}_gene_annotations".format(self.sample))

        command_3 = "annotate_variation.pl -geneanno --separate -dbtype refGene -buildver hg38 {} -out {} /risapps/rhel7/annovar/2019.10.24/humandb_hg38/".format(
            avinput_file, output_annotated_gene_file )
        try:
            subprocess.check_call(command_3, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running annovar gene annotation!")

    def compare_annotation_w_repeat_regions(self):
        input_annotated_file = os.path.join(self.output, "annotation_annovar", "{}_annotations.hg38_multianno.csv".format(self.sample))
        #remove header for mthe annotation file
        input_annotated_file_no_header = os.path.join(self.output, "annotation_annovar", "{}_annotations.hg38_multianno.csv.no_header.bed".format(self.sample))

        input_annotated_file_df = pd.read_csv(input_annotated_file)
        col_names = input_annotated_file_df.columns
        print(col_names)
        input_annotated_file_df.to_csv(input_annotated_file_no_header, header = None, index = False, sep ="\t")
        ##compare with repeats from tandem repeats and repeat masker data bases (hg38)
        tandem_repeat_database_hg38 = "/rsrch4/scratch/lym_myl_rsch/vravanmehr/duplex_pipeline/data/Allogene_715/duplex_cons_annotations/repeats/hg38.trf.bed" #change the path!
        repeat_masker_database_hg38 = "/rsrch4/scratch/lym_myl_rsch/vravanmehr/duplex_pipeline/data/Allogene_715/duplex_cons_annotations/repeats/repeatMasker.bed" #change the path!

        coverage_with_repeat_masker = os.path.join(self.output, "annotation_annovar", "{}_annotations_cov_masker.bed".format(self.sample))
        repeat_masker_filtered = os.path.join(self.output, "annotation_annovar", "{}_annotations_filtered_repeat_masker.bed".format(self.sample))
        command_1 = "bedtools coverage -a {} -b {} > {}".format(input_annotated_file_no_header, repeat_masker_database_hg38, coverage_with_repeat_masker)
        try:
            subprocess.check_call(command_1, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running bedtools coverage!")

        cov_with_masker = pd.read_csv(coverage_with_repeat_masker, sep = "\t", header =None)
        repeat_masker_filter = cov_with_masker[cov_with_masker.iloc[:,57] == 0]
        repeat_masker_filter.to_csv(repeat_masker_filtered, sep ="\t", index = False, header = None)
        coverage_with_tandem = os.path.join(self.output, "annotation_annovar", "{}_annotations_cov_with_tandem.bed".format(self.sample))
        command_2 = "bedtools coverage -a {} -b {} > {}".format(repeat_masker_filtered, tandem_repeat_database_hg38, coverage_with_tandem)
        try:
            subprocess.check_call(command_2, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running bedtools coverage!")
        repeat_filtered = os.path.join(self.output, "annotation_annovar", "{}_annotations_filtered_repeats".format(self.sample))
        cov_with_tandem = pd.read_csv(coverage_with_tandem, sep = "\t", header =None)
        tandem_repeat_filter = cov_with_tandem[cov_with_tandem.iloc[:,61] == 0]
        tandem_repeat_filter_drop_columns = tandem_repeat_filter.drop([54,55,56,57,58,59,60,61],axis=1)
        print(tandem_repeat_filter_drop_columns.shape)
        #final_filtered_repeats = pd.DataFrame(tandem_repeat_filter_drop_columns, columns= col_names)
        tandem_repeat_filter_drop_columns.columns = col_names
        tandem_repeat_filter_drop_columns.to_csv(repeat_filtered, sep ="\t", index = False)
        
        os.remove(input_annotated_file_no_header)
        os.remove(coverage_with_repeat_masker)
        os.remove(coverage_with_tandem)
        os.remove(repeat_masker_filtered)



    def annotate_variants_vep(self):
        """
        Annotate variants with VEP:
        """
        input_vcf_file = os.path.join(self.output, self.sample, "combined_variants_varscan",
                                      "{}_combined_varscan_variants.vcf".format(self.sample))

        output_annotated_path = os.path.join(self.output, "annotation_vep")
        if not os.path.exists(output_annotated_path):
            os.makedirs(output_annotated_path)

        output_annotated_file = os.path.join(output_annotated_path, "{}_annotations_vep.txt".format(self.sample))

        command = "vep -i {} -o {} --dir_cache /rsrch3/scratch/reflib/REFLIB_data/ensembl-vep-v101/  \
        --assembly GRCh38 --offline --fork 10 --species homo_sapiens --everything --force_overwrite".format(input_vcf_file, output_annotated_file)
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as error:
            sys.exit("Error in running vep!")
	


   