import argparse
#import ATAC_preprocessing
#import MINT_preprocessing
#import RNA_preprocessing
import WES_preprocessing


def parse_args():
    '''
    Parses the PDX analysis arguments.
    '''
    parser = argparse.ArgumentParser(description="Run PDX preprocessing.")

    parser.add_argument('--sample', nargs='?', help="Sample name")

    parser.add_argument('--output', nargs='?', default="output",
                        help='path to the output directory')

    parser.add_argument('--input', nargs='?', default="input",
                        help='path to the input files such as genome references, black list regions, capture targets'
                             ' and bait interval files, etc')

    parser.add_argument('--threads', type=int,  default=12,
                        help='number of threads')

    parser.add_argument('--sequence', nargs='?', default="RNA",
                        help='The type of sequencing data. It should be either "RNA", "MINT-CHIP", "ATAC", "WES" or "lpWGS" ')

    parser.add_argument('--input1', nargs='?',
                        help='The path to the first fastq file (RNA and WES).')

    parser.add_argument('--input2', nargs='?',
                        help='The path to the second fastq file (RNA or WES).')

    parser.add_argument('--atac_1_r1', nargs='?',
                        help='The path to the first ATAC-seq fastq file from the first library.')

    parser.add_argument('--atac_1_r2', nargs='?',
                        help='The path to the second ATAC-seq fastq file from the first library.')
    parser.add_argument('--atac_2_r1', nargs='?',
                        help='The path to the first ATAC-seq fastq file from the second library.')

    parser.add_argument('--atac_2_r2', nargs='?',
                        help='The path to the second ATAC-seq fastq file from the second library.')

    parser.add_argument('--mint_r1', nargs='?',
                        help='The path to the first mint-chip fastq file.')

    parser.add_argument('--mint_r2', nargs='?',
                        help='The path to the second mint-chip fastq file.')

    parser.add_argument('--mint_type', nargs='?', default="H3K27Ac",
                        help=' It should be "H3K27Ac" or "input" , "H3K4Me3", ...')

    parser.add_argument('--wes_type',  nargs='?', default="Nimblegen",
                        help='The type of WES sequencing data. It should be either "Nimblegen" or "IDT" ')

    parser.add_argument('--wes_model', nargs='?', default="PDX",
                        help='The model of WES data. It should be either "PDX" or "PBMC" ')
    return parser.parse_args()


def main(args):
    '''
    Pipeline for PDX preprocessing 
    '''
    if args.sequence == "RNA":

        preprocess = RNA_preprocessing.RNAPreprocess(args.input1, args.input2, args.sample, args.output)
        preprocess.fastqc(args.input1, args.output)
        preprocess.fastqc(args.input2, args.output)
        preprocess.trimmomatic(args.threads, args.input)
        if args.model == "PDX":
            preprocess.star_mapping_to_human_mouse(args.threads, args.genome_dir_human, args.genome_dir_mouse,
                                                   args.star_mode)
            preprocess.xenofilter(args.threads)
            preprocess.sort_bam()

        elif args.model == "human":
            preprocess.star_mapping_to_human(args.threads, args.genome_dir_human, args.star_mode)
        elif args.model == "mouse":
            preprocess.star_mapping_to_mouse(args.threads, args.genome_dir_mouse, args.star_mode)
        preprocess.htseq(args.input, args.model)

    elif args.sequence == "ATAC":
        preprocess = ATAC_preprocessing.ATACPreprocess(args.atac_1_r1, args.atac_1_r2, args.atac_2_r1, args.atac_2_r2, args.sample, args.output, args.input)
        preprocess.fastqc(args.atac_1_r1, args.output)
        preprocess.fastqc(args.atac_1_r2, args.output)
        preprocess.fastqc(args.atac_2_r1, args.output)
        preprocess.fastqc(args.atac_2_r2, args.output)
        preprocess.trimmomatic(args.threads,  "ATAC_1")
        preprocess.trimmomatic(args.threads, "ATAC_2")
        preprocess.bwa_mem_mapping_to_human_mouse(args.threads)
        preprocess.xenofilter(args.threads, "ATAC_1")
        preprocess.xenofilter(args.threads,  "ATAC_2")
        preprocess.mark_duplicate("ATAC_1")
        preprocess.mark_duplicate("ATAC_2")
        preprocess.merge_human_filtered_alignments()
        preprocess.sort_bam()
        preprocess.mark_duplicate_merged_file()
        preprocess.remove_chr_x_y_m_un()
        preprocess.remove_blacklist()
        preprocess.remove_duplicates()
        preprocess.remove_non_primary_alignment_reads()
        preprocess.remove_unmapped_reads()
        preprocess.remove_non_unique_reads()
        preprocess.remove_more_than_4_mismatches_reads()
        preprocess.remove_soft_clipp_reads()
        preprocess.remove_insert_size_greater_2kb_reads()
        preprocess.remove_reads_mapped_to_diff_chromosomes()
        preprocess.remove_not_fr_orientation()
        preprocess.alignment_QC()
        preprocess.estimate_library_complexity()
        preprocess.chrom_size()
        preprocess.normalize()
        preprocess.to_bigwig()
        preprocess.peak_calling()
        preprocess.annotate_peaks()
        preprocess.atac_qc_visualization()
        preprocess.deeptools_compute_matrix()
        preprocess.deeptools_plot_heatmap()
        preprocess.deeptools_plot_profile()

    elif args.sequence == "MINT-CHIP":
        preprocess = MINT_preprocessing.MINTPreprocess(args.mint_r1, args.mint_r2, args.sample, args.output, args.input,
                                                   args.sequence, args.mint_type)
        preprocess.fastqc(args.mint_r1, args.output)
        preprocess.fastqc(args.mint_r2, args.output)
        preprocess.trimmomatic(args.threads)
        preprocess.bwa_mem_mapping_to_human_mouse(args.threads)
        preprocess.xenofilter(args.threads)
        preprocess.sort_bam()
        preprocess.mark_duplicate()
        preprocess.remove_chr_x_y_un()
        preprocess.remove_blacklist()
        preprocess.remove_duplicates()
        preprocess.remove_non_primary_alignment_reads()
        preprocess.remove_unmapped_reads()
        preprocess.remove_non_unique_reads()
        preprocess.remove_more_than_4_mismatches_reads()
        preprocess.remove_insert_size_greater_2kb_reads()
        preprocess.remove_reads_mapped_to_diff_chromosomes()
        preprocess.remove_not_fr_orientation()
        preprocess.fix_reads()
        preprocess.alignment_QC()
        preprocess.estimate_library_complexity()
        # # preprocess.chrom_size()
        preprocess.normalize()
        preprocess.to_bigwig()
        preprocess.peak_calling_sicer(args.threads)
        preprocess.phantom()

    elif args.sequence == "WES":
        preprocess = WES_preprocessing.WESPreprocess(args.input1, args.input2, args.sample, args.output, args.input,
                                                  args.sequence, args.wes_type, args.wes_model)
        preprocess.fastqc(args.input1, args.output)
        preprocess.fastqc(args.input2, args.output)
        preprocess.trimmomatic(args.threads)
        preprocess.create_read_group_info()
        preprocess.bwa_mem_mapping_to_human_mouse(args.threads) #maybe this bwa works?

        #don't run these 2 lines (PDX only)
        preprocess.xenofilter(args.threads)
        preprocess.sort_bam()

        #these 5 lines are part of the run
        preprocess.remove_duplicates()
        preprocess.add_replace_groups()
        preprocess.index_bam_file()
        preprocess.target_create()
        preprocess.indel_realigner()

        # #don't run these four lines
        # preprocess.bedtointerval()
        # preprocess.bedtointerval_3()
        # preprocess.bedtointerval_4()
        # preprocess.bedtointerval_5()

        #everything below here is part of the run
        preprocess.base_quality_score_recalibration()
        preprocess.apply_bqsr()
        preprocess.metrics()
        preprocess.create_mpileup()
        preprocess.snp_variant_calling_varscan()
        preprocess.indel_variant_calling_varscan()
        preprocess.snp_varscan_update()
        preprocess.indel_varscan_update()
        preprocess.update_seq_dict_vcf()
        preprocess.merge_variants()
        preprocess.annotate_variants_annovar()
        preprocess.compare_annotation_w_repeat_regions()


if __name__ == "__main__":
    args = parse_args()
    main(args)