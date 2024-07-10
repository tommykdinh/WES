import argparse
import sys

from genome_index_generating import GenomeIndexGenerator


def parse_args():
    '''
    Parses the PDX analysis arguments.
    '''
    parser = argparse.ArgumentParser(description="Generate human and mouse index files.")

    parser.add_argument('--threads', type=int, default=12,
                        help='number of threads')

    parser.add_argument('--input', nargs='?', default="input",
                        help='input path')

    parser.add_argument('--genome_dir_human', nargs='?', default="genome_dir_human",
                        help='path to the human genome directory')

    parser.add_argument('--genome_dir_mouse', nargs='?', default="genome_dir_mouse",
                        help='path to the mouse genome directory')

    parser.add_argument('--sequence', nargs='?', default="RNA",
                        help='The sequence type. It should be either RNA, ATAC, MINT-CHIP, WES or lpWGS')

    parser.add_argument('--model', nargs='?', default="PDX",
                        help='It should be either PDX, human or mouse.')
    return parser.parse_args()


def genome_index_generate(args):
    index_generator = GenomeIndexGenerator(args.threads, args.input, args.genome_dir_human, args.genome_dir_mouse)
    if args.sequence == "RNA":
        if args.model == "PDX":
            index_generator.star_indexing_human()
            index_generator.star_indexing_mouse()
        elif args.model == "human":
            index_generator.star_indexing_human()
        elif args.model == "mouse":
            index_generator.star_indexing_mouse()
        else:
            sys.exit("Error! The model should be PDX, human or mouse.")


    else:# for MINT-CHIP, ATAC, WGS, lpWGS, use BWA for indexing
        index_generator.bwa_indexing_human()
        index_generator.bwa_indexing_mouse()


if __name__ == '__main__':
    args = parse_args()
    genome_index_generate(args)