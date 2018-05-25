import argparse
import pickle
import sys

from colorama import Fore, Style  
from file_io.read import *
from process.process import *

def main(args):
    sys.stdout.write(Fore.GREEN + '[MAIN PROGRAM]\n' + Style.RESET_ALL)
    start_time = time.time()
    sys.stdout.write(Fore.GREEN + '[PROCESS]' + Style.RESET_ALL + ' argument validation\n')
    input_file, args = process_and_validate_argument(args)
    if args.extract:
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' extract repeat header\n')
        if 'repeat_list_pickle' in input_file:
            with open(input_file['output_dir'] + '/repeat_list.pickle', 'rb') as in_pickle:
                repeat_list = pickle.load(in_pickle)
        else:
            repeat_list = read_table(input_file)
            with open(input_file['output_dir'] + '/repeat_list.pickle', 'wb') as out_pickle:
                pickle.dump(repeat_list, out_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' extract genome sequence\n')
        if 'genome_seq_pickle' in input_file:
            with open(input_file['output_dir'] + '/genome_seq.pickle', 'rb') as in_pickle:
                genome_seq = pickle.load(in_pickle)
        else:
            genome_seq = read_fasta(input_file['genome_fasta'])
            with open(input_file['output_dir'] + '/genome_seq.pickle', 'wb') as out_pickle:
                pickle.dump(genome_seq, out_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' extract repeat sequence\n')
        extract_repeat_seq(repeat_list, genome_seq, input_file['TE_fasta'])
    if args.kallisto:
        asynchronous_quant(input_file)
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' merge xprs result\n')
        merge_xprs_result(input_file)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', action='store_true', help='extract TE sequence')
    parser.add_argument('--kallisto', action='store_true', help='perform kallisto quantification')
    
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', type=str, help='input file in JSON format', required=True)
    
    args = parser.parse_args()
    
    main(args)