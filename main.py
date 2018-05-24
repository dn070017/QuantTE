import argparse
import pickle
import shutil
import sys

from colorama import Fore, Style  
from file_io.read import *
from process.process import *

#with open('merge.fa','wb') as merge_file:
#    shutil.copyfileobj(open('../NGS1050418/QuantTE/refMrna.fa', 'rb'), merge_file)
#    shutil.copyfileobj(open('TE.fasta', 'rb'), merge_file)

def main(args):
    sys.stdout.write(Fore.GREEN + '[MAIN PROGRAM]\n' + Style.RESET_ALL)
    start_time = time.time()
    sys.stdout.write(Fore.GREEN + '[PROCESS]' + Style.RESET_ALL + ' argument validation\n')
    input_file, args = process_and_validate_argument(args)
    if args.extract:
        #sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + ' extract repeat header\n')
        repeat_list = read_table(input_file)
        with open(input_file['output_dir'] + '/repeat_list.pickle', 'wb') as out_pickle:
            pickle.dump(repeat_list, out_pickle, protocol=pickle.HIGHEST_PROTOCOL)
    if args.kallisto:
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', action='store_true', help='extract TE sequence')
    parser.add_argument('--kallisto', action='store_true', help='perform kallisto quantification')
    
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', type=str, help='input file in JSON format', required=True)
    
    args = parser.parse_args()
    
    main(args)