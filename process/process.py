import json
import os
import re
import sys
import time

from Bio.Seq import Seq
from colorama import Fore, Style  
from subprocess import Popen

def process_and_validate_argument(args):
    terminate = False
    if os.path.exists(args.input):
        json_file = os.path.abspath(args.input)
        with open(json_file) as input_json:
            try:
                input_file = json.load(input_json)
            except:
                sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' failed to load JSON format\n')
                terminate = True
            else:
                required_file = list()
            
                if 'output_dir' not in input_file:
                    input_file['output_dir'] = os.path.abspath('./')
                    sys.stderr.write(Fore.RED + '[WARNING]' + Style.RESET_ALL + ' set output_dir to ./\n')
                    
                if args.kallisto:
                    required_file.extend(['TE_fasta', 'transcript_fasta', 'read_fastq'])
                if args.extract:
                    if 'genome_pickle' in input_file:
                        required_file.append('genome_pickle')
                    else:
                        required_file.append('genome_fasta')
                    if 'repeat_masker_pickle' in input_file:
                        required_file.append('repeat_masker_pickle')
                    elif 'repeat_masker_dir' in input_file:
                        required_file.append('repeat_masker_dir')
                    else:
                        required_file.append('repeat_masker_table')
                    if 'TE_fasta' not in input_file:
                        input_file['TE_fasta'] = os.path.abspath(input_file['output_dir'] + '/TE.fasta')
                        sys.stderr.write(Fore.RED + '[WARNING]' + Style.RESET_ALL + ' set TE_fasta to ' + \
                                         input_file['output_dir'] + '/TE.fasta\n')
                
                for required in required_file:
                    if required not in input_file:
                        sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' required')
                        sys.stderr.write(Fore.RED + ' ' + required + Style.RESET_ALL + ' in input file\n')
                        terminate = True

                for var in input_file.keys():
                    file = input_file[var]
                    if isinstance(file, list):
                        if var == 'read_label':
                            for label in file:
                                sys.stderr.write(Fore.RED + '[WARNING] ' + Style.RESET_ALL + 'create directory for read_label\n')
                                os.makedirs(input_file['output_dir'] + '/' + label, exist_ok=True)
                        continue
                    elif var == 'output_dir':
                        path = os.path.abspath(file)
                        if not os.path.exists(path):
                            sys.stderr.write(Fore.RED + '[WARNING] ' + Style.RESET_ALL + 'create directory for outdir\n')
                            os.makedirs(path, exist_ok=True)
                            input_file[var] = os.path.abspath(file)
                    elif var not in ['TE_fasta']:
                        if os.path.exists(file):
                            input_file[var] = os.path.abspath(file)
                        else:
                            sys.stderr.write(Fore.RED + '[ERROR] ' + file + Style.RESET_ALL + ' does not exist\n') 
                            terminate = True
    else:
        sys.stderr.write(Fore.RED + '[ERROR]' + Style.RESET_ALL + ' input JSON does not exsist\n')
        terminate = True

    if terminate:
        sys.exit(1)
    else:
        return input_file, args

def extract_repeat_seq(repeat_list, genome_seq, output_file):
    with open(output_file, 'w') as output:
        for repeat_name in repeat_list:
            data = repeat_name.split(':')
            repeat_class = data[0]
            repeat = data[1]
            chromosome = data[2]
            strand = data[3]
            start = int(data[4].split('-')[0]) - 1
            end = int(data[4].split('-')[1]) - 1

            print('>' + repeat_name, file=output)
            if chromosome not in genome_seq:
                continue
            if strand == '+':
                seq = genome_seq[chromosome][start:end+1] 
                print('\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)), file=output)
            elif strand == 'C' or strand == '-':
                seq = str(Seq(genome_seq[chromosome][start:end+1]).reverse_complement())
                print('\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)), file=output)

    return

def merge_xprs_result(xprs_list):
    xprs_dict = defaultdict(dict)
    with open(sys.argv[1], 'r') as xprs:
        for i, xprs_line in enumerate(xprs):
            if i == 0:
                continue
            xprs_data = xprs_line.split('\t')
            target_name = xprs_data[0]
            length = int(xprs_data[1])
            count = float(xprs_data[3])
            tpm = float(xprs_data[4])
            regex = re.match('(\S+?):(\S+?):', target_name)
            if regex:
                target_name = regex.group(1) + ':' + regex.group(2)
            try:
                xprs_dict[target_name]['seq'] += 1
                xprs_dict[target_name]['length'] += length
                xprs_dict[target_name]['count'] += count
                xprs_dict[target_name]['tpm'] += tpm
            except:
                xprs_dict[target_name]['seq'] = 1
                xprs_dict[target_name]['length'] = length
                xprs_dict[target_name]['count'] = count
                xprs_dict[target_name]['tpm'] = tpm

    for name in sorted(xprs_dict.keys()):
        print(name, round(xprs_dict[name]['length']/xprs_dict[name]['seq'], 3),
              round(xprs_dict[name]['count'], 3), round(xprs_dict[name]['tpm'], 3), sep='\t')
    return

def asynchronous_quant():
    command_list = ['kallisto quant -t 32 -o ./kallisto/A -i ./kallisto/kallisto_index --single -l 80 -s 15 ../../NGS1050418/A.fastq.gz',
                    'kallisto quant -t 32 -o ./kallisto/B -i ./kallisto/kallisto_index --single -l 80 -s 15 ../../NGS1050418/B.fastq.gz']
    queue = [Popen(command.split()) for command in command_list]
    
    while True:
        if len(queue) == 0:
            break
        for process in queue:
            retcode = process.poll()
            if retcode is not None:
                queue.remove(process)
            else: # No process is done, wait a bit and check again.
                time.sleep(.1)
                continue

#os.system('kallisto index -i kallisto/kallisto_index ../merge.fa')
#asynchronous_quant()