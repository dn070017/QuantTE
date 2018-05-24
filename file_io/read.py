import os
import pandas as pd
import re
import sys

from colorama import Fore, Style 

def read_table(input_file, filter_class=['LINE', 'LTR', 'DNA', 'Satellite', 'SINE'], length_threshold=75):
    if 'repeat_masker_table' in input_file:
        table_list = [input_file['repeat_masker_table']]
    elif 'repeat_masker_dir' in input_file:
        dir = input_file['repeat_masker_dir']
        table_list = os.listdir(dir)
        table_list = [dir + '/' + table for table in table_list]
    
    regex = '|'.join(filter_class)
    repeat_list = list()
    for input_table in table_list:
        sys.stdout.write(Fore.CYAN + '[PROCESS]' + Style.RESET_ALL + \
                        ' extract repeat header from ' + input_table + '\n')
        table = pd.read_csv(input_table, sep='\s+', skiprows=[0, 1, 2], header=None)
        table.columns = ['SW score', 'perc div.', 'perc del.', 'perc ins.', 'query seq', 
                        'query begin', 'query end', 'query (left)', 'strand', 'repeat',
                        'repeat class', 'repeat begin', 'repeat end', 'repeat (left)', 'ID']
        for i, data in table.iterrows():
            if (data['query end'] - data['query begin'] + 1) < length_threshold:
                continue
            if not re.match(regex, data['repeat class']):
                continue
            repeat_name = data['repeat class'] + ':' + data['repeat'] + ':' + \
                          data['query seq'] + ':' + str(data['strand']) + ':' + \
                          str(data['query begin']) + '-' + str(data['query end']) 
            repeat_list.append(repeat_name)

    return repeat_list

def read_fasta(fasta):
    sequence = ''
    genome_seq = dict()
    
    with open(fasta, 'r') as genome_file:
        for i, genome_line in enumerate(genome_file):
            if genome_line[0] == '>':
                if sequence != '':
                    genome_seq[genome_name] = sequence
                    sequence = ''
                regex = re.match('>(\S+)', genome_line)
                genome_name = regex.group(1)
            else:
                sequence += genome_line.rstrip()
        genome_seq[genome_name] = sequence

    return genome_seq