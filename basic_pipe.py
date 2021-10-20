#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:21:56 2018

@author: omer
"""
from pathlib import Path
import re
import subprocess
from run_qiime_single_ITSx import run_qiime
import pandas as pd
#import pdb
import numpy as np
import psutil
from sys import argv



def merge(reads_1, reads_2, paths):
    out_file = paths['merged'] / reads_1.name
    log_file = (paths['logs'] / reads_1.name.replace('fastq.gz', 'txt'))
    #subprocess.run('module load pear', shell=True)
    
    #Write you PEAR path here
    pear_path = 'pear '
    
    #for optimizing memory usage
    #free_mem = np.floor(psutil.virtual_memory().available / 10**9)-1
    #print(free_mem)
    #cmd = '%s -y %dG -f %s -r %s -o %s' %(pear_path, 8, reads_1, reads_2, out_file)
    
    cmd = '%s -f %s -r %s -o %s' %(pear_path, reads_1, reads_2, out_file)
    #print(cmd)
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(reads_1)
    log_file.write_text(proc.stdout.decode())
    stats = re.findall('\nAssembled reads \.*: (.*)', proc.stdout.decode())[0]
    stats = [float(x.replace(',','')) for x in re.findall('\d+[\.|,]?\d+' ,stats)]
    
    return stats

def trim(in_file, file_type, paths, primers):
    if not in_file.exists():
        print('fnf')
        return
        
    if file_type == 'merged':
        out_file = paths['trimmed'] / re.sub('.assembled.fastq','',in_file.name)
        cmd = 'cutadapt --minimum-length 80 --discard-untrimmed --cores 0 -g %s...%s -o %s %s' \
        %(primers['F'], reverse_complement(primers['R']), out_file, in_file)
    else:
        out_file = paths['trimmed'] / in_file.name
        if file_type == '1':
            cmd = 'cutadapt --minimum-length 50 --cores 0 -g %s -a %s -o %s %s' \
            %(primers['F'], reverse_complement(primers['R']), out_file, in_file)
            
        elif file_type == '2':
            cmd = 'cutadapt --minimum-length 50 --cores 0 -g %s -a %s -o %s %s' \
            %(primers['R'], reverse_complement(primers['F']), out_file, in_file)
    
    
    log_file = (paths['logs'] / in_file.name.split('.')[0]).with_suffix('.txt')
    
    #print(cmd)
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(proc.stderr.decode())
    log_file.write_text(proc.stdout.decode())
    
    stats = re.findall('Total reads processed:\s*(.*)',proc.stdout.decode())
    s = re.findall('Reads with adapters:\s*(.*)',proc.stdout.decode())
    stats = stats + re.findall('\d+[,|\.]?\d+', s[0])
    s = re.findall('Reads written \(passing filters\):\s*(.*)',proc.stdout.decode())
    stats = stats + re.findall('\d+[,|\.]?\d+', s[0])
    
    return stats

def trim_ITS(in_file, paths): # APPLIES ONLY TO FILES THAT WERE ALREADY MERGED
    out_file = paths['qiime_ready'] / re.sub('.assembled.fastq','.fastq.gz',in_file.name)
    cmd = 'itsxpress --fastq  %s --single_end \
    --region ITS2 --taxa Fungi --outfile %s --threads 4' %(in_file, out_file)
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    
def reverse_complement(seq):
    bases = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
    return "".join(reversed(list(map(lambda x: bases[x], seq))))
    
def define_params(library):
    paths = dict()
    #Define where output should be generated:
    paths['home'] = Path('<your-main-directory>')
    paths['samples'] = paths['home'] / 'samples'
    paths['fastq'] = paths['samples'] / library / 'fastq'
    paths['merged'] = paths['samples'] / library / 'merged'
    paths['trimmed'] = paths['samples'] / library / 'trimmed'
    paths['qiime_ready'] =  paths['samples'] / library / 'qiime_ready'
    paths['logs'] = paths['samples'] / library / 'logs'
    paths['metadata'] = paths['samples'] / library / 'metadata'
    paths['metadata_file'] = paths['samples'] / library / 'metadata' / (library + '.tsv')
    paths['unite'] = paths['home'] / 'unite'
    paths['itsxpress'] = paths['samples'] / library / 'itsxpress'
    
    paths['merged'].mkdir(exist_ok=True)
    paths['trimmed'].mkdir(exist_ok=True)
    paths['qiime_ready'].mkdir(exist_ok=True)
    (paths['qiime_ready'] / 'fastq').mkdir(exist_ok=True)
    paths['logs'].mkdir(exist_ok=True)
    paths['itsxpress'].mkdir(exist_ok=True)
    
    
    primers = {'F' : 'GTGAATCATCGAATCTTTGAA', 'R' : 'TCCTCCGCTTATTGATATGC'}
    primers['RC_R'] = reverse_complement(primers['R'])
    
    params = dict()
    params['merge'] = True
    params['trim'] = True
    params['itsxpress_standalone'] = False
    
    return paths, primers, params


library = argv[1]
paths, primers, params = define_params(library)

    
if params['merge']:
    print('merging')
    merged_df = pd.DataFrame(columns=('File', 'n_input', 'n_merged', 'percentage'))
    fastq_R_files = paths['fastq'].glob('*_R[2-9]_*.fastq.gz')
    file_n = 1
    
    for r_file in fastq_R_files:
        s = "file %d: %s" %(file_n, r_file)
        file_n +=1
        f_file = r_file.parent / re.sub('_R\d_', '_R1_', r_file.name)
        if not f_file.exists():
            print('not here :( ' + f_file)
            continue
        stats = merge(f_file, r_file, paths)
            
        row = [f_file.name.replace('.fasta.gz', '')] + stats
        merged_df.loc[len(merged_df)] = row
    
    merged_df.to_csv(paths['logs'] / 'merging_table.csv')
        

if params['trim']:
    print('trimming')
    trimmed_df = pd.DataFrame(columns=('File', 'n_input', 'n_trimmed', '%_trimmed', 'n_written', '%_written'))
    for f in paths['merged'].glob('*.assembled.fastq'):
        stats = trim(f, 'merged', paths, primers)
        row = [f.name] + [float(x.replace(',', '')) for x in stats]
        trimmed_df.loc[len(trimmed_df)] = row
    trimmed_df.to_csv(paths['logs'] / 'trimmed_table.csv')
        
if params['itsxpress_standalone']:
    file_n = 0
    for f in paths['merged'].glob('*.assembled.fastq'):
        file_n +=1
        print('file # ' + str(file_n))
        trim_ITS(f, paths)
    
 

run_qiime(paths, primers)
    
