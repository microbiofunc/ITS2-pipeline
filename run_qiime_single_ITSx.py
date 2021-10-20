#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 14:41:24 2018

@author: omer
"""

import subprocess
import multiprocessing
import os
from Bio import SeqIO
import pandas as pd

def modify_rec(rec):
    rec.id = rec.id.split('|')[0]
    rec.description = ''
    return rec

def run_qiime(paths, primers):
        
    n_cores = multiprocessing.cpu_count()
    
    print('start')
    c = 1
    cmd = "qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path %s \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path %s/demux-seqs.qza" %(paths['trimmed'] , paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    
    cmd = "qiime demux summarize \
    --i-data %s/demux-seqs.qza \
    --o-visualization %s/demux-seqs.qzv" %(paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime dada2 denoise-single \
      --i-demultiplexed-seqs %s/demux-seqs.qza \
      --p-trim-left 0 \
      --p-trunc-len 0 \
      --p-n-threads 0 \
      --o-representative-sequences %s/rep-seqs-untrimmed.qza \
      --o-table %s/table.qza \
      --o-denoising-stats %s/stats.qza" %(paths['qiime_ready'], paths['qiime_ready'], paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return

    cmd = "qiime metadata tabulate \
    --m-input-file %s/stats.qza \
    --o-visualization %s/stats.qzv" %(paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime tools export \
    --input-path %s/rep-seqs-untrimmed.qza \
    --output-path %s/"  %(paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "ITSx -i %s/dna-sequences.fasta -o %s/its --cpu %i  --save_regions ITS2" %(paths['qiime_ready'], paths['qiime_ready'], n_cores)
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
    
    
    
    c = 1
    its_path = '%s/its.ITS2.fasta' %paths['qiime_ready']
    
    with open(its_path) as f:
        reads = [r for r in SeqIO.parse(f, 'fasta')]
        
    itsx_df = pd.DataFrame()
    itsx_df['feature'],itsx_df['taxa'],_ = zip(*list(map(lambda xx: xx.id.split('|'), reads)))
    itsx_df.set_index('feature', drop=True, inplace=True)
    
    rep_seqs_path = '%s/rep-seqs.qza' %paths['qiime_ready']
    rep_seqs_fasta = rep_seqs_path.replace('qza', 'fasta')
    
    reads_out = list(map(modify_rec, reads))
    
    with open(rep_seqs_fasta, 'wt') as f:
        SeqIO.write(reads_out, f, 'fasta')
        
    feature_df = pd.DataFrame(list(map(lambda x: x.id, reads_out)), columns = ['feature id'])
    feature_df.to_csv(paths['qiime_ready'] / 'features.csv', sep='\t', index=False)
    
    cmd = "qiime feature-table filter-features \
    --i-table %s/table.qza \
    --m-metadata-file %s/features.csv \
    --o-filtered-table %s/filtered_table.qza" %(paths['qiime_ready'], paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    c+=1

    if proc.stderr:
        print(c)
        print(proc.stderr)
        return

    cmd = "qiime tools import \
    --input-path %s \
    --output-path %s \
    --type 'FeatureData[Sequence]'" %(rep_seqs_fasta, rep_seqs_path)
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    c+=1

    if proc.stderr:
        print(c)
        print(proc.stderr)
        return

    
    
    cmd = "qiime feature-table summarize \
    --i-table %s/filtered_table.qza \
    --o-visualization %s/filtered_table.qzv \
    --m-sample-metadata-file %s"%(paths['qiime_ready'], paths['qiime_ready'], paths['metadata_file'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime feature-table tabulate-seqs \
    --i-data %s/rep-seqs.qza \
    --o-visualization %s/rep-seqs.qzv" %(paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    #%% Train classifier
    
    cmd = "qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path %s/sh_qiime_release_04.02.2020/sh_refs_qiime_ver8_dynamic_04.02.2020.fasta \
    --output-path %s/unite_dev_dynamic_otus.qza" %(paths['unite'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path %s/sh_qiime_release_04.02.2020/sh_taxonomy_qiime_ver8_dynamic_04.02.2020.txt \
    --output-path %s/ref-taxonomy_unite_dev_dynamic.qza" %(paths['unite'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime feature-classifier extract-reads \
    --i-sequences %s/unite_dev_dynamic_otus.qza \
    --p-f-primer %s \
    --p-r-primer %s \
    --o-reads %s/ref-seqs.qza" %(paths['qiime_ready'], primers['F'], primers['RC_R'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        return
    
    cmd = "qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads %s/unite_dev_dynamic_otus.qza \
    --i-reference-taxonomy %s/ref-taxonomy_unite_dev_dynamic.qza \
    --o-classifier %s/classifier.qza"  %(paths['qiime_ready'], paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    #%% back to qiime pipeline
    
    cmd = "qiime feature-classifier classify-sklearn \
    --i-classifier %s/classifier.qza \
    --i-reads %s/rep-seqs.qza \
    --o-classification %s/taxonomy.qza" %(paths['qiime_ready'], paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime metadata tabulate \
    --m-input-file %s/taxonomy.qza \
    --o-visualization %s/taxonomy.qzv" %(paths['qiime_ready'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return
    
    cmd = "qiime taxa barplot \
    --i-table %s/filtered_table.qza \
    --i-taxonomy %s/taxonomy.qza \
    --m-metadata-file %s \
    --o-visualization %s/taxa-bar-plots.qzv" %(paths['qiime_ready'], paths['qiime_ready'], paths['metadata_file'], paths['qiime_ready'])
    
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    c+=1
    print(c)
    if proc.stderr:
        print(c)
        print(proc.stderr)
        return

    cmd = 'qiime tools export --input-path %s/taxonomy.qza --output-path %s/' %(paths['qiime_ready'],paths['qiime_ready'])
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmd = 'qiime tools export --input-path %s/filtered_table.qza --output-path %s/' %(paths['qiime_ready'],paths['qiime_ready'])
    subprocess.run(cmd, shell=True)
    cmd = 'qiime tools export --input-path %s/rep-seqs.qza --output-path %s/' %(paths['qiime_ready'],paths['qiime_ready'])
    subprocess.run(cmd, shell=True)
    
    cmd = 'biom convert -i %s/feature-table.biom -o %s/table.from_biom.txt --to-tsv' %(paths['qiime_ready'],paths['qiime_ready'])
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    
    
    with open (str(paths['qiime_ready'] / 'dna-sequences.fasta')) as f:
        reads = {r.id : str(r.seq) for r in SeqIO.parse(f, 'fasta')}
        
    df = pd.read_csv(paths['qiime_ready'] / 'taxonomy.tsv', index_col=0, delimiter='\t')
    df = pd.concat([df, itsx_df], axis=1)
    df['seq'] = df.index.map(lambda x: reads[x])
    
    biom_df = pd.read_csv(paths['qiime_ready'] / 'table.from_biom.txt', skiprows=1, delimiter='\t', index_col=0)
    
    df = pd.concat([df, biom_df], axis=1)
    df.to_csv(paths['qiime_ready'] / 'summary.csv')
