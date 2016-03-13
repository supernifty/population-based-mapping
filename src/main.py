#!/usr/bin/env python

import argparse
import collections
import datetime
import os
import sys

def write_log(log_fh, msg):
    '''
        debugging
    '''
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_fh.write('{0}: {1}\n'.format(now, msg))

def run(cmd, log_fh):
    write_log(log_fh, 'executing {0}'.format(cmd))
    os.system(cmd)

def write_read(target_fh, sequence, sequence_id):
    '''
        write fastq record
    '''
    sequence_quality = '~' * len(sequence)
    target_fh.write('@{0}\n{1}\n{2}\n{3}\n'.format(sequence_id, sequence, '+', sequence_quality))

def generate_reads(source_fh, kmer, resolution, target_fh, chromosomes, log_fh):
    '''
        fasta -> fastq
    '''
    write_log(log_fh, 'generating reads with length {0} and resolution {1}...'.format(kmer, resolution))
    current = ''
    current_pos = 0
    chromosome_id = 0
    read_start = 0
    for idx, line in enumerate(source_fh):
        if line.startswith('>'):
            chromosome_id += 1
            chromosomes.append(line[1:].split(' ')[0].strip())
        else:
            current += line.strip()
        # any sequence to write?
        while current_pos + len(current) >= read_start + kmer:
            write_read(target_fh, current[read_start - current_pos: read_start - current_pos + kmer], '{0}_{1}'.format(chromosome_id, read_start))
            read_start += resolution
            current = current[read_start - current_pos:]
            current_pos = read_start
    # write any remaining sequence
    while current_pos + len(current) > read_start:
        write_read(target_fh, current[read_start - current_pos:], '{0}_{1}'.format(chromosome_id, read_start))
        read_start += resolution
        current = current[read_start - current_pos:]
        current_pos = read_start
    write_log(log_fh, 'generating reads: done')

def index_reference(fasta_file, log_fh):
    run('bwa index %s' % fasta_file, log_fh)

def align_reference(fastq, reference, sam, log_fh):
    run('bwa mem -t 8 %s %s > %s' % (reference, fastq, sam), log_fh)

def analyze_sam(sam_fh, unique, multi, other, log_fh):
    for line in sam_fh:
        if line.startswith('@'):
            continue
        fields = line.strip('\n').split('\t')
        pos = fields[0].split('_')[1] # ignore chromosome
        flag = int(fields[1])
        mapq = int(fields[4])
        if flag == 0 and mapq > 0: # mapped uniquely
            unique[pos] += 1
            multi[pos] += 0 # ensure all values in each dict
            other[pos] += 0
        elif flag & 0x04 != 0 or mapq == 0: # multi map
            unique[pos] += 0
            multi[pos] += 1
            other[pos] += 0
        else:
            unique[pos] += 0
            multi[pos] += 0
            other[pos] += 1
    
def write_bedgraph_header(output):
    #output.write('track type=bedGraph name=track_label description=center_label visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16\n')
    output.write('track type=bedGraph\n')

def analyze(source, references, kmer, resolution, output, log_fh):
    '''
        align source to targets
        evaluate and generate data
    '''
    # generate reads
    fastq = '{0}.fastq'.format(source)
    chromosomes = []
    generate_reads(open(source, 'r'), kmer, resolution, open(fastq, 'w'), chromosomes, log_fh)
    
    # index align analyze
    sams = []
    unique = collections.defaultdict(int)
    multi = collections.defaultdict(int)
    other = collections.defaultdict(int)
    for idx, reference in enumerate(references):
        index_reference(reference, log_fh) # index
        sam = '{0}_{1}.sam'.format(source, idx)
        align_reference(fastq, reference, sam, log_fh) # align
        sams.append(sam)
        write_log(log_fh, 'analyzing sam file...')
        analyze_sam(open(sam, 'r'), unique, multi, other, log_fh) # analyze
        write_log(log_fh, 'analyzing sam file: done')

    # write bedgraph result
    write_log(log_fh, 'writing bedgraph...')
    write_bedgraph_header(output)
    for key, value in unique.items():
        pos = int(key)
        score = 1. * value / len(references)
        output.write('{0} {1} {2} {3}\n'.format(chromosomes[0], pos, pos + resolution, score))
    write_log(log_fh, 'writing bedgraph: done')

def main():
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='Calculate population wide reference bias')
    parser.add_argument('--source', required=True, help='source reference file')
    parser.add_argument('--targets', required=True, help='target reference files')
    parser.add_argument('--kmer', type=int, default=100, help='size of kmers')
    parser.add_argument('--resolution', type=int, default=10, help='resolution of analysis (1 means analyze every base)')
    args = parser.parse_args()
    analyze(args.source, args.targets.split(','), args.kmer, args.resolution, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()

