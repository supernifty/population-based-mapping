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

def update_range(start, count, deltas, dicts):
    '''
        read not mapped, update all positions with the unmapped status
    '''
    for x in xrange(start, start+count):
        for delta, target in zip(deltas, dicts):
            target[x] += delta

def handle_mapped(start, count, mapq, seen_once, seen_multiple, unique, multi, unmapped):
    '''
        read was mapped. iterate over each position in the read and update depending on if it's been seen before:
        - if the position has been seen at least once before, 
          - if it's been seen exactly once before: take it off the unique list, add to the multi list
          - else it's been seen more than once already, do nothing
        - otherwise, this is the first time for this position...
          - if mapq is acceptable, add to the unique list
          - else consider it a multimap
    '''
    if start in seen_once:
        if start not in seen_multiple: # seen exactly once before now
            seen_multiple.add(start)
            update_range(start, count, (-1, 1, 0), (unique, multi, unmapped))
        else:
            pass # seen more than once - do nothing
    elif mapq > 0: # position specified, not seen before, positive mapq
        seen_once.add(start)
        update_range(start, count, (1, 0, 0), (unique, multi, unmapped))
    else: # position specified, not seen before, mapq = 0 -> assume it's a multi map 
        seen_once.add(start)
        seen_multiple.add(start)
        update_range(start, count, (0, 1, 0), (unique, multi, unmapped))
 
def analyze_sam(sam_fh, unique, multi, unmapped, log_fh):
    seen_once = set() # seen once
    seen_multiple = set() # seen more than once
    for line in sam_fh:
        if line.startswith('@'): # skip comments
            continue
        fields = line.strip('\n').split('\t')
        original_pos = int(fields[0].split('_')[1]) # ignore chromosome
        #flag = int(fields[1])
        new_pos = int(fields[3])
        mapq = int(fields[4])
        if new_pos == 0: # unmapped
            update_range(original_pos, len(fields[9]), (0, 0, 1), (unique, multi, unmapped))
        else:
            handle_mapped(original_pos, len(fields[9]), mapq, seen_once, seen_multiple, unique, multi, unmapped)
       
def write_bedgraph_header(output):
    #output.write('track type=bedGraph name=track_label description=center_label visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16\n')
    output.write('track type=bedGraph\n')

def mean(values):
    result = 0.
    for value in values:
        result += value
    return 1. * result / len(values)

def write_stats(out, unique, multi, unmapped, total):
    out.write('{0: >10}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('Type', 'Count', 'Max', 'Min', 'Mean', 'Mean %'))
    out.write('{0: >10}\t{1}\t{2}\t{3}\t{4:.2f}\t{5:.2f}\n'.format('Unique', len(unique), max(unique.values()), min(unique.values()), mean(unique.values()), 100. * mean(unique.values()) / total))
    out.write('{0: >10}\t{1}\t{2}\t{3}\t{4:.2f}\t{5:.2f}\n'.format('Multi', len(multi), max(multi.values()), min(multi.values()), mean(multi.values()), 100. * mean(multi.values()) / total))
    out.write('{0: >10}\t{1}\t{2}\t{3}\t{4:.2f}\t{5:.2f}\n'.format('Unmapped', len(unmapped), max(unmapped.values()), min(unmapped.values()), mean(unmapped.values()), 100. * mean(unmapped.values()) / total))

def write_tsv(out, unique, multi, unmapped):
    out.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format('Position', 'Unique', 'Multi', 'Unmapped', 'Total'))
    all_keys = set(unique.keys() + multi.keys() + unmapped.keys())
    for key in sorted(all_keys):
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(key, unique[key], multi[key], unmapped[key], unique[key] + multi[key] + unmapped[key]))

def analyze(source, references, kmer, resolution, output, tsv_out, log_fh, stage=0):
    '''
        align source to targets
        evaluate and generate data
    '''
    progress = 0
    # stage 0 generate reads for single source
    chromosomes = []
    if progress >= stage:
        fastq = '{0}.fastq'.format(source)
        generate_reads(open(source, 'r'), kmer, resolution, open(fastq, 'w'), chromosomes, log_fh)
    else:
        with open(source, 'r') as fasta_in:
            line = fasta_in.readline()
            chromosomes.append(line[1:].split(' ')[0].strip())
    write_log(log_fh, 'stage {0} complete'.format(progress))
    progress += 1
    
    # stage 1: index reference and align source to it
    if progress >= stage:
        for idx, reference in enumerate(references): # each reference genome
            index_reference(reference, log_fh) # index
            sam = '{0}_{1}.sam'.format(source, idx)
            align_reference(fastq, reference, sam, log_fh) # align source to this reference
    write_log(log_fh, 'stage {0} complete'.format(progress))
    progress += 1

    # stage 2: analyze
    sams = []
    unique = collections.defaultdict(int) # uniquely mapped to reference
    multi = collections.defaultdict(int) # mapped to multiple locations
    unmapped = collections.defaultdict(int) # didn't map to reference
    write_log(log_fh, 'references: {0}'.format(references))
    overlaps = kmer / resolution
    for idx, reference in enumerate(references): # each reference
        sam = '{0}_{1}.sam'.format(source, idx)
        sams.append(sam)
        write_log(log_fh, 'analyzing sam file {0}: {1}...'.format(idx, sam))
        analyze_sam(open(sam, 'r'), unique, multi, unmapped, log_fh) # analyze: update unique/multi/unmapped counts
        write_stats(log_fh, unique, multi, unmapped, overlaps * len(sams)) # general stats
        write_log(log_fh, 'analyzing sam file: done')

    # write bedgraph result
    write_log(log_fh, 'writing bedgraph...')
    write_bedgraph_header(output)
    for key, value in unique.items():
        pos = int(key)
        score = 1. * overlaps * value / len(references) # normalize
        output.write('{0} {1} {2} {3:.3f}\n'.format(chromosomes[0], pos, pos + 1, score))
    write_log(log_fh, 'writing bedgraph: done')

    write_log(log_fh, 'stats...')
    write_stats(log_fh, unique, multi, unmapped, overlaps * len(sams))

    if tsv_out:
        write_log(log_fh, 'tsv...')
        write_tsv(open(tsv_out, 'w'), unique, multi, unmapped)


def main():
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='Calculate population wide reference bias')
    parser.add_argument('--source', required=True, help='source reference file')
    parser.add_argument('--targets', required=True, nargs='+', help='target reference files')
    parser.add_argument('--kmer', type=int, default=100, help='size of kmers')
    parser.add_argument('--resolution', type=int, default=10, help='resolution of analysis (1 means analyze every base)')
    parser.add_argument('--out', required=False, help='tsv file of combined results')
    parser.add_argument('--stage', required=False, default=0, type=int, help='stage to start from (0=generate reads,1=index and align,2=analyze)')
    args = parser.parse_args()
    analyze(args.source, args.targets, args.kmer, args.resolution, sys.stdout, args.out, sys.stderr, args.stage)

if __name__ == '__main__':
    main()

