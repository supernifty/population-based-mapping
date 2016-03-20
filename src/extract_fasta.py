#!/usr/bin/env python
'''
    simple script to extract combined fasta references out into multiple files
    write n.fasta in the current directory for each sequence in sys.stdin
'''

import sys

count = 0
for line in sys.stdin:
    if line.startswith('>'):
        count += 1
        target = open('{0}.fasta'.format(count), 'w')
    target.write(line)

