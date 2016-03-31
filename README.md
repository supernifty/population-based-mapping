
# Population-based mapping to measure reference bias

## Synopsis

Using a population of references (target references) we map each k-mer from one reference (source reference) to each alternative reference (including itself).

We then look at the resulting set of how mapping went for each k-mer. Possible outcomes for each reference:
* didn't map
* mapped uniquely
* mapped to multiple locations

Then for each location (k-mer) on the source reference, can visualize its mappability based on the proportion of alignments that map uniquely (high confidence), not at all, or to multiple locations.

We can also generate a population mappability track.

## Installation

## Testing
To run unit tests:
```
./test.sh
```

## Usage

Generate statistics by running the following:
```
python src/main.py --source source.fasta  --targets target1.fasta target2.fasta ... targetn.fasta --kmer 100 --resolution 100 --out all.tsv > map.bedgraph
```

Command used for e-coli analysis:
```
python src/main.py --source tests/e-coli-mg1655.fasta --targets data/ecoli/*.fasta --kmer 100 --resolution 100 --out e-coli-62.tsv > e-coli-62.bedgraph
python src/main.py --source tests/e-coli-mg1655.fasta --targets data/ecoli/*.fasta --kmer 100 --resolution 10 --out e-coli-62.tsv > e-coli-62.bedgraph
```

Command used for staph analysis:
```
python src/main.py --source data/staph/1.fasta --targets data/staph/*.fasta --kmer 100 --resolution 10 --out staph.160331.tsv > staph.160331.bedgraph
```

Visualization by running:
```
python src/plot.py --start 200000 --finish 300000 --out zoomed.pdf < e-coli-62.tsv
```

## Previous work

## Future work

This implementation uses BWA. We could compare to other aligners.

Look for correlation with this track to 
* exon/intron
* repeat tracks

Pairwise comparisons using Mauve?

## Links and references

* bedgraph format: https://genome.ucsc.edu/goldenPath/help/bedgraph.html
