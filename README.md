
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

## Usage

```
python population_mappability.py --source source.fasta  --targets target1.fasta,target2.fasta,...,targetn.fasta --kmer 100 --resolution 100
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
