#!/usr/bin/env python3

"""
Sums the bitscores of identical query/target sequence pairs in BLAST results.
"""

import sys

# Parse DIAMOND output generated with --outfmt 6 qseqid sseqid bitscore
with open(sys.argv[1]) as fp:
    x = {}
    
    sources = set()
    targets = set()
    for line in fp:
        source, target, score = line.strip().split()
        score = float(score)
        
        sources.add(source)
        targets.add(target)

        # Save summed bitscores for each pair in nested dictionary
        if source in x:
            if target in x[source]:
                x[source][target] += score
            else:
                x[source][target] = score
        else:
            x[source] = {}
            x[source][target] = score

with open(sys.argv[2], "w") as fp:
    for source, targets in x.items():
        for target, score in targets.items():
            fp.write(f"{source},{target},{score}\n")
