#!/usr/bin/env python
# Modified by Chang Xu, on top of the original version.
# Python version of the perl script originally developed by Markus Rasmussen, Sebastian DiLorenzo (mostly borrowed @ Seqanswers)
# samtools mpileup -f <reference>.fa <tumorfile>.bam > <output>.pileup
# samtools mpileup -uf <reference>.fa <tumorfile>.bam | bcftools view -vcg - > <output>.vcf
# usage: this.py <vcf> <hg_build>

import sys
import gzip

with gzip.open(sys.argv[1]) as vcf:
    try:
        for l in vcf:
            if l[0] == '#':
                continue
            sp = l.strip().split('\t')
            ref = sp[3].upper()
            alt = sp[4].upper()
            if len(alt) > 1 or len(ref) > 1 or alt == 'N':
                continue
            for tag in sp[7].split(';'):
                if tag.startswith('DP='):
                    depth = int(tag[3:])
                if tag.startswith('DP4='):
                    dp4_alt = sum(map(int, tag[4:].split(',')[-2:]))
            pct = dp4_alt * 1.0 / depth
            print '{}\t{}\t{} > {}\t{}\t{}\t{}\t{:.2f}'.format(sp[0], sp[1], 
                                                               ref, alt, 
                                                               sp[5], depth, dp4_alt, pct)
    except StopIteration:
        pass
