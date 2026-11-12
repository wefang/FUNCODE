#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import base_module as bm
from collections import defaultdict

inputfile = './metadata/ENCODE_DHS/ENCFF503GCK.tsv'
outputfile = './processed_data/indexed_dhs_hg38_susScr11_mapped/hg38tosusScr11_base_alignment.txt'
o = open(outputfile, 'w')
o.write('identifier\thg38_summit\tsusScr11_summit\thg38_extend_region\tsusScr11_extend_region\taxtnet_id\taxtnet_direction\tgap_status\n')
error = []

# coordiantes transformation if alignment on negative strand
ref_spec = 'hg38'
cmp_spec = 'susScr11'
#load chromosome information
ref_chrom_dict = bm.chrom_size('./1_map_regions/chrom_size/' + ref_spec + '.chrom.sizes.txt')
cmp_chrom_dict = bm.chrom_size('./1_map_regions/chrom_size/' + cmp_spec + '.chrom.sizes.txt')
# genome-wide axt net file, blocks split by hg38 chromosome
axt_fp = './raw_data/hg38_susScr11/' + ref_spec + '.' + cmp_spec + '.net.axt'
chrom_block_dict = defaultdict(dict)
for each in bm.readblock(axt_fp):
    block_dict = chrom_block_dict[each[0].split()[1]]
    block_dict[len(block_dict)] = each[:]
# extent region length
setlen = 200

with open(inputfile,'r') as i:
    next(i)
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[3]
            chrom = linelist[0]
            # no chrY, and summit - 1 is the center of the 200 bp window, as in dhs_summit_extended of human_dhs_noY.rda
            if chrom == 'chrY':
                continue
            base = int(linelist[6]) - 1
            block_dict = chrom_block_dict[chrom]
            posi = bm.binarySearch(0, len(block_dict)-1, base, block_dict)
            if posi != 'no match':
                out = bm.outputline(posi, base, setlen, block_dict, ref_chrom_dict, cmp_chrom_dict)
                final_line = index + '\t' + out
                o.write(final_line)
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
