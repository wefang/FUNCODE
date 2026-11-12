#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import base_module as bm
from collections import defaultdict

inputfile = './processed_data/pig_atac_halstead/pig_dhs_200bp.bed'
outputfile = './processed_data/indexed_dhs_hg38_susScr11_mapped/susScr11tohg38_base_alignment.txt'
o = open(outputfile, 'w')
o.write('identifier\tsusScr11_summit\thg38_summit\tsusScr11_extend_region\thg38_extend_region\taxtnet_id\taxtnet_direction\tgap_status\n')
error = []

# bm.outputline without the coordinate transformation on negative strand
def outputline_noflip(posi, base, setlen, blockdict, ref_chrom_dict, cmp_chrom_dict):
    theline = blockdict[posi][0].split()
    align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]),theline[7]
    seq_dict = bm.seqdict(blockdict[posi][1],blockdict[posi][2])
    ref_chrom_len = ref_chrom_dict[ref_chrom]
    cmp_chrom_len = cmp_chrom_dict[cmp_chrom]
    for key,value in seq_dict.items():
        if ref_start + value['ref_trueindex'] == base:
            out_cmp = cmp_start + value['cmp_trueindex']
            if value['cmp_base'] != '-':
                status = 'non-gap'
            else:
                status = 'gap'
                closest, thekey = bm.searchclosest(key, seq_dict)
                if closest == 'after':
                    out_cmp = cmp_start + value['cmp_trueindex'] + 1
            break
    my_base = ref_chrom + ':' + str(base)
    out_posi = cmp_chrom + ':' + str(out_cmp)
    my_region = ref_chrom + ':' + '-'.join(bm.extend_region(base, setlen, ref_chrom_len))
    out_region = cmp_chrom + ':' + '-'.join(bm.extend_region(out_cmp, setlen, cmp_chrom_len))
    return my_base + '\t' + out_posi + '\t' + my_region + '\t' + out_region + '\t' + align_id + '\t' + align_dire + '\t' + status + '\n'

ref_spec = 'susScr11'
cmp_spec = 'hg38'
#load chromosome information
ref_chrom_dict = bm.chrom_size('./1_map_regions/chrom_size/' + ref_spec + '.chrom.sizes.txt')
cmp_chrom_dict = bm.chrom_size('./1_map_regions/chrom_size/' + cmp_spec + '.chrom.sizes.txt')
# hg38 to susScr11 axt net file searched on the susScr11 side:
# swap the two species in each block and sort blocks by susScr11 start
axt_fp = './raw_data/hg38_susScr11/' + cmp_spec + '.' + ref_spec + '.net.axt'
chrom_block_list = defaultdict(list)
for each in bm.readblock(axt_fp):
    h = each[0].split()
    chrom_block_list[h[4]].append([' '.join([h[0]] + h[4:7] + h[1:4] + h[7:]), each[2], each[1]])
chrom_block_dict = defaultdict(dict)
for chrom, blocks in chrom_block_list.items():
    blocks.sort(key=lambda b: int(b[0].split()[2]))
    chrom_block_dict[chrom] = dict(enumerate(blocks))
# extent region length
setlen = 200

with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[3]
            chrom = linelist[0]
            base = (int(linelist[1]) + int(linelist[2])) // 2 + 1
            block_dict = chrom_block_dict[chrom]
            posi = bm.binarySearch(0, len(block_dict)-1, base, block_dict)
            if posi != 'no match':
                out = outputline_noflip(posi, base, setlen, block_dict, ref_chrom_dict, cmp_chrom_dict)
                final_line = index + '\t' + out
                o.write(final_line)
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
