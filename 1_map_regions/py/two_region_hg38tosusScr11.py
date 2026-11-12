#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import region_module as rm
from collections import defaultdict

inputfile = './processed_data/indexed_dhs_hg38_susScr11_mapped/hg38_susScr11_core_set.txt'
outputfile = './intermediate_data/hg38_susScr11_percentage.txt'
o = open(outputfile,'w')
o.write('identifier\tref_region\tcmp_region\taxt_id\tdirection\tmatched_bases\tpercentage\n')
error = []

# coordiantes transformation if alignment on negative strand
rev_spec = 'susScr11'
#load chromosome information of aligned species
dict_fp = './1_map_regions/chrom_size/' + rev_spec + '.chrom.sizes.txt'
rev_chrom_dict = rm.chrom_size(dict_fp)
# genome-wide axt net file, blocks split by hg38 chromosome
axt_fp = './raw_data/hg38_susScr11/hg38.susScr11.net.axt'
chrom_block_dict = defaultdict(dict)
for each in rm.readblock(axt_fp):
    block_dict = chrom_block_dict[each[0].split()[1]]
    block_dict[len(block_dict)] = each[:]

with open(inputfile,'r') as i:
    next(i)
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[0]
            ref_region = linelist[3]
            chrom = ref_region.split(':')[0]
            start = int(ref_region.split(':')[1].split('-')[0])
            end = int(ref_region.split(':')[1].split('-')[1])
            mylist = [int(start),int(end)]
            block_dict = chrom_block_dict[chrom]

            con,posi = rm.binarySearch(0, len(block_dict)-1, mylist, block_dict)
            if con == 'in':
                out = rm.outputline_ref_percentage([posi], mylist, block_dict, rev_chrom_dict)
                if out != 'no_match\n':
                    final_line = index + '\t' + out
                    o.write(final_line)

            if con in ['ovlow','ovhigh','out']:
                posilist = rm.contSearch(con, posi, mylist, block_dict)
                out = rm.outputline_ref_percentage(posilist, mylist, block_dict, rev_chrom_dict)
                if out != 'no_match\n':
                    final_line = ''
                    for eachline in out.split('\n'):
                        if eachline != '':
                            final_line = final_line + index + '\t' + eachline + '\n'
                    o.write(final_line)

        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
