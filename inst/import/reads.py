#!/usr/bin/env python3

import pysam
import pandas as pd
import numpy as np
import re
import numbers
from collections import defaultdict

########################################################################
class Read:
    def __init__(self, 
                 read1, 
                 read2, 
                 umi1_idx=None, 
                 umi2_idx=None, 
                 min_mapq=30, 
                 max_insert=1000, 
                 max_no_conv=3,
                 infer_umi_field=False,
                 template_conv_tag='XT',
                 meth_calls_tag='XP',
                 bismark_flags=False):

        self.pairing = self._get_pairing(read1, read2, min_mapq, max_insert)
        self.id = read1.query_name
        
        self.template_conv_tag = template_conv_tag
        self.meth_calls_tag = meth_calls_tag 
        self.bismark_flags = bismark_flags
        
        if self.pairing in ['good', 'single_R1', 'single_R2']:
            if self.pairing == 'single_R2':
                read1, read2 = read2, read1

            self.umi1, self.umi2 = self._get_umis(read1.query_name, umi1_idx, umi2_idx, infer_umi_field)
            self.chrom = read1.reference_name            
            self.conv = read1.get_tag(template_conv_tag)

            self.start1 = read1.pos
            self.qual1 = read1.qual
            self.patt1 = read1.get_tag(self.meth_calls_tag)
            
            self.H = 0 
            self.h = 0
            self.x = 0
            self.X = 0
            self.z = 0
            self.Z = 0           
            
            self._get_conv(self.patt1)

            self.insert_len = read1.tlen

            if self.single:
                self._process_single(read1)
            else:
                self._process_double(read1, read2)

            if self.conv == 'GA':
                self.start1 -= 1
                if not self.single:
                    self.start2 -= 1
            
            if self.H + self.X > max_no_conv:
                self.pairing = 'no_conv'
    
    def _get_conv(self, patt):
        self.H += patt.count('H')
        self.h += patt.count('h')
        self.x += patt.count('x')
        self.X += patt.count('X')
        self.z += patt.count('z')
        self.Z += patt.count('Z')
        
    def _process_single(self, read1):
        self.umi2 = '-'
        self.start = read1.pos
        self.end = '-'
        if read1.is_reverse:
            self.strand = '-'
            if not self.bismark_flags:
                self.patt1 = self.patt1[::-1]
                self.qual1 = self.qual1[::-1]
        else:
            self.strand = '+'

    def _process_double(self, read1, read2):
        self.start2 = read2.pos
        self.qual2 = read2.qual
        self.patt2 = read2.get_tag(self.meth_calls_tag)
        self._get_conv(self.patt2)
        if read1.is_reverse:
            self.start = read2.pos
            self.end = read1.pos + read1.alen
            self.strand = '-'

            [self.start1, self.start2] = [self.start2, self.start1]
            [self.qual1, self.qual2] = [self.qual2, self.qual1]
            [self.patt1, self.patt2] = [self.patt2, self.patt1]

        else:
            self.start = read1.pos
            self.end = read2.pos + read2.alen
            self.strand = '+'
        
        if not self.bismark_flags:
            self.patt2 = self.patt2[::-1]
            self.qual2 = self.qual2[::-1]

            if read1.is_reverse and read2.is_reverse:
                self.patt2 = self.patt2[::-1]
                self.qual2 = self.qual2[::-1]

    def _get_umis(self, read_id, umi1_idx=None, umi2_idx=None, infer_umi_field=False):
        read_id = read_id.split(':')
        if umi1_idx is None:
            if infer_umi_field:
                umi1 = [x for x in read_id if re.search(r'^umi=', x)]
                if len(umi1) == 0:
                    umi1 = '-'
                else:
                    umi1 = umi1[0].translate(None, 'umi=')
            else:
                umi1 = '-'
        else:
            umi1 = read_id[umi1_idx].translate(None, 'umi=')

        if umi2_idx is None:
            if infer_umi_field:
                umi2 = [x for x in read_id if re.search(r'^umi1=', x)]
                if len(umi2) == 0:
                    umi2 = '-'
                else:
                    umi2 = umi2[0].translate(None, 'umi1=')
            else:
                umi2 = '-'
        else:
            umi2 = read_id[umi2_idx].translate(None, 'umi1=')
        return [umi1, umi2]

    def _valid_cigar(self, cigar):
        pattern = re.compile("^[0-9]+M$")
        return pattern.match(cigar)

    def _get_pairing(self, read1, read2, min_mapq=30, max_insert=1000):
        self.single = False

        # single reads were marked in the past by 512 flag / read2 starting with '-'. 
        # Since we currently to not know how to deal with 2 reads mapped to the same strand we treat it also as single
        if read2 is not None:
            if (read2.flag & 512) or (read1.query_name != read2.query_name) or (read1.is_reverse == read2.is_reverse):
                read2 = None

        if read2 is None:
            if read1.mapping_quality < min_mapq:
                return 'unmapped'
            if not self._valid_cigar(read1.cigarstring):
                return 'bad_cigar'
            self.single = True
            return 'single_R1'
        if read1 is None:
            if read2.mapping_quality < min_mapq:
                return 'unmapped'
            if not self._valid_cigar(read2.cigarstring):
                return 'bad_cigar'
            self.single = True
            return 'single_R2'

        pairing1 = read1.get_tag('YT')
        pairing2 = read2.get_tag('YT')
        if pairing1 != pairing2:
            raise Exception('Unpaired paring:' + pairing1 + ' ' + pairing2)

        if pairing1 == 'DP':
            return 'discordant'

        valid_R1 = read1.mapping_quality >= min_mapq and self._valid_cigar(read1.cigarstring)
        valid_R2 = read2.mapping_quality >= min_mapq and self._valid_cigar(read2.cigarstring)

        if valid_R1 and valid_R2:
            if read1.reference_name != read2.reference_name or abs(read1.pos - read2.pos) > max_insert:
                return 'discordant'
            return 'good'
        if valid_R1 and not valid_R2:
            self.single = True
            return 'single_R1'
        if valid_R2 and not valid_R1:
            self.single = True
            return 'single_R2'
        if read1.mapping_quality < min_mapq and read2.mapping_quality < min_mapq:
            return 'unmapped'

        return 'bad_cigar'

    def in_range(self, genomic_range):
        start_in_range = isinstance(self.start, numbers.Number) and (genomic_range[0] <= self.start <= genomic_range[1])
        end_in_range = isinstance(self.end, numbers.Number) and (genomic_range[0] <= self.end <= genomic_range[1])
        return start_in_range or end_in_range

########################################################################
def bam_iter(bam, 
             umi1_idx=None, 
             umi2_idx=None, 
             infer_umi_field=False, 
             single_end=False, 
             min_mapq=30, 
             max_no_conv=3,
             template_conv_tag='XT', 
             meth_calls_tag='XP', 
             bismark_flags=False):    
    if single_end:
        for read in bam:
            yield Read(read, None, umi1_idx, umi2_idx, infer_umi_field=infer_umi_field, min_mapq=min_mapq, max_no_conv=max_no_conv, template_conv_tag=template_conv_tag, meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags)
    else:
        for read in bam:
            if read.is_read1:
                read1 = read
            if read.is_read2 or read.flag & 772:
                yield Read(read1, read, umi1_idx, umi2_idx, infer_umi_field=infer_umi_field, min_mapq=min_mapq, max_no_conv=max_no_conv, template_conv_tag=template_conv_tag, meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags)
