#!/usr/bin/env python2
from __future__ import print_function
import argparse
from textwrap import dedent
import sys
import pysam
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from tqdm import tqdm

#TODO: REMEMBER TO REMOVE 772/512 FLAGS BEFORE PRODUCTION

########################################################################
PROG = 'tidy_cpgs.py'

########################################################################
class Read:
    def __init__(self, 
                 read1, 
                 read2, 
                 umi1_idx=None, 
                 umi2_idx=None, 
                 min_mapq=30, 
                 max_insert=1000, 
                 infer_umi_field=False,
                 template_conv_tag='XT',
                 meth_calls_tag='XP',
                 bismark_flags=False):

        self.pairing = self._get_pairing(read1, read2, min_mapq, max_insert)
        
        self.template_conv_tag = template_conv_tag
        self.meth_calls_tag = meth_calls_tag 
        self.bismark_flags = bismark_flags
        
        if self.pairing in ['good', 'single_R1', 'single_R2']:
            if self.pairing == 'single_R2':
                read1, read2 = read2, read1

            self.umi1, self.umi2 = self._get_umis(read1.query_name, umi1_idx, umi2_idx, infer_umi_field)
            self.chrom = read1.reference_name
            self.id = read1.query_name
            self.conv = read1.get_tag(template_conv_tag)

            self.start1 = read1.pos
            self.qual1 = read1.qual
            self.patt1 = read1.get_tag(self.meth_calls_tag)

            self.insert_len = read1.tlen

            if self.single:
                self._process_single(read1)
            else:
                self._process_double(read1, read2)

            if self.conv == 'GA':
                self.start1 -= 1
                if not self.single:
                    self.start2 -= 1

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
                umi1 = filter(lambda x: re.search(r'^umi=', x), read_id)
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
                umi2 = filter(lambda x: re.search(r'^umi1=', x), read_id)
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


########################################################################
def bam_iter(bam, 
             umi1_idx=None, 
             umi2_idx=None, 
             infer_umi_field=False, 
             single_end=False, 
             min_mapq=30, 
             template_conv_tag='XT', 
             meth_calls_tag='XP', 
             bismark_flags=False):
    if single_end:
        for read in bam:
            yield Read(read, None, umi1_idx, umi2_idx, infer_umi_field=infer_umi_field, min_mapq=min_mapq, template_conv_tag=template_conv_tag, meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags)
    else:
        for read in bam:
            if read.is_read1:
                read1 = read
            if read.is_read2 or read.flag & 772:
                yield Read(read1, read, umi1_idx, umi2_idx, infer_umi_field=infer_umi_field, min_mapq=min_mapq, template_conv_tag=template_conv_tag, meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags)

########################################################################
def bam_reader(bam, chunk_size, stats, umi1_idx=None, umi2_idx=None, infer_umi_field=False, chrom=None,
               genomic_range=None, single_end=False, show_progress=True, min_mapq=30, 
               template_conv_tag='XT', meth_calls_tag='XP', add_chr_prefix=False,
               bismark_flags=False):
    while True:
        ids = []
        chroms = []
        starts = []
        ends = []
        strands = []
        umis1 = []
        umis2 = []
        insert_lengths = []
        read_starts = []
        patts = []
        quals = []

        for read in tqdm(
                bam_iter(bam, umi1_idx=umi1_idx, umi2_idx=umi2_idx, infer_umi_field=infer_umi_field,
                         single_end=single_end, min_mapq=min_mapq, template_conv_tag=template_conv_tag, 
                         meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags),
                unit='reads',
                unit_scale=True, disable=not show_progress):

            stats[read.pairing] += 1            

            if read.pairing in ['good', 'single_R1', 'single_R2'] and \
                    ((chrom is None) or \
                     (not add_chr_prefix and chrom == read.chrom) or \
                     (add_chr_prefix and chrom == 'chr' + read.chrom)):
                    
                if (genomic_range is None) \
                        or (genomic_range[0] <= read.start <= genomic_range[1]) \
                        or (genomic_range[0] <= read.end <= genomic_range[1]):

                    ids.append(read.id)
                    chroms.append(read.chrom)
                    starts.append(read.start)
                    ends.append(read.end)
                    strands.append(read.strand)
                    umis1.append(read.umi1)
                    umis2.append(read.umi2)
                    insert_lengths.append(read.insert_len)
                    read_starts.append(read.start1)
                    patts.append(read.patt1)
                    quals.append(read.qual1)

                    if not read.single:
                        ids.append(read.id)
                        chroms.append(read.chrom)
                        starts.append(read.start)
                        ends.append(read.end)
                        strands.append(read.strand)
                        umis1.append(read.umi1)
                        umis2.append(read.umi2)
                        insert_lengths.append(read.insert_len)
                        read_starts.append(read.start2)
                        patts.append(read.patt2)
                        quals.append(read.qual2)

                    if len(ids) >= chunk_size:
                        break

        if not chroms:
            break

        ids = np.array(ids, dtype=np.str)
        chroms = np.array(chroms, dtype=np.str)
        starts = np.array(starts, dtype=np.int)
        ends = np.array(ends, dtype=np.str)
        strands = np.array(strands, dtype=np.str)
        umis1 = np.array(umis1, dtype=np.str)
        umis2 = np.array(umis2, dtype=np.str)
        insert_lengths = np.array(insert_lengths, dtype=np.int)
        read_starts = np.array(read_starts, dtype=np.int)

        # reshape patts and quals to nXseq_len array, where n=number of patterns 
        # and seq_len is the longest pattern (shorter patterns / quals would be padded)
        seq_len = len(max(patts, key=len))
        patts = np.array(patts)
        patts = patts.view('S1').reshape((-1, seq_len))
        quals = np.array(quals).view(np.uint8).reshape((-1, seq_len))

        ids_df = pd.DataFrame({
            'read_id': ids,
            'chrom': chroms,
            'start': starts,
            'end': ends,
            'strand': strands,
            'umi1': umis1,
            'umi2': umis2,
            'insert_len': insert_lengths,
            'read_start': read_starts,
            'index': range(0, len(chroms), 1)})

        yield ids_df, patts, quals


########################################################################
def extract_cgs(ids_df, patts, quals, columns, min_qual=10, genomic_range=None, cgs_mask=None, trim=None, add_chr_prefix=False):
    pats_df = pd.DataFrame(patts, index=ids_df.index).reset_index()
    pats_df = pd.melt(pats_df, var_name='pos', value_name='meth', id_vars=['index'])
    pats_df['qual'] = quals.ravel()

    pats_df = pats_df[((pats_df.meth == 'z') | (pats_df.meth == 'Z')) & (pats_df.qual >= min_qual)]
    df = pd.merge(pats_df, ids_df, how='left', on='index')

    df['cg_pos'] = df.read_start + df.pos

    if genomic_range is not None:
        df = df[(df['cg_pos'] >= genomic_range[0]) & (df['cg_pos'] <= genomic_range[1])]

    if trim is not None:
        df = df[abs(df.start - df.cg_pos) > trim]
        df = df[(abs(pd.to_numeric(df.end, errors='coerce') - df.cg_pos) > trim) | (df.end == '-')]

    if add_chr_prefix:
        df['chrom'] = 'chr' + df['chrom'].astype(str)
        
    if cgs_mask is not None:
        df = pd.merge(df, cgs_mask, how='left', copy=False)
        df = df[pd.isnull(df.f)]         

    df = df[columns]

    return df


########################################################################
def read_cgs_mask(file):
    cgs_mask = pd.read_csv(file, usecols=['chrom', 'start']).rename(columns={"start": "cg_pos"})
    cgs_mask['f'] = True
    return cgs_mask


########################################################################
def main(argv):
    try:
        args = parse_args(argv[1:])
    except ValueError as err:
        eprint("Error parsing arguments: %s" % err)
        return 2
    
    if args.bismark:
        args.template_conv_tag = 'XG'
        args.meth_calls_tag = 'XM'
        args.bismark_flags = True
    
    if args.cgs_mask is not None:
        cgs_mask = read_cgs_mask(args.cgs_mask)
    else:
        cgs_mask = None
        
    stats = defaultdict(int)

    columns = ['read_id', 'chrom', 'start', 'end', 'strand', 'umi1', 'umi2', 'insert_len', 'cg_pos', 'meth', 'qual']
    with pysam.AlignmentFile(args.input, "rb") as in_bam, \
            open(args.output, 'w') if args.output is not '-' else sys.stdout as out_file:
        stats_out = open(args.stats, 'w') if args.stats is not '-' else sys.stderr

        out_file.write(','.join(columns) + '\n')

        for ids_df, patts, quals in bam_reader(in_bam,
                                               chunk_size=args.chunk_size,
                                               stats=stats,
                                               umi1_idx=args.umi1_idx,
                                               umi2_idx=args.umi2_idx,
                                               infer_umi_field=args.infer_umi_field,
                                               chrom=args.chrom,
                                               genomic_range=args.genomic_range,
                                               single_end=args.single_end,
                                               show_progress=not args.no_progress,
                                               min_mapq=args.min_mapq,
                                               template_conv_tag=args.template_conv_tag,
                                               meth_calls_tag=args.meth_calls_tag,
                                               add_chr_prefix=args.add_chr_prefix,
                                               bismark_flags=args.bismark_flags):
            cpgs = extract_cgs(ids_df, patts, quals, columns=columns,
                               min_qual=args.min_qual,
                               genomic_range=args.genomic_range,
                               cgs_mask=cgs_mask,
                               trim=args.trim,
                               add_chr_prefix=args.add_chr_prefix)

            if args.sort_chunk:
                cpgs.sort_values(by=['chrom', 'start', 'end', 'strand', 'umi1', 'umi2'], inplace=True)
            cpgs.to_csv(out_file, header=False, index=False, float_format='%.0f', mode='a')

        stats_out.write(','.join(list(stats.keys())) + '\n')
        stats_out.write(','.join(str(x) for x in list(stats.viewvalues())))

    return 0


########################################################################
def parse_args(argv):    
    parser = argparse.ArgumentParser(prog=PROG,
                                     description = 'Generates a tidy csv file from bissli bam.',
                                     epilog=dedent('''
                                     Output file has the following fields for each CpG:      
                                     
                                       read_id       id of the read of the CpG
                                       chrom         chromosome of the read
                                       start         start position of the read
                                       end           end position of the read ('-' for single end/singly mapped reads)
                                       strand        strand position of the read
                                       umi1          umi sequence on read1 ('-' when there is no umi1 field)
                                       umi2          umi sequence on read2 ('-' when there is no umi2 field)
                                       insert_len    insert length of the read
                                       cg_pos        genomic position of the CpG
                                       meth          methylation call of the CpG 'Z' for methylated 'z' for unmethylated
                                       qual          quality value of the base of the methylation call
                                                     ('C' base for CT conversion, 'G' base for GA conversion)  
                                         '''),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                    usage='%(prog)s [options] -i <bam> -o <output> -s <stats>')
    
    io = parser.add_argument_group('io')
    parser.add_argument('-i', '--input', type=str, default='-', metavar='<bam>',
                        help="Name of input bam/sam file. a value of '-' indicates that the file would be read from STDIN. (default: %(default)s)")
    parser.add_argument('-o', '--output', type=str, default='-', metavar='<output>',
                        help="Name of output file. a value of '-' indicates that the output would be written to STDOUT. (default: %(default)s)")
    parser.add_argument('-s', '--stats', type=str, default='-', metavar='<stats>',
                        help="Name of stats file. a value of '-' indicates that the output would be written to STDERR. (default: %(default)s)")
    
    umi = parser.add_argument_group('umi')
    umi.add_argument('--umi1-idx', metavar='UMI1_INDEX_FIELD', type=int, default=None,
                        help="position of umi1 in index (0 based), when delimiter is ':'. e.g. in  HISEQ:164:C6H4FANXX:1:1113:2797:umi=ACTAGGA:26364_1:N:0:CGATGT, use --umi1-idx 6")
    umi.add_argument('--umi2-idx', metavar='UMI2_INDEX_FIELD', type=int, default=None,
                        help="position of umi2 in index (0 based). e.g. in  HISEQ:164:C6H4FANXX:1:1113:2797:umi=ACTAGGA:umi1=GGATATA:26364_1:N:0:CGATGT, use --umi2-idx 7")
    umi.add_argument('--infer-umi-field', action='store_true',
                        help="infer position of umi in index. Position is inferred by looking for a field that starts with 'umi=' for umi1 and 'umi1=' for umi2. Note that this may slow the process.")   
      
    genomic_limits = parser.add_argument_group('genomic_limits', 'limit the genomic scope (e.g. for parallelization)')
    genomic_limits.add_argument('--chrom', type=str, default=None,
                        help="process only reads from a given chromosome")
    genomic_limits.add_argument('--genomic-range', nargs=2, type=int, default=None, 
                                metavar=('START', 'END'), 
                                help='process only reads from genomic coordinates in the chromosome specified in --chrom')
    
    protocol = parser.add_argument_group('protocol', 'protocol specific options')
    protocol.add_argument('--cgs-mask', type=str, default=None, metavar='FILE',
                        help="comma separated file with positions of cpgs to mask (e.g. MSP1 sticky ends). File needs to have 'chrom' and 'start' fields with the chromosome and position of 'C' in the cpgs to mask")    
    protocol.add_argument('--trim', metavar='BASE_PAIRS', type=int, default=None,
                        help="trim cpgs that are --trim basepairs from the beginning/end of the read (e.g. first bp in RRBS data)")
    protocol.add_argument('--single-end', action='store_true', help='single end reads')
    
    
    mapper = parser.add_argument_group('mapper', 'mapper output options (bismark/bissli)')
    mapper.add_argument('--template-conv-tag', type=str, default='XT',
                        help="sam tag for template conversion (XG for bismark, XT for bissli). (default: %(default)s)")
    mapper.add_argument('--meth-calls-tag', type=str, default='XP',
                        help="sam tag for methylation calls (XM for bismark, XP for bissli). (default: %(default)s)") 
    mapper.add_argument('--bismark-flags', action='store_true',
                        help='flag field was generated by bismark ("tries to take the strand a bisulfite read originated from into account"). Reverse strands would be treated the same as forward strand. In bissli mode, the methylation string of reverse strand is reversed')
    mapper.add_argument('--bismark', action='store_true', help='default parameters for bismark alignments (same as --template_conv_tag XG --meth_calls_tag XM --bismark-flags)')    
    mapper.add_argument('--add-chr-prefix', action='store_true', help='add "chr" prefix for chromosomes (e.g. in order to import to misha)')   
    
    other = parser.add_argument_group('other')
    other.add_argument('--no-progress', action='store_true', help='do not display progress bar')   
    other.add_argument('--min-mapq', metavar='MIN_MAPQ', type=int, default=30,
                        help="minimum mapping quality. (default: %(default)s)")
    other.add_argument('--min-qual', metavar='MIN_QUAL', type=int, default=10,
                        help="minimum phred quality for CpG. (default: %(default)s)")
    other.add_argument('--chunk-size', metavar='CHUNK_SIZE', type=int, default=1000000,
                        help="chunk size. (default: %(default)s)")
    other.add_argument('--sort-chunk', action='store_true', help='sort each chunk')
    
    

    args, unknown = parser.parse_known_args(argv)

    return args

########################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

########################################################################
if __name__ == '__main__':
    sys.exit(main(sys.argv))
