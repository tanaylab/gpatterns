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
from reads import Read, bam_iter

#TODO: REMEMBER TO REMOVE 772/512 FLAGS BEFORE PRODUCTION

########################################################################
PROG = 'tidy_cpgs.py'


########################################################################
def bam_reader(bam, 
                chunk_size,
                stats,
                conv_stats,
                umi1_idx=None,
                umi2_idx=None,
                infer_umi_field=False,
                chrom=None,
                genomic_range=None,
                single_end=False,
                show_progress=True,
                min_mapq=30,
                max_no_conv=3,
                template_conv_tag='XT',
                meth_calls_tag='XP',
                add_chr_prefix=False,
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
        hs = []
        Hs = []
        xs = []
        Xs = []

        prev_id = None

        for read in tqdm(
                bam_iter(bam, umi1_idx=umi1_idx, umi2_idx=umi2_idx, infer_umi_field=infer_umi_field,
                         single_end=single_end, min_mapq=min_mapq, max_no_conv=max_no_conv, template_conv_tag=template_conv_tag, 
                         meth_calls_tag=meth_calls_tag, bismark_flags=bismark_flags),
                unit='reads',
                unit_scale=True, disable=not show_progress):

            if read.id != prev_id:
                prev_id = read.id
                stats[read.pairing] += 1                            
            
            if read.pairing in ['good', 'single_R1', 'single_R2'] and \
                    ((chrom is None) or \
                     (not add_chr_prefix and chrom == read.chrom) or \
                     (add_chr_prefix and chrom == 'chr' + read.chrom)):
                    
                if (genomic_range is None) \
                        or (genomic_range[0] <= read.start <= genomic_range[1]) \
                        or (genomic_range[0] <= read.end <= genomic_range[1]):

                    conv_stats['h'] += read.h
                    conv_stats['H'] += read.H
                    conv_stats['x'] += read.x
                    conv_stats['x'] += read.X
                    conv_stats['z'] += read.z
                    conv_stats['Z'] += read.Z
                                
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

                    hs.append(read.h)
                    Hs.append(read.H)
                    xs.append(read.x)
                    Xs.append(read.X)

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

                        hs.append(read.h)
                        Hs.append(read.H)
                        xs.append(read.x)
                        Xs.append(read.X)

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

        hs = np.array(hs, dtype=np.int)
        Hs = np.array(Hs, dtype=np.int)
        xs = np.array(xs, dtype=np.int)
        Xs = np.array(Xs, dtype=np.int)

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
            'h': hs,
            'H': Hs,
            'x': xs,
            'X': Xs,
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
    conv_stats = defaultdict(int)

    all_columns = ['read_id', 'chrom', 'start', 'end', 'strand', 'umi1', 'umi2', 'insert_len', 'cg_pos', 'meth', 'qual', 'h', 'H', 'x', 'X']
#    out_columns = all_columns[0:11]
    out_columns = all_columns
    reads_columns = ['read_id', 'chrom', 'start', 'end', 'strand', 'umi1', 'umi2', 'insert_len', 'h', 'H', 'x', 'X']
    
    with pysam.AlignmentFile(args.input, "rb") as in_bam, \
            open(args.output, 'w') if args.output is not '-' else sys.stdout as out_file, \
            open(args.stats, 'w') if args.stats is not '-' else sys.stderr as stats_out:


        reads_out = open(args.reads, 'w') if args.reads is not None else None  

        out_file.write(','.join(out_columns) + '\n')
        if reads_out is not None:
            reads_out.write(','.join(reads_columns) + '\n')   
        
        for ids_df, patts, quals in bam_reader(in_bam,
                                               chunk_size=args.chunk_size,
                                               stats=stats,
                                               conv_stats=conv_stats,
                                               umi1_idx=args.umi1_idx,
                                               umi2_idx=args.umi2_idx,
                                               infer_umi_field=args.infer_umi_field,
                                               chrom=args.chrom,
                                               genomic_range=args.genomic_range,
                                               single_end=args.single_end,
                                               show_progress=not args.no_progress,
                                               min_mapq=args.min_mapq,
                                               max_no_conv=args.max_no_conv,
                                               template_conv_tag=args.template_conv_tag,
                                               meth_calls_tag=args.meth_calls_tag,
                                               add_chr_prefix=args.add_chr_prefix,
                                               bismark_flags=args.bismark_flags):
            cpgs = extract_cgs(ids_df, patts, quals, columns=all_columns,
                               min_qual=args.min_qual,
                               genomic_range=args.genomic_range,
                               cgs_mask=cgs_mask,
                               trim=args.trim,
                               add_chr_prefix=args.add_chr_prefix)

            if args.sort_chunk:
                cpgs.sort_values(by=['chrom', 'start', 'end', 'strand', 'umi1', 'umi2'], inplace=True)
            cpgs.to_csv(out_file, header=False, index=False, float_format='%.0f', mode='a')

            if reads_out is not None:
                ids_df[reads_columns].drop_duplicates(['read_id']).to_csv(reads_out, header=False, index=False, float_format='%.0f', mode='a')
        
        np.seterr(divide='ignore')
        
        stats['CHH'] = np.float64(conv_stats['H']) / np.float64(conv_stats['H'] + conv_stats['h'])
        stats['CHG'] = np.float64(conv_stats['X']) / np.float64(conv_stats['X'] + conv_stats['x'])
        stats['CpG'] = np.float64(conv_stats['Z']) / np.float64(conv_stats['Z'] + conv_stats['z'])

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
    parser.add_argument('-r', '--reads', type=str, default=None, metavar='<reads>',
                        help="Name of reads file. Outputs statistics for every read (read_id,chrom,start,end,strand,umi1,umi2,insert_len)")
    
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
    other.add_argument('--max-no-conv', metavar='MAX_NO_CONV', type=int, default=3,
                        help="maximal number of non converted C (or G in GA mode) not in CpG context. (default: %(default)s)")
    
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
