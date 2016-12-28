#!/usr/bin/env python2
import argparse
import sys
#from tqdm import tqdm
import itertools
import csv
import gzip
import subprocess as sp
from os.path import splitext

########################################################################
PROG = 'filter_dups_cpgs.py'

########################################################################
def hamming(str1, str2):
    return(sum(itertools.imap(str.__ne__, str1, str2)))

########################################################################
def umi_distinct(umi1, umi2, min_hamming):
    return(hamming(umi1, umi2) >= min_hamming)

########################################################################
class Read:
    def __init__(self, read_id, chrom, start, end, strand, umi1, umi2, insert_len, cg_pos, meth, qual):
        self.read_id = read_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.umi1 = umi1
        self.umi2 = umi2
        self.insert_len = insert_len
        self.cg_pos = cg_pos
        self.meth = meth
        self.qual = qual
        
        if (umi1 is not None) and (umi2 is not None):
            self.umi = ''.join([umi1, umi2])
        elif umi1 is not None:
            self.umi = self.umi1
        elif umi2 is not None:
            self.umi = self.umi2
        
        self.key = "%s:%s-%s" % (self.chrom, self.start, self.strand)
        self.key2 = self.end    
        
    def write(self, out_file, n_umis):
        for cg_pos, meth, qual in zip(self.cg_pos, self.meth, self.qual):
            out_file.write(','.join([self.read_id, 
                                     self.chrom, 
                                     self.start, 
                                     self.end, 
                                     self.strand, 
                                     self.umi1, 
                                     self.umi2, 
                                     self.insert_len,
                                     str(n_umis),
                                     cg_pos, 
                                     meth,
                                     qual]) + '\n') 
            
    def is_duplicate(self, other_read, two_sides=True, use_seq=False, only_seq=False, min_hamming=4):
        if (not only_seq):
            if self.key != other_read.key:      
                return False
            if (two_sides and (self.key2 != other_read.key2) and (other_read.key2 != '-') and (self.key2 != '-')):
                return False            
        if (use_seq or only_seq) and umi_distinct(self.umi, other_read.umi, min_hamming):
            return False
        
        return True

########################################################################
def cgs_reader(csv_file):    
    first_line = True
    cg_pos = []
    meth = []
    qual = []
    last_line = None
    n = 0
    
    try:
        for line in csv.DictReader(csv_file):                
            n += 1
            if first_line:            
                last_read_id = line['read_id']
                first_line = False
            
            if last_read_id != line['read_id']:
                yield Read(last_line['read_id'], 
                           last_line['chrom'], 
                           last_line['start'], 
                           last_line['end'], 
                           last_line['strand'], 
                           last_line['umi1'], 
                           last_line['umi2'], 
                           last_line['insert_len'],
                           cg_pos, 
                           meth,
                           qual)         
                
                cg_pos = [line['cg_pos']]
                meth = [line['meth']]
                qual = [line['qual']]         
                last_read_id = line['read_id']
                last_line = line
                
            else:
                #if we have the same position+read twice            
                if len(meth) > 0 and len(cg_pos) > 0 and cg_pos[-1] == line['cg_pos']:                
                    #if meth calls differ discard the call
                    if meth[-1] != line['meth']:                    
                        cg_pos.pop()
                        meth.pop()
                    continue
                
                
                cg_pos.append(line['cg_pos'])
                meth.append(line['meth'])
                qual.append(line['qual'])
                last_line = line
    except csv.Error:
        # ignore malformed csv's
        pass
    
    if n > 0:
        yield Read(last_line['read_id'], 
                       last_line['chrom'], 
                       last_line['start'], 
                       last_line['end'], 
                       last_line['strand'], 
                       last_line['umi1'], 
                       last_line['umi2'], 
                       last_line['insert_len'],
                       cg_pos, 
                       meth,
                       qual)         

########################################################################
def get_sort_fields(only_seq):
    if only_seq:
        return('6,7')
    return('2,7')

########################################################################
def main(argv):
    try:
        args = parse_args(argv[1:])
    except ValueError as err:
        eprint("Error parsing arguments: %s" % err)
        return 2
    
    stats = {'total_reads':0,
             'uniq_reads':0}
    
    with open(args.input, "rb") if args.input is not '-' else sys.stdin as in_file, \
        open(args.output, 'w') if args.output is not '-' else sys.stdout as out_file:
        
        out_file.write(','.join(['read_id', 'chrom', 'start', 'end', 'strand', 'umi1', 'umi2',  'insert_len', 'num', 'cg_pos', 'meth', 'qual']) + '\n')        
        
        cmd_prefix = ''
        if args.input is not '-':
            ext = splitext(args.input)[1]            
            if ext == '.gz':
                cmd_prefix = 'gzip -d -c | '
        
        if args.sorted:
            in_file_sorted = in_file            
        else:
            p = sp.Popen(cmd_prefix + 'awk \'NR==1; NR > 1 {print $0 | "sort --field-separator=, -k' + str(get_sort_fields(args.only_seq)) + ' -k1 -k9"}\'', stdin = in_file, stdout=sp.PIPE, shell=True, universal_newlines=True)            
            in_file_sorted = p.communicate()[0].decode(args.encoding).splitlines()        
        
        first_read = True
        n_umis = 0
        reads = []
        
        for read in cgs_reader(in_file_sorted):  
            stats['total_reads'] += 1            
            if (first_read):             
                reads.append(read)
                last_read = read
                first_read = False                

            if (not read.is_duplicate(last_read, two_sides=not args.only_R1, use_seq=args.use_seq, only_seq=args.only_seq, min_hamming=args.min_hamming)):
                for r in reads:
                    r.write(out_file, n_umis)
                stats['uniq_reads'] += 1                
                reads = [read]                
                last_read = read
                n_umis = 1
            else:
                n_umis += 1
        
        for r in reads:
            r.write(out_file, n_umis)
    
    with open(args.stats, 'w') if args.stats is not '-' else sys.stderr as out_stats:
        out_stats.write('\t'.join(['total_reads', 'uniq_reads', 'uniq_frac']) + '\n')
        
        try:
            uniq_frac = str(float(stats['uniq_reads'])/float(stats['total_reads']))
        except ZeroDivisionError: 
            uniq_frac = 'NA'
        out_stats.write('\t'.join(
                    [str(stats['total_reads']), 
                     str(stats['uniq_reads']), 
                     uniq_frac
                    ]) + '\n')
    
    print >> sys.stderr, 'Unique mols: %d/%d' % (stats['uniq_reads'], stats['total_reads'])
    return 0

########################################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(prog=PROG,
                                     description='Filters tidy_cpgs files from duplicates.',                                     
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-i', '--input', type=str, default='-', metavar='<input>',
                        help="Name of input tidy cpgs file (can be gzipped). a value of '-' indicates that the input would be read from STDIN (and then it should no be gzipped).")   
    parser.add_argument('-o', '--output', type=str, default='-', metavar='<output>',
                        help="Name of output file. a value of '-' indicates that the output would be written to STDOUT")    
    parser.add_argument('-s', '--stats', type=str, default='stats', metavar='<stats>',
                        help="Name of stats file. a value of '-' indicates that the output would be written to STDERR")
    parser.add_argument("--only_R1", help="use only read1",
                    action="store_true")
    parser.add_argument("--use-seq", help="use UMI sequence (not only position)",
                    action="store_true")
    parser.add_argument("--only-seq", help="use only UMI sequence (without positions)",
                    action="store_true")
    parser.add_argument("--sorted", help="file is already sorted",
                    action="store_true")
    parser.add_argument('--min-hamming', metavar='NUM', type=int, default=4,
                        help="Minimal hamming distance to consider UMIs different (applies only when using UMI sequence)")
    parser.add_argument('--encoding', type=str, default='ascii', metavar='<output>',
                        help="encoding of the csv file. default: 'ascii'")                            
    args, unknown = parser.parse_known_args(argv)
  
    return args

        

        
########################################################################
if (__name__ == '__main__'):
    sys.exit(main(sys.argv))
    
