#!/usr/bin/env python2
import argparse
import pysam
import sys
import collections
from tqdm import tqdm
import numpy

########################################################################
PROG = 'conv_stats.py'

########################################################################
def conversion_stats(bam, min_mapq):
    char_counts = collections.Counter()    
    for read in tqdm(bam):           
        if (read.mapping_quality >= min_mapq):
            char_counts += collections.Counter(read.get_tag('XP'))            
            
    stats = dict()
    stats['CHH'] = numpy.float64(char_counts['H']) / numpy.float64(char_counts['H'] + char_counts['h'])
    stats['CHG'] = numpy.float64(char_counts['X']) / numpy.float64(char_counts['X'] + char_counts['x'])
    stats['CpG'] = numpy.float64(char_counts['Z']) / numpy.float64(char_counts['Z'] + char_counts['z'])
    return(stats)


########################################################################
def main(argv):
    try:
        args = parse_args(argv[1:])
    except ValueError as err:
        eprint("Error parsing arguments: %s" % err)
        return 2
    
    with pysam.AlignmentFile(args.input, "rb") as in_bam:        
        stats = conversion_stats(in_bam, args.min_qual)    
        sys.stdout.write('\t'.join(['CHH', 'CHG', 'CpG']) + '\n')
        sys.stdout.write('\t'.join([str(stats['CHH']), str(stats['CHG']), str(stats['CpG'])]) + '\n')
                
    return 0


########################################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(prog=PROG,
                                     description='Returns conversion (%CHH, %CHG, %CpG) stats from bissli aligned bam file',                                     
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-i', '--input', type=str, default='-', metavar='<input_bam>',
                        help="Name of input bam/sam file. a value of '-' indicates that the file would be read from STDIN") 
    parser.add_argument('--min-qual', metavar='NUM', type=int, default=30,
                        help="minimal mapping quality")
    
    args = parser.parse_args(argv)
  
    return args

        

        
########################################################################
if (__name__ == '__main__'):
    sys.exit(main(sys.argv))