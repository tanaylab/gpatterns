#!/usr/bin/env python3
from __future__ import print_function
import sys
import argparse
import collections
import re
import subprocess
import os
import numpy
import numpy.core.defchararray as npchar
import pandas
import pysam
from functools import reduce
from tqdm import tqdm

import time

########################################################################
PROG = 'map_fastq'
VERSION = '0.1'

CHUNK_SIZE = 1000000


########################################################################
def main(argv):
    try:
        args = parse_args(argv[1:])
    except ValueError as err:
        eprint("Error parsing arguments: %s" % err)
        return 2
    
    index, read_regions, umi_regions, idx_regions, target_names = read_index(args.index)
    
    writer = TargetWriter(read_regions, target_names, max_reads=args.reads_per_file)

    r1 = fastq_reader(args.read1, CHUNK_SIZE, not args.no_progress)
    r2 = None
    if (not args.read2 is None):
        r2 = fastq_reader(args.read2, CHUNK_SIZE, not args.no_progress)

    ids2, seqs2, quals2 = None, None, None
    for ids1, seqs1, quals1 in r1:
        if (not r2 is None):
            try:
                ids2, seqs2, quals2 = r2.__next__()
                if (len(ids1) != len(ids2)):
                    raise StopIteration
            except StopIteration:
                raise ValueError("FASTQs holding read1 and read2 are of unequal length")


        ids = [ids1, ids2]
        seqs = [seqs1, seqs2]
        quals = [quals1, quals2]

        mapping = extract_seqs(seqs, idx_regions, umi_regions, quals, args.min_qual, df_colname='barcode', keep_components=True)        
        mapping['idx'] = numpy.arange(len(mapping))
        mapping = mapping.merge(index, on='barcode', how='left').sort_values(by='idx')
        
        writer.write(ids, seqs, quals, mapping)

    writer.close()

    return 0


########################################################################
Region = collections.namedtuple('Region', ('name', 'read', 'start', 'end'))


########################################################################
class TargetWriter:
    def __init__(self, regions, target_names, max_forks=None, max_reads=None):
        self.regions = regions
        self.reads = [False, False]
        self.targets = []
        self.max_forks = max_forks if (not max_forks is None) else len(target_names)
        self.procs = []
        self.pids = []
        for region in self.regions:
            self.reads[region.read] = True
        for name in target_names:
            sub_targets = []
            for region in self.regions:
                if max_reads == None:
                    print('gzip -c > %s_%s.fastq.gz' % (name, region.name))
                    proc = subprocess.Popen('gzip -c > %s_%s.fastq.gz' % (name, region.name), shell=True, stdin=subprocess.PIPE)
                else:
                    proc = subprocess.Popen('%s/split_fastq.py -p %s_%s -s fastq -r %s ' % (os.path.dirname(os.path.realpath(__file__)), name, region.name, max_reads), shell=True, stdin=subprocess.PIPE)
                self.procs.append(proc)
                sub_targets.append(proc.stdin)                
            self.targets.append(sub_targets)
            

    def wait(self):
        while (self.pids):
            pid = 0
            while (pid != self.pids[0]):
                pid, _ = os.waitpid(self.pids[0], 0)
            del self.pids[0]

        for pid in self.pids:
            rc = 0
            while (rc != pid):
                pid, _ = os.waitpid(pid, 0)

    
    def close(self):
        self.wait()
        for proc in self.procs:
            proc.communicate()


    def write(self, ids, seqs, quals, mapping):
        self.wait()

        for i, qual in enumerate(quals):
            if (not qual is None):
                quals[i] = qual.view('U1')

        for target_idx, sub_targets in enumerate(self.targets):
            pid = os.fork()
            if (pid != 0):
                self.pids.append(pid)
                continue
            
            mask = (mapping['target_idx'] == target_idx).values
            mapping = mapping[mask].drop(columns=['idx', 'target_idx', 'barcode'])
            for i, read in enumerate(self.reads):
                if (read):                     
                    ids[i] = ids[i][mask]
                    seqs[i] = seqs[i][mask,:]
                    quals[i] = quals[i][mask,:]             
            for region, target in zip(self.regions, sub_targets):
                cut_ids = ids[region.read]
                cut_seqs = extract_seqs(seqs, [region])
                cut_quals = extract_seqs(quals, [region])                
                if (mapping.shape[1]):
                    for id, seq, qual, named in zip(cut_ids, cut_seqs, cut_quals, mapping.itertuples(index=False)):
                        id = ['@', id]
                        for name, value in zip(named._fields, named):
                            id += [':', name, '=', value]
                        l = id + ['\n', seq, '\n+\n', qual, '\n']
                        l = [x.encode() for x in l]                        
                        target.writelines(l)
                else:
                    for id, seq, qual in zip(cut_ids, cut_seqs, cut_quals):    
                        l = ['@', id , '\n', seq, '\n+\n', qual, '\n']      
                        l = [x.encode() for x in l]
                        target.writelines(l)
                target.flush()
            sys.exit(0)                


########################################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(prog=PROG,
                                     description='Split FASTQs based on internal barcodes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r1', '--read1', type=str, default='-', metavar='<read1.fastq>',
                        help="Name of input R1 file. A value of '-' indicates that read 1 is read from stdin.")
    
    parser.add_argument('-r2', '--read2', type=str, default=None, metavar='<read2.fastq>',
                        help="Name of input R2 file.")
    
    parser.add_argument('-i', '--index', type=str, default='-', metavar='<index.tab>',
                        help="File holding the index definitions. A value of '-' indicates that the index should be read from stdin.")
    
    parser.add_argument('-q', '--min-qual', type=int, default=0, metavar='<qual>',
                        help="Minimal base quality for index. Bases with a lower quality are masked.")
                        
    parser.add_argument('-r', '--reads-per-file', type=int, default=None, metavar='<reads_per_file>',
                        help="Number of reads per file.")

    parser.add_argument('--no-progress', action='store_true', help='do not display progress bar')

    args = parser.parse_args(argv)
    
    if ((args.read1 == '-') and (args.index == '-')):
        raise ValueError("It is impossible to read both read 1 and the index from stdin.")
    
    return args


########################################################################
def read_index(index_fn):
    if (index_fn == '-'):
        fh = sys.stdin
    else:
        fh = open(index_fn, 'r')

    index_tab = pandas.read_table(fh)
    index_cols = []
    other_cols = []
    read_regions = []
    umi_regions = []
    idx_regions = []
    
    umi_counter = 0
    
    for column in index_tab.columns:
        if (column in ('target', 'counter')):
            other_cols.append(column)
            continue
        
        match = re.match('(idx|umi|read)(?::([a-zA-Z_][a-zA-Z0-9_]*))?\|(?:r([12]):)?([0-9]+)-([0-9]+)', column)
        if (match is None):
            raise ValueError("Failed to parse index title: %s" % column)
        mode = match.group(1)
        name = match.group(2)
        read = match.group(3)
        start = int(match.group(4)) - 1
        end = int(match.group(5))
        if (read is None):
            read = 0
        else:
            read = int(read) - 1

        if (mode == 'read'):
            read_regions.append(Region(name, read, start, end))
            continue
        
        if (mode == 'umi'):
            if (name is None):
                name = 'umi' if (umi_counter == 0) else 'umi%d' % umi_counter
                umi_counter += 1
            umi_regions.append(Region(name, read, start, end))
            continue
        
        if (mode == 'idx'):
            idx_regions.append(Region(name, read, start, end))            
            index = numpy.array(index_tab[column]).astype('U')            
            if (len(index[0]) != (end-start)):
                raise ValueError('Number of bases in index %s does not match its region' % column)
            index_cols.append(index)

    if (not index_cols):
        raise ValueError("Index table does not contain any indices")

    targets = index_tab['target'].unique()
    target_tab = pandas.DataFrame({'target': targets,
                                  'target_idx': numpy.arange(len(targets), dtype=numpy.uint32)})
    index_tab = index_tab.merge(target_tab, how='left')

    index_tab['barcode'] = reduce(npchar.add, index_cols)
    
    index_tab = index_tab[['barcode', 'target_idx']]
    
    return index_tab, read_regions, umi_regions, idx_regions, targets


########################################################################
def fastq_reader(filename, chunk_size, show_progress=True):
    fastq = pysam.FastxFile(filename, persist=False)

    while True:
        ids = []
        seqs = []
        quals = []

        for read in tqdm(fastq, unit = 'reads', unit_scale=True, disable=not show_progress):            
            ids.append(read.name)
            seqs.append(read.sequence)
            quals.append(read.quality)
            if (len(ids) >= chunk_size):
                break
        
        if (not ids):
            break

        ids = numpy.array(ids, dtype=object)        
        seqs = numpy.array(seqs)                
        seqs = seqs.view('U1').reshape((seqs.size, -1))
        quals = numpy.array(quals)
        quals = quals.view(numpy.uint8).reshape((quals.size, -1))

        yield ids, seqs, quals


########################################################################
def extract_seqs(seqs, idx_regions, umi_regions=[], quals=None, min_qual=0, df_colname=None, keep_components=False):
    return_array = False
    if ((df_colname is None) and (not keep_components)):
        return_array = True

    width = 0
    index = []
    components = {}
    for region in idx_regions:
        name = region.name
        read = region.read
        start = region.start
        end = region.end
        width += end-start
        bases = seqs[read][:, start:end]
        if ((not quals is None) and (min_qual > 0)):
            mask  = quals[read][:, start:end] < min_qual + 33
            bases = bases.ravel()
            bases[mask.ravel()] = 'N'
            bases.reshape((len(seqs[read]), -1))
        else:
            bases = numpy.copy(bases)
        index.append(bases)
        if (keep_components and name):
            components[name] = bases.view('U%d'%(end-start)).ravel()
            
    if (keep_components):
        for region in umi_regions:
            name = region.name
            read = region.read
            start = region.start
            end = region.end
            components[name] = numpy.copy(seqs[read][:, start:end]).view('U%d'%(end-start)).ravel()

    if ((not df_colname is None) or return_array):
        if (len(index) == 1):
            merged = index[0]
        else:
            merged = numpy.hstack(index)
        merged = merged.view('U%d'%width).ravel()
        
    if (df_colname is None):
        return merged
        
    if (not df_colname is None):
        components[df_colname] = merged

    return pandas.DataFrame(components)


########################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


########################################################################
if (__name__ == '__main__'):
    sys.exit(main(sys.argv))
    
