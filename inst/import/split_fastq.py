#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess

########################################################################
PROG = 'split_fastq.py'
VERSION = '0.1'

CHUNK_SIZE = (2**30)

file_counter = 0
children = set()


########################################################################
def main(argv):
    args = parse_args(argv)

    saved = []
    saved_len = 0
    while True:
        # Read lines from input
        lines = sys.stdin.buffer.readlines(CHUNK_SIZE)        
        
        # Release any zombies
        while (len(children) > 0):
            pid, _ = os.waitpid(-1, os.WNOHANG)
            if (pid == 0):
                break
            children.remove(pid)

        # Check for EOF
        if (len(lines) == 0):        
            write(args.prefix, args.suffix, args.counter_width, args.reads, saved, [], saved_len, flush=True)
            break

        # Write what we have
        saved, saved_len = write(args.prefix, args.suffix, args.counter_width, args.reads, saved, lines, saved_len)
    
    # Wait for any remaining children
    while (len(children) > 0):
        pid, _ = os.waitpid(-1, 0)
        children.remove(pid)


########################################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(prog=PROG,
                                     description='Split lines from stdin into multiple gzipped files.',
                                     epilog='Output filenames follow the pattern <prefix>.<counter>.<suffix>.gz',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--prefix', required=True, type=str, metavar='<prefix>',
                        help="Prefix of output file names")
    parser.add_argument('-s', '--suffix', default='fastq', type=str, metavar='<suffix>',
                        help="Suffix of output file names")
    parser.add_argument('-c', '--counter-width', type=int, default=2, metavar='<num>',
                        help="Number of digits in the couter part of the output file names.")
    parser.add_argument('-r', '--reads', type=int, default=4000000, metavar='<num>',
                        help="Number of reads to put in each output file.")

    args = parser.parse_args(argv[1:])

    return args
    

########################################################################
def write(prefix, suffix, counter_width, split_len, saved, lines, saved_len, flush=False):
    global file_counter
    global children
    
    split_len *= 4 
    
    if (flush):
        file_counter += 1
        file_name =  '{prefix}.{counter:0{width}}.{suffix}.gz'.format(prefix=prefix, suffix=suffix, counter=file_counter, width=counter_width)
        gzip = subprocess.Popen('gzip -c > %s'%file_name, shell=True, stdin=subprocess.PIPE)
        for block in saved:
            gzip.stdin.writelines(block)
        gzip.communicate()
        return [], 0

    total = len(lines) + saved_len 
    
    if (total < split_len):
        saved.append(lines)
        saved_len += len(lines)
        return saved, saved_len

    from_line = 0
    to_line = split_len - saved_len
    while(to_line < len(lines)):
        file_counter += 1
        pid = os.fork()
        if (pid == 0):
            file_name =  '{prefix}.{counter:0{width}}.{suffix}.gz'.format(prefix=prefix, suffix=suffix, counter=file_counter, width=counter_width)
            gzip = subprocess.Popen('gzip -c > %s'%file_name, shell=True, stdin=subprocess.PIPE)
            for block in saved:
                gzip.stdin.writelines(block)
            gzip.stdin.writelines(lines[from_line:to_line])
            gzip.communicate()
            sys.exit()

        children.add(pid)
        saved = []
        from_line = to_line
        to_line += split_len
    
    leftover = len(lines) - from_line    

    return [lines[-leftover:]], leftover


########################################################################
if (__name__ == '__main__'):
    sys.exit(main(sys.argv))