#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Temp;
use IO::Uncompress::Gunzip;
use Bio::Cigar;

my $PROGRAM = 'bissli2';
my $VERSION = '0.2.4';

my @HARDCODED_BOWTIE2_OPTIONS = ('-q', '--quiet', '--reorder');

########################################################################

my $progname = $0;
$progname =~ s|^.*/||;
my $cmdline =  join(' ', $progname, @ARGV);

my ($genome_dir, $index, $m1, $m2, $output_fname, $bowtie2, $tmp_dir, $keep_tmp, $ct, $ga, @bowtie2_options) = parse_args();

my $tmp_dirobj = File::Temp->newdir("$tmp_dir/bissli_XXXXXX", CLEANUP => (1-$keep_tmp));
$tmp_dir       = $tmp_dirobj->dirname;
if ($keep_tmp) {
    print STDERR "Temporary files are written to: $tmp_dir\n";
}


my ($m1_orig_fname, $m1_ct_fname, $m1_ga_fname, $m2_orig_fname, $m2_ct_fname, $m2_ga_fname) = convert_input($m1, $m2, $tmp_dir, $ct, $ga);
my @fnames =  ($m1_orig_fname, $m1_ct_fname, $m1_ga_fname, $m2_orig_fname, $m2_ct_fname, $m2_ga_fname);
my ($sam_ct_fname, $sam_ga_fname) = run_bowtie($bowtie2, $index, $m1_ct_fname, $m1_ga_fname, $m2_ct_fname, $m2_ga_fname, $tmp_dir, $ct, $ga, @bowtie2_options);

my $genome = construct_genome($genome_dir);
merge_sams($output_fname, $sam_ct_fname, $sam_ga_fname, $m1_orig_fname, $m2_orig_fname, $genome, $cmdline);


########################################################################
sub parse_args
{
    my $help         = 0;
    my $m1           = undef;
    my $m2           = undef;
    my $bowtie2      = 'bowtie2';
    my $c            = 0;
    my $ct           = 0;
    my $f            = 0;
    my $genome_dir   = undef;
    my $ga           = 0;
    my $keep_tmp     = 0;
    my $qseq         = 0;
    my $r            = 0;
    my $output_fname = '-';
    my $tmp_dir      = '/tmp';
    my $U            = undef;
    my $index        = undef;
    my @bowtie2_options;

    GetOptions( '1=s'       => \$m1,
                '2=s'       => \$m2,
                'bowtie2=s' => \$bowtie2,
                'c'         => \$c,
                'ct'        => \$ct,
                'f'         => \$f,
                'g=s'       => \$genome_dir,
                'ga'        => \$ga,
                'keep-tmp'  => \$keep_tmp,
                'qseq'      => \$qseq,
                'r'         => \$r,
                'S=s'       => \$output_fname,
                'tmp-dir=s' => \$tmp_dir,
                'U=s'       => \$U,
                'x=s'       => \$index,
                'h|help'    => \$help)  or pod2usage(-exitval => 2, -verbose => 0);

    pod2usage(-exitval => 0, -verbose => 1) if ($help);
    die("The flag -U cannot be used along either -1 or -2") if ((defined $U) and ((defined $m1) or (defined $m2)));
    die("The flags -1 and -2 must always be used together") if ((defined $m1) ^ (defined $m2));
    die("No input files") unless ((defined $U) or (defined $m1));
    die("Genome (-g) must be specified.") unless (defined $genome_dir);
    dir("Index (-x) must be specified.") unless (defined $index);
    die("Only FASTQ input file format is supported") if ($c or $f or $qseq or $r);
#    die("--ct and --ga flags are mutually exclusive") if ($ct and $ga);
    die("Currently, either --ct or --ga must be specified") unless ($ct or $ga);
    
    if (defined $U) {
        $m1 = $U;
    }
    unless ($ct or $ga) {
        $ct = 1;
        $ga = 1;
    }
    
    @bowtie2_options = @ARGV;
    return ($genome_dir, $index, $m1, $m2, $output_fname, $bowtie2, $tmp_dir, $keep_tmp, $ct, $ga, @bowtie2_options);
}


########################################################################
sub convert_input
{
    my ($m1, $m2, $tmp_dir, $ct, $ga) = @_;
    
    my ($m1_orig_fname, $m1_ct_fname, $m1_ga_fname);
    my ($m2_orig_fname, $m2_ct_fname, $m2_ga_fname);
    my ($m1_orig, $m1_ct, $m1_ga);
    my ($m2_orig, $m2_ct, $m2_ga);
    
    $m1_orig_fname = "$tmp_dir/m1_orig.fastq";
    open($m1_orig, '>', $m1_orig_fname) || die("Cannot create temporary FASTQ file '$m1_orig_fname'");
    if ($ct) {
        $m1_ct_fname = "$tmp_dir/m1_ct.fastq";
        open($m1_ct,   '>', $m1_ct_fname) || die("Cannot create temporary FASTQ file '$m1_ct_fname'");
    }
    if ($ga) {
        $m1_ga_fname = "$tmp_dir/m1_ga.fastq";
        open($m1_ga,   '>', $m1_ga_fname) || die("Cannot create temporary FASTQ file '$m1_ga_fname'");
    }

    if (defined $m2) {
        # Note that read 2 CT and GA version will be created based on whether we perform a GA or CT alignment respectively 
        $m2_orig_fname = "$tmp_dir/m2_orig.fastq";
        open($m2_orig, '>', $m2_orig_fname) || die("Cannot create temporary FASTQ file '$m2_orig_fname'");
        if ($ga) {
            $m2_ct_fname = "$tmp_dir/m2_ct.fastq";
            open($m2_ct,   '>', $m2_ct_fname) || die("Cannot create temporary FASTQ file '$m2_ct_fname'");
        }
        if ($ct) {
            $m2_ga_fname = "$tmp_dir/m2_ga.fastq";
            open($m2_ga,   '>', $m2_ga_fname) || die("Cannot create temporary FASTQ file '$m2_ga_fname'");
        }        
    }
    
    my (@m1_fnames, @m2_fnames); 
    @m1_fnames = split(',', $m1);
    if (defined $m2) {
        @m2_fnames = split(',', $m2);
        die("Input file lists <m1> and <m2> are of unequal length.") if (scalar @m1_fnames != scalar @m2_fnames);
    }
    for my $fname (@m1_fnames) {
        convert_fastq($m1_orig, $m1_ct, $m1_ga, $fname);
    }
    if (defined $m2) {
        for my $fname (@m2_fnames) {
            convert_fastq($m2_orig, $m2_ct, $m2_ga, $fname);
        }
    }

    close($m1_orig) if (defined $m1_orig);
    close($m1_ct)   if (defined $m1_ct);
    close($m1_ga)   if (defined $m1_ga);
    close($m2_orig) if (defined $m2_orig);
    close($m2_ct)   if (defined $m2_ct);
    close($m2_ga)   if (defined $m2_ga);
    
    return ($m1_orig_fname, $m1_ct_fname, $m1_ga_fname, $m2_orig_fname, $m2_ct_fname, $m2_ga_fname);
}

########################################################################
sub convert_fastq
{
    my ($orig_fastq, $ct_fastq, $ga_fastq, $input_fname) = @_;

    my $input = new IO::Uncompress::Gunzip($input_fname);
    
    while (1) {
        my ($header, $seq, $qual) = read_fastq_read($input, $input_fname);
        last unless (defined $header);
        
        write_fastq_read($orig_fastq, $header, $seq, $qual);
        
        if (defined $ct_fastq) {
            my $ct_seq = $seq;
            $ct_seq =~ tr/cC/tT/;
            write_fastq_read($ct_fastq, $header, $ct_seq, $qual);
        }

        if ($ga_fastq) {
            my $ga_seq = $seq;
            $ga_seq =~ tr/gG/aA/;
            write_fastq_read($ga_fastq, $header, $ga_seq, $qual);
        }
    }
}


########################################################################
sub run_bowtie
{
    my ($bowtie2, $index, $m1_ct_fname, $m1_ga_fname, $m2_ct_fname, $m2_ga_fname, $tmp_dir, $ct, $ga, @bowtie2_options) = @_;
    my $sam_ct_fname;
    my $sam_ga_fname;

    if ($ct) {
        $sam_ct_fname = "$tmp_dir/align_ct.sam";
        my @cmd;
        if (defined $m2_ga_fname) {
            @cmd = ($bowtie2, '-1', $m1_ct_fname, '-2', $m2_ga_fname);
        }
        else {
            @cmd = ($bowtie2, '-U', $m1_ct_fname);        
        }
        @cmd = (@cmd, '-x', $index, '-S', $sam_ct_fname, '--norc', @HARDCODED_BOWTIE2_OPTIONS, @bowtie2_options);
        print STDERR "Running: @cmd\n";
        my $rc = system(@cmd);
        die("Failed running CT bowtie") unless ($rc == 0);
    }
    if ($ga) {
        $sam_ga_fname = "$tmp_dir/align_ga.sam";
        my @cmd;
        if (defined $m2_ct_fname) {
            @cmd = ($bowtie2, '-1', $m1_ga_fname, '-2', $m2_ct_fname);
        }
        else {
            @cmd = ($bowtie2, '-U', $m1_ga_fname);        
        }
        @cmd = (@cmd, '-x', $index, '-S', $sam_ga_fname, '--nofw', @HARDCODED_BOWTIE2_OPTIONS, @bowtie2_options);
        print STDERR "Running: @cmd\n";
        my $rc = system(@cmd);
        die("Failed running GA bowtie") unless ($rc == 0);    
    }
    
    return ($sam_ct_fname, $sam_ga_fname);
}


########################################################################
sub merge_sams
{
    my ($output_fname, $sam_ct_fname, $sam_ga_fname, $m1_orig_fname, $m2_orig_fname, $genome, $cmdline) = @_;
    
    my ($sam_ct,  $sam_ga);
    my ($m1_orig, $m2_orig);
    if (defined $sam_ct_fname) {
        open($sam_ct, '<', $sam_ct_fname) || die("Cannot open temporary SAM '$sam_ct_fname'");
    }
    if (defined $sam_ga_fname) {
        open($sam_ga, '<', $sam_ga_fname) || die("Cannot open temporary SAM '$sam_ga_fname'");
    }
    if (defined $m1_orig_fname) {
        open($m1_orig, '<', $m1_orig_fname) || die ("Cannot open temporary FASTQ '$m1_orig_fname'");
    }
    if (defined $m2_orig_fname) {
        open($m2_orig, '<', $m2_orig_fname) || die ("Cannot open temporary FASTQ '$m2_orig_fname'");
    }

    my $output;
    if ($output_fname eq '-') {
        $output = *STDOUT;
    }
    else {
        open($output, '>', $output_fname) || die("Cannot open output file '$output_fname'");
    }
    
    my $headers_ct;
    my $headers_ga;
    $headers_ct = read_headers($sam_ct) if (defined $sam_ct);
    $headers_ga = read_headers($sam_ga) if (defined $sam_ga);
    die("Bowtie2 generated different SAM headers between CT and GA runs.") if (diff_headers($headers_ct, $headers_ga));
    
    my ($headers, $chrom_lens) = process_headers($headers_ct, $headers_ga, $cmdline);
    for my $header (@{$headers}) {
        print $output $header, "\n";
    }
    
    while (1) {
        my ($align1, $align2) = read_alignment($m1_orig, $m2_orig, $sam_ct, $sam_ga, $chrom_lens);
        last unless (defined $align1);
        
        calc_meth_pattern($align1, $genome);
        calc_meth_pattern($align2, $genome) if defined ($align2);
        
        write_sam_alignment($output, $align1);
        write_sam_alignment($output, $align2) if defined ($align2);
    }

    close ($output) unless ($output_fname eq '-');
    close ($m2_orig) if defined($m2_orig);
    close ($m1_orig) if defined($m1_orig);
    close ($sam_ga)  if defined($sam_ga);
    close ($sam_ct)  if defined($sam_ct);
}


########################################################################
sub read_headers
{
    my ($file) = @_;
    
    my @headers;
    while (my $line = <$file>) {
        if (substr($line, 0, 1) ne '@') {
            seek ($file, -length($line), SEEK_CUR);
            last;
        }
        chomp($line);
        push(@headers, $line) unless (substr($line, 1, 2) eq 'PG');
    }

    return \@headers;
}


########################################################################
sub diff_headers
{
    my ($headers_ct, $headers_ga) = @_;
    
    return 0 unless((defined $headers_ct) and (defined $headers_ga));
    
    if (scalar @{$headers_ct} != scalar @{$headers_ga}) {
        return 1;
    }
    
    for (my $i=0; $i<@{$headers_ct}; ++$i) {
        if ($headers_ct->[$i] ne $headers_ga->[$i]) {
            return 1 if (substr($headers_ct->[$i], 1, 2) ne 'PG');
            return 1 if (substr($headers_ga->[$i], 1, 2) ne 'PG');
        }
    }
    
    return 0;
}


########################################################################
sub process_headers
{
    my ($headers_ct, $headers_ga, $cmdline) = @_;

    my $old_headers;
    $old_headers = $headers_ct;
    $old_headers = $headers_ga unless (defined $headers_ct);

    my %chrom_lens;
    my @new_headers;
    HEADER: for my $header (@{$old_headers}) {
        my $type = substr($header, 1, 2);
        
        if ($type eq "PG") {
            # TBD: Create a proper PG chain using the old PG headers!
            next HEADER;
        }
        
        if ($type eq "SQ") {
            my @fields = split("\t", $header);
            my $chr;
            for my $field (@fields) {
                if ($field =~ m/^SN:([a-zA-Z0-9_]+)_(fwd|rev)$/) {
                    $chr = $1;
                    $field = "SN:$1"; 
                }
                elsif ($field =~ m/^LN:([0-9]+)/) {
                    # We assume that the LN tag appears after the SN tag.
                    # It is unclear if this is a requirement by the SAM format.
                    my $len = $chrom_lens{$chr};
                    if (defined $len) {
                        die ("Conflicting length for chromosome '$chr'") if ($len != $1);
                        next HEADER;
                    }
                    $chrom_lens{$chr} = $1;
                }
            }
            push (@new_headers, join("\t", @fields));
            next;
        }
        
        push (@new_headers, $header);
    }

    push (@new_headers, join("\t", '@PG', "ID:$PROGRAM", "PN:$PROGRAM", "VN:$VERSION", "CL:\"$cmdline\""));

    return (\@new_headers, \%chrom_lens);
}


########################################################################
sub read_alignment
{
    my ($m1_orig, $m2_orig, $sam_ct, $sam_ga, $chrom_lens) = @_;
    
    my ($m1_header, $m2_header);
    my ($m1_seq,    $m2_seq);
    ($m1_header, $m1_seq) = read_fastq_read($m1_orig, "Temporary read 1 orig");
    return (undef, undef) unless (defined $m1_header);
    my $id = substr((split(/\s/, $m1_header))[0], 1);
    $m1_seq  = '*' if ($m1_seq eq '');
    
    if (defined $m2_orig) {
        ($m2_header, $m2_seq) = read_fastq_read($m2_orig, "Temporary read 2 orig");
        die("Temporary read 2 orig FASTQ is shorter than read1") unless (defined $m2_header);
        my $m2_id = substr((split(/\s/, $m2_header))[0], 1);
        die ("Mismatched read ID between read1 and read2 (FASTQ): '$id'!='$m2_id'") if ($m2_id ne $id);
        $m2_seq = '*' if ($m2_seq eq '');
    }

    my ($ct1, $ct2);
    my ($ga1, $ga2);
    if (defined $sam_ct) {
        ($ct1, $ct2) = read_sam_alignment($sam_ct, "Temporary CT SAM", (defined $m2_orig), 'CT');
        die ("Missing reads in temporary CT SAM") unless (defined $ct1);
    }
    if (defined $sam_ga) {
        ($ga1, $ga2) = read_sam_alignment($sam_ga, "Temporary GA SAM", (defined $m2_orig), 'GA');
        die ("Missing reads in temporary GA SAM") unless (defined $ga1);
    }

    my ($align1, $align2, $ct_choice) = choose_alignemt($ct1, $ct2, $ga1, $ga2);
    
    if (defined $align2) {
        # Correct swapped bowtie2 alignment
        my ($fastq_seq, $sam_seq, $paired_seq); 
        if ($ct_choice) {
            $fastq_seq  = $m1_seq;
            $sam_seq    = $align1->[0]->[9];
            $paired_seq = $align2->[0]->[9];
        }
        else {
            $fastq_seq  = $m2_seq;
            $sam_seq    = $align2->[0]->[9];            
            $paired_seq = $align1->[0]->[9];
        }
        $fastq_seq  =~ tr/cC/tT/;            
        if ($fastq_seq ne $sam_seq) {
            die unless ($fastq_seq = $paired_seq);
            ($align1, $align2) = ($align2, $align1);
        }

        $align1->[0]->[9] = $m1_seq;
        $align2->[0]->[9] = $m2_seq;
        fix_alignment($align1, length($m2_seq), $chrom_lens);
        fix_alignment($align2, length($m1_seq), $chrom_lens);        
    }
    else {
        $align1->[0]->[9] = $m1_seq;
        fix_alignment($align1, 0, $chrom_lens);
    }
    
    return ($align1, $align2);
}


########################################################################
sub choose_alignemt
{
    my ($ct1, $ct2, $ga1, $ga2) = @_;
    
    if ((defined $ct1) and (defined $ga1)){
        my $ct_score = $ct1->[0]->[4] + $ct2->[0]->[4];
        my $ga_score = $ga1->[0]->[4] + $ga2->[0]->[4];        
        if ($ct_score > $ga_score){
            return($ct1, $ct2, 1);
        } else {
            return($ga1, $ga2, 0);
        }
    }
    
    #die("Logic for merging SAM file is yet to be implemented...") if ((defined $ct1) and (defined $ga1));
    return ($ct1, $ct2, 1) if (defined $ct1);
    
    return ($ga1, $ga2, 0);
}


########################################################################
sub fix_alignment
{
    my ($align, $paired_length, $chrom_lens) = @_;

    my ($chrom, $paired_chrom);
    my ($reversed, $paired_reversed) = (0, 0);
    if ($align->[0]->[2] =~ m/^([a-zA-Z0-9_]+)_(fwd|rev)$/) {
        $chrom = $1;
        $align->[0]->[2] = $1;
        $reversed = 1 if ($2 eq 'rev');
    }
    elsif ($align->[0]->[2] eq '*') {
        $align->[2] = undef;
    }
    if ($align->[0]->[6] eq '=') {
        $paired_reversed = $reversed;
        $paired_chrom = $chrom;
    }
    elsif ($align->[0]->[6] =~ m/^([a-zA-Z0-9_]+)_(fwd|rev)$/) {
        $paired_chrom = $1;
        $align->[0]->[6] = $1;
        $paired_reversed = ($2 eq 'rev');
    }

    if ($reversed) {
        $align->[0]->[1] ^= 0x10;
        if ($align->[0]->[3] != 0) {
            my $ref_length =  length($align->[0]->[9]);
            if ($align->[0]->[5] ne '*'){
                my $cigar = Bio::Cigar->new($align->[0]->[5]);
                $ref_length = $cigar->reference_length;
            }
            $align->[0]->[3] = $chrom_lens->{$chrom} - $align->[0]->[3] - $ref_length + 2;
        }
    }
    if ($align->[0]->[1] & 0x10) {
        $align->[2] = reverse_conversion($align->[2]);        
    }

    if ($paired_reversed) {
        $align->[0]->[1] ^= 0x20;
        if ($align->[0]->[7] != 0) {
            $align->[0]->[7] = $chrom_lens->{$paired_chrom} - $align->[0]->[7] - $paired_length + 2;
        }        
    }
}


########################################################################
sub calc_meth_pattern
{
    my ($align, $genome) = @_;
    
    if ($align->[0]->[1] & 0x4) {
        $align->[2] = undef;
        return;
    }

    my $flags = $align->[0]->[1];
    my $chrom = $align->[0]->[2];
    my $start = $align->[0]->[3];
    my $seq   = $align->[0]->[9];
    $seq =~ s/^\s+|\s+$//g;
    my $revcomp = 0;
    my $cigar = Bio::Cigar->new($align->[0]->[5]);
        
    if ($align->[1] eq 'GA') {
        $seq = reverse_complement($seq);
        $revcomp = 1;
    }
    
    my $crick = $revcomp;
    if ($flags & 0x10) {
        $crick = 1 - $crick;
    }

    my $ref;
    if ($crick) {
        $ref     = get_reference($genome, $chrom, $start-2, $cigar->reference_length+2);
        $ref     = reverse_complement($ref);
    }
    else {
        $ref = get_reference($genome, $chrom, $start, $cigar->reference_length+2);
    }
    
    $seq = substr($seq, ($cigar->rpos_to_qpos(1) - 1), $cigar->reference_length);
    
    $align->[3] = calc_ct_meth_pattern($seq, $ref, $revcomp);

}

########################################################################
sub calc_ct_meth_pattern
{
    my ($seq, $ref, $revcomp) = @_;
    
    my $next_next = chop($ref);
    my $next      = chop($ref);
    my @pattern;
    
    while ($seq) {
        my $before = chop($ref);
        my $after  = chop($seq);
        my $call;
        if ($before ne 'C') {
            push(@pattern, '.');
        }
        else {
            if ($next eq 'G') {
                $call = 'z';
            }
            elsif ($next eq 'N') {
                $call = 'u';
            }
            elsif ($next_next eq 'G') {
                $call = 'x';
            }
            elsif ($next_next eq 'N') {
                $call = 'u'; 
            }
            else {
                $call = 'h';
            }
            
            if ($after eq 'T') {
                push(@pattern, $call);
            }
            elsif ($after eq 'C') {
                push(@pattern, uc($call));
            }
            else {
                push(@pattern, '.');                
            }
        }
        $next_next = $next;
        $next      = $before;
    }

    if ($revcomp) {
        return join('', @pattern);        
    }    

    return reverse(join('', @pattern));
}

########################################################################
sub read_fastq_read
{
    my ($fastq, $fastq_name) = @_;
    
    my $header = <$fastq> || return (undef, undef, undef);
    my $seq    = <$fastq> || die("FASTQ file '$fastq_name' terminated unexpectadly");
    my $plus   = <$fastq> || die("FASTQ file '$fastq_name' terminated unexpectadly");
    my $qual   = <$fastq> || die("FASTQ file '$fastq_name' terminated unexpectadly");

    chomp($header);
    chomp($seq);
    chomp($qual);
    
    return ($header, $seq, $qual);
}


########################################################################
sub write_fastq_read
{
    my ($fastq, $header, $seq, $qual) = @_;
    
    print $fastq $header, "\n";
    print $fastq $seq, "\n";
    print $fastq "+\n";
    print $fastq $qual, "\n";
}


########################################################################
sub read_sam_alignment
{
    my ($sam, $sam_name, $paired_end, $conversion) = @_;
    
    # Format of an alignmet:
    # Vector of:
    # [0] Reference to a vector holding SAM fields
    # [1] Template conversion
    # [2] Genome conversion   
    # [3] Methylation pattern
    
    my ($align1, $align2);
    
    my $line = <$sam> || return (undef, undef);
    chomp($line);
    my @fields1 = split("\t", $line);
    $align1 = [\@fields1, $conversion, $conversion, undef];
    
    if ($paired_end) {
        my $line = <$sam> || die("SAM file '$sam_name' terminated unexpectadly");
        chomp($line);
        my @fields2 = split("\t", $line);
        $conversion = reverse_conversion($conversion);
        $align2 = [\@fields2, $conversion, $conversion, undef];
        die ("Mismatched read ID between read1 and read2 (SAM): '$fields1[0]'!='$fields2[0]'") if ($fields1[0] ne $fields2[0]);
    }
    
    return ($align1, $align2);
}


########################################################################
sub write_sam_alignment
{
    my ($sam, $align) = @_;
    
    print $sam join("\t", @{$align->[0]});
    print $sam "\tXP:Z:", $align->[3] if (defined $align->[3]);
    print $sam "\tXR:Z:", $align->[1] if (defined $align->[1]);
    print $sam "\tXT:Z:", $align->[2] if (defined $align->[2]);
    print $sam "\n";
}


########################################################################
sub reverse_conversion
{
    my ($conversion) = @_;
    
    return undef unless (defined $conversion);
    
    return 'GA' if ($conversion eq 'CT');
    return 'CT' if ($conversion eq 'GA');
    
    die ("Unknown conversion encountered: '$conversion'");
}


########################################################################
sub print_stats
{
    my ($read_count, $align_count, $multiple_count, $conv_counts) = @_;
    
    my @conversions = ("CT/CT", "CT/GA", "GA/CT", "GA/GA");
    my $failed_count = $read_count - $align_count - $multiple_count;
    
    print  STDERR "# reads processed: $read_count\n";
    printf STDERR "# reads with at least one reported alignment: %d (%.2f%%)\n", $align_count, 100.*$align_count/$read_count;
    printf STDERR "# reads that failed to align: %d (%.2f%%)\n", $failed_count, 100.*$failed_count/$read_count;
    printf STDERR "# reads with alignments suppressed due to -m: %d (%.2f%%)\n", $multiple_count, 100.*$multiple_count/$read_count;
    print  STDERR "# read/genome conversion counts:\n";
    for my $conversion (@conversions) {
        my $count = $conv_counts->{$conversion};
        $count = 0 unless (defined $count);
        print STDERR "#     $conversion: $count\n";        
    }
        
    
}


########################################################################
sub construct_genome
{
    my ($genome_dir) = @_;
    
    return [ $genome_dir, { } ];
}


########################################################################
sub get_reference
{
    my ($genome, $chr, $start, $len) = @_;
    my $prepad = '';

    --$start;
    if ($start < 0) {
        $prepad = 'N' x (-$start);
        $start = 0;
    }

    $start = 0 if ($start < 0);

    my $ref = $genome->[1]->{$chr};
    unless (defined $ref) {
        $ref = read_fasta("$genome->[0]/$chr.fa");
        $genome->[1]->{$chr} = $ref;
    }
    
    my $segment = substr(${$ref}, $start, $len);
    if (length($segment) < $len) {
        $segment .= 'N' x ($len - length($segment));
    }
    
    return $segment;
}


########################################################################
sub read_fasta
{
    my ($fname) = @_;
    
    open(my $fasta, '<', $fname) or die("Cannot open FASTA file '$fname'");
    my @lines = <$fasta>;
    close($fasta);
    
    shift(@lines);
    chomp(@lines);
    my $seq  = join('', @lines);
    $seq =~ tr|actg|ACTG|;
    
    return \$seq;
}


########################################################################
sub reverse_complement
{
    my ($seq) = @_;
    
    $seq = reverse($seq);
    $seq =~ tr/ACGT/TGCA/;
    
    return $seq;
}

########################################################################


__END__


=head1 NAME

bissli2-align.pl - Align bisulfite converted reads to a reference genome using bowtie2.


=head1 SYNOPSIS

 bissli2-align.pl  [OPTIONS]  -g <genome>  -x <index>  {-1 <m1> -2 <m2> | -U <r>}
 

=head1 ARGUMENTS

=over 4

=item B<ebwt>

The basename of the index to be searched. The index is created using bissli2-build. 

=item B<input-file>

Input FASTQ file holding reads. Note that other input file formats are unsupported.

=item B<output-file>

Output SAM file that will hold aligned reads. Note that output file must always be in SAM format.

=back


=head1 OPTIONS

=over 4

=item B<-h, --help>

Show this help message and exit.

=item B<-1 <m1>> B<-2 <m2>>

Files containing paired-end reads. <m1> and <m2> are comma-separated lists of filenames containing read 1 and read 2 respectively. 

=item B<--bowtie2 <path>>

Path to the bowtie2 executable. [bowtie2]

=item B<--ct>

Assume that the reads to be aligned underwent C->T conversion.
If paired-end reads are aligned, then assume read 1 underwent C->T conversion while read 2 underwent G->A conversion. 
The --ct and --ga options are mutually exclusive.

=item B<-g <dir>>

A directory holding the reference genome, one FASTA file per chromosome.

=item B<--ga>

Assume that the reads to be aligned underwent G->A conversion.
If paired-end reads are aligned, then assume read 1 underwent G->A conversion while read 2 underwent C->T conversion. 
The --ct and --ga options are mutually exclusive.

=item B<--keep-tmp>

Don't remove temporary directory once the run complete. For debug purposes only.

=item B<-S <filename>>

Name of output SAM file. If this option is not set, output will be written to stdout in SAM format.

=item B<--tmp-dir <dir>>

Directorhy for storing temporary files. If this option is not set, temporary files will be stored under /tmp. 

=item B<-U <r>>

Files containing single-end reads. <r> is a comma-separated list of filenames. 

=item B<bowtie2 Options>

Any unknown options are passed transparently to bowtie2. Note that the following bowtie2 options are always applied: -q --quiet --reorder.

=back


=head1 OUTPUT FORMAT

Output is always in SAM format. All optional fields generated by bowtie2 are maintained in the output. In addition, the following optional fields are added:

=over 4

=item B<XP>

Methylation call string, using Bismark's conventions. This tag is identical to Bismark's XM tag.

=item B<XR>

Read conversion state for the alignment. This tag is identical to Bismark's XR tag.

= item B<XT>

Genome conversion state for the alignment. This tag is identical to Bismark's XG tag.

=back



=back
