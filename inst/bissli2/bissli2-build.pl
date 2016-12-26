#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
#use Getopt::ArgParse;

########################################################################

my ($bowtie2_build, $reference, $idx_base, @bowtie2_options) = parse_args();

# Split list of reference files
my @reference = split(',', $reference);
my @refs;

# Create a temporary directory
my $tmpdir_obj = File::Temp->newdir('bissli2-build.XXXXXXXX', TMPDIR => 1, CLEANUP => 1) or die("Could not create temporary directory.");
my $tmpdir     = $tmpdir_obj->dirname;

# Create a forward and reverse CT converted version for each input file
for my $input_fname (@reference) {
    my $basename = $input_fname;
    $basename =~ s|^.*/||;
    $basename =~ s|\..*$||;

    #print STDERR "Processing $basename\n";
    
    my $output_fwd_fname = "$tmpdir/${basename}_fwd.fa";
    my $output_rev_fname = "$tmpdir/${basename}_rev.fa";
        
    open(my $input, '<', $input_fname) or die("Cannot open input file '$input_fname'");
    open(my $output_fwd, '>', $output_fwd_fname) or die("Cannot open output file '$output_fwd_fname'");
    open(my $output_rev, '>', $output_rev_fname) or die("Cannot open output file '$output_rev_fname'");

    my $header;    
    my @seqbuf;
    while (my $line = <$input>) {
        chomp($line);
        
        if ((length($line) > 0) and (substr($line, 0, 1) eq '>')) {
            write_converted($output_fwd, $output_rev, $header, \@seqbuf);
            @seqbuf = ();
            $header = $line;
            next;
        }
        
        push(@seqbuf, $line)
    }
    
    write_converted($output_fwd, $output_rev, $header, \@seqbuf);    
    
    close($output_rev);
    close($output_fwd);
    close($input);
    
    push (@refs, $output_fwd_fname);
    push (@refs, $output_rev_fname);    
}

my @cmd      = ($bowtie2_build, '-f', @bowtie2_options, join(',',@refs), $idx_base);

print STDERR "Running command:\n";
print STDERR "@cmd\n";
my $rc = system(@cmd);

die("Failed running bowtie-build (rc: $rc)") unless ($rc == 0);

print STDERR "\nbissli2-build completed successfully.\n";


########################################################################
sub parse_args
{
    my $help           = 0;
    my $reference      = undef;
    my $idx_base       = undef;
    my $bowtie2_build  = 'bowtie2-build';
    my @bowtie2_options;
    
    GetOptions( 'bowtie2-build=s' => \$bowtie2_build ) or pod2usage(-exitval => 2, -verbose => 0);
    
    pod2usage(-exitval => 0, -verbose => 1) if ($help);
    
    $idx_base  = pop(@ARGV);
    $reference = pop(@ARGV);
    
    die("The reference FASTAs and the index's basename must be defined.") unless ((defined $idx_base) and ($reference));
    
    @bowtie2_options = @ARGV;
    
    return ($bowtie2_build, $reference, $idx_base, @bowtie2_options);
}


########################################################################
sub write_converted
{
    my ($output_fwd, $output_rev, $header, $seqbuf) = @_;
    
    return unless ((defined $header) and (@{$seqbuf} > 0));
    
    my $header_fwd = $header;
    $header_fwd =~ s/^([^\s]*)/$1_fwd/;
    print $output_fwd "$header_fwd\n";

    for my $line (@{$seqbuf}) {
        my $converted = $line;
        $converted =~ tr/cC/tT/;
        print $output_fwd $converted, "\n";
    }
    
    my $header_rev = $header;
    $header_rev =~ s/^([^\s]*)/$1_rev/;
    print $output_rev "$header_rev\n";
    
    @{$seqbuf} = reverse(@{$seqbuf}); 

    for my $line (@{$seqbuf}) {
        my $converted = reverse($line);
        $converted =~ tr/aAcCgGtT/tTgGtTaA/;
        print $output_rev $converted, "\n";
    }
}


#######################################################################


__END__


=head1 NAME

bissli2-build.pl - Create an bowtie2 index for a bisulfite converted genome.


=head1 SYNOPSIS

 bissli2-align.pl  [OPTIONS]  <reference>  <idx_base>
 

=head1 ARGUMENTS

=over 4

=item B<reference>

A comma-separated list of FASTA files holding the reference genome. 

=item B<idx_base>

The basename of the index files that will be created.

=back


=head1 OPTIONS

=over 4

=item B<-h, --help>

Show this help message and exit.

=item B<--bowtie2-build <PATH>>

The path of the bowtie2-build executable. [bowtie2-build]

=item B<bowtie2 options>

Any unknown options are passed transparently to bowtie2-build.
Note that bissli2-build can only handle FASTA files as input. Therefore the '-f' option is always passed to bowtie2-build.

=back

