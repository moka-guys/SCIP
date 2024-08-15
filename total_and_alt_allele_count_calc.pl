#!/usr/bin/perl

use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use Getopt::Std;
use List::Util;
use List::Util qw(any);
use List::Util qw(sum);
use Statistics::Basic qw(:all nofill);


# get options from commandline
my $usage = "usage: dedup_test.pl [options] <arguments...>
options:
	-d <input>
	-i <input>
	-o <output>";

die $usage unless @ARGV;

my %opt;
getopts('d:i:o:',\%opt);

# store commandline options in variables
my $depth = $opt{d};
my $mpileUpFile = $opt{i};
my $output = $opt{o};


# create read / write streams
open(my $r1, "<", $mpileUpFile) or die "Cannot open $mpileUpFile $!";
open(my $fw2, ">", $output) or die "Cannot create $output $!";

my $row;
my @fetalFraction;
my $informativeSNPs = 0;

#print $fw2 "chr\tstart\tend\tnum_reads\tA\tA_Frac\tG\tG_Frac\tC\tC_Frac\tT\tT_Frac.\n";
print $fw2 "chr\tstart\tend\tnum_reads\tA\tG\tC\tT\n"; # setting up header of output txt file

while(my $row = <$r1>){
	
	chomp $row;

	my @columns = split("\t",$row);
	my $chr = $columns[0];
	my $start = $columns[1];
    my $end = $start + 1;
    my $num_reads = $columns[3];
    my $calls = $columns[4];

################################################
    # Ensure num_reads is defined and numeric
    next unless defined $num_reads && $num_reads =~ /^\d+$/;

    # Convert dots and commas to the actual bases
    $calls =~ tr/.,/Aa/ if $columns[2] eq 'A';
    $calls =~ tr/.,/Cc/ if $columns[2] eq 'C';
    $calls =~ tr/.,/Gg/ if $columns[2] eq 'G';
    $calls =~ tr/.,/Tt/ if $columns[2] eq 'T';
###################################################
    # User specifies depth to get accurate fetal fraction
    if($num_reads < $depth){ 
                next;
        }

    my $num_A = 0;
    my $num_G = 0;
    my $num_C = 0;
    my $num_T = 0;

    # Count individual base counts
    while ($calls =~ /[aA]/g) { $num_A++ }
    while ($calls =~ /[gG]/g) { $num_G++ }
    while ($calls =~ /[cC]/g) { $num_C++ }
    while ($calls =~ /[tT]/g) { $num_T++ }

    #print $fw2 $chr."\t".$start."\t".$end."\t".$num_reads."\t".$num_A."\t".$num_G."\t".$num_C."\t".$num_T."\n";
    print $fw2 $chr."\t".$start."\t".$end."\t".$num_reads."\t".$num_A."\t".$num_G."\t".$num_C."\t".$num_T."\n"; # prints to output txt
    print $chr."\t".$start."\t".$end."\t".$num_reads."\t".$num_A."\t".$num_G."\t".$num_C."\t".$num_T."\n"; # prints to console
   }

