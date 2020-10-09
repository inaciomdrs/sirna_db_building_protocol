#!/usr/bin/perl -w
#Written by Malachi Griffith
#This script gets Tm scores for probe sequences using the nearest neighbour method as described at 
#http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo
#It can be used to get the Tm for a single probe or alternatively you can provide a tab-delimited input file 
#containing probe sequences as long as you specify the number of the column containing these sequences
#The results will be printed to output with the Tm appended onto the input data.

use strict;
use Data::Dumper;
use Getopt::Std;

getopts("s:f:c:");
use vars qw($opt_s $opt_f $opt_c);

#Usage instructions 
unless ($opt_s || $opt_f){
  print "\nThis function calculates the Nearest Neighbor Tm of an oligo";
  print "\nSequence must be at least 8 oligos long and contain at least one G/C";
  print "\nTo test with a single probe sequence try the following";
  print "\nUsage: get_NearNeighbor_Tm.pl -s ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA\n";
  print "\nTo generate for a large number of probes, provide a tab-delimited input file containing probe sequences";
  print "\nSpecify which column contains the probe sequences, the file will be regenerated with the Tm's appended";
  print "\nUsage: get_NearNeighbor_Tm.pl -c 3 -f ~/sequences.tsv\n\n";
  exit();
}

#Allow test with a single probe 
if ($opt_s){
  my $probe_seq = $opt_s;
  my $temp_k = tmCalc ('-sequence'=>$probe_seq, '-silent'=>0);
  my $Tm_celsius = tmConverter ('-tm'=>$temp_k, '-scale'=>'Kelvin');
  $Tm_celsius = sprintf("%.1f", $Tm_celsius);
  print "\tTm_C=\t$Tm_celsius\n";
  exit();
}

unless ($opt_f && $opt_c){
  print "\nPlease specify an input file and the column containing the probe sequences";
  exit();
}

my $probe_column = $opt_c-1;

my $infile = $opt_f;

open (INFILE, "$infile") || die "\nCould not open $infile";

#Get all data from the probe input file
my $probe_count = 0;
while (<INFILE>){

  #Deal with the header line
  if ($probe_count == 0){
    $probe_count++;
    chomp($_);
    print "$_\tTm_Celsius\n";
    next;
  }
  $probe_count++;

  my @arr = split("\t", $_);
  my $probe_seq = $arr[$probe_column];

  my $temp_k = tmCalc ('-sequence'=>$probe_seq, '-silent'=>1);
  my $Tm_celsius = tmConverter ('-tm'=>$temp_k, '-scale'=>'Kelvin');

  chomp($_);
  my $Tm_rounded = sprintf("%.3f", $Tm_celsius);
  print "$_\t$Tm_rounded\n";

}

close INFILE;

exit();


#######################
#tmCalc               #
#######################
sub tmCalc{
  my %args = @_;

  my $seq = $args{'-sequence'};
  my $silent = $args{'-silent'};

  #Create a hash to store the deltaH and deltaS values for every possible 2 bp combination (AA/TT, AT/TA, etc.)
  my %thermo;

  #DeltaH values are calories/mol - Taken from Sugimoto, et al., NAR, 1996 - As done in Pan et al., Molecular Cell 2004.
  $thermo{deltaH}{aa} = -8000;
  $thermo{deltaH}{tt} = -8000;
  $thermo{deltaH}{at} = -5600;
  $thermo{deltaH}{ta} = -6600;
  $thermo{deltaH}{ca} = -8200;  $thermo{deltaH}{tg} = -8200;
  $thermo{deltaH}{ct} = -6600;  $thermo{deltaH}{ag} = -6600;
  $thermo{deltaH}{ga} = -8800;  $thermo{deltaH}{tc} = -8800;
  $thermo{deltaH}{gt} = -9400;  $thermo{deltaH}{ac} = -9400;
  $thermo{deltaH}{cg} = -11800;
  $thermo{deltaH}{gc} = -10500;
  $thermo{deltaH}{gg} = -10900;
  $thermo{deltaH}{cc} = -10900;

  #DeltaS values are calories/mol/K
  $thermo{deltaS}{aa} = -21.9;
  $thermo{deltaS}{tt} = -21.9;
  $thermo{deltaS}{at} = -15.2;
  $thermo{deltaS}{ta} = -18.4;
  $thermo{deltaS}{ca} = -21.0;  $thermo{deltaS}{tg} = -21.0;
  $thermo{deltaS}{ct} = -16.4;  $thermo{deltaS}{ag} = -16.4;
  $thermo{deltaS}{ga} = -23.5;  $thermo{deltaS}{tc} = -23.5;
  $thermo{deltaS}{gt} = -25.5;  $thermo{deltaS}{ac} = -25.5;
  $thermo{deltaS}{cg} = -29.0;
  $thermo{deltaS}{gc} = -26.4;
  $thermo{deltaS}{gg} = -28.4;
  $thermo{deltaS}{cc} = -28.4;

  #Initialize counters for each of the possible adjacent nucleotide pairs
  my %neighbour_counts;
  $neighbour_counts{aa} = 0;
  $neighbour_counts{tt} = 0;
  $neighbour_counts{at} = 0;
  $neighbour_counts{ta} = 0;
  $neighbour_counts{ca} = 0;  $neighbour_counts{tg} = 0;
  $neighbour_counts{ct} = 0;  $neighbour_counts{ag} = 0;
  $neighbour_counts{ga} = 0;  $neighbour_counts{tc} = 0;
  $neighbour_counts{gt} = 0;  $neighbour_counts{ac} = 0;
  $neighbour_counts{cg} = 0;
  $neighbour_counts{gc} = 0;
  $neighbour_counts{gg} = 0;
  $neighbour_counts{cc} = 0;

  my $seq_length = length($seq);

  unless($seq_length > 1){
    print "\nProbe sequence has length of 1 or less\n\n";
    exit();
  }

  #Go through the probe sequence and get all the nucleotide neighbours 
  for (my $i = 0; $i <= $seq_length-2; $i++){

    my $basepair1 = substr($seq, $i, 1);
    my $basepair2 = substr($seq, $i+1, 1);
    my $couple = "$basepair1"."$basepair2";
    my $lower_couple = lc($couple);

    #Count the neighbors by incrementing the appropriate counter
    $neighbour_counts{$lower_couple}++;
  }

  my $h_nuc = 5000; #Calories per mol
  $h_nuc = 3400; #According to Sugimoto - at least thats what this website says 

  #Reproduce the calculation at http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo
  #NOTE:  ** This calculation will give you the Nearest Neighbor Value shown on this website.**

  #The nearest neighbor and thermodynamic calculations are done essentially as described by Breslauer et al., (1986) Proc. Nat. Acad. Sci. 83:3746-50
  #but using the values published by Sugimoto et al., (1996) Nucl. Acids Res. 24:4501-4505 (Abstract). 
  #assumes that the sequences are not symmetric and contain at least one G or C. The minimum length for the query sequence is 8.

  #To find the deltaH multiply the deltaH value for each couple by the number of occurences of each pair, add all these up and add the initiation value 'h_nuc'
  #in Calories/mol
  my $deltaH = - (($neighbour_counts{aa} * $thermo{deltaH}{aa}) + ($neighbour_counts{tt} * $thermo{deltaH}{tt}) +
		  ($neighbour_counts{at} * $thermo{deltaH}{at}) + ($neighbour_counts{ta} * $thermo{deltaH}{ta}) +
		  ($neighbour_counts{ca} * $thermo{deltaH}{ca}) + ($neighbour_counts{tg} * $thermo{deltaH}{tg}) +
		  ($neighbour_counts{ct} * $thermo{deltaH}{ct}) + ($neighbour_counts{ag} * $thermo{deltaH}{ag}) +
		  ($neighbour_counts{ga} * $thermo{deltaH}{ga}) + ($neighbour_counts{tc} * $thermo{deltaH}{tc}) +
		  ($neighbour_counts{gt} * $thermo{deltaH}{gt}) + ($neighbour_counts{ac} * $thermo{deltaH}{ac}) +
		  ($neighbour_counts{cg} * $thermo{deltaH}{cg}) + ($neighbour_counts{gc} * $thermo{deltaH}{gc}) +
		  ($neighbour_counts{gg} * $thermo{deltaH}{gg}) + ($neighbour_counts{cc} * $thermo{deltaH}{cc}));



  #in calories/mol/K
  my $deltaS = - (($neighbour_counts{aa} * $thermo{deltaS}{aa}) + ($neighbour_counts{tt} * $thermo{deltaS}{tt}) +
		  ($neighbour_counts{at} * $thermo{deltaS}{at}) + ($neighbour_counts{ta} * $thermo{deltaS}{ta}) +
		  ($neighbour_counts{ca} * $thermo{deltaS}{ca}) + ($neighbour_counts{tg} * $thermo{deltaS}{tg}) +
		  ($neighbour_counts{ct} * $thermo{deltaS}{ct}) + ($neighbour_counts{ag} * $thermo{deltaS}{ag}) +
		  ($neighbour_counts{ga} * $thermo{deltaS}{ga}) + ($neighbour_counts{tc} * $thermo{deltaS}{tc}) +
		  ($neighbour_counts{gt} * $thermo{deltaS}{gt}) + ($neighbour_counts{ac} * $thermo{deltaS}{ac}) +
		  ($neighbour_counts{cg} * $thermo{deltaS}{cg}) + ($neighbour_counts{gc} * $thermo{deltaS}{gc}) +
		  ($neighbour_counts{gg} * $thermo{deltaS}{gg}) + ($neighbour_counts{cc} * $thermo{deltaS}{cc}));

  #This formula corresponds to that found at the website:
  #http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo

  #Assume a salt concentration of 50 mM (milliMolar)
  my $na_salt_conc = (50*1e-3);

  #Assume a primer concentration of 50 nM (nanoMolar)
  my $primer_conc = (50*1e-9);

  #The Gas constant R = 1.987 cal/mol * K
  #Assumes annealing occurs at pH of 7.0
  #Tm assumes the sequences are not symmetric and contain at least one G or C
  #The oligo should be at least 8 nucleotides long for accurate calculations
  #Calculations use 50nM primer conc, and 50 mM salt conc.

  
  my $tm_k = ($deltaH - $h_nuc)/($deltaS + (1.987 * log(1/$primer_conc))) + (16.6*(log($na_salt_conc)/log(10)));

  unless ($silent == 1){
      $a= $deltaH/1000;
      print "DeltaH=\t$a\tDeltaS=\t$deltaS\tTm_K=\t$tm_k";
  }

  return($tm_k);

}

############################
#tmConverter               #
############################
sub tmConverter{

  my %args = @_;
  my $tm = $args{'-tm'};
  my $scale = $args{'-scale'};

  my $tm_convert;

  unless ($tm){
    print "\nRequired parameter for tmConverter missing\n\n";
    exit();
  }

  unless ($scale eq "Celcius" || $scale eq "Kelvin"){
    print "\nTemperature scale specification not understood by tmConverter\n\n";
    exit();
  }

  #If the temperature provide was in degrees Kelvin:
  if ($scale eq "Kelvin"){
    $tm_convert = $tm-273.15;
  }

  #If the temperature provided was in degrees celcius:
  if ($scale eq "Celcius"){
    $tm_convert = 273.15+$tm;
  }
  return($tm_convert)
}

