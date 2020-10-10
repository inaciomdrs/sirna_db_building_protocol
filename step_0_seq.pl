#!/usr/bin/perl

$len = $ARGV[0];
$len += 2;

# Replace path/to/SARS_CoV_2/Wuhan_strain/ASM985889v3.fa
# with the path to fasta file containing the sequence of
# SARS_CoV_2 Wuhan strain at NCBI Assembly (code ASM985889v3)
$a = qx(more path/to/SARS_CoV_2/Wuhan_strain/ASM985889v3.fa | grep -v "^>");
$a =~ s/\n//eg;

@COVID = split(//,$a);

for ($k=0; $k<=$#COVID-$len; $k++){

    print "Seq_$k\t";
    for ($i=$k; $i<$k+$len; $i++){
	print "$COVID[$i]";
    }print "\n";
}
