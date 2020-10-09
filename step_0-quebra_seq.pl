#!/usr/bin/perl

$len = $ARGV[0];

$a = qx(more ~jorge/ref/sars-cov-2.ref/ASM985889v3.fa | grep -v "^>");
$a =~ s/\n//eg;

@COVID = split(//,$a);

for ($k=0; $k<=$#COVID-$len; $k++){

    print "Seq_$k\t";
    for ($i=$k; $i<$k+$len; $i++){
	print "$COVID[$i]";
    }print "\n";
}
