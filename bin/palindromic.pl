#!/usr/bin/perl

$seq = $ARGV[0];
$seq =~ s/U/T/eg;

chomp($seq);
$min = $ARGV[1];
$max = $ARGV[2];
$a="";
$a =  &find_palindromic_seqs($seq,$min,$max);

if ( $a eq "" ){$a = "0000\n"}
print "$a"; 


sub find_palindromic_seqs {    
    my $seq = $_[0];
    my $min = $_[1];
    my $max = $_[2];

    my $result="";
    my $seq_len=length($seq);

    for($i=0;$i<$seq_len-$min+1;$i++){
	$j=$min;
	while($j<$max+1 and ($i+$j)<=$seq_len){
	    $sub_seq=substr($seq,$i,$j);
		
	    my $aux = &RevComp($sub_seq);
	    my $p   = $i+1;
	    if ( $sub_seq eq  $aux  ){
		$result .= "$p $sub_seq\n";
	    }
	    $j++;
	}
	
    }
    
    return $result;
}

sub RevComp {
    my $seq = $_[0];
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return reverse($seq);
}
