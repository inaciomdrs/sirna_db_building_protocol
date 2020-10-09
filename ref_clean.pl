#!/usr/bin/perl

open(A,"$ARGV[0]");

$f=0;
while(<A>){
    chomp();
    $_ =~ s/\r/""/eg;
    
    if ($_ =~ /^>/){
	
	if ($f == 1 ){
	    #print "$print\n";
	    
	    $len = length($seq);
	    $N = $seq =~ y/N/N/;
	    $qtd_n =$N ;
	    
	    $PP = $qtd_n/29903;
	    if ( $PP <= 0.10 ){
		print "$print\n$seq\n";
	    }
	    $seq   ="";
	    $print ="";
	    $f=0;
	}
	if ($f == 0 ){
	    $print = "$_";
	    $f = 1;
	} 
    } else {
	$seq  .= uc("$_");
	$f = 1;
    }

}


$len = length($seq);
$N = $seq =~ y/N/N/;
$qtd_n =$N ;

$PP = $qtd_n/29903;
if ( $PP <= 0.10 ){
    print "$print\n$seq\n";
}
