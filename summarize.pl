#!/usr/bin/perl

#cut -f 1,2,3,14 antisense.fa.SARS2010.sam
%olig=();
%olig_hits=();
%hits=();

open(A,"$ARGV[0]");
while(<A>){
    chomp();
    @t =split(/\t/, $_);
    $t[13]=~ s/NM:i:/""/eg;

    if ($t[13] ne "" ){
	if( $olig{$t[0]} eq ""    ){ $olig{$t[0]} = 10000; };
	if( $olig{$t[0]} > $t[13] ){ $olig{$t[0]} = $t[13] };
	#print "$t[1]\t$t[9]\t$t[13]\n";
    }
    if ( $t[1] != 4 ){

	if ( $olig_hits{"$t[0] $t[2]"} != 1 ){
	    $hits{$t[0]} += 1;
	    $olig_hits{"$t[0] $t[2]"} = 1;
	}
	
    }
}
	

while (($key, $value) = each (%olig)) {

    print "$key\t$value\t$hits{$key}\n";
}
