#!/usr/bin/perl

$set2 = $ARGV[0];

&Fheader;

sub Fheader {

    @h1 = ("Seq_id","Target_region","","Nat_Sense","oNat_Sense","","Sin_Sense","oSin_Sense","","AntiSense","oAntiSense","oAntiSense","","Annot","","sAA","sTT","sGG","sCC","","fTT","fAA","fCC","fGG","","","A","U","G","C","GC","AU","UUUU","GCCA","QtdPenta80%","","","A","U","G","C","GC","AU","UUUU","GCCA","QtdPenta80%","","","A","U","G","C","GC","AU","UUUU","GCCA","QtdPenta80%","","hepta_Sense","hepta_AS","","Palindromic_N","Palindromic_S","Palindromic_AS","","");

     @h3 = ("Tm","TmSalt","TmNN","RlogK","deltaG","deltaH","deltaS","Hairpin","SelfAnnealing","3´comp","","","Tm","TmSalt","TmNN","RlogK","deltaG","deltaH","deltaS","Hairpin","SelfAnnealing","3´comp","","","Tm","TmSalt","TmNN","RlogK","deltaG","deltaH","deltaS","Hairpin","SelfAnnealing","3´comp","","","Predicted Eficacy","#GG","dG(-1)","dG(-2..-7)","dG(2.6.13)","dG_Best","dG target","dG duplex","dG(18)","dG_self","","","Predicted Eficacy","#GG","dG(-1)","dG(-2..-7)","dG(2.6.13)","dG_Best","dG target","dG duplex","dG(18)","dG_self","","","Predicted Eficacy","#GG","dG(-1)","dG(-2..-7)","dG(2.6.13)","dG_Best","dG target","dG duplex","dG(18)","dG_self","","GOOD","","DDG","DGss","Dgem","DG","","","DDG","DGss","Dgem","DG","","","DDG","DGss","Dgem","DG","","","GOOD","DG","DDG","","","GOOD","DG","DDG","","","GOOD","DG","DDG");


    $print1="";
    for(my $i=0; $i <=$#h1; $i++){
	$print1 .= "$h1[$i]\t";
    }

    $print3="";
    for(my $i=0; $i <=$#h3; $i++){
	$print3 .= "$h3[$i]\t";
    }
    chop($print3);
    
    $hh2=0;
    &makeH2("NAT");
    &makeH2("SIN");
    &makeH2("ANT");

    $print2="";
    for(my $i=0; $i <=$#h2; $i++){
	$print2 .= "$h2[$i]\t";
    }

    #print "$print2\n";
    qx( echo "$print1$print2$print3"  > head.txt );

}

sub makeH2 {
    $h2[$hh2++] = "hs";
    $h2[$hh2++] = "hs_cds";
    $h2[$hh2++] = "hs_ncrna";

    $h2[$hh2++] = "mers";
    $h2[$hh2++] = "sars";
    $h2[$hh2++] = "h1n1";
    $h2[$hh2++] = "";

    my @tmp = qx(ls STS/ | grep "^$_[0]" | grep ".$set2.sts"); chomp(@tmp);
    for(my $i=0; $i <= $#tmp; $i++){
	my @a = split(/\./, $tmp[$i]);
	if (($a[1] !~ /_cds/) and ($a[1] !~ /genome/) and ($a[1] !~ /MERS/) and ($a[1] !~ /SARS2010/) and ($a[1] !~ /H1N1/)) {
	    #print "$a[1]\t"; #pais
	    $h2[$hh2++] = "$a[1]";
	}
    }
    $h2[$hh2++] = "";
    $h2[$hh2++] = "";
}
