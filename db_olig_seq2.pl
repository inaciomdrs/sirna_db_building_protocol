#!/usr/bin/perl

$dep0 = "NC_045512.bed";
%ANN=();
%GH =();

$name = $ARGV[0];
$set  = $name;
$set  =~ s/input./""/eg;
$set2 = $ARGV[1];

&G_hit;
&annot;

open(A,"$name");
while(<A>){
    chomp;

    @aux = split(/\t/,$ _);

    $t_len = length($aux[1]);
    $o_len = $t_len-2;
    
    ###Target Seq        -----------------------------
    print "$aux[0]\t";
    print "5'-$aux[1]-3'\t\t";

    
    ###Nat-Sense Olig    -----------------------------
    $Nat_Sense  =  substr($aux[1],2,$o_len);
    $oNat_Sense =  &olig_seq($Nat_Sense);
 
    print "5'-$Nat_Sense-3'\t";
    print "5'-$oNat_Sense-3'\t\t";

    
    ###Sin-Sense Olig    -----------------------------
    $Sin_Sense  =  substr($aux[1],2,$o_len-2)."TT";
    $oSin_Sense =  &olig_seq($Sin_Sense);
    
    print "5'-$Sin_Sense-3'\t";
    print "5'-$oSin_Sense-3'\t\t";

    
    ###AntiSense Olig    -----------------------------
    $AntiSense  =  &Comp(substr($aux[1],0,$o_len));
    $oAntiSense =  &olig_seq($AntiSense);    

    print "3'-$AntiSense-5'\t";
    print "3'-$oAntiSense-5'\t";

    $aux = reverse($oAntiSense);
    print "5'-$aux-3'\t\t";

    
    ###Annot region Genes -----------------------------
    $ss   = $aux[0];   $ss =~ s/Seq_/""/eg;
    $ann1 = $ANN{$ss};
    $ann2 = $ANN{$ss+$t_len};
    $ann  = &mAnnt($ann1,$ann2);

    if ($ann eq ''){ $ann = "--" };
    print "$ann\t\t";

    
    ###Olig composition -----------------------------

    my $iAA = 0;  if($aux[1] =~ /^AA/){ $iAA = 1;}
    my $iGG = 0;  if($aux[1] =~ /^GG/){ $iGG = 1;}
    my $iCC = 0;  if($aux[1] =~ /^CC/){ $iCC = 1;}
    my $iTT = 0;  if($aux[1] =~ /^TT/){ $iTT = 1;}

    my $fTT = 0;  if($aux[1] =~ /TT$/){ $fTT = 1;}
    my $fCC = 0;  if($aux[1] =~ /CC$/){ $fCC = 1;}
    my $fGG = 0;  if($aux[1] =~ /GG$/){ $fGG = 1;}
    my $fAA = 0;  if($aux[1] =~ /AA$/){ $fAA = 1;}

    print "$iAA\t";
    print "$iTT\t";
    print "$iGG\t";
    print "$iCC\t\t";

    print "$fTT\t";
    print "$fAA\t";
    print "$fCC\t";
    print "$fGG\t\t";

    print "oNat_Sense:\t";
    &olig_composition($oNat_Sense);
    print "\toSin_Sense:\t";
    &olig_composition($oSin_Sense);
    print "\toAntiSense:\t";
    $aux = reverse($oAntiSense);
    &olig_composition($aux);


    ###Olig start hepta composition -----------------------------

    $sub_seqS=substr($oNat_Sense,0,7);
    $sub_seqF=substr(reverse($oAntiSense),0,7);

    $sA =  $sub_seqS =~ y/A/A/; $sT =  $sub_seqS =~ y/U/U/; $pATs = sprintf("%.2f", ($sA+$sT)/7);
    $fG =  $sub_seqF =~ y/G/G/; $fC =  $sub_seqF =~ y/C/C/; $pGCf = sprintf("%.2f", ($fG+$fC)/7);

    print "\t$pATs\t$pGCf\t\t";
    
    
    ###Palindromic sequences  -----------------------------
    
    $palo =qx(./bin/palindromic.pl  $Nat_Sense 6 $o_len | head -n 1 ); chomp($palo); if ($palo eq '0000' ){$palo=0;} else {$palo=1;}
    print "$palo\t";

    $palo =qx(./bin/palindromic.pl  $Sin_Sense 6 $o_len | head -n 1 ); chomp($palo); if ($palo eq '0000' ){$palo=0;} else {$palo=1;}
    print "$palo\t";
    
    $aux = reverse($AntiSense);
    $palo =qx(./bin/palindromic.pl  $aux 6 $o_len | head -n 1 ); chomp($palo); if ($palo eq '0000' ){$palo=0;} else {$palo=1;}
    print "$palo\t\t";
    

    
    ###Genomic, Transcriptome and Covid matches  -----------------------------
    print "Nat_Sense:\t"; &alin_seq2($aux[0],"NAT");
    print "\tSin_Sense:\t"; &alin_seq2($aux[0],"SIN");
    print "\tAntiSense:\t"; &alin_seq2($aux[0],"ANT");
   
            
    
    ###Melting ----------------------------
    $TmI = qx( ./bin/OligoCalc.py $oNat_Sense ); chomp($TmI);
    $hI  = qx( ./bin/HairpinCalc.py        $oNat_Sense 5 4 ); chomp($hI);
    $sI  = qx( ./bin/SelfAnnealingSites.py $oNat_Sense 5 4 ); chomp($sI);
    $tI  = qx( ./bin/complementarity3.py   $oNat_Sense 5 4 ); chomp($tI);
							   
    print "Olig_Calc\toNat_Sense:\t$TmI\t";
    print "$hI\t$sI\t$tI\t";

    $TmI = qx( ./bin/OligoCalc.py $oSin_Sense ); chomp($TmI);
    $hI  = qx( ./bin/HairpinCalc.py        $oSin_Sense 5 4 ); chomp($hI);
    $sI  = qx( ./bin/SelfAnnealingSites.py $oSin_Sense 5 4 ); chomp($sI);
    $tI  = qx( ./bin/complementarity3.py   $oSin_Sense 5 4 ); chomp($tI);
    
    print "\toSin_Sense:\t$TmI\t";
    print "$hI\t$sI\t$tI\t";
    
    $aux = reverse($oAntiSense);
    $TmI = qx( ./bin/OligoCalc.py $aux ); chomp($TmI);
    $hI  = qx( ./bin/HairpinCalc.py        $aux 5 4 ); chomp($hI);
    $sI  = qx( ./bin/SelfAnnealingSites.py $aux 5 4 ); chomp($sI);
    $tI  = qx( ./bin/complementarity3.py   $aux 5 4 ); chomp($tI);
    
    print "\toAntiSense:\t$TmI\t";
    print "$hI\t$sI\t$tI\t";

    
    ###Delta ThermoComposition21  --------------------------------------------
    print "ThermoComposition21\tNat_Sense:\t";
    $aux2 = &RevComp($Nat_Sense); &delta_Composition21($aux2);

    print "\tSin_Sense:\t";
    $aux2 = &RevComp($Sin_Sense); &delta_Composition21($aux2);

    print "\tAntiSense:\t";
    $aux2 = reverse($AntiSense);
    $aux2 = &RevComp($aux2); &delta_Composition21($aux2);

    ###SSD software  --------------------------------------------
    $ssd = qx(./bin/deltacalculator.py $oNat_Sense); chomp($ssd);

    @aux2  = split(/\s+/, $ssd);
    
    $good = 0; if($aux2[1] >= 0 ){ $good =1; }
    print "SSD\t$good\t";
    print "Nat_Sense:\t$ssd\t";
    
    $ssd = qx(./bin/deltacalculator.py $oSin_Sense); chomp($ssd);
    print "\tSin_Sense:\t$ssd\t";

    $aux2 = reverse($oAntiSense);
    $ssd = qx(./bin/deltacalculator.py $aux2); chomp($ssd);
    print "\tAntiSense:\t$ssd\t";

    ###si_shRNA_selector software  --------------------------------------------
    print "si_shRNA_selector\t";
    &si_shRNA_selector($aux[1]);

    ###End                         --------------------------------------------
    print "\n";
}

sub si_shRNA_selector {

    qx(printf ">teste\n$_[0]\n" > $set..1);
    qx(./bin/si_shRNA_selector -len $set2 -i $set..1 -o $set..2 );
    
    my $natural = qx(tail -n 1 $set..2); chomp($natural);
    qx(rm $set..*);
    
    
    my $aux  =  substr($_[0],0,$o_len)."TT";
    qx(printf ">teste\n$aux\n" > $set..1);
    qx(./bin/si_shRNA_selector -len $set2 -i $set..1 -o $set..2 );
    
    my $sin = qx(tail -n 1  $set..2); chomp($sin);
    qx(rm $set..*);
    
    my $aux = &RevComp($_[0]);
    qx(printf ">teste\n$aux\n" > $set..1);
    qx(./bin/si_shRNA_selector -len $set2 -i $set..1 -o $set..2 );
    
    my $ant = qx(tail -n 1  $set..2); chomp($ant);
    qx(rm $set..*);

    my $eff_N = 0;
    my $eff_S = 0;
    my $eff_A = 0;
    if ($natural =~ /efficient/){ $eff_N = 1; }
    if ($sin     =~ /efficient/){ $eff_S = 1; }
    if ($ant     =~ /efficient/){ $eff_A = 1; }

    my @aux = split(/\s+/, $natural);
    print "Nat_Sense:\t$eff_N\t$aux[3]\t$aux[4]\t";

    my @aux = split(/\s+/, $sin);
    print "\tSin_Sense:\t$eff_S\t$aux[3]\t$aux[4]\t";

    my @aux = split(/\s+/, $ant);
    print "\tAntiSense:\t$eff_A\t$aux[3]\t$aux[4]\t";

}

sub delta_Composition21 {

    if ($set2 == 21 ) {
	qx(printf ">teste\n$_[0]\n" > $set.fa);
	qx(./bin/ThermoComposition21 $set.fa > /dev/null);
    
	my $delta0 = qx(tail -n1 $set.csv | cut -d',' -f 4,8,11,12,13,14,16,17,18,19; rm $set.* );
    
	chomp($delta0);
	$delta0 =~ s/\s+/""/eg;
	$delta0 =~ s/,/"\t"/eg;
	print "$delta0\t";
    } else {
	print "-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t";

    }
    
}

sub alin_seq2 {
    
    $genome ="genome";           &genomic_match2($_[0],$genome,$_[1],"hs");
    $genome ="GRCh38_cds";       &genomic_match2($_[0],$genome,$_[1],"hs_cds");
    $genome ="GRCh38_cds_ncrna"; &genomic_match2($_[0],$genome,$_[1],"hs_ncrna");
    
    $genome ="MERS";             &genomic_match2($_[0],$genome,$_[1],"mers");
    $genome ="SARS2010";         &genomic_match2($_[0],$genome,$_[1],"sars");
    $genome ="H1N1";             &genomic_match2($_[0],$genome,$_[1],"h1n1");
    print "\t";
    #------------------------------------------------
    $genome ="Brazil";           &covid_match2($_[0],$genome,$_[1],"Brazil");
    $genome ="Wuhan";            &covid_match2($_[0],$genome,$_[1],"Wuhan");
    $genome ="China";            &covid_match2($_[0],$genome,$_[1],"China");
    $genome ="England";          &covid_match2($_[0],$genome,$_[1],"England");
    $genome ="Germany";          &covid_match2($_[0],$genome,$_[1],"Germany");
    $genome ="Italy";            &covid_match2($_[0],$genome,$_[1],"Italy");
    $genome ="Russia";           &covid_match2($_[0],$genome,$_[1],"Russia");
    $genome ="Spain";            &covid_match2($_[0],$genome,$_[1],"Spain");
    $genome ="USA";              &covid_match2($_[0],$genome,$_[1],"USA");
}

sub G_hit {
    my @tmp = qx(ls STS/ | grep ".$set2.sts"); chomp(@tmp);

    for(my $i=0; $i <= $#tmp; $i++){

	my @a = split(/\./, $tmp[$i]);
	open (K, "STS/$tmp[$i]");
	while (<K>) {
	    chomp;
	    my @b = split(/\t/, $_);
	    $GH{"$a[0] $a[1] $b[0]"} = "$b[1]\t$b[2]";
	}
	close(K);
    }
}


sub genomic_match2 {
    ($nm,) = split(/\t/, $GH{"$_[2] $_[1] $_[0]"} );
    if ($nm eq "" ) {$nm = "-1";}
    print "$nm\t";
}

sub covid_match2 {
    ($l,$nm) = split(/\t/, $GH{"$_[2] $_[1] $_[0]"} );  
    if ($nm eq "" ) {$nm = "0";}
    print "$nm\t";
}

sub olig_composition {

    my $A = $_[0] =~ y/A/A/;
    my $U = $_[0] =~ y/U/U/;
    my $G = $_[0] =~ y/G/G/;
    my $C = $_[0] =~ y/C/C/;
    
    my $gc = sprintf("%.2f", ($G+$C)/$o_len);
    my $au = sprintf("%.2f", ($A+$U)/$o_len);

    print "$A\t";
    print "$U\t";
    print "$G\t";
    print "$C\t";

    print "$gc\t";
    print "$au\t";

    my $UUUU = scalar($_[0] =~ m/(UUUU)/g); if ( $UUUU > 0 ) {$UUUU = 1 } else {$UUUU = 0};
    my $GCCA = scalar($_[0] =~ m/(GCCA)/g); if ( $GCCA > 0 ) {$GCCA = 1 } else {$GCCA = 0};

    print "$UUUU\t";
    print "$GCCA\t";

    my $qtdPenta = &penta($_[0]);
    print "$qtdPenta\t";
}

sub olig_seq {
    my $o = $_[0];
    $o =~ s/T/U/eg;
    $o;
}

sub penta {
    my $seq = $_[0];
    my $qtd = 0;

    for($i=0;$i<= length($seq)-5 ;$i++){
	my $sub_seq=substr($seq,$i,5);
	my $A = $sub_seq =~ y/A/A/;
	my $T = $sub_seq =~ y/U/U/;

	my $AU = ($A+$T)/5;
	if ($AU >= 0.8 ){
	    $qtd++;
	}
    }
    return $qtd ;
}

sub mAnnt {
    my $a1 = $_[0];
    my $a2 = $_[1];

    my @aux = split(/\,/,$a2);
    for (my $i=0; $i<=$#aux; $i++){
	if ($a1 !~ /$aux[$i],/ ) {
	    $a1 .= "$aux[$i],";
	}
    }
    chop($a1);
    return $a1;
}


sub annot {
    open(B, $dep0);
    while (<B>){
	chomp;
	my @aux = split(/\t/,$_);

	for(my $i=$aux[1]; $i<= $aux[2]; $i++){
	    if ($ANN{$i} !~ /$aux[3],/ ) {
		$ANN{$i} .= "$aux[3],";
	    }
	}      
    }
    close(B);
}

sub RevComp {
    my $seq = $_[0];
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return reverse($seq);
}

sub Comp {
    my $seq = $_[0];
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return $seq;
}
