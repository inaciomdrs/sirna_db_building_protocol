#!/usr/bin/perl

open(A,$ARGV[0]);

open(NAT,">NAT.$ARGV[1]");
open(SIN,">SIN.$ARGV[1]");
open(AS,">ANT.$ARGV[1]");
while(<A>){
    chomp;

    @aux = split(/\t/,$ _);

    $t_len = length($aux[1]);
    $o_len = $t_len-2;
    
    ###Nat-Sense Olig    -----------------------------
    $Nat_Sense  =  substr($aux[1],2,$o_len);
    print NAT ">$aux[0]\n$Nat_Sense\n";
        
    ###Sin-Sense Olig    -----------------------------
    $Sin_Sense  =  substr($aux[1],2,$o_len-2)."TT";
    print SIN ">$aux[0]\n$Sin_Sense\n";
        
    ###AntiSense Olig    -----------------------------
    $AntiSense  =  &Comp(substr($aux[1],0,$o_len));
    $AntiSense  =  reverse($AntiSense);
    print AS ">$aux[0]\n$AntiSense\n";
}
close(NAT);
close(SIN);
close(AS);
close(A);


&alin_seq("NAT.$ARGV[1]");
&alin_seq("SIN.$ARGV[1]");
&alin_seq("ANT.$ARGV[1]");


sub alin_seq {

    $genome ="~jorge/ref/hg_genome/hg19/Sequence/BowtieIndex/genome";
    &genomic_match($_[0],$genome,150);
    
    $genome ="~jorge/ref/hs_cDNA/GRCh38_cds";
    &genomic_match($_[0],$genome,220);
    
    $genome ="~jorge/ref/hs_cDNA/GRCh38_cds_ncrna";
    &genomic_match($_[0],$genome,220);
    
    $genome ="~jorge/ref/genomes_v_set/MERS";
    &genomic_match($_[0],$genome,220);
    
    $genome ="~jorge/ref/genomes_v_set/SARS2010";
    &genomic_match($_[0],$genome,220);
    
    $genome ="~jorge/ref/genomes_v_set/H1N1";
    &genomic_match($_[0],$genome,220);
    
    ##########
    
    $genome ="~jorge/ref/gisaid/Brazil";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/Wuhan";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/China";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/England";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/Germany";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/Italy";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/Russia";
    &covid_match($_[0],$genome,10);

    $genome ="~jorge/ref/gisaid/Spain";
    &covid_match($_[0],$genome,10);
    
    $genome ="~jorge/ref/gisaid/USA";
    &covid_match($_[0],$genome,10);
}



sub genomic_match {   
    my @tmp =split(/\//, $_[1]);
    $_[0]=~ s/\.$o_len/""/eg;
    
    print "bowtie -S -a --pairtries 4 -p 30 -n 3 -e $_[2] -l 7 $_[1] -f $_[0].$o_len | ./summarize.pl - > $_[0].$tmp[$#tmp].$o_len.sts\n";
}

sub covid_match {
    my @tmp =split(/\//, $_[1]);
    $_[0]=~ s/\.$o_len/""/eg;
    
    print "bowtie -S -a --pairtries 4 -p 30 -n 3 -e $_[2] -l 7 $_[1] -f $_[0].$o_len | ./summarize.pl - > $_[0].$tmp[$#tmp].$o_len.sts\n";
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
