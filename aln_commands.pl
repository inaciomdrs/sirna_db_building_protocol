#!/usr/bin/perl

open(A,$ARGV[0]);

open(NAT,">NAT.$ARGV[1]");
open(SIN,">SIN.$ARGV[1]");
open(AS,">ANT.$ARGV[1]");
while(<A>){
    chomp;

    @aux = split(/\t/,$_);

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

    # If you do not wish to perform alignment steps of sequences
    # against human genome, human coding transcriptome, human
    # non-coding transcriptome, MERS, SARS, and H1N1 genomes,
    # you can download the results files at 
    # https://osf.io/a3x8m/, and delete or comment lines
    # 54, 55, 62, 63, 70, 71, 78, 79, 86, 87, 94 and 95

    # Replace path/to/hg_genome/hg19/Sequence/BowtieIndex/genome
    # with the path to fasta file containing the sequence of
    # Human (GRCh37) genome.
    # IMPORTANT NOTE: do not alter the name "genome", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/hg_genome/hg19/Sequence/BowtieIndex/genome";
    &genomic_match($_[0],$genome,150);
    
    # Replace path/to/hs_cDNA/GRCh38_cds
    # with the path to fasta file containing the sequence of
    # Human (GRCh38) coding transcriptome.
    # IMPORTANT NOTE: do not alter the name "GRCh38_cds", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/hs_cDNA/GRCh38_cds";
    &genomic_match($_[0],$genome,220);
    
    # Replace path/to/hs_cDNA/GRCh38_cds_ncrna
    # with the path to fasta file containing the sequence of
    # Human (GRCh38) non-coding transcriptome.
    # IMPORTANT NOTE: do not alter the name "GRCh38_cds_ncrna", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/hs_cDNA/GRCh38_cds_ncrna";
    &genomic_match($_[0],$genome,220);
    
    # Replace path/to/genomes_v_set/MERS
    # with the path to fasta file containing the sequence of
    # MERS genome.
    # IMPORTANT NOTE: do not alter the name "MERS", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/genomes_v_set/MERS";
    &genomic_match($_[0],$genome,220);
    
    # Replace path/to/genomes_v_set/SARS2010
    # with the path to fasta file containing the sequence of
    # SARS genome.
    # IMPORTANT NOTE: do not alter the name "SARS2010", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/genomes_v_set/SARS2010";
    &genomic_match($_[0],$genome,220);
    
    # Replace path/to/genomes_v_set/H1N1
    # with the path to fasta file containing the sequence of
    # H1N1 genome.
    # IMPORTANT NOTE: do not alter the name "H1N1", otherwise
    # subsequent programs won't be able to find it.
    $genome ="path/to/genomes_v_set/H1N1";
    &genomic_match($_[0],$genome,220);
    
    ##########################################################

    # Here you can use any country name, so as to add or remove 
    # countries. Add or remove as many countries as you wish or
    # may be necessary. 
    # 
    # IMPORTANT NOTE: keep the pattern 
    # `path/to/gisaid/CountryName`, where "CountryName" may not
    # have spaces or special characters. 

    # Replace path/to/gisaid/Brazil
    # with the path to fasta file containing the sequence of
    # Brazil SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Brazil";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/Wuhan
    # with the path to fasta file containing the sequence of
    # China (Wuhan-region only) SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Wuhan";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/China
    # with the path to fasta file containing the sequence of
    # China (without Wuhan-region) SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/China";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/England
    # with the path to fasta file containing the sequence of
    # England SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/England";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/Germany
    # with the path to fasta file containing the sequence of
    # Germany SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Germany";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/Italy
    # with the path to fasta file containing the sequence of
    # Italy SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Italy";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/Russia
    # with the path to fasta file containing the sequence of
    # Russia SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Russia";
    &covid_match($_[0],$genome,10);
    
    # Replace e ="path/to/gisaid/Spain
    # with the path to fasta file containing the sequence of
    # Spain SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/Spain";
    &covid_match($_[0],$genome,10);
    
    # Replace path/to/gisaid/USA
    # with the path to fasta file containing the sequence of
    # USA SARS-CoV-2 strains genomes at GISAID 
    $genome ="path/to/gisaid/USA";
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
