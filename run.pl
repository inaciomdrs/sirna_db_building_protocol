#!/usr/bin/perl

$sLen = $ARGV[0];

if ($sLen ne "" ){
    
    qx(split -d -a 4  -l 748 input input. );
    
    @roda = qx(ls input.*); chomp(@roda);

    qx(./makeH.pl $sLen);
    qx(cat head.txt >  input.0000.res );
    
    for ($i=0; $i <= $#roda; $i++){
	qx(./db_olig_seq2.pl $roda[$i] $sLen >> $roda[$i].res &);
    }
    
} else {
    print "Missing the size of the sirna (between 18 and 21)\n";
}


#Do this after finishing the background processes;
#qx(cat head *.res > db.txt );
#qx(rm input.*);
