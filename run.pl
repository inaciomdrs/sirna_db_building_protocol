#!/usr/bin/perl

$sLen = $ARGV[0];

if ($sLen ne "" ){
    
    qx(split -d -a 4  -l 748 input input. );
    @roda = qx(ls input.*); chomp(@roda);

    for ($i=0; $i <= $#roda; $i++){

	qx(./db_olig_seq2.pl $roda[$i] $sLen > $roda[$i].res &);
	
    }
} else {
    print "Faltou o Tamanho do sirna (entre 18 e 21)\n";
}
#qx(cat *.res > db.txt );
#qx(rm input.*);
