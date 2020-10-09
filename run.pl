#!/usr/bin/perl

$sLen = $ARGV[0];

if ($sLen ne "" ){
    $sLen += 2;
    qx(./step_0-quebra_seq.pl $sLen  > input );
    #qx(head -n 1 input > t; mv t input);
    qx(split -d -a 4  -l 500 input input. );

    $sLen -= 2;
    @roda = qx(ls input.*); chomp(@roda);

    for ($i=0; $i <= $#roda; $i++){

	qx(./db_olig_seq2.pl $roda[$i] $sLen > $roda[$i].res &);
	
    }
} else {
    print "Faltou o Tamanho do sirna (entre 18 e 21)\n";
}
#qx(cat *.res > resultado.txt );
#qx(rm input.*);
