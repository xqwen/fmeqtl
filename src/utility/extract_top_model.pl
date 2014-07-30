@files = <output/*.rst>;
foreach $f (@files){
    
    $f =~/(ENSG\d+)/;
    $g = $1;
    open FILE, "$f";
    while(<FILE>){
	next if $_ !~/\[/;
	print "$g $_" if $_ !~ /NULL/;
	last;
    }


    chomp $out;



}
