#!/usr/bin/perl

while(<>){
    
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    if($data[0]=~/size/){
	shift @data;
	@size = @data;
    }

    if($data[0]=~/corr/){
	shift @data;
	@corr = @data;
    }

}

#print "@size\n@corr\n";

foreach $c (@corr){
    foreach $s (@size){
	
	$y = $s*$c**0.5;
	$x = $s*(1-$c)**0.5;
	
	printf "%.5f  %.5f\n", $x, $y;
    }
}
	    
    
