# get ngt R input
while(<>){
    
    next if $_ !~ /^\s*(ENSG\d+)\s+.*\d\s+(\[.*\])\s*$/;
    
    $g = $1;
    $sig = $2;
    
    run_check($g, $sig);

}

sub run_check{
    
    my ($gene, $sig) = @_;

    `rm _ngt.dat` if -e "_ngt.dat";
    `rm _t2d.dat` if -e "_t2d.dat";
    `rm ngt.dat` if -e "ngt.dat";
    `rm t2d.dat` if -e "t2d.dat";
    

    $sig=~s/\[/ /g;
    $sig =~s/\]/ /g;
    my @list = split /\s+/,$sig;
    shift @list until $list[0]=~/^\S/;
    my @ngt;
    my @t2d;
    foreach $id (@list){
	$id =~/^(\S+)\:(\d)$/;
	push @ngt, $1;
	push @t2d, $1; 
    }

    
    `grep ngt eqtl_data/$gene.dat | head -1 > ngt.dat`;
    foreach $snp (@ngt){
	`grep $snp eqtl_data/$gene.dat >> ngt.dat`;
    }

    `grep t2d eqtl_data/$gene.dat | head -1 > t2d.dat`;
    foreach $snp (@t2d){
	`grep $snp eqtl_data/$gene.dat >> t2d.dat`;
    }
    
    `grep ngt ngt.dat | trans_file - | grep -v ngt | grep -v geno > _ngt.dat`;
    `grep t2d t2d.dat | trans_file - | grep -v t2d | grep -v geno  > _t2d.dat`;
    `rm ngt.dat t2d.dat`;

    $out = `Rscript reg.R`;
    #`rm _ngt.dat _t2d.dat`;
    
    my @list = split /\n/, $out;
    
    print "\n=========================== $gene $sig ============================\n";

    foreach $l (@list){
	
	if($l =~ /NGT/){
	    $nv = join " + ", @ngt; 
	    print "\n\t\tNGT: $gene \~ $nv\n";
	    print "\nAllele frequency\n";
	    $flag =1;
	    next;
	}
	
	
	if($l =~/T2D/){
	    $nv = join " + ", @t2d; 
	    print "\n\t\tT2D: $gene \~ $nv\n";
	    print "\nAllele frequency\n";
	    $flag =1;
	    next;
	}

	if($l=~/COR/){
	    print "\nGenotype Correlation\n";
	    next;
	}

	if($l=~/REG/){
	    $flag = 0;
	    
	    next;
	}
	
	if($l=~/Estimate/){
	    print "\n";
	    $flag = 1;
	}
	if($l=~/^\s*---/||$l=~/^\s*Residual/){
	    $flag = 0;
	}
	
	next if $flag==0;
	print "$l\n";

    }
        

}


