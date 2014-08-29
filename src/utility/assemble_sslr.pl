$subgrp_file = $ARGV[0];
$cis_file = $ARGV[1];
$expr_file = $ARGV[2];
$geno_file = $ARGV[3];
$out_dir = $ARGV[4];


################ read subgroup defn file ################

open SUBDEF, "$subgrp_file";
while(<SUBDEF>){
    chomp;
    next if $_ !~ /\S/;
    /^\s*(\S+)\s+(\S+)/;
    my $id = $1;
    my $grp = $2;
    
    if(!defined($subs{$grp})){
	$subs{$grp} = [];
    }
    push @{$subs{$grp}}, $id;
    
}

# remove single key entry
foreach $k (keys %subs){
    #printf "$k  %d\n",scalar(@{$subs{$k}});
    if(scalar(@{$subs{$k}})==1){
	delete $subs{$k};
    }
}

@grps = sort (keys %subs);



############# read cis defn file ##################


open CISDEF, "$cis_file";
while(<CISDEF>){
    
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    if(!defined($rcd{$data[0]})){
	$rcd{$data[0]} = [];
	push @glist, $data[0];
    }
    $grcd{$data[1]}=1;
    push  @{$rcd{$data[0]}}, $data[1];
}


###################### read expr file ####################


my @ex_files;

if($ARGV[2]=~/\#/){
    foreach $gp (@grps){
	my $f = $ARGV[2];
	$f =~s/\#/$gp/;
	push @ex_files, $f;
    }
}else{
    push @ex_files, $ARGV[2];
}

#print "$ARGV[2]\n@ex_files\n";
foreach $exf (@ex_files){
    
    if($exf =~/\.gz$/){
	open EXPR, "zcat $exf |";
    }else{
	open EXPR, "$exf";
    }

    my $counter = 1;
    my %ex_index;
    my @ngrps;

    while(<EXPR>){

	next if $_ !~ /\S/;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
	
	if ($counter==1){
	    my %map;
	    @map{@data} = (0..$#data);
	    foreach $g (@grps){
		my @iset;
		foreach $id (@{$subs{$g}}){
		    if(defined($map{$id})){
			push @iset, $map{$id};
		    }
		}
		if($#iset >=0){
		    push @ngrps,$g;
		    $ex_index{$g} = \@iset;
		}
	    }

	    #foreach $g (@ngrps){
	    #    printf "$g  %d\n", scalar(@{$ex_index{$g}});
	    #}
	

	    $counter++;
	    next;
	
	}
    
	if(defined($rcd{$data[0]})){
	    
	    foreach $g (@ngrps){
		$gene{$data[0]} .= sprintf "pheno $data[0] $g  @data[@{$ex_index{$g}}]\n"; 
	    
	    }
	}
    
    }
}


###################### read geno file ####################



my @ge_files;

if($ARGV[3]=~/\#/){
    foreach $gp (@grps){
        my $f = $ARGV[3];
        $f =~s/\#/$gp/;
        push @ge_files, $f;
    }
}else{
    push @ge_files, $ARGV[3];
}


foreach $gef (@ge_files){
    
    if($gef =~/\.gz$/){
	open GENO, "zcat $gef |";
    }else{
	open GENO, "$gef";
    }

    $counter = 1;
    my %g_index;
    my @ngrps;
    
    while(<GENO>){

	next if $_ !~ /\S/;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
    
	if ($counter==1){
	    my %map;
	    @map{@data} = (0..$#data);
	    
	    foreach $g (@grps){
		my @iset;
		foreach $id (@{$subs{$g}}){
		    if(defined($map{$id})){
			push @iset, $map{$id};
		    }
		}
		if($#iset >=0){
		    push @ngrps,$g;
		    $g_index{$g} = \@iset;
		}
	   }	
	
	    $counter++;
	    next;
	}
    
	if(defined($grcd{$data[0]})){
	
	    foreach $g (@ngrps){
		$snp{$data[0]} .= sprintf "geno $data[0] $g  @data[@{$g_index{$g}}]\n";
	    }
	}
    }
}


###################### assemble everything ##################

foreach $g (@glist){

    open OUT, ">$out_dir/$g.dat";
    print OUT "$gene{$g}";
    foreach $s ( @{$rcd{$g}}){
	print OUT "$snp{$s}";
    }
    close OUT;
}
