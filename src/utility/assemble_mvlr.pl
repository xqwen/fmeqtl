$subgrp_file = $ARGV[0];
$cis_file = $ARGV[1];
$expr_file = $ARGV[2];
$geno_file = $ARGV[3];
$out_dir = $ARGV[4];


################ read subgroup defn file ################

open SUBDEF, "$subgrp_file";
while(<SUBDEF>){
    chomp;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    $data[0]=~tr/a-z/A-Z/;
    if($data[0] eq "GRP" || $data[0] eq "GROUP"){
	shift @data;
	@grps = @data;
    }
    if($data[0] eq "ID"){
	shift @data;
	@ids = @data;
    }
    
}



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


foreach $gp (@grps){
    my $f = $ARGV[2];
    $f =~s/\#/$gp/;
    push @ex_files, $f;
}


for $i (0..$#grps){

    $exf = $ex_files[$i];
    $cgp = $grps[$i];

    if($exf =~/\.gz$/){
	open EXPR, "zcat $exf |";
    }else{
	open EXPR, "$exf";
    }

    my $counter = 1;
    my @ex_ids;

    while(<EXPR>){

	next if $_ !~ /\S/;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
	
	if ($counter==1){
	    my %map;
	    @map{@data} = (0..$#data);
	    
	    
	    foreach $id (@ids){
		if(defined($map{$id})){
		    push @ex_ids, $map{$id};
		}
	    }

	    
	    $counter++;
	    next;
	    
	}
    
	if(defined($rcd{$data[0]})){
	    $gene{$data[0]} .= sprintf "pheno $cgp  @data[@ex_ids]\n"; 	    

	}
    
    }
}


###################### read geno file ####################




$gef =  $ARGV[3];


if($gef =~/\.gz$/){
    open GENO, "zcat $gef |";
}else{
    open GENO, "$gef";
}

$counter = 1;
my @g_ids;

    
while(<GENO>){
  
    next if $_ !~ /\S/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    if ($counter==1){
	my %map;
	@map{@data} = (0..$#data);
	
	    
	
	foreach $id (@ids){
	    if(defined($map{$id})){
		push @g_ids, $map{$id};
	    }
		
	    $counter++;
	    next;
	}
    }
    if(defined($grcd{$data[0]})){
	
	$snp{$data[0]} .= sprintf "geno $data[0]  @data[@g_ids]\n";
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
