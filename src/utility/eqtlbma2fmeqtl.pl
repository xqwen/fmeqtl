#!/usr/bin/perl
# usage: eqtlbma2fmeqtl --geno list_genotypes.txt --scoord snp_coords.bed.gz --exp list_explevels.txt --gcoord gene_coords.bed.gz  --cis 1000 --outdir eqtl_dat

parse_args();
cis_snp();
assemble();
## parse parameters

sub parse_args{


    for $i (0..$#ARGV){
	
	if($ARGV[$i] eq "--help"){
	    print STDERR "usage: perl eqtlbma2fmeqtl.pl --geno list_genotypes.txt --scoord snp_coords.bed.gz --exp list_explevels.txt --gcoord gene_coords.bed.gz  --cis 1000 --outdir eqtl_dat\n";
	    exit(0);
	}


	if($ARGV[$i] eq "--geno"){
	    $geno_info = $ARGV[++$i];
	    next;
	}

	if($ARGV[$i] eq "--exp"){
	    $exp_info = $ARGV[++$i];
	    next;
	}
	
	if($ARGV[$i] eq "--scoord"){
	    $snp_map = $ARGV[++$i];
	    next;
	}
	
	if($ARGV[$i] eq "--gcoord"){
	    $gene_map = $ARGV[++$i];
	    next;
	}

	if($ARGV[$i] eq "--cis"){
	    $radius = 1000*$ARGV[++$i];
	    next;
	}

	if($ARGV[$i] eq "--outdir"){
	    $outdir = $ARGV[++$i];
	    next;
	}
	

	#print STDERR "Error: unknown command line option \"$ARGV[0]\"\n";
	#print STDERR "usage: perl eqtlbma2fmeqtl.pl --geno list_genotypes.txt --scoord snp_coords.bed.gz --exp list_explevels.txt --gcoord gene_coords.bed.gz  --cis 1000 --outdir eqtl_dat\n";
	#exit(0);


    }
    
}




sub cis_snp {
    
    
    if($gene_map =~/\.gz$/){
	open FILE, "zcat $gene_map | ";
    }else{
	open FILE, "$gene_map";
    }
    
    while(<FILE>){
	chomp;
	next if $_ !~ /\d/;

	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
	$rcd{$data[0]}->{$data[3]} = $data[1];
    }

    if($snp_map =~/\.gz$/){
	open FILE, "zcat $snp_map | ";
    }else{
	open FILE, "$snp_map";
    }

    foreach $chr (keys %rcd){
	my @list = sort {$rcd{$chr}->{$a} <=> $rcd{$chr}->{$b}} keys %{$rcd{$chr}};
	$glist{$chr} = \@list;
    }
	


    while(<FILE>){
	next if $_ !~ /\d/;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;

	my $chr = $data[0];
	my $pos = $data[1];
	my $snp = $data[3];

	foreach $g (@{$glist{$chr}}){
	    my $diff = $rcd{$chr}->{$g} - $pos;
	    last if($diff > $radius);
	    if(abs($diff)<=$radius){
		$gmap{$snp}->{$g} = 1;
		#print "$chr $g $snp      $pos  $rcd{$chr}->{$g}  $diff\n";
	    }
	}
    }
}



sub assemble {


    open FILE, "$geno_info";
     
    my %geno_rcd;
    while(<FILE>){
	next if $_!~ /\s*(\S+)\s+(\S+)\s*/;
	$geno_data{$1} = $2;
	push @grps, $1;
	$geno_rcd{$2}++;
    }
    if(scalar(keys %geno_rcd)>1){
	print STDERR "Data format SSLR\n";
	assemble_sslr();
	
    }else{
	print STDERR "Data format MVLR\n";
	assemble_mvlr();
    }
    
    print STDERR "Assembling ...\n";

    `mkdir -p $outdir`;
    foreach $chr (keys %rcd){
	foreach $g (@{$glist{$chr}}){
	    open OUT, ">$outdir/$g.dat";
	    print OUT "$con{$g}";
	}
    }


}



sub assemble_sslr{ 

    open FILE, "$exp_info";
    while(<FILE>){
	ext if $_!~ /\s*(\S+)\s+(\S+)\s*/;
	$exp_data{$1} = $2;
    }


    my %index;

    foreach $s (@grps){
	
	if($exp_data{$s} =~/\.gz$/){
	    open FILE, "zcat $exp_data{$s} | ";
	}else{
	    open FILE, "$exp_data{$s}";
	}

	$start = 0;
	while(<FILE>){
	    
	    next if $_ !~ /\S/;
	    chomp;
	    my @data = split /\s+/, $_;
	    shift @data until $data[0]=~/^\S/;
	    
	    if($start==0){
		$index{$s} = \@data;
		$start = 1;
		next;
	    }

	    my $g = $data[0];
	    shift @data;
	    $con{$g} .= sprintf "pheno $g $s @data\n";
	}
	    
	
    }


    foreach $s (@grps){

        if($geno_data{$s} =~/\.gz$/){
            open FILE, "zcat $geno_data{$s} | ";
        }else{
	    open FILE, "$geno_data{$s}";
        }

        $start = 0;
        my @ind_index;
        while(<FILE>){

            next if $_ !~ /\S/;
            chomp;
            my @data = split /\s+/, $_;
            shift @data until $data[0]=~/^\S/;
	    if($start==0){
		my %hash;
	        @hash{@data} = (0..$#data);
		@ind_index = @hash{@{$index{$s}}};
		$start = 1;
		next;
	    }

            my $snp = $data[0];
            shift @data;
	    foreach $g (keys %{$gmap{$snp}}){
		$con{$g} .= sprintf "geno $snp $s @data[@ind_index]\n";
	    }
	}


    }

}




sub assemble_mvlr{ 

    open FILE, "$exp_info";
    while(<FILE>){
	ext if $_!~ /\s*(\S+)\s+(\S+)\s*/;
	$exp_data{$1} = $2;
    }

    
    my @index_v;


    foreach $s (@grps){
	
	if($exp_data{$s} =~/\.gz$/){
	    open FILE, "zcat $exp_data{$s} | ";
	}else{
	    open FILE, "$exp_data{$s}";
	}

	$start = 0;
	while(<FILE>){
	    
	    next if $_ !~ /\S/;
	    chomp;
	    my @data = split /\s+/, $_;
	    shift @data until $data[0]=~/^\S/;
	    
	    if($start==0){
		if(scalar(@index_v)==0){
		    @index_v = @data;
		}
		$start = 1;
		next;
	    }

	    my $g = $data[0];
	    shift @data;
	    $con{$g} .= sprintf "pheno $s @data\n";
	}
	    
	
    }


    $s = $grps[0];

    if($geno_data{$s} =~/\.gz$/){
	open FILE, "zcat $geno_data{$s} | ";
    }else{
	open FILE, "$geno_data{$s}";
    }

    my @ind_index;
    $start = 0;

    while(<FILE>){
	next if $_ !~ /\S/;
	chomp;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
	if($start==0){
            my %hash;
	    @hash{@data} = (0..$#data);
	    @ind_index = @hash{@index_v};
	    $start = 1;
	    next;
	}

	my $snp = $data[0];
	shift @data;
	foreach $g (keys %{$gmap{$snp}}){
	    $con{$g} .= sprintf "geno $snp @data[@ind_index]\n";
	}
    }

}
