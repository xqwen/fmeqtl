
$radius = $ARGV[2];
$gene_map = $ARGV[0];
$snp_map = $ARGV[1];
$out_file = $ARGV[3];


my %rcd;
my @list;

if($gene_map =~/\.gz$/){
    open FILE, "zcat $gene_map | ";
}else{
    open FILE, "$gene_map";
}

while(<FILE>){
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    $rcd{$data[0]} = {tss => $data[-1]};
    push @list, $data[0];
}

@list = sort {$rcd{$a}->{tss} <=> $rcd{$b}->{tss}} keys %rcd;


if($snp_map =~/\.gz$/){
    open FILE, "zcat $snp_map | ";
}else{
    open FILE, "$snp_map";
}

open OUT, ">$out_file";

while(<FILE>){
    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    my $pos = $data[-1];
    foreach $g (@list){
	my $diff = $rcd{$g}->{tss}-$pos;
	last if($diff > $radius);	    
	if(abs($diff)<=$radius){
	    print OUT "$g $data[0]\n"; #$rcd{$g}->{tss} $pos $diff \n";
	}
    }
}



