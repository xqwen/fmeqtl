open FILE, "fusion.gene.list";
@list = <FILE>;
chomp @list;
@rcd{@list} = @list;


open FILE, "gencode.v17.tss.tab";
$chr = "chr1";
open OUT, ">chr.gene.list";
while(<FILE>){
    next if $_ !~ /(ENSG\d+)/;
    next if ! defined($rcd{$1});
    $g = $1;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    if($data[1] ne $chr){
	close OUT;	
	`sort -nk3 chr.gene.list > $chr.gene.list`; 
	$chr = $data[1];	
	open OUT, ">chr.gene.list";
    }
    
    print OUT "$g $chr:$data[2]  $data[2]\n";
}

`sort -nk3 chr.gene.list > $chr.gene.list`;
