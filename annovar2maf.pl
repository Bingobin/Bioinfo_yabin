#!/usr/bin/perl

use warnings;
use strict;

my $file = shift or die $!;
my $samid = shift or die $!;

my $title = "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tTumor_Sample_Barcode\tMutation_Status\tAAChange\tTranscript_Id\tTxChange\tExon_Number\tDetailGene\tRMSK\tGeneHancer\tCOSMIC\tt_depth\tt_ref_count\tt_alt_count\tTumorVAF\tnum_soft\tsoft_name";
print "$title\n";

my $entrez_gene_id = "NA";
#my $center = "WL1109";
my $center = "TCGA";
my $ncbi_build = "hg38";
my $mutation_status = "Somatic";
my $strand = "+";

my %variant = (
		"frameshift deletion" => "Frame_Shift_Del",
		"frameshift insertion" => "Frame_Shift_Ins",
		"nonframeshift deletion" => "In_Frame_Del",
		"nonframeshift insertion" => "In_Frame_Ins",
		"nonsynonymous SNV" => "Missense_Mutation",
		"stopgain" => "Nonsense_Mutation",
		"stoploss" => "Nonstop_Mutation",
		"synonymous SNV" => "Silent",
		"unknown" => "UNKNOWN",
		"UNKNOWN" => "UNKNOWN",
		"NA" => "UNKNOWN",
		"downstream" => "3'Flank",
		"upstream" => "5'Flank",
		"splicing" => "Splice_Site",
		"UTR3" => "3'UTR",
		"UTR5" => "5'UTR",
		"intergenic" => "IGR",
		"intronic" => "Intron",
		"ncRNA_exonic" => "RNA",
		"ncRNA_intronic" => "RNA",
		"ncRNA_splicing" => "RNA",
		"ncRNA_UTR5" => "RNA",
		"ncRNA_UTR3" => "RNA",
		"ncRNA" => "RNA"
);

open IN, "$file" or die $!;
while(<IN>){
	chomp;
	next if (/^Chr/);
	my @tmp = split /\t/;
	for (my $i=0;$i<@tmp;$i++){
		if($tmp[$i] eq "."){
			$tmp[$i] = "NA";
		}
	}
	my $hugo_symbol = $tmp[6];
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $end = $tmp[2];
	my $ref_allele = $tmp[3];
	my $var_allele1 = $tmp[3];
	my $var_allele2 = $tmp[4];
	my $snp_rs = $tmp[13];
	my $genedetail = $tmp[7];
	my $rmsk = $tmp[10];
	my $genehancer = $tmp[11];
	my $cosmic = $tmp[12];
	my $ref_r = $tmp[14];
	my $alt_r = $tmp[15];
	my $dep = $ref_r + $alt_r;
	my $t_vaf = $tmp[16];
	my $soft_num = $tmp[17];
	my $soft_nam = $tmp[18];
	my $variant_class;
	my ($aa_change, $tr_id, $tx_change, $exon_num) = ("NA", "NA", "NA", "NA");
	my @func  = split(/;/, $tmp[5]);
	if(@func > 1){
		my @gene = split(/;/, $tmp[6]);
		$hugo_symbol = $gene[0];
	}
	if ($func[0] eq "exonic"){
		my @exonic = split(/;/, $tmp[8]);
		$variant_class = $variant{$exonic[0]};
		$genedetail = $tmp[9];
		$tmp[9] = ":NA:NA:NA:NA" if ($tmp[9] eq "UNKNOWN");
		my @change = split(/:/, $tmp[9]);
		if ($tmp[9] =~ /whole/){
			push @change, "NA";
			push @change, "NA";
		}
		$aa_change = (split(/,/, $change[4]))[0];
		$tr_id = $change[1];
		$tx_change = $change[3];
		$exon_num = $change[2];
	}else{
		$variant_class = $variant{$func[0]};
	}
	my $variant_type;
	if ($tmp[3] eq "-"){
		$variant_type = "INS";
	}elsif($tmp[4] eq "-"){
		$variant_type = "DEL";
	}else{
		$variant_type = "SNP";
	}
	print "$hugo_symbol\t$entrez_gene_id\t$center\t$ncbi_build\t$chr\t$start\t$end\t$strand\t$variant_class\t$variant_type\t$ref_allele\t$var_allele1\t$var_allele2\t$snp_rs\t$samid\t$mutation_status\t$aa_change\t$tr_id\t$tx_change\t$exon_num\t$genedetail\t$rmsk\t$genehancer\t$cosmic\t$dep\t$ref_r\t$alt_r\t$t_vaf\t$soft_num\t$soft_nam\n";
}
close IN;
