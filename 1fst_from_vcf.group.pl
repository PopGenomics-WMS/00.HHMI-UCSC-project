#!/usr/bin/perl

open result, ">$ARGV[0].1genotype0.txt";


open sample, "<group.Horse.list";
while (<sample>) {
	chomp $_;
#	($population, $id)= split (/\s+/, $_);
        ($id,$population)= split (/\s+/, $_);
	push @ {$pop_sample_ind {$population}}, $id;
}

@pop_id = sort(keys%pop_sample_ind);

print result "chrom\tposition\t";

for ($i=0;$i<=@pop_id-1;$i++) {
	print result "$pop_id[$i]_0\t$pop_id[$i]_1\t$pop_id[$i]\t";
}

print result "\n";






open vcf, "<$ARGV[0]";
while (<vcf>) {
	my $snp_chr=0;
	chomp $_;
	if ($_ =~ /^\#\#/) {
		next;
	}

	if ($_=~ /^\#CHROM/) {
			@each_feature_id =split (/\t/, $_);
		    for ($j=0;$j<=@each_feature_id-1;$j++) {
				$feature_list {$each_feature_id[$j]} = $j;
			}

	}


	else {


	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	
		@each_snp_line=split (/\t/, $_);
		$snp_chr = $each_snp_line[0];
		$snp_position = $each_snp_line[1];
		print result "$snp_chr\t$snp_position\t";



        for ($i=0;$i<=@pop_id-1;$i++) {
			$allele1_count{0}=0;
			$allele1_count{1}=0;

			foreach $ind_name ( @ {$pop_sample_ind{$pop_id[$i]}} ) {
				#print "$feature_list{$ind_name}\n";

				$genotype_information = $each_snp_line [$feature_list{$ind_name}];
				$allele1= substr ($genotype_information, 0, 1);
				$allele2= substr ($genotype_information, 2, 1);
				$allele1_count {$allele1} ++;
				$allele1_count {$allele2} ++;
			}

			$pop_allele0_count = $allele1_count{0};   # 0 is ref, 1 is alter
			$pop_allele1_count = $allele1_count{1};
			$pop_genotyped_ind_number = $allele1_count{0}+ $allele1_count{1};
			print result "$pop_allele0_count\t$pop_allele1_count\t$pop_genotyped_ind_number\t";
		}
			print result "\n";

	}

}


