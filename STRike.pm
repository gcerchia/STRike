=head1 LICENSE

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

     	http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.

=head1 CONTACT

	Giovanni Cerchia <gcerchia@engenome.com>

=cut

=head1 NAME

	STRike.1.0

=head1 SYNOPSIS

	mv STRike.pm ~/.vep/Plugins
	./vep -i variants.vcf --plugin STRike

=head1 DESCRIPTION

	A VEP plugin for the Ensembl Variant Effect Predictor (VEP) that returns
	HGVS notations for Short Tandem Repeats (STR).

=cut

package STRike;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
	return ['Transcript'];
}

sub get_header_info {
	return {
		HGVS_G_STR => "HGVS_g notation for STR",
		HGVS_C_N_STR => "HGVS_c (or n) notation for STR",
	};
}

# Returns the minimum repeated motif and the relative counts
sub get_minimum_motif_and_counts {
	my $sequence = shift();
	my $max_i = length($sequence);
	my @a = (1..$max_i);

	my $max_kmer;
	my $max_counts;

	# Checks for each k-mer starting from the first character of the string
	for my $i (@a) {
		my $kmer = substr($sequence, 0, $i);
		my @occurrences = $sequence =~ /$kmer/g;
		my $counts = @occurrences;
		my $total_nts = $counts*length($kmer);
		if($total_nts == length($sequence)) {
			$max_kmer = $kmer;
			$max_counts = $counts;
			last;
		}
	}
	my $min_repeat_unit = substr($sequence, 0, length($max_kmer));

	return ($min_repeat_unit, $max_counts);
}

# Returns the sequence of the transcript with the mutation, including 5' and 3' UTRs.
sub get_variation_cdna_seq{
    my $tva = $_[0];

    # Transcript sequence contains 5' and 3' regions.
    my $seq = $tva->transcript->seq->seq;
    if (!defined($tva->transcript_variation->cdna_start) || !defined($tva->transcript_variation->cdna_end)){
    	my $pos = $tva->transcript_variation->variation_feature->start;
	my $alleles = $tva->allele_string;
	my $transcript_id = $tva->transcript->display_id;
        print("[WARNING] cDNA coordinates not defined - ".$transcript_id." ".$pos." ".$alleles." (INTRON VARIANT)\n");
        return $seq;
    }

    # Variation position counting UTR regions.
    my $variation_start = $tva->transcript_variation->cdna_start - 1;
    my $variation_end = $tva->transcript_variation->cdna_end - 1;

    # If is a deletion, feature_seq is '-', so we will use '' instead to build the final sequence.
    my $feature_seq = $tva->feature_seq eq "-" ? "" : $tva->feature_seq;
    substr($seq, $variation_start, $variation_end - $variation_start + 1) = $feature_seq;

    return $seq;
}

# Returns the number of repeats of a specific motif
sub get_repeats_count_on_seq {
	my $sequence = shift();
	my $motif = shift();
	
	my $flag = 1;
	my $start = 0;
	my $counts = 0;
	while($flag == 1) {
		my $s = substr($sequence, $start, length($motif));
		if($s eq $motif) {
			$counts += 1;
			$start += length($motif);
		} else {
			$flag = 0;
		}
	}
	return $counts;
}

# Returns the number of repeats of a specific motif (starting from the end of the sequence)
sub get_repeats_count_on_seq_reverse {
	my $sequence = shift();
	my $motif = shift();

	my $flag = 1;
	my $start = length($sequence) - length($motif);
	my $counts = 0;
	while($flag == 1) {
		my $s = substr($sequence, $start, length($motif));
		if($s eq $motif) {
			$counts += 1;
			$start -= length($motif);
		} else {
			$flag = 0;
		}
	}
	return $counts;
}

# Returns the fixed start position after the application of the 3' rule and the new motif
sub get_most_3_prime_pos() {
	my $sequence = shift();
	my $motif = shift();
	my $hgvs_c_start = shift();

	my $offset = 0;
	if(!(length($motif) % 3)) {		
		my $flag = 1;
		my $start_pos = 0;
		my $rep = 0;
		
		# Checks repeats count for the given motif
		while($flag) {
			my $tmp = substr($sequence, $start_pos, length($motif));
			if ($tmp eq $motif) {
				$rep += 1;
				$start_pos += length($motif);
			} else {
				$flag = 0;
			}
		}
		my $end_pos = $rep*length($motif);
		my $after_variation_seq = substr($sequence, $end_pos, length($motif));

		my @a = (1..length($motif)-1);
		my $sub_motif = "";

		# Checks if the sequence after the variation site partially contains the motif and calculates the offset
		for my $i (@a) {
			my $tmp = substr($after_variation_seq, 0, $i);
			my $tmp_motif = substr($motif, 0, $i);
			if($tmp eq $tmp_motif and length($tmp)>length($sub_motif)) {
				$sub_motif = $tmp;
				$offset = length($tmp);
			}
		}

		if($offset != 0) {
			my $new_motif = substr($motif, $offset)."".substr($motif, 0, $offset);
			return ($new_motif, $hgvs_c_start + $offset);
		}
	}
	return ($motif, $hgvs_c_start);
}

sub run {
	my ($self, $transcript_variation_allele) = @_;

	# Gets the initial HGVS transcript notation
	my $c = $transcript_variation_allele->hgvs_transcript;
	my @hgvs_c_arr = split(':', $c);
	my $reference_seq = $hgvs_c_arr[0];
	my $orig_hgvs_c = $hgvs_c_arr[1];
	my $hgvs_c = $c;
	my $final_repeats = 0;

	if(!defined $hgvs_c) {
		my $pos = $transcript_variation_allele->transcript_variation->variation_feature->start;
		my $alleles = $transcript_variation_allele->allele_string;
		print("[WARNING] HGVS transcript not defined (".$pos." ".$alleles.")\n");
	}
	elsif(index($orig_hgvs_c, 'ins') != -1 || index($orig_hgvs_c, 'dup') != -1 || index($orig_hgvs_c, 'del') != -1) {
		(my $hgvs_c_pos, my $hgvs_c_seq) = split(m[(?:ins|dup|del)+], $orig_hgvs_c);
		if(!defined $hgvs_c_seq) {
			$hgvs_c_seq = "";
		}

		my $hgvs_c_type = substr($orig_hgvs_c, 0, 1);
		my $hgvs_c_var = "";
		if(index($orig_hgvs_c, 'ins') != -1) {
			my $idx = index($orig_hgvs_c, 'ins');
			$hgvs_c_var = substr($orig_hgvs_c, $idx, 3);
		} elsif(index($orig_hgvs_c, 'dup') != -1) {
			my $idx = index($orig_hgvs_c, 'dup');
			$hgvs_c_var = substr($orig_hgvs_c, $idx, 3);
		} elsif(index($orig_hgvs_c, 'del') != -1) {
			my $idx = index($orig_hgvs_c, 'del');
			$hgvs_c_var = substr($orig_hgvs_c, $idx, 3);
		}

		my $hgvs_c_motif = "";
		my $hgvs_c_counts = 0;
		my $hgvs_c_start = 0;
		my $hgvs_c_end = 0;
		my $total_repeats = 0;
		my $seq_to_check = "";

		if(index($orig_hgvs_c, 'del') != -1) {
			($hgvs_c_motif, $hgvs_c_counts) = &get_minimum_motif_and_counts((split("/", $transcript_variation_allele->transcript_variation->cdna_allele_string))[0]);
		} else {
			($hgvs_c_motif, $hgvs_c_counts) = &get_minimum_motif_and_counts($transcript_variation_allele->feature_seq);
		}

		my $alt_cdna_seq = &get_variation_cdna_seq($transcript_variation_allele);
		my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($transcript_variation_allele->transcript);
		if(index($hgvs_c_pos, "-") != -1 && !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
			# Case 5' UTR variants
			my $sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_alt_cdna_seq = substr($alt_cdna_seq, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_alt_cdna_seq2 = substr($alt_cdna_seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			$seq_to_check = $transcript_variation_allele->transcript_variation->transcript->five_prime_utr->seq;
			my $g_start = $transcript_variation_allele->transcript_variation->variation_feature->start;
			my $g_end = $transcript_variation_allele->transcript_variation->variation_feature->end;
			my $g_strand = $transcript_variation_allele->transcript_variation->variation_feature->strand;
			my @cdna_coords = $trmapper->genomic2cdna($g_start, $g_end, $g_strand);
			$hgvs_c_start = $cdna_coords[0]->start;
			# Applies 3' rule
			($hgvs_c_motif, $hgvs_c_start) = &get_most_3_prime_pos($sub_alt_cdna_seq, $hgvs_c_motif, $hgvs_c_start);
			$sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $hgvs_c_start - 1);
			$sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $hgvs_c_start - 1);
			my $ref_repeats = &get_repeats_count_on_seq_reverse($sub_cdna_seq2, $hgvs_c_motif) + &get_repeats_count_on_seq($sub_cdna_seq, $hgvs_c_motif);

			if($transcript_variation_allele->transcript->strand == -1) {
				if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
					$total_repeats =  $ref_repeats - $hgvs_c_counts;
					$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$total_repeats - length($seq_to_check) - 1;
				} else {
					$total_repeats = $hgvs_c_counts;
					$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$ref_repeats - length($seq_to_check) - 1;
				}
			} else {
				if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
					$total_repeats =  $ref_repeats - $hgvs_c_counts;
				} else {
					$total_repeats = $hgvs_c_counts + $ref_repeats;
				}
				$hgvs_c_start = $hgvs_c_start - length($seq_to_check) - 1;
			}
			$hgvs_c_end = $hgvs_c_start + length($hgvs_c_motif)*$total_repeats - 1;

			if($hgvs_c_end >= 0) {
				$hgvs_c_end = "";
			}
		} elsif(index($hgvs_c_pos, "*") != -1 && !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
			# Case 3' UTR variants
			my $sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_alt_cdna_seq = substr($alt_cdna_seq, $transcript_variation_allele->transcript_variation->cdna_start);
			my $sub_alt_cdna_seq2 = substr($alt_cdna_seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			$seq_to_check = $transcript_variation_allele->transcript_variation->transcript->three_prime_utr->seq;
			my $g_start = $transcript_variation_allele->transcript_variation->variation_feature->start;
			my $g_end = $transcript_variation_allele->transcript_variation->variation_feature->end;
			my $g_strand = $transcript_variation_allele->transcript_variation->variation_feature->strand;
			my @cdna_coords = $trmapper->genomic2cdna($g_start, $g_end, $g_strand);
			my $three_utr_start = index($transcript_variation_allele->transcript->seq->seq, $seq_to_check) + 1;
			my $ref_repeats = &get_repeats_count_on_seq_reverse($sub_cdna_seq2, $hgvs_c_motif) + &get_repeats_count_on_seq($sub_cdna_seq, $hgvs_c_motif);
			$total_repeats = $hgvs_c_counts + $ref_repeats;
			$hgvs_c_start = $cdna_coords[0]->start;
			# Applies 3' rule
			(my $hgvs_c_motif_fixed, my $hgvs_c_start_fixed) = &get_most_3_prime_pos($sub_alt_cdna_seq, $hgvs_c_motif, $hgvs_c_start);
			if(!$hgvs_c_motif_fixed eq $hgvs_c_motif) {
				$sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $hgvs_c_start - 1);
				$sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $hgvs_c_start - 1);
				my $ref_repeats = &get_repeats_count_on_seq_reverse($sub_cdna_seq2, $hgvs_c_motif) + &get_repeats_count_on_seq($sub_cdna_seq, $hgvs_c_motif);
				$hgvs_c_start = $hgvs_c_start_fixed;
			}

			if($transcript_variation_allele->transcript->strand == -1) {
				if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
					$total_repeats =  $ref_repeats - $hgvs_c_counts;
					$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$total_repeats - $three_utr_start + 1;
					$hgvs_c_end = $hgvs_c_start + length($hgvs_c_motif)*$total_repeats - 1;
					$hgvs_c_start = "*".$hgvs_c_start;
					$hgvs_c_end = "*".$hgvs_c_end;
				} else {
					$total_repeats =  $ref_repeats + $hgvs_c_counts;
					$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$ref_repeats - $three_utr_start + 1;
					$hgvs_c_end = $hgvs_c_start + length($hgvs_c_motif)*$total_repeats - 1;
					$hgvs_c_start = "*".$hgvs_c_start;
					$hgvs_c_end = "*".$hgvs_c_end;
				}
			} else {
				if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
					$total_repeats =  $ref_repeats - $hgvs_c_counts;
				} else {
					$total_repeats =  $ref_repeats + $hgvs_c_counts;
				}
				$hgvs_c_start = "*".($hgvs_c_start - $three_utr_start + 1);
				$hgvs_c_end = "*".($hgvs_c_start + length($hgvs_c_motif)*$total_repeats - 1);
			}
		} elsif(!(defined $transcript_variation_allele->transcript_variation->intron_number)) {
			# Case exon variants
			my $sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_alt_cdna_seq = substr($alt_cdna_seq, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $sub_alt_cdna_seq2 = substr($alt_cdna_seq, 0, $transcript_variation_allele->transcript_variation->cdna_start - 1);
			my $g_start = $transcript_variation_allele->transcript_variation->variation_feature->start;
			my $g_end = $transcript_variation_allele->transcript_variation->variation_feature->end;
			my $g_strand = $transcript_variation_allele->transcript_variation->variation_feature->strand;
			my @cdna_coords = $trmapper->genomic2cdna($g_start, $g_end, $g_strand);
			my $five_prime_utr = $transcript_variation_allele->transcript_variation->transcript->five_prime_utr;
			my $ref_repeats = &get_repeats_count_on_seq_reverse($sub_cdna_seq2, $hgvs_c_motif) + &get_repeats_count_on_seq($sub_cdna_seq, $hgvs_c_motif);
			$total_repeats = $hgvs_c_counts + $ref_repeats;
			$hgvs_c_start = $cdna_coords[0]->start;
			$hgvs_c_end = $cdna_coords[0]->start + length($hgvs_c_motif)*$total_repeats - 1;
			if(defined $five_prime_utr) {
				# Applies 3' rule
				(my $hgvs_c_motif_fixed, my $hgvs_c_start_fixed) = &get_most_3_prime_pos($sub_alt_cdna_seq, $hgvs_c_motif, $hgvs_c_start);
				if(!($hgvs_c_motif_fixed eq $hgvs_c_motif)) {
					$sub_cdna_seq = substr($transcript_variation_allele->transcript->seq->seq, $hgvs_c_start - 1);
					$sub_cdna_seq2 = substr($transcript_variation_allele->transcript->seq->seq, 0, $hgvs_c_start - 1);
					my $sub_alt_cdna_seq = substr($alt_cdna_seq, $hgvs_c_start);
					my $sub_alt_cdna_seq2 = substr($alt_cdna_seq, 0, $hgvs_c_start - 1);
					my $ref_repeats = &get_repeats_count_on_seq_reverse($sub_cdna_seq2, $hgvs_c_motif) + &get_repeats_count_on_seq($sub_cdna_seq, $hgvs_c_motif);
					$total_repeats = $hgvs_c_counts + $ref_repeats;
					$hgvs_c_start = $hgvs_c_start_fixed;
					$hgvs_c_motif = $hgvs_c_motif_fixed;
				}
				
				if($transcript_variation_allele->transcript->strand == -1) {
					if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
						$total_repeats =  $ref_repeats - $hgvs_c_counts;
						$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$total_repeats - length($five_prime_utr->seq);
					} else {
						$hgvs_c_start = $hgvs_c_start - length($hgvs_c_motif)*$ref_repeats - length($five_prime_utr->seq);
					}
				} else {
					if($transcript_variation_allele->variation_feature->allele_string =~ /-$/) {
						$total_repeats =  $ref_repeats - $hgvs_c_counts;
					}
					$hgvs_c_start = $hgvs_c_start - length($five_prime_utr->seq);
				}
				$hgvs_c_end = $hgvs_c_start + length($hgvs_c_motif)*$total_repeats - 1;
				
				if($total_repeats == 1) {
					if(index($orig_hgvs_c, 'del') != -1) {
						$hgvs_c_end = $hgvs_c_start + length($hgvs_c_motif) - 1;
					} else {
						$hgvs_c_start -= 1;
						$hgvs_c_end = $hgvs_c_start + 1;
					}
				}
			} else {
				print("[WARNING] 5' UTR sequence not defined for ".$transcript_variation_allele->transcript->display_id."\n");
			}
		} else {
			# Case intron variants
			my @all_introns = $transcript_variation_allele->transcript_variation->transcript->get_all_Introns;
			my $intron_number = (split("/",$transcript_variation_allele->transcript_variation->intron_number))[0];
			my $intron_idx = $intron_number - 1;
			my $intron_seq = $all_introns[0][$intron_idx]->seq;
			my $g_start = $transcript_variation_allele->transcript_variation->variation_feature->start;
			my $near_exon;
			my $offset_start;
			my $offset_end;
			my $intron_strand = $all_introns[0][$intron_idx]->seq_region_strand;
			$total_repeats = $hgvs_c_counts;
			my $ref_repeats_len;
			# Applies 3' rule
			if(index($orig_hgvs_c, "+") != -1) {
				if($intron_strand == 1) {
					$offset_start = $g_start - $all_introns[0][$intron_idx]->prev_Exon->end;
					my $sub_intron_seq = substr($intron_seq, $g_start - $all_introns[0][$intron_idx]->prev_Exon->end - 1);
					$ref_repeats_len = &get_repeats_count_on_seq($sub_intron_seq, $hgvs_c_motif)*length($hgvs_c_motif);
					my @tmp = split("_", $orig_hgvs_c);
					$near_exon = substr($tmp[1], 0, index($tmp[1], '+'));
				} elsif($intron_strand == -1) {
					my $sub_intron_seq = substr($intron_seq, 0, $all_introns[0][$intron_idx]->prev_Exon->start - $g_start);
					$ref_repeats_len = &get_repeats_count_on_seq_reverse($sub_intron_seq, $hgvs_c_motif)*length($hgvs_c_motif);
					$offset_start = $all_introns[0][$intron_idx]->prev_Exon->start - $g_start - $ref_repeats_len + 1;
					my @tmp = split("_", $orig_hgvs_c);
					$near_exon = substr($tmp[1], 0, index($tmp[1], '+'));
				}
				$total_repeats += ($ref_repeats_len/length($hgvs_c_motif));
				$offset_end = $offset_start + length($hgvs_c_motif)*$total_repeats - 1;
				if( !(length($hgvs_c_motif) % 3)) {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$near_exon."+".$offset_start."_".$near_exon."+".$offset_end.$hgvs_c_motif."[".$total_repeats."]";
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$near_exon."+".$offset_start."_".$near_exon."+".$offset_end.$hgvs_c_var.$hgvs_c_seq;
				}
			} elsif(index($orig_hgvs_c, "-") != -1) {
				if($intron_strand == 1) {
					$offset_start = abs($g_start - $all_introns[0][$intron_idx]->next_Exon->start);
					my $sub_intron_seq = substr($intron_seq, $g_start - $all_introns[0][$intron_idx]->prev_Exon->end - 1);
					$ref_repeats_len = &get_repeats_count_on_seq($sub_intron_seq, $hgvs_c_motif)*length($hgvs_c_motif);
					my @tmp = split("_", $orig_hgvs_c);
					if(substr($tmp[1], 0, 1) eq "-") {
						$near_exon = substr($tmp[1], 0, index($tmp[1], '-'));
					} else {
						$near_exon = substr($tmp[1], 0, index($tmp[1], '-'));
					}
				} elsif($intron_strand == -1) {
					my $sub_intron_seq = substr($intron_seq, 0, $all_introns[0][$intron_idx]->prev_Exon->start - $g_start);
					$ref_repeats_len = &get_repeats_count_on_seq_reverse($sub_intron_seq, $hgvs_c_motif)*length($hgvs_c_motif);
					$offset_start = abs($g_start - $all_introns[0][$intron_idx]->next_Exon->end);
					my @tmp = split("_", $orig_hgvs_c);
					if(substr($tmp[1], 0, 1) eq "-") {
						$near_exon = substr($tmp[1], 0, rindex($tmp[1], '-'));
					} else {
						$near_exon = substr($tmp[1], 0, index($tmp[1], '-'));
					}
				}
				$total_repeats += ($ref_repeats_len/length($hgvs_c_motif));
				$offset_end = abs($offset_start - length($hgvs_c_motif)*$total_repeats - 1);
				if(!(length($hgvs_c_motif) % 3)) {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$near_exon."-".$offset_start."_".$near_exon."-".$offset_end.$hgvs_c_motif."[".$total_repeats."]";
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$near_exon."-".$offset_start."_".$near_exon."-".$offset_end.$hgvs_c_var.$hgvs_c_seq;
				}
			} else {
				print("[ERROR] Wrong hgvs_c format: ".$reference_seq.":".$orig_hgvs_c." (".$g_start." ".$hgvs_c_motif.")\n");
			}
		}

		# Final HGVS notation formatting
		if(index($orig_hgvs_c, 'del') != -1) {
			if(!(length($hgvs_c_motif) % 3) and $total_repeats >= 1 and !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
				if($hgvs_c_end eq "") {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start.$hgvs_c_motif."[".$total_repeats."]";
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start."_".$hgvs_c_end.$hgvs_c_motif."[".$total_repeats."]";
				}
			} elsif(!(length($hgvs_c_motif) % 3) and !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
				if($hgvs_c_end eq "") {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start.$hgvs_c_var.$hgvs_c_seq;
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start."_".$hgvs_c_end.$hgvs_c_var.$hgvs_c_seq;
				}
			}
		} else {
			if(!(length($hgvs_c_motif) % 3) and $total_repeats > 1 and !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
				if($hgvs_c_end eq "") {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start.$hgvs_c_motif."[".$total_repeats."]";
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start."_".$hgvs_c_end.$hgvs_c_motif."[".$total_repeats."]";
				}
			} elsif(!(length($hgvs_c_motif) % 3) and !(defined $transcript_variation_allele->transcript_variation->intron_number)) {
				if($hgvs_c_end eq "") {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start.$hgvs_c_var.$hgvs_c_seq;
				} else {
					$hgvs_c = $reference_seq.":".$hgvs_c_type.".".$hgvs_c_start."_".$hgvs_c_end.$hgvs_c_var.$hgvs_c_seq;
				}
			}
		}

		$final_repeats = $total_repeats;
	}

	# Generates genomic HGVS notation
	my $hgvs_g = "";
	my @alleles = split('/', $transcript_variation_allele->allele_string);
	if(!defined $hgvs_c and length($alleles[0]) <= length(length($alleles[1]))) {
		(my $motif, my $counts) = &get_minimum_motif_and_counts($alleles[1]);
		my $start = $transcript_variation_allele->base_variation_feature->seq_region_start;
		my $end = $start + ($counts*length($motif)) - 1;
		$hgvs_g = 'g.'.$start.'_'.$end.$motif.'['.$counts.']';
	} elsif ((!defined $hgvs_c and length($alleles[0]) > length(length($alleles[1]))) or index($orig_hgvs_c, 'del') != -1) {
		(my $motif, my $counts) = &get_minimum_motif_and_counts($alleles[0]);
		my $start = $transcript_variation_allele->base_variation_feature->seq_region_start;
		my $end = $start + ($counts*length($motif)) - 1;
		$hgvs_g = 'g.'.$start.'_'.$end.'del'.$motif.'['.$counts.']';
	} else {
		(my $motif, my $counts) = &get_minimum_motif_and_counts($alleles[1]);
		my $start = $transcript_variation_allele->base_variation_feature->seq_region_start;
		my $end = $start + ($final_repeats*length($motif)) - 1;
		$hgvs_g = 'g.'.$start.'_'.$end.$motif.'['.$final_repeats.']';
	}

	return {
		HGVS_G_STR => $hgvs_g,
		HGVS_C_N_STR => $hgvs_c,
	};

}

1;
