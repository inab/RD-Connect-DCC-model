#!/usr/bin/perl -w

use v5.12;
no warnings qw(experimental);
use strict;

package Heap::Elem::VCFElem;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use Heap::Elem;

require Exporter;

@ISA = qw(Exporter Heap::Elem);

# No names exported.
@EXPORT = ( );

# Available for export: VCFElem (to allocate a new RDConnect::VCFElem value)
@EXPORT_OK = qw( compareVCFlines compareVCFcols genPriorityString VCFElem);

$VERSION = '0.01';

use constant {
	VCF_SORTED_FILE	=>	0,
	VCF_SAMPLES	=>	1,
	VCF_HANDLER	=>	2,
	VCF_LINE	=>	3,
	VCF_LINE_KEYS	=>	4,
	VCF_METADATA	=>	5,
	VCF_FILE	=>	6,
	VCF_JOIN_POS	=>	7,
};

use constant {
	VCF_CHROM	=>	0,
	VCF_POS		=>	1,
	VCF_REF		=>	2,
	VCF_ALT		=>	3,
};

sub compareVCFlines($$);
sub compareVCFcols($$);
sub genPriorityString($);

sub compareVCFlines($$) {
	my($nleft,$nright)=@_;
	
	return compareVCFcols($nleft->[VCF_LINE_KEYS],$nright->[VCF_LINE_KEYS]);
}

sub compareVCFcols($$) {
	my($left,$right)=@_;
	
	if($left ~~ $right) {
		return 0;
	} elsif($left->[VCF_CHROM] lt $right->[VCF_CHROM] || (
			$left->[VCF_CHROM] eq $right->[VCF_CHROM] && (
				$left->[VCF_POS] < $right->[VCF_POS] || (
					$left->[VCF_POS] == $right->[VCF_POS] && (
						$left->[VCF_REF] lt $right->[VCF_REF] || (
							$left->[VCF_REF] eq $right->[VCF_REF] && $left->[VCF_ALT] lt $right->[VCF_ALT]
						)
					)
				)
			)
		)
	) {
		# It's over!!!!!
		return -1;
	} else {
		return 1;
	}
}

sub genPriorityString($) {
	my $vcfObj = shift;
	
	my $line = $vcfObj->[VCF_LINE_KEYS];
	
	return join(',',$line->[VCF_CHROM],sprintf('%.8x',$line->[VCF_POS]),$line->[VCF_REF],$line->[VCF_ALT]);
}

sub VCFElem {	# exportable synonym for new
	__PACKAGE__->new(@_);
}

sub priority {
	return genPriorityString($_[0]->val);
}

# compare two VCF elems
sub cmp {
	return compareVCFlines($_[0][0], $_[1][0]);
}

1;
