#!/usr/bin/perl -w

use v5.12;
no warnings qw(experimental);
use strict;

package Heap::Elem::VCFElemRev;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use Heap::Elem::VCFElem;

require Exporter;

@ISA = qw(Exporter Heap::Elem::VCFElem);

# No names exported.
@EXPORT = ( );

# Available for export: VCFElem (to allocate a new RDConnect::VCFElemRev value)
@EXPORT_OK = qw( VCFElemRev );

$VERSION = '0.01';

sub VCFElemRev {	# exportable synonym for new
	__PACKAGE__->new(@_);
}

# compare two VCF elems
sub cmp {
	return Heap::Elem::VCFElem::compareVCFlines($_[1][0], $_[0][0]);
}

1;
