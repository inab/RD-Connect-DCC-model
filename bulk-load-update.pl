#!/usr/bin/perl -w

use strict;

use diagnostics;
use FindBin;
use lib $FindBin::Bin."/schema+tools/lib";

use Carp;
use Config::IniFiles;
use IO::File;
use File::Basename;
use File::Spec;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Elasticsearch;


use Search::Elasticsearch 1.12;
use JSON;

use constant {
	UPDATESIZE	=>	1000,
};

sub compareVCFlines($$);
sub compareVCFcols($$);

sub compareVCFlines($$) {
	my($nleft,$nright)=@_;
	
	return compareVCFcols($nleft->[4],$nright->[4]);
}

sub compareVCFcols($$) {
	my($left,$right)=@_;
	
	if($left ~~ $right) {
		return 0;
	} elsif($left->[0] lt $right->[0] || ($left->[0] eq $right->[0] && ($left->[1] < $right->[1] || ($left->[1] == $right->[1] && ($left->[2] lt $right->[2] || ($left->[2] eq $right->[2] && $left->[3] lt $right->[3])))))) {
		# It's over!!!!!
		return -1;
	} else {
		return 1;
	}
}

my %IDENTfields = (
	'RS'	=>	'RScode',
	'dbSNPBuildID'	=>	'dbSNPBuildID'
);

my $j = JSON->new->pretty;

if(scalar(@ARGV)>=2) {
	my $iniFile = shift(@ARGV);
	my $workingDir = shift(@ARGV);
	
	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
		
	# And create the working directory
	File::Path::make_path($workingDir);
	
	# Let's parse the model
	my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
	# Setting up the right path on relative cases
	$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

	print "Parsing model $modelFile...\n";
	my $model = undef;
	eval {
		$model = BP::Model->new($modelFile);
	};
	
	if($@) {
		Carp::croak('ERROR: Model parsing and validation failed. Reason: '.$@);
	}
	print "\tDONE!\n";
	
	my %storageModels = ();
	
	# Is there any file whose data has to be mapped?
	if(scalar(@ARGV)>0) {
		# Setting up the loader storage model(s)
		Carp::croak('ERROR: undefined destination storage model')  unless($ini->exists($BP::Loader::Mapper::SECTION,'loaders'));
		my $loadModelNames = $ini->val($BP::Loader::Mapper::SECTION,'loaders');
		
		my @loadModels = ();
		foreach my $loadModelName (split(/,/,$loadModelNames)) {
			unless(exists($storageModels{$loadModelName})) {
				$storageModels{$loadModelName} = BP::Loader::Mapper->newInstance($loadModelName,$model,$ini);
				push(@loadModels,$loadModelName);
			}
		}
		# Now, do we need to push the metadata there?
		if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
			foreach my $mapper (@storageModels{@loadModels}) {
				$mapper->storeNativeModel();
			}
		}
		
		# Two hacks in a row... Yuck!
		my $mapper = $storageModels{'elasticsearch'};

		my $concept = $model->getConceptDomain('ssm')->conceptHash->{'p'};
		my $corrConcept = BP::Loader::CorrelatableConcept->new($concept);
		
		$mapper->setDestination($corrConcept);
		
		my $destination = $mapper->getInternalDestination();
		
		my @sortedFiles = ();
		my $numsorted = 1;
		# Let's pre-process the data files
		foreach my $vcffile (@ARGV) {
			# First, its sections from the report, with byte offsets
			if(open(my $VCF,BP::Loader::CorrelatableConcept::GUNZIP.' -c '.$vcffile.' | grep -nF "#" |')) {
				print STDERR "Reading header from $vcffile\n";
				my($baselineno,$p_samples);
				while(my $line=<$VCF>) {
					my($lineno,$header) = split(/:/,$line,2);
					
					next  unless(substr($header,0,1) eq '#');
					
					chomp($header);
					
					# Let's only read the headers, not the metadata
					if(substr($header,1,1) ne '#') {
						$baselineno = $lineno+1;
						my @columns = split(/\t/,$header);
						$p_samples = [ @columns[9..$#columns] ];
						
						last;
					}
				}
				close($VCF);
				
				if(defined($baselineno)) {
					print STDERR "Sorting $vcffile\n";
					# Now we can resave and sort the file by chromosome, coordinates and mutation
					my $sortedFile = File::Spec->catfile($workingDir,$numsorted.'_'.File::Basename::basename($vcffile,'.vcf.gz').'.tab.gz');
					
					system(BP::Loader::CorrelatableConcept::GUNZIP.' -c '.$vcffile.' | tail -n +'.$baselineno.' | '.BP::Loader::CorrelatableConcept::SORT.' --parallel='.$BP::Loader::CorrelatableConcept::NUMCPUS.' -S 50% '.' -k1,1 -k2,2n -k4,4 -k5,5 | '.BP::Loader::CorrelatableConcept::GZIP.' -9c > '.$sortedFile);
					
					# We record the sorted file and the sample names
					push(@sortedFiles,[$sortedFile,$p_samples]);
					$numsorted++;
				}
			} else {
				Carp::carp("Unable to open file $vcffile. Skipping...");
			}
		}
			
		# Second pass, reopen all the files and read each first line
		my @joiningFiles = @sortedFiles;
		foreach my $p_sorted (@joiningFiles) {
			if(open(my $VCF,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$p_sorted->[0])) {
				my $first = <$VCF>;
				chomp($first);
				my @colvalues = split(/\t/,$first);
				$p_sorted->[2] = $VCF;
				$p_sorted->[3] = \@colvalues;
				$p_sorted->[4] = [@colvalues[0,1,3,4]];
			} else {
				Carp::croak("Unable to process file $p_sorted->[0]. Dying...");
			}
		}
		
		@joiningFiles = sort(compareVCFlines @joiningFiles);
		
		# Third, read all the existing id entries, and open the file
		print STDERR "Fetching existing entry ids\n";
		my $existingFile = File::Spec->catfile($workingDir,'000_existing.tab.gz');
		my $numEntries = $mapper->existingEntries(['chromosome','chromosome_start','mutated_from_allele','mutated_to_allele'],$existingFile);
		my $existingId = undef;
		my @existingCols = ();
		my $p_existingCols = \@existingCols;
		my $EXISTING = undef;
		if($numEntries > 0) {
			if(open($EXISTING,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$existingFile)) {
				my $first = <$EXISTING>;
				if(defined($first)) {
					chomp($first);
					($existingId,@existingCols) = split(/\t/,$first);
				} else {
					$EXISTING->close();
					$EXISTING = undef;
				}
			} else {
				Carp::croak("Unable to process file $existingFile. Dying...");
			}
		}
		
		# And now, coordinate all of them to push the data!!!!!!
		print STDERR "Bulk upsert starts (in updates of ",UPDATESIZE,")\n";
		
		# my @bulkEntries = ();
		my $nBulk = 0;
		my $totalBulk = 0;
		my $hasOneFile = scalar(@joiningFiles)>0;
		my $hasManyFiles = scalar(@joiningFiles)>1;
		my $numIns = 0;
		my $numUpd = 0;
		while($hasOneFile) {
			# Choose the least value from the candidates
			my $reprvalues = $joiningFiles[0]->[3];
			
			# Process the winners
			# unlinking from the original value
			my $colvalrepr = [@{$joiningFiles[0]->[4]}];
			
			# The data entry
			my $mutation_type=undef;
			
			my @mutationData = ();
			my $expressed_allele = $reprvalues->[4];
			my $expressed_allele0 = undef;
			
			my $chromosome_start = $reprvalues->[1];
			my $chromosome_end = undef;
			
			if(index($expressed_allele,',')!=-1) {
				$expressed_allele = [ split(/,/,$expressed_allele) ];
				$expressed_allele0 = $expressed_allele->[0];
			} else {
				$expressed_allele0 = $expressed_allele;
			}

			if(length($reprvalues->[3])==1 && length($expressed_allele0)==1) {
				$mutation_type = ($reprvalues->[3] eq '-')?'I':(($expressed_allele0 eq '-')?'D':'S');
			} else {
				$mutation_type = 'M';
			}
			
			$chromosome_end = int($chromosome_start);
			$chromosome_end += length($reprvalues->[3]) - 1  if($reprvalues->[3] ne '-');
			

			my %entry =(
				'chromosome'	=> $reprvalues->[0],
				'chromosome_start'	=> $chromosome_start,
				'chromosome_end'	=> $chromosome_end,
				'mutation_type'	=> $mutation_type,
				'reference_genome_allele'	=> $reprvalues->[3],
				'expressed_allele'	=> $expressed_allele,
				'control_genotype'	=> $reprvalues->[3].'/'.$reprvalues->[3],
				'sample_genotype'	=> $reprvalues->[3].'/'.$expressed_allele,
				'mutated_from_allele'	=> $reprvalues->[3],
				'mutated_to_allele'	=> $reprvalues->[4],
				#'probability'	=>,
				'data' => \@mutationData,
			);
			
			# The identified mutation
			if($reprvalues->[2] ne '.') {
				my %mutationIdent =(
					'other'	=> {},
				);
				
				foreach my $keyval (split(/;/,$reprvalues->[2])) {
					my $eqpos = index($keyval,'=');
					my $key;
					my $value=undef;
					my $p_data = undef;
					if($eqpos!=-1) {
						$key = substr($keyval,0,$eqpos);
						$value = substr($keyval,$eqpos+1);
						my @kvdata = split(/,/,$value);
						$p_data = (scalar(@kvdata)>1)?\@kvdata:$value;
					} else {
						$key = $keyval;
					}
					if(exists($IDENTfields{$key})) {
						$mutationIdent{$IDENTfields{$key}} = $p_data;
					} else {
						$mutationIdent{'other'}{$key} = $p_data;
					}
				}
				$entry{'mutation_ident'} = \%mutationIdent
			}
			
			my $equalPos = 0;
			my @leftFiles = ();
			my $erased = undef;
			foreach my $chosen (@joiningFiles) {
				last  unless($chosen->[4] ~~ $colvalrepr);
				$equalPos++;
				
				do {
					my $tabvalues = $chosen->[3];
					
					my @sampleFields = split(/:/,$tabvalues->[8]);
					# Consolidating the entry
					my %sampleMutInfo = ();
					my $total_read_count = -1;
					my $mutant_allele_read_count = undef;
					
					foreach my $keyval (split(/;/,$tabvalues->[7])) {
						my $eqpos = index($keyval,'=');
						my $key;
						my $value=undef;
						my $p_data = undef;
						if($eqpos!=-1) {
							$key = substr($keyval,0,$eqpos);
							$value = substr($keyval,$eqpos+1);
							my @kvdata = split(/,/,$value);
							$p_data = (scalar(@kvdata)>1)?\@kvdata:$value;
						} else {
							$key = $keyval;
						}
						$sampleMutInfo{$key} = $p_data;
						if($key eq 'DP') {
							$total_read_count = $value;
						} elsif($key eq 'IS') {
							($mutant_allele_read_count) = split(/,/,$value);
						}
					}
					
					my $samplePos = 9;
					foreach my $sample (@{$chosen->[1]}) {
						my @sampleValues = split(/:/,$tabvalues->[$samplePos]);
						my %sampleMutData = map { $sampleFields[$_] => $sampleValues[$_] } (0..$#sampleFields);
						my %dataEntry = (
							'analysis_id'	=> $sample,
							'analyzed_sample_id'	=> $sample,
							'total_read_count'	=> $total_read_count,
							'quality_score'	=> $tabvalues->[5],
							'sample_mut_info'	=> \%sampleMutInfo,
							'mut_data'	=> \%sampleMutData
						);
						$dataEntry{'mutant_allele_read_count'} = $mutant_allele_read_count  if(defined($mutant_allele_read_count));
						
						push(@mutationData,\%dataEntry);
						$samplePos++;
					}
					
					# Preparing next read
					if($chosen->[2]->eof) {
						$chosen->[2]->close();
						$chosen->[2] = undef;
					} else {
						my $line = $chosen->[2]->getline();
						chomp($line);
						my @tokens = split(/\t/,$line);
						@{$chosen->[3]} = @tokens;
						@{$chosen->[4]} = @tokens[0,1,3,4];
					}
					
				} while(defined($chosen->[2]) && $colvalrepr ~~ $chosen->[4]);
				
				# Save for next round, but only if it can give us more data!
				if(defined($chosen->[2])) {
					push(@leftFiles,$chosen);
				} else {
					$erased = 1;
				}
			}
			
			# Now, is it an insert or an update?
			if(defined($EXISTING)) {
				while(defined($EXISTING) && compareVCFcols($p_existingCols,$colvalrepr)<0) {
					if($EXISTING->eof()) {
						$EXISTING->close();
						$EXISTING = undef;
					} else {
						my $first = <$EXISTING>;
						chomp($first);
						($existingId,@existingCols) = split(/\t/,$first);
					}
				}
			}
			
			# Validation and default values filling
			my $entorp = $mapper->validateAndEnactEntry(\%entry);
			
			# Pushing the compound entry to Elasticsearch
			my $bulkDef = undef;
			if(defined($EXISTING) && $p_existingCols ~~ $colvalrepr) {
				#$bes->update({
				#	id => $existingId,
				#	lang => 'mvel',
				#	script => 'ctx._source.data += newdoc',
				#	params => {
				#		newdoc => \@mutationData
				#	}
				#});
				$mapper->_incrementalUpdate($existingId,['data',\@mutationData]);
				$numUpd++;
			} else {
				#print $j->encode(\%entry),"\n";
				# push(@bulkEntries,\%entry);
				$mapper->_bulkInsert($destination,$entorp);
				$numIns++;
			}
			
			$nBulk ++;
			#if($nBulk >= UPDATESIZE) {
			#	$totalBulk += $nBulk;
			#	# $mapper->bulkInsert($destination,\@bulkEntries);
			#	# @bulkEntries = ();
			#	$nBulk = 0;
			#	print STDERR "INFO: Upserted $totalBulk entries...\n";
			#}
			
			# Next round!
			if($hasManyFiles) {
				unless($erased) {
					@joiningFiles = sort compareVCFlines @joiningFiles;
				} else {
					@joiningFiles = sort compareVCFlines (@joiningFiles[$equalPos..$#joiningFiles],@leftFiles);
					$hasOneFile = scalar(@joiningFiles)>0;
					$hasManyFiles = scalar(@joiningFiles)>1;
				}
			} elsif($erased) {
				$hasOneFile = undef;
			}
		}
		if($totalBulk > 0 || $nBulk > 0) {
			#$mapper->bulkInsert($destination,\@bulkEntries);
			$totalBulk += $nBulk;
			$mapper->freeDestination();
			print STDERR "INFO: Upserted $totalBulk entries ($numIns insertions, $numUpd updates)...\n";
		}
	}
}
