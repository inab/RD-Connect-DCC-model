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
use Sys::CPU;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Elasticsearch;


use Search::Elasticsearch 1.12;
use JSON;

use constant {
	GZIP	=>	'pigz',
	GUNZIP	=>	'unpigz',
	UPDATESIZE	=>	1000,
};

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
	
	my $NUMCPUS = Sys::CPU::cpu_count();

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
		
		my $bes = $mapper->getInternalDestination();
		
		my @sortedFiles = ();
		my $numsorted = 0;
		# Let's pre-process the data files
		foreach my $vcffile (@ARGV) {
			# First, its sections from the report, with byte offsets
			if(open(my $VCF,GUNZIP.' -c '.$vcffile.' | grep -nF "#" |')) {
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
					
					system(GUNZIP.' -c '.$vcffile.' | tail -n +'.$baselineno.' | sort --parallel='.$NUMCPUS.' -S 50% '.' -k1,1 -k2,2n -k4,4 -k5,5 | '.GZIP.' -9c > '.$sortedFile);
					
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
			if(open(my $VCF,'-|',GUNZIP,'-c',$p_sorted->[0])) {
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
		
		# And now, coordinate all of them to push the data!!!!!!
		print STDERR "Bulk insert starts (in updates of ",UPDATESIZE,")\n";
		# my @bulkEntries = ();
		my $nBulk = 0;
		my $totalBulk = 0;
		while(scalar(@joiningFiles)>0) {
			# Choose the least value from the candidates
			my @chosenFiles = ();
			my @leftFiles = ();
			my $reprvalues = undef;
			my $colvalrepr = undef;
			
			foreach my $joining (@joiningFiles) {
				if(defined($colvalrepr) && $colvalrepr ~~ $joining->[4]) {
					push(@chosenFiles,$joining);
				} elsif(defined($colvalrepr) && ($colvalrepr->[0] lt $joining->[4][0] || ($colvalrepr->[0] eq $joining->[4][0] && ($colvalrepr->[1] < $joining->[4][1] || ($colvalrepr->[1] == $joining->[4][1] && ($colvalrepr->[2] lt $joining->[4][2] || ($colvalrepr->[2] eq $joining->[4][2] && $colvalrepr->[3] lt $joining->[4][3]))))))) {
					# It's over!!!!!
					push(@leftFiles,$joining);
				} else {
					push(@leftFiles,@chosenFiles);
					$reprvalues = $joining->[3];
					$colvalrepr = $joining->[4];
					@chosenFiles = ($joining);
				}
			}
			
			# Process the winners
			# unlinking from the original value
			$colvalrepr = [@{$colvalrepr}];
			
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
					'other'	=> [],
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
						push(@{$mutationIdent{'other'}},{'key'=>$key,'value'=>$p_data});
					}
				}
				$entry{'mutation_ident'} = \%mutationIdent
			}
			
			foreach my $chosen (@chosenFiles) {
				do {
					my $tabvalues = $chosen->[3];
					
					my @sampleFields = split(/:/,$tabvalues->[8]);
					# Consolidating the entry
					my @sampleMutInfo = ();
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
						push(@sampleMutInfo,{'key'=>$key,'value'=>$p_data});
						if($key eq 'DP') {
							$total_read_count = $value;
						} elsif($key eq 'IS') {
							($mutant_allele_read_count) = split(/,/,$value);
						}
					}
					
					my $samplePos = 9;
					foreach my $sample (@{$chosen->[1]}) {
						my @sampleValues = split(/:/,$tabvalues->[$samplePos]);
						my @sampleMutData = map { {'key' => $sampleFields[$_], 'value' => $sampleValues[$_] } } (0..$#sampleFields);
						my %dataEntry = (
							'analysis_id'	=> $sample,
							'analyzed_sample_id'	=> $sample,
							'total_read_count'	=> $total_read_count,
							'quality_score'	=> $tabvalues->[5],
							'sample_mut_info'	=> \@sampleMutInfo,
							'mut_data'	=> \@sampleMutData
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
				push(@leftFiles,$chosen)  if(defined($chosen->[2]));
			}
			
			# Pushing the compound entry to Elasticsearch
			#print $j->encode(\%entry),"\n";
			# push(@bulkEntries,\%entry);
			$bes->index({source=>\%entry});
			$nBulk ++;
			if($nBulk >= UPDATESIZE) {
				$totalBulk += $nBulk;
				# $mapper->bulkInsert($destination,\@bulkEntries);
				# @bulkEntries = ();
				$nBulk = 0;
				print STDERR "INFO: Inserted $totalBulk entries...\n";
			}
			
			# Next round!
			@joiningFiles = @leftFiles;
		}
		if($totalBulk > 0) {
			#$mapper->bulkInsert($destination,\@bulkEntries);
			$mapper->freeDestination();
			print STDERR "INFO: Inserted $totalBulk entries...\n";
		}
	}
}
