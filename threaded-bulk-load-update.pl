#!/usr/bin/perl -w

use v5.12;
use strict;
use warnings qw(all);
no warnings qw(experimental);
use diagnostics -warntrace;

use threads;
use threads::shared;
use Thread::Queue;

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
use BP::Loader::Mapper::MongoDB;

use constant MAX_WORKERS => 4;

use constant {
	VCF_SORTED_FILE	=>	0,
	VCF_SAMPLES	=>	1,
	VCF_HANDLER	=>	2,
	VCF_LINE	=>	3,
	VCF_LINE_KEYS	=>	4,
	VCF_METADATA	=>	5,
	VCF_FILE	=>	6,
};

use constant {
	VCF_CHROM	=>	0,
	VCF_POS		=>	1,
	VCF_REF		=>	2,
	VCF_ALT		=>	3,
};

use constant {
	VCF_CHROM_COL	=>	0,
	VCF_POS_COL	=>	1,
	VCF_ID_COL	=>	2,
	VCF_REF_COL	=>	3,
	VCF_ALT_COL	=>	4,
	VCF_QUAL_COL	=>	5,
	VCF_FILTER_COL	=>	6,
	VCF_INFO_COL	=>	7,
	VCF_FORMAT_COL	=>	8,
	VCF_FIRST_SAMPLE_COL	=>	9,
};

my @VCF_COL_KEYS = (VCF_CHROM_COL,VCF_POS_COL,VCF_REF_COL,VCF_ALT_COL);

sub compareVCFlines($$);
sub compareVCFcols($$);

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

# Shared by all the instances
my $numsorted :shared;

sub sortWorker($$) {
	my($workingDir,$queue)=@_;
	
	my @sortedFiles = ();
	while(defined(my $vcffile = $queue->dequeue())) {
		# First, its sections from the report, with byte offsets
		if(open(my $VCF,BP::Loader::CorrelatableConcept::GUNZIP.' -c '.$vcffile.' | grep -nF "#" |')) {
			print STDERR "Reading header from $vcffile\n";
			
			my %metadata = ();
			
			my($baselineno,$p_samples);
			while(my $line=<$VCF>) {
				my($lineno,$header) = split(/:/,$line,2);
				
				next  unless(substr($header,0,1) eq '#');
				
				chomp($header);
				
				# Let's only read the headers, not the metadata
				if(substr($header,1,1) ne '#') {
					$baselineno = $lineno+1;
					my @columns = split(/\t/,$header,-1);
					$p_samples = [ @columns[9..$#columns] ];
					
					last;
				} else {
					my($metakey,$metaval) = split(/=/,substr($header,2),2);
					
					if($metaval =~ /^<([^>]+)>/) {
						my %metameta = ();
						
						my @metacomps = split(/,/,$1,-1);
						foreach my $metacomp (@metacomps) {
							my($subkey,$subval) = split(/=/,$metacomp,2);
							$metameta{$subkey} = $subval;
						}
						
						$metadata{$metakey}{$metameta{'ID'}} = \%metameta;
					} else {
						$metadata{$metakey} = $metaval;
					}
				}
			}
			close($VCF);
			
			if(defined($baselineno)) {
				print STDERR "Sorting $vcffile\n";
				# Now we can resave and sort the file by chromosome, coordinates and mutation
				my $localnumsorted = undef;
				{
					lock($numsorted);
					$localnumsorted = $numsorted;
					$numsorted++;
				}
				
				my $sortedFile = File::Spec->catfile($workingDir,$localnumsorted.'_'.File::Basename::basename($vcffile,'.vcf.gz').'.tab.gz');
				
				system(BP::Loader::CorrelatableConcept::GUNZIP.' -c '.$vcffile.' | tail -n +'.$baselineno.' | '.BP::Loader::CorrelatableConcept::SORT.' --parallel='.$BP::Loader::CorrelatableConcept::NUMCPUS.' -S 50% '.' -k1,1 -k2,2n -k4,4 -k5,5 | '.BP::Loader::CorrelatableConcept::GZIP.' -9c > '.$sortedFile);
				
				# We record the sorted file and the sample names
				my $p_sorted = [];
				@{$p_sorted}[VCF_FILE,VCF_SORTED_FILE,VCF_SAMPLES,VCF_METADATA] = ($vcffile,$sortedFile,$p_samples,\%metadata);
				
				push(@sortedFiles,$p_sorted);
			}
		} else {
			Carp::carp("Unable to open file $vcffile. Skipping...");
		}
	}
	return @sortedFiles;
}


my %IDENTfields = (
	'RS'	=>	'RScode',
	'dbSNPBuildID'	=>	'dbSNPBuildID'
);

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
		foreach my $loadModelName (split(/,/,$loadModelNames,-1)) {
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
		
		my $concept = $model->getConceptDomain('ssm')->conceptHash->{'p'};
		my $corrConcept = BP::Loader::CorrelatableConcept->new($concept);
		
		my $mapper = $storageModels{$loadModels[0]};
		$mapper->setDestination($corrConcept);
		
		my @sortedFiles = ();
		$numsorted = 1;
		# Let's pre-process the data files
		{
			my $queue = Thread::Queue->new();
			
			my @workers = map { threads->create(\&sortWorker,$workingDir,$queue) } 1..MAX_WORKERS;
			
			$queue->enqueue(@ARGV);
			$queue->end();
			foreach my $worker (@workers) {
				push(@sortedFiles,$worker->join());
			}
		}
		
		# Second pass, reopen all the files and read each first line
		my @joiningFiles = @sortedFiles;
		foreach my $p_sorted (@joiningFiles) {
			if(open(my $VCF,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$p_sorted->[VCF_SORTED_FILE])) {
				my $first = <$VCF>;
				chomp($first);
				my @colvalues = split(/\t/,$first,-1);
				$p_sorted->[VCF_HANDLER] = $VCF;
				$p_sorted->[VCF_LINE] = \@colvalues;
				$p_sorted->[VCF_LINE_KEYS] = [@colvalues[@VCF_COL_KEYS]];
			} else {
				Carp::croak("Unable to process file $p_sorted->[VCF_SORTED_FILE]. Dying...");
			}
		}
		
		@joiningFiles = sort(compareVCFlines @joiningFiles);
		
		# Third, read all the existing id entries, and open the file
		print STDERR "Fetching existing entry ids\n";
		my $existingFile = File::Spec->catfile($workingDir,'000_existing.tab.gz');
		my $numEntries = $mapper->existingEntries($existingFile);
		my $existingId = undef;
		my @existingCols = ();
		my $p_existingCols = \@existingCols;
		my $EXISTING = undef;
		if($numEntries > 0) {
			if(open($EXISTING,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$existingFile)) {
				my $first = <$EXISTING>;
				if(defined($first)) {
					chomp($first);
					($existingId,@existingCols) = split(/\t/,$first,-1);
				} else {
					$EXISTING->close();
					$EXISTING = undef;
				}
			} else {
				Carp::croak("Unable to process file $existingFile. Dying...");
			}
		}
		
		# And now, coordinate all of them to push the data!!!!!!
		print STDERR "Bulk upsert starts (in updates of ",$mapper->bulkBatchSize,")\n";
		
		# my @bulkEntries = ();
		my $nBulk = 0;
		my $totalBulk = 0;
		my $hasOneFile = scalar(@joiningFiles)>0;
		my $hasManyFiles = scalar(@joiningFiles)>1;
		my $numIns = 0;
		my $numUpd = 0;
		while($hasOneFile) {
			# Choose the least value from the candidates
			my $reprvalues = $joiningFiles[0]->[VCF_LINE];
			my $p_metadata = $joiningFiles[0]->[VCF_METADATA];
			
			# Process the winners
			# unlinking from the original value
			my $colvalrepr = [@{$joiningFiles[0]->[VCF_LINE_KEYS]}];
			
			# The data entry
			my $mutation_type=undef;
			
			my @mutationData = ();
			my $expressed_allele = $reprvalues->[VCF_ALT_COL];
			my $expressed_allele0 = undef;
			
			my $chromosome_start = $reprvalues->[VCF_POS_COL];
			my $chromosome_end = undef;
			
			if(index($expressed_allele,',')!=-1) {
				$expressed_allele = [ split(/,/,$expressed_allele,-1) ];
				$expressed_allele0 = $expressed_allele->[0];
			} else {
				$expressed_allele0 = $expressed_allele;
			}

			if(length($reprvalues->[VCF_REF_COL])==1 && length($expressed_allele0)==1) {
				$mutation_type = ($reprvalues->[VCF_REF_COL] eq '-')?'I':(($expressed_allele0 eq '-')?'D':'S');
			} else {
				$mutation_type = 'M';
			}
			
			$chromosome_end = int($chromosome_start);
			$chromosome_end += length($reprvalues->[VCF_REF_COL]) - 1  if($reprvalues->[VCF_REF_COL] ne '-');
			

			my %entry =(
				'chromosome'	=> $reprvalues->[VCF_CHROM_COL],
				'chromosome_start'	=> $chromosome_start,
				'chromosome_end'	=> $chromosome_end,
				'mutation_type'	=> $mutation_type,
				'reference_genome_allele'	=> $reprvalues->[VCF_REF_COL],
				'expressed_allele'	=> $expressed_allele,
				'control_genotype'	=> $reprvalues->[VCF_REF_COL].'/'.$reprvalues->[VCF_REF_COL],
				'sample_genotype'	=> $reprvalues->[VCF_REF_COL].'/'.$expressed_allele,
				'mutated_from_allele'	=> $reprvalues->[VCF_REF_COL],
				'mutated_to_allele'	=> $reprvalues->[VCF_ALT_COL],
				#'probability'	=>,
				'data' => \@mutationData,
			);
			
			# The identified mutation
			if($reprvalues->[VCF_ID_COL] ne '.') {
				my %mutationIdent =(
					'other'	=> {},
				);
				
				my $prevkey = undef;
				foreach my $keyval (split(/;/,$reprvalues->[VCF_ID_COL],-1)) {
					my $eqpos = index($keyval,'=');
					my $key = undef;
					my $value=undef;
					my $p_data = undef;
					if($eqpos!=-1) {
						$key = substr($keyval,0,$eqpos);
						$value = substr($keyval,$eqpos+1);
						my @kvdata = split(/,/,$value,-1);
						$p_data = (scalar(@kvdata)>1)?\@kvdata:$value;
					} elsif(exists($p_metadata->{'ID'}{$keyval})) {
						$key = $keyval;
					} else {
						$value = $keyval;
					}
					
					if(defined($key)) {
						if(exists($IDENTfields{$key})) {
							unless(exists($mutationIdent{$IDENTfields{$key}})) {
								$mutationIdent{$IDENTfields{$key}} = $p_data;
							} elsif(defined($value) && $value ne '.') {
								# We don't know what to do here!
								my $oldval = $mutationIdent{$IDENTfields{$key}};
								if(ref($oldval)) {
									push(@{$oldval},(ref($p_data)?@{$p_data}:($p_data)));
								} elsif($oldval ne $value) {
									$mutationIdent{$IDENTfields{$key}} = [$oldval,(ref($p_data)?@{$p_data}:($p_data))];
								}
								#Carp::carp('Illformed VCF file '.$joiningFiles[0]->[VCF_FILE]." on ID column, repeated token $key ($oldval <=> $value): ".join(' ',
								#		'CHROM' => $reprvalues->[VCF_CHROM_COL],
								#		'POS' => $chromosome_start,
								#		'ID' => $reprvalues->[VCF_ID_COL],
								#	)
								#)  if($oldval ne $value);
							}
						} elsif(exists($mutationIdent{'other'}{$key})) {
							if(defined($value) && $value ne '.') {
								# We don't know what to do here!
								my $oldval = $mutationIdent{'other'}{$key};
								if(ref($oldval)) {
									push(@{$oldval},(ref($p_data)?@{$p_data}:($p_data)));
								} elsif($oldval ne $value) {
									$mutationIdent{'other'}{$key} = [$oldval,(ref($p_data)?@{$p_data}:($p_data))];
								}
								#Carp::carp('Illformed VCF file '.$joiningFiles[0]->[VCF_FILE]." on ID column, repeated token $key ($oldval <=> $value): ".join(' ',
								#		'CHROM' => $reprvalues->[VCF_CHROM_COL],
								#		'POS' => $chromosome_start,
								#		'ID' => $reprvalues->[VCF_ID_COL],
								#	)
								#)  if($oldval ne $value);
							}
						} else {
							$mutationIdent{'other'}{$key} = $p_data;
						}
						$prevkey = $key;
					} else {
						#Carp::carp('Illformed VCF file '.$joiningFiles[0]->[VCF_FILE]." on ID column, after token $prevkey: ".join(' ',
						#		'CHROM' => $reprvalues->[VCF_CHROM_COL],
						#		'POS' => $chromosome_start,
						#		'ID' => $reprvalues->[VCF_ID_COL],
						#	)
						#);
						if(exists($IDENTfields{$prevkey})) {
							if(ref($IDENTfields{$prevkey})) {
								push(@{$mutationIdent{$IDENTfields{$prevkey}}},(ref($p_data)?@{$p_data}:($p_data)));
							} else {
								$mutationIdent{$IDENTfields{$prevkey}} = [$mutationIdent{$IDENTfields{$prevkey}},(ref($p_data)?@{$p_data}:($p_data))];
							}
						} elsif(ref($mutationIdent{'other'}{$prevkey})) {
							push(@{$mutationIdent{'other'}{$prevkey}},(ref($p_data)?@{$p_data}:($p_data)));
						} else {
							$mutationIdent{'other'}{$prevkey} = [$mutationIdent{'other'}{$prevkey},(ref($p_data)?@{$p_data}:($p_data))];
						}
					}
				}
				$entry{'mutation_ident'} = \%mutationIdent
			}
			
			my $equalPos = 0;
			my @leftFiles = ();
			my $erased = undef;
			foreach my $chosen (@joiningFiles) {
				last  unless($chosen->[VCF_LINE_KEYS] ~~ $colvalrepr);
				$equalPos++;
				
				do {
					my $tabvalues = $chosen->[VCF_LINE];
					my $p_metadata = $chosen->[VCF_METADATA];
					
					my @sampleFields = split(/:/,$tabvalues->[VCF_FORMAT_COL],-1);
					# Consolidating the entry
					my %sampleMutInfo = ();
					my $total_read_count = -1;
					my $mutant_allele_read_count = undef;
					
					my $prevkey = undef;
					foreach my $keyval (split(/;/,$tabvalues->[VCF_INFO_COL],-1)) {
						my $eqpos = index($keyval,'=');
						my $key = undef;
						my $value=undef;
						my $p_data = undef;
						if($eqpos!=-1) {
							$key = substr($keyval,0,$eqpos);
							$value = substr($keyval,$eqpos+1);
							my @kvdata = split(/,/,$value,-1);
							$p_data = (scalar(@kvdata)>1)?\@kvdata:$value;
						} elsif(exists($p_metadata->{'INFO'}{$keyval})) {
							$key = $keyval;
						} else {
							$value = $keyval;
						}
						
						if(defined($key)) {
							unless(exists($sampleMutInfo{$key})) {
								$sampleMutInfo{$key} = $p_data;
								
								if($key eq 'DP') {
									$total_read_count = $value;
								} elsif($key eq 'IS') {
									$mutant_allele_read_count = $p_data->[0];
								}
							} elsif(defined($value) && $value ne '.') {
								# We don't know what to do here!
								my $oldval = $sampleMutInfo{$key};
								if(ref($oldval)) {
									push(@{$oldval},(ref($p_data)?@{$p_data}:($p_data)));
								} elsif($oldval ne $value) {
									$sampleMutInfo{$key} = [$oldval,(ref($p_data)?@{$p_data}:($p_data))];
								}
								#Carp::carp('Illformed VCF file '.$joiningFiles[0]->[VCF_FILE]." on INFO column, repeated token $key ($oldval <=> $value): ".join(' ',
								#		'CHROM' => $reprvalues->[VCF_CHROM_COL],
								#		'POS' => $chromosome_start,
								#		'INFO' => $reprvalues->[VCF_INFO_COL],
								#	)
								#)  if($oldval ne $value);
							}
							$prevkey = $key;
						} else {
							#Carp::carp('Illformed VCF file '.$joiningFiles[0]->[VCF_FILE]." on INFO column, after token $prevkey: ".join(' ',
							#		'CHROM' => $reprvalues->[VCF_CHROM_COL],
							#		'POS' => $chromosome_start,
							#		'INFO' => $reprvalues->[VCF_INFO_COL],
							#	)
							#);
							if(ref($sampleMutInfo{$prevkey})) {
								push(@{$sampleMutInfo{$prevkey}},(ref($p_data)?@{$p_data}:($p_data)));
							} else {
								$sampleMutInfo{$prevkey} = [$sampleMutInfo{$prevkey},(ref($p_data)?@{$p_data}:($p_data))];
							}
						}
					}
					
					my $samplePos = VCF_FIRST_SAMPLE_COL;
					foreach my $sample (@{$chosen->[VCF_SAMPLES]}) {
						my @sampleValues = split(/:/,$tabvalues->[$samplePos],-1);
						my %sampleMutData = map { $sampleFields[$_] => $sampleValues[$_] } (0..$#sampleFields);
						my %dataEntry = (
							'analysis_id'	=> $sample,
							'analyzed_sample_id'	=> $sample,
							'total_read_count'	=> $total_read_count,
							'quality_score'	=> $tabvalues->[VCF_QUAL_COL],
							'sample_mut_info'	=> \%sampleMutInfo,
							'mut_data'	=> \%sampleMutData
						);
						$dataEntry{'mutant_allele_read_count'} = $mutant_allele_read_count  if(defined($mutant_allele_read_count));
						
						push(@mutationData,\%dataEntry);
						$samplePos++;
					}
					
					# Preparing next read
					if($chosen->[VCF_HANDLER]->eof) {
						$chosen->[VCF_HANDLER]->close();
						$chosen->[VCF_HANDLER] = undef;
					} else {
						my $line = $chosen->[VCF_HANDLER]->getline();
						chomp($line);
						my @tokens = split(/\t/,$line,-1);
						@{$chosen->[VCF_LINE]} = @tokens;
						@{$chosen->[VCF_LINE_KEYS]} = @tokens[@VCF_COL_KEYS];
					}
					
				} while(defined($chosen->[VCF_HANDLER]) && $colvalrepr ~~ $chosen->[VCF_LINE_KEYS]);
				
				# Save for next round, but only if it can give us more data!
				if(defined($chosen->[VCF_HANDLER])) {
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
						($existingId,@existingCols) = split(/\t/,$first,-1);
					}
				}
			}
			
			my $isUpdate = defined($EXISTING) && $p_existingCols ~~ $colvalrepr;
			$entry{BP::Loader::Mapper::COL_INCREMENTAL_UPDATE_ID} = $existingId  if($isUpdate);
			
			# Validation and default values filling
			# Pushing the compound entry to Elasticsearch
			eval {
				$mapper->bulkInsert(\%entry);
			};
			
			if($@) {
				use Data::Dumper;
				
				print STDERR "ERROR: $@\n\n",Dumper(\%entry),"\n";
				
				exit 1;
			}
			if($isUpdate) {
				$numUpd++;
			} else {
				$numIns++;
			}
			
			$nBulk ++;
			
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
