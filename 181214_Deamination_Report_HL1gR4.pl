#!/usr/bin/perl -w
##########################################################################################
# date:   		2018-12-14																 
# author:  		CRISPRengineer																 
# description: 																			
# usage: 		
##########################################################################################

		


use strict;
use Cwd;

#---------------------------------------------------------------------------------------#
my $start=(times)[0];			# save the start time
#---------------------------------------------------------------------------------------#

our($opt_i, $opt_b, $opt_l); # set up global variables for Getopt
get_args();			# get the command line arguments 

#print "\$opt_i:\t\t$opt_i\n\n";

#my @inputFiles = @{get_inputFiles($opt_i)};

###################################
# get input files #
###################################

my ($barcode,$reference) = get_inputFiles($opt_i);
my @fastq = @{get_fastq($opt_i)};
my @bam = @{get_bam($opt_i)};
my @sam = @{get_sam($opt_i)};
my @fasta = @{get_fasta($opt_i)}; 
my $fasta;

open(my $FA, "$opt_i/$reference");

while(<$FA>){
	chomp;
	if ($_=~/^>/) {}
	else {
		$fasta.=$_;
	}
}

close($FA);

#print "\n";

#print "$fasta\n";

my @bpfasta=split //, $fasta;

#print join "@bpfasta\n";

#print "$bpfasta[269]\n";
#print @fastq;


###################################
# Output data file #
###################################

my $data="$opt_i/SNV_report.txt";
my $metadata="$opt_i/SNV_mutations_report.txt";
open(DATA, ">$data");
open(METADATA, ">$metadata");

print "#data file:\t\t\t$data\n\n";

my $indir="$opt_i/";
my $outdir="$opt_i/";

#print join ("\n" , @sam);
#print join ("\n" , @bam);

my $mpileupcmd;
my $mpileupcapture;
my @capturearray;
my ($Tcount,$Acount,$Gcount,$Ccount);
my $ref="HL1gR4";
my $StartPos="3407";
my $EndPos="3453";


my $c;
my $adjBase;
my $altNucCount;

my $totalReads=0;
my $totalAltReads=0;

my	$Cfreq;
my	$Tfreq;
my	$Afreq;
my	$Gfreq;


print DATA "sample\tpos\ttotal reads\tG\tA\tT\tC\trefBase\tAltBaseCount\tGfreq\tAfreq\tTfreq\tCfreq\n";

print METADATA "sample\tStartPos\tEndPos\tTotalReads\tTotalAltReads\tFreqAltReads\n";

	

foreach(@bam){	
	if ($_=~/(\S*)\.srt.bam/){	
		print METADATA "$1\t$StartPos\t$EndPos\t";		
		for ($c=$StartPos;$c<=$EndPos;$c++){	
			$mpileupcapture=`samtools mpileup -d 10000000 -r $ref:$c-$c $1.srt.bam`;
		
			print "samtools mpileup -r $ref:$c-$c $1.srt.bam\n";
		
			@capturearray=split /\t/, $mpileupcapture;
		
			#my @lines = split /\n/, $str;
		
			print DATA "$1\t$capturearray[1]\t$capturearray[3]\t";
		
			
			$Gcount = ($capturearray[4] =~ tr/Gg//);
			print DATA "$Gcount\t";
			$Acount = ($capturearray[4] =~ tr/Aa//);
			print DATA "$Acount\t";
			$Tcount = ($capturearray[4] =~ tr/Tt//);
			print DATA "$Tcount\t";
			$Ccount = ($capturearray[4] =~ tr/Cc//);
			print DATA "$Ccount\t";
			
			$adjBase=$c-1;
			print DATA "$bpfasta[$adjBase]\t";
			if ($bpfasta[$adjBase]=~/G/){$altNucCount=$Acount+$Tcount+$Ccount}
			elsif ($bpfasta[$adjBase]=~/A/){$altNucCount=$Gcount+$Tcount+$Ccount}
			elsif ($bpfasta[$adjBase]=~/T/){$altNucCount=$Acount+$Gcount+$Ccount}
			elsif ($bpfasta[$adjBase]=~/C/){$altNucCount=$Acount+$Tcount+$Gcount}
			print DATA "$altNucCount\t";
			
			
			$Cfreq=100*$Ccount/$capturearray[3];
			$Tfreq=100*$Tcount/$capturearray[3];
			$Afreq=100*$Acount/$capturearray[3];
			$Gfreq=100*$Gcount/$capturearray[3];
			
			print DATA "$Gfreq\t$Afreq\t$Tfreq\t$Cfreq\n";
			
			#$string =~ tr/A-Z//;
			#print "$mpileupcapture\n";
			#system($mpileupcmd);
			#$result = `command arg1 arg2`;
			#print "samtools mpileup -r Ug2:90-90 -b /Users/cory/Church/NGS/2016-8-8_293T_AS_CORB_fastq/alignment/Ug2/bamfiles.txt\n";
			
			$totalReads=$totalReads+$capturearray[3];
			
			$totalAltReads=$totalAltReads+$altNucCount;
	
			$altNucCount=0;
			
		}
		print METADATA "$totalReads\t";
		print METADATA "$totalAltReads\t".$totalAltReads/$totalReads."\n";
		
		$totalReads=0;
		$totalAltReads=0;
	}
	
}

close(DATA);

close(METADATA);

print "#data file:\n$data\n\n";
my $display =`cat $data`;
print "$display\n";


#=====================================================#
#================== END Main Script ==================#
#=====================================================#


###############################################################
# get_args parses the commandline arguments using Getopt::Std #
###############################################################


#---------------------------------------------------------------------------------------#
sub get_args{
	use Getopt::Std;
	
	# set in main script:
	# our($opt_i); # set up global variables for Getopt
	# -i = input directory

	getopt('ibl');		# return variables following -i -b -l
	
	if (!defined $opt_i) {$opt_i=getcwd()}
	if (!defined $opt_l) {$opt_l=20}	
}
#---------------------------------------------------------------------------------------#

##################################################
# Indir parsing to obtain fastq and barcodes.fil #  
##################################################

#---------------------------------------------------------------------------------------#
sub get_inputFiles {
	my @inputFiles;
	my $dir=shift @_;
	my $reference;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.fil$/) {
    		push(@inputFiles, $file);
    		$barcode=$file;
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fa$/) {
    		
    		$reference=$file;
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fasta$/) {
    		
    		$reference=$file;
    	}
	}
	#return \@inputFiles;
	return ($barcode, $reference);
}
#---------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
sub get_fastq {
	my @inputFastq;
	my $dir=shift @_;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.fastq.gz$/) {
    		push(@inputFastq, $file);
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fastq$/) {
    		push(@inputFastq, $file);
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fq$/) {
    		push(@inputFastq, $file);
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fq.gz$/) {
    		push(@inputFastq, $file);
    	}
	}
	return \@inputFastq;
	
}
#---------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
sub get_bam {
	my @inputBam;
	my $dir=shift @_;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.bam$/) {
    		push(@inputBam, $file);
    	}
	}
	return \@inputBam;
	
}
#---------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
sub get_sam {
	my @inputsam;
	my $dir=shift @_;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.sam$/) {
    		push(@inputsam, $file);
    	}
	}
	return \@inputsam;
}
#---------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
sub get_fasta {
	my @inputFasta;
	my $dir=shift @_;
	opendir(DIR, $dir) or die $!;
	while (my $file = readdir(DIR)) {
		next unless (-f "$dir/$file");
		if ($file !~ /^\./ & $file =~ /\.fa$/) {
    		push(@inputFasta, $file);
    	}
    	elsif ($file !~ /^\./ & $file =~ /\.fasta$/) {
    		push(@inputFasta, $file);
		}
	}	
	return \@inputFasta;
}
#---------------------------------------------------------------------------------------#


################################################
# Calculate and print execution time of script #
################################################

my $end=(times)[0];			# save the end time
my $dt =$end-$start;		# difference is the execution time
print STDERR "Execution time = $dt seconds\n"; # outputs to STDERR