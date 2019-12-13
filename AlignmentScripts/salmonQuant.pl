#!/usr/bin/perl
#
#
#use warnings;
use strict;

if ($#ARGV != 5) {
    print "usage: Enter full path to target genome directory, output directory, directory with R1 reads, directory with R2 reads, and directory with contaminant genome \ne.g. countRNASeqReads.pl /home/genome /home/out /home/reads/1 /home/reads/2 /home/other_genome";
    exit;
}
my $GenomeDir = $ARGV[0];
my $out = $ARGV[1];
my $pair1 = $ARGV[2];
my $pair2 = $ARGV[3];
my $Hs = $ARGV[4];
my $GTF = $ARGV[5];


chomp $GenomeDir;
if ( -d $GenomeDir )
{
    print "\nGenome Directory: $GenomeDir\n";
}
else { die "\nGenome Directory does not exist\n"; }

getReference("$GenomeDir");

chomp $out;
makeOutDirs("$out");

chomp $pair1;
chomp $pair2;
chomp $Hs;
chomp $GTF;

my @PairedReadfiles1 = ();
my @PairedReadfiles2 = ();
#my $trim1;
#my $trim2;
#my @trimReads1 = ();
#my @trimReads2 = ();

chomp $pair1;

chomp $pair2;

@PairedReadfiles1 = getReads("$pair1");
@PairedReadfiles2 = getReads("$pair2");

#QCreads("$pair1", "$pair2", "$out", \@PairedReadfiles1, \@PairedReadfiles2, "fastq");

#($trim1, $trim2) = trimReads("$pair1", "$pair2", "$out", \@PairedReadfiles1, \@PairedReadfiles2);

#@trimReads1 = getReads("$trim1");
#@trimReads2 = getReads("$trim2");

align($GenomeDir, "$pair1", "$pair2", "$out", "$Hs", "$GTF", \@PairedReadfiles1, \@PairedReadfiles2);

############################################################
sub getReference
{
    my $GenomeDir = shift();
    my @Files = ();
    my $file;
    my @Ext = ();
    my $fasta;
    my $ext;
    my $i = 0;
    my $found = 0;
    my $GTFfile;
    my $Index;
    if (-d $GenomeDir)
    {
        opendir (GENOMEDIR, "$GenomeDir");
        @Files = readdir GENOMEDIR;
        closedir GENOMEDIR;
        splice (@Files, 0, 2);
        foreach $file (@Files)
        {
            ($Ext[$i]) = $file =~ /(\.[^.]+)$/;
		if ($Ext[$i] eq "")
		{
			$i++;
		}
            	elsif ((($Ext[$i] eq ".fasta") || ($Ext[$i] eq ".fa")) && ($found == 0))
            	{
                	$fasta = $file;
                	$found = 1;
                	print "\n\tUsing genome FASTA file $GenomeDir/$fasta\n";
            	}
            	else{$i++};
        }
    }
    else {die "\n\t**Reference directory not found**\n";}

}
############################################################
sub makeOutDirs
{
    my $topOut = shift();
    unless (-d $topOut)
    {system "mkdir $topOut";}
    if ( -d $topOut)
    {
        unless (-d "$topOut/quant")
        {system "mkdir $topOut/quant";}
        print "\tSalmon quant files stored in $topOut/quant\n";

    }
}
############################################################
sub getReads
{
    my $readsDir = shift();
    my @Reads = ();
    my $elem;
    if ( -d $readsDir)
    {
        opendir (READS, "$readsDir");
        @Reads = readdir READS;
        closedir READS;
        splice (@Reads, 0, 2);
    }
    else { die "\n\t**Reads not found**\n"; }

    my @sortReads = sort(@Reads);
    return @sortReads;
    }
############################################################
sub QCreads
{
    my $read1Dir = shift;
    my $read2Dir = shift;
    my $outDir = shift;
    my $pair1Reads = shift;
    my $pair2Reads = shift;
    my $format = shift;
    my @Reads1 = @{$pair1Reads};
    my @Reads2 = @{$pair2Reads};
    my $i;
    my $size = @Reads1;
    my $size2 = @Reads2;
    if ($size != $size2) {die "Paired end read files unequal";}

    for ($i = 0; $i < $size; $i++)
    {
        system "fastqc -o $outDir/fastqc_in -f $format $read1Dir/$Reads1[$i] $read2Dir/$Reads2[$i]";
    }
}
############################################################
sub trimReads
{
    my $read1Dir = shift;
    my $read2Dir = shift;
    my $outDir = shift;
    my $pair1Reads = shift;
    my $pair2Reads = shift;
    my @Reads1 = @{$pair1Reads};
    my @Reads2 = @{$pair2Reads};
    my $i;
    my $size = @Reads1;
    my $size2 = @Reads2;
    my @PRE;
    my $prefix;
    if ($size != $size2) {die "\n\t***Paired-end read file names unequal***\n";}

    for ($i = 0; $i < $size; $i++)
    {
        print "\n$read1Dir/$Reads1[$i] $read2Dir/$Reads2[$i]\n";
        @PRE = split('_', $Reads1[$i]);
        $prefix = $PRE[0];
        print "\n\tTrimming $prefix reads\n";
        my @PRE2 = split('_', $Reads2[$i]);
        my $pre2 = $PRE2[0];
        if($prefix eq $pre2)
        {
           system "java -jar /gpfs/fs1/home/jws48/Software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -summary $out/trim/Summary/$prefix.summary $read1Dir/$Reads1[$i] $read2Dir/$Reads2[$i] $outDir/trim/1/$prefix.1.fastq.gz $outDir/trim/singleton/1/$prefix.1_unpaired.fq.gz $outDir/trim/2/$prefix.2.fastq.gz $outDir/trim/singleton/2/$prefix.2_unpaired.fq.gz ILLUMINACLIP:/home/jws48/Software/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50";

           system "fastqc -o $outDir/fastqc_trim $outDir/trim/1/$prefix.1.fastq $outDir/trim/2/$prefix.2.fastq";
        }
    }

    my $trimReads1 = "$outDir/trim/1";
    my $trimReads2 = "$outDir/trim/2";

    return ($trimReads1, $trimReads2);
}
############################################################
sub align
{
    my $GenomeDir = shift;
    my $read1Dir = shift;
    my $read2Dir = shift;
    my $out = shift;
    my $Hs = shift;
    my $GTF = shift;
    my $pair1Reads = shift;
    my $pair2Reads = shift;
    my @Reads1 = @{$pair1Reads};
    my @Reads2 = @{$pair2Reads};
    my $i;
    my $size = @Reads1;
    my $size2 = @Reads2;
    my @PRE;
    my $prefix;

    if ($size != $size2) {die "Paired-end read file names unequal";}

    for ($i = 0; $i < $size; $i++)
    {
        @PRE = split(/\./, $Reads1[$i]);
        $prefix = $PRE[0];

        print "\n$prefix\n";

        my $read1 = "$read1Dir/$prefix" . ".1.fastq.gz";
        my $read2 = "$read2Dir/$prefix" . ".2.fastq.gz";

		system "mkdir $out/quant/$prefix";

		system "salmon quant -i $Hs -l A -1 $read1 -2 $read2 -p 8 --validateMappings -o $out/quant/$prefix.Hs";

		system "salmon quant -i $GenomeDir -l A -1 $read1 -2 $read2 -p 8 -g $GTF --validateMappings -o $out/quant/$prefix";
	}
}
