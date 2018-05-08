#!/usr/bin/perl

# parseGenBankFile
# parse RNA/Protein GenBank file
# extracts symbol and RefSeq IDs information
#
# Input
#   -gbk  : GenBank compressed file as downloaded for NCBI
#   -help : prints a help message
#
# Output
#   tab-delimited table of RefSeq IDs, Official Symbols and
#   records description
#
# Version 1.0
# Date Sep 2016
# Georgi Tushev
# Scientific Computing Facility
# Max-Planck Institute for Brain Research
# send bug reports to sciclist@brain.mpg.de

use warnings;
use strict;
use Getopt::Long();

sub parseGBKRecord($$);
sub printInfo($);
sub printCodons($$);
sub usage($);

# main body
MAIN:{
    
    # query file name
    my $gbk_file;
    my $gbk_prot;
    my $gbk_rna;
    my $help;
    
    # set-up input parameters
    Getopt::Long::GetOptions(
    "gbk=s" => \$gbk_file,
    "help" => \$help
    ) or usage("Error :: invalid command line options");
    
    
    # default help output
    usage("version 1.0, Sep 2016") if($help);
    
    # check if genbank file is provided
    usage("Error :: gz compressed GenBank is required") unless defined($gbk_file);
    
    # record counter
    my $counter = 0;
    my %codons = ();
    
    my @codon_tbl = (
    "ATT", "ATC", "ATA", # I isoleucine
    "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", # L leucine
    "GTT", "GTC", "GTA", "GTG", # V valine
    "TTT", "TTC", # F phenylalanine
    "ATG", # M methionine
    "TGT", "TGC", # C cysteine
    "GCT", "GCC", "GCA", "GCG", # A alanine
    "GGT", "GGC", "GGA", "GGG", # G glycine
    "CCT", "CCC", "CCA", "CCG", # P proline
    "ACT", "ACC", "ACA", "ACG", # T threonine
    "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", # S serine
    "TAT", "TAC", # Y tyrosine
    "TGG", # W tryptophan
    "CAA", "CAG", # Q glutamine
    "AAT", "AAC", # N asparagine
    "CAT", "CAC", # H histidine
    "GAA", "GAG", # E glutamic acid
    "GAT", "GAC", # D aspartic acid
    "AAA", "AAG", # K lysine
    "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", # R arginine
    "TAA", "TAG", "TGA" # stop codons
    );
    
    
    # open file to read
    open(my $fh, "gzcat $gbk_file |") or die("Can't open $gbk_file to read!\n");
    
    # record separator
    local $/ = "//\n";
    
    # print header line
    #print "#Symbol\tRefSeqID\tRefSeqRef\tEntrezID\tGI\tLength\tOrganism\tDescription\n";
    print "#Symbol\tRefSeqID\t",join("\t",@codon_tbl),"\n";
    
    # loop through each record
    while(my $record = <$fh>)
    {
        # remove new line
        chomp($record);
        
        # parse GBK record
        my $info = parseGBKRecord($record, \%codons);
        
        # print info
        #printInfo($info);
        
        # increment record counter
        $counter++;
        
        #last;
    }
    close($fh);
    
    #print condon table
    printCodons(\%codons, \@codon_tbl);
    
    #print $record,"\n";
    print STDERR $counter,"\n";
    
}

### --- SUBROUTINES --- ###
sub parseGBKRecord($$)
{
    # get locally a record
    my $record = $_[0];
    my $codons_ref = $_[1];
    
    # parse features
    my $symbol = ($record =~ m/\/gene=\"(.+)\"/) ? $1 : "<symbol>";
    
    my $refseq_id = ($record =~ m/LOCUS\s+([ANXY][MPR]\_[0-9]+)/) ? $1 : "<refseq_id>";
    
    my $entrez_id = ($record =~ m/\/db_xref=\"GeneID:([0-9]+)\"/) ? $1 : "<entrez_id>";
    
    my $refseq_ref = ($record =~ m/\/protein_id=\"(?:.*)([ANXY][MRP]\_[0-9]+)/) ? $1 : "<refseq_ref>";
    if ($refseq_ref eq "<refseq_ref>")
    {
        $refseq_ref = ($record =~ m/\/coded_by=\"(?:.*)([ANXY][MRP]\_[0-9]+)/) ? $1 : "<refseq_ref>";
    }
    
    my $gi_prot = ($record =~ m/GI\:([0-9]+)/) ? $1 : "<gi_prot>";
    
    my $organism = ($record =~ m/ORGANISM\s+(.+)\n/) ? $1 : "<organism>";
    
    my $length = ($record =~ m/LOCUS\s+[ANXY][MPR]\_[0-9]+\s+([0-9]+)\s/) ? $1 : "<length>";
    
    my $description = ($record =~ m/DEFINITION\s+(.+)ACCESSION/s) ? $1 : "<description>";
    
    # clean description
    $description =~ s/\n//g;
    $description =~ s/ +/ /g;
    $description =~ s/ \[$organism\]\.//g;
    $description =~ s/\s?$organism//g;
    
    # get CDS positions
    my $cds_start = ($record =~ m/CDS\s+(\d+)/) ? $1 : "<unknown>";
    my $cds_end = ($record =~ m/CDS\s+\d+\.\.(\d+)/) ? $1 : "<unknown>";
    
    # get sequence
    my $seq_start = index($record, "ORIGIN");
    my $seq_record = substr($record, $seq_start, length($record) - $seq_start);
    $seq_record =~ s/[^acgtn]//g;
    $seq_record =~ tr/[acgtn]/[ACGTN]/;
    
    # count codons
    if (($cds_start ne "<unknown>") && ($cds_end ne "<unknown>"))
    {
        for (my $k = $cds_start - 1; $k < $cds_end; $k+=3)
        {
            my $codon = substr($seq_record, $k, 3);
            $codons_ref->{$symbol}{$refseq_id}{$codon}++;
        }
    }
    
    # Return result
    return [$symbol, $refseq_id, $refseq_ref, $entrez_id, $gi_prot, $length, $organism, $description];
    
}


sub printInfo($)
{
    my $info_ref = shift;
    my $info_siz = scalar(@{$info_ref});
    for(my $k = 0; $k < $info_siz; $k++)
    {
        print $info_ref->[$k],"\t" if ($k < ($info_siz - 1));
        print $info_ref->[$k],"\n" if ($k == ($info_siz - 1));
    }
    
}

sub printCodons($$)
{
    my $codon_ref = $_[0];
    my $codon_tbl = $_[1];
    
    foreach my $gene (sort keys %{$codon_ref})
    {
        foreach my $transcript (sort keys %{$codon_ref->{$gene}})
        {
            my @count = ();
            foreach my $codon (@{$codon_tbl})
            {
                my $value = exists($codon_ref->{$gene}{$transcript}{$codon}) ? $codon_ref->{$gene}{$transcript}{$codon} : 0;
                push(@count, $value);
            }
            
            print $gene,"\t",$transcript,"\t",join("\t", @count),"\n";
        }
    }
}


sub usage($)
{
    my $message = $_[0];
    if (defined $message && length $message)
    {
        $message .= "\n" unless($message =~ /\n$/);
    }
    
    my $command = $0;
    $command =~ s#^.*/##;
    
    print STDERR (
    $message,
    "usage: $command -gbk genebank_file.gbk.gz\n" .
    "description: parse RNA/Protein GenBank file and extracts symbol and RefSeq IDs information\n" .
    "parameters:\n" .
    "-gbk\n" .
    "\tGZIP compressed GenBank file as downloaded from NCBI ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/RNA/\n" .
    "-help\n" .
    "\tdefine usage\n"
    );
    
    die("\n");
}
