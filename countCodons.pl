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

sub parseGenBankFile($$);
sub printCodonTable($$$);
sub usage($);

# main body
MAIN:{
    
    # query file name
    my $gbk_file;
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
    my %codon_table = ();
    
    my %genetic_code = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
    
    my %aminoacid_map = (
    'A' => 'Ala', 'R' => 'Arg', 'N' => 'Asn', 'D' => 'Asp',
    'C' => 'Cys', 'Q' => 'Gln', 'E' => 'Glu', 'G' => 'Gly',
    'H' => 'His', 'I' => 'Ile', 'L' => 'Leu', 'K' => 'Lys',
    'M' => 'Met', 'F' => 'Phe', 'P' => 'Pro', 'S' => 'Ser',
    'T' => 'Thr', 'W' => 'Trp', 'Y' => 'Tyr', 'V' => 'Val',
    'B' => 'Asx', 'Z' => 'glx', 'X' => 'Xaa', '*' => 'Stop'
    );
    
    # parse GenBank file
    parseGenBankFile($gbk_file, \%codon_table);
    
    # print codon table
    #printCodonTable(\%codon_table, \%genetic_code, \%aminoacid_map);
}

### --- SUBROUTINES --- ###
sub printCodonTable($$$)
{
    my $codon_table_ref = $_[0];
    my $genetic_code = $_[1];
    my $aminoacid_map = $_[2];
    
    # print header
    print "# Gene.Symbol\t";
    print "Transcript\t";
    my @list_codons = sort keys %{$genetic_code};
    my @list_aa = ();
    my @list_abbreviation = ();
    
    foreach (@list_codons)
    {
        my $aa = $genetic_code->{$_};
        push(@list_aa, $aa);
        push(@list_abbreviation, $aminoacid_map->{$aa});
    }
    
    # sort based on AminoAcid letter
    my @index_sort = sort {$list_aa[$a] cmp $list_aa[$b]} 0..$#list_aa;
    @list_codons = @list_codons[@index_sort];
    @list_aa = @list_aa[@index_sort];
    @list_abbreviation = @list_abbreviation[@index_sort];
    
    
    print join("\t",@list_codons),"\n";
    print "#\t#\t",join("\t", @list_aa),"\n";
    print "#\t#\t",join("\t", @list_abbreviation),"\n";
    
    # print counts
    foreach my $gene (sort keys %{$codon_table_ref})
    {
        foreach my $transcript (sort keys %{$codon_table_ref->{$gene}})
        {
            my @count = ();
            foreach my $codon (@list_codons)
            {
                my $value = exists($codon_table_ref->{$gene}{$transcript}{$codon}) ? $codon_table_ref->{$gene}{$transcript}{$codon} : 0;
                push(@count, $value);
            }
            
            print $gene,"\t",$transcript,"\t",join("\t", @count),"\n";
        }
    }
}



sub parseGenBankFile($$)
{
    my $gbk_file = $_[0];
    my $codon_table_ref = $_[1];
    my $count_records = 0;
    my $count_protein_coding = 0;
    
    # open file to read
    open(my $fh, "gunzip -c $gbk_file |") or die("Can't open $gbk_file to read!\n");
    
    # record separator
    local $/ = "//\n";
    
    # loop through each record
    while(my $record = <$fh>)
    {
        # remove new line
        chomp($record);
        $count_records++;
        next if ($record !~ m/\/protein_id=/);
        
        # parse features
        my $symbol = ($record =~ m/\/gene=\"(.+)\"/) ? $1 : "<symbol>";
        my $refseq_id = ($record =~ m/LOCUS\s+([ANXY][MPR]\_[0-9]+)/) ? $1 : "<refseq_id>";
        my $cds_start = ($record =~ m/CDS\s+(\d+)/) ? $1 : -1;
        my $cds_end = ($record =~ m/CDS\s+\d+\.\.(\d+)/) ? $1 : -1;
        next if (($cds_start == -1) || ($cds_end == -1));
        
        # parse sequence
        my $seq_start = index($record, "ORIGIN");
        my $seq_record = substr($record, $seq_start, length($record) - $seq_start);
        $seq_record =~ s/[^acgtn]//g;
        $seq_record =~ tr/[acgtn]/[ACGTN]/;
        
        if (($cds_end - $cds_start) > 0)
        {
            print $refseq_id,";",$symbol,"\t",substr($seq_record, $cds_start - 1, $cds_end - $cds_start + 1),"\n";
        }
        
        
        # count codons
=head
        for (my $k = $cds_start - 1; $k < $cds_end; $k+=3)
        {
            my $codon = substr($seq_record, $k, 3);
            $codon_table_ref->{$symbol}{$refseq_id}{$codon}++;
        }
=cut
        $count_protein_coding++;
    }
    close($fh);
    
    print STDERR "GenBank records: ", $count_records, "\n";
    print STDERR "GenBank protein coding: ", $count_protein_coding, "\n";
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
