#!/usr/bin/env perl

use strict;
use warnings;

#
# imports and whatnot
#

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use WWW::Mechanize;
use Bio::Tools::CodonTable;

use lib dirname(__FILE__) . '/lib';
use ExonFrames ':all';

my $dir = dirname(__FILE__);

#
# get options and args and throw errors if they're wrong
#

my $help        = '';
my $file_prefix = '';

GetOptions(
    'help' => \$help,
    'd=s'  => \$file_prefix
) || pod2usage( -verbose => 2 );

my @files = @ARGV;

if ($help) {
    pod2usage( -verbose => 2 );
}
elsif ( scalar(@files) < 2 ) {
    die "Clustal requires at least 2 gene inputs for an alignment. Run with "
      . "--help option for usage information.\n";
}
elsif ( not @files ) {
    die "No files were provided. Run with --help option for usage "
      . "information\n";
}

$file_prefix = $file_prefix || "aln";

#
# make the output directory if possible, die if not; define file path
#

unless ( -e uc($file_prefix) || mkdir( uc($file_prefix) ) ) {
    die "Could not create directory: " . uc($file_prefix);
}

my $aln_path = uc($file_prefix) . "/$file_prefix";

#
# get the AA seqs, framelines, and names for each gene
#

my @aa_seqs;
my @framelines;
my @gene_names;
my $invalid;

my $codon_cable = Bio::Tools::CodonTable->new();

foreach my $file (@files) {
    my $gene_name;
    my $ex_seqs;

    while ( not $gene_name ) {
        print "\nEnter a gene name for file $file: ";
        $gene_name = <STDIN>;
        $gene_name =~ s/\s//g;
        print "  ERROR: no name given - try again\n" unless $gene_name;
    }

    open( my $infile, '<', $file ) || die "Can't open $file: $!\n";
    my @lines = <$infile>;
    close($infile);

    ( $invalid, $ex_seqs ) = exon_split(@lines);

    my @seqs = @{$ex_seqs};

    if ($invalid) {
        print "  WARNING: invalid characters in exon sequence were removed\n";
    }

    my $dna_seq = join( '', @seqs );

    my $aa_seq = $codon_cable->translate($dna_seq);
    $aa_seq =~ s/\*.*$//g;

    my $frameline = frame_line( $dna_seq, $aa_seq, @seqs );

    push( @aa_seqs,    $aa_seq );
    push( @framelines, $frameline );
    push( @gene_names, $gene_name );
}

#
# generate a FASTA file of all protein seqs, send it to clustal
#

print "\n";

my @fasta_output = fasta_format( \@aa_seqs, \@gene_names );

open( my $output, '>', "$aln_path.sequence.txt" )
  || die "Can't open $aln_path.sequence.txt: $!\n";
print $output @fasta_output;
close($output);

my $email;

while ( not $email ) {
    print "\nEnter an email address for Clustal Omega (you will receive your "
      . "results immediately, but Clustal requires an email address be given): ";
    $email = <STDIN>;
    $email =~ s/\s//g;
    print "  ERROR: no email given - try again\n" unless $email;
}

#
# runs Clustal twice because BOXSHADE needs the no numbers version, but people
# like the version WITH numbers. so get both.
#

print "\n==== RUNNING CLUSTAL ====\n";

system(
    "perl",                   "$dir/clustalo_soaplite.pl",
    "$aln_path.sequence.txt", "--outfile",
    $aln_path,                "--email",
    $email,                   "--outfmt",
    "clustal",                "--quiet",
    "--noguidetreeout",       "--nodismatout"
);

system(
    "perl",                   "$dir/clustalo_soaplite.pl",
    "$aln_path.sequence.txt", "--outfile",
    $aln_path,                "--email",
    $email,                   "--outfmt",
    "clustal_num",            "--noguidetreeout",
    "--nodismatout"
);

print "Creating result file: $aln_path.aln-clustal.clustal\n";

open( my $input, '<', "$aln_path.aln-clustal_num.clustal_num" )
  || die "Can't open $aln_path.aln-clustal_num.clustal_num: $!\n";
my @clustal_in = <$input>;
close($input);

print "Creating output file: $aln_path" . "_frames.txt\n";

my ( $clustal, $boxframes ) = clustal_frames( \@clustal_in, \@framelines );

my @clustal_out = @{$clustal};
my @box_frames  = @{$boxframes};

open( $output, '>', $aln_path . '_frames.txt' )
  || die "Can't open file $aln_path" . "_frames.txt: $!\n";
print $output @clustal_out;
close($output);

print "===== CLUSTAL DONE ======\n\n";

#
# now we take the outputs from clustal_frames() and send them to BOXSHADE
#

my $box_url = 'http://ch.embnet.org/software/BOX_form.html';

my $box_mech = WWW::Mechanize->new();
$box_mech->agent_alias('Windows Mozilla');

print "=== RUNNING BOXSHADE ====\n";
print "FETCHING PAGE\n";

$box_mech->get($box_url);

print "SETTING OPTIONS\n";

open( $input, '<', "$aln_path.aln-clustal.clustal" )
  || die "Can't open $aln_path.aln-clustal.clustal: $!\n";
my @clustaln = <$input>;
close($input);

my $clust_box = join( '', @clustaln );

$box_mech->select( 'outp',   'RTF_new' );
$box_mech->select( 'inform', 'ALN' );
$box_mech->field( 'seq', $clust_box );

print "RUNNING\n";

$box_mech->click_button( value => "Run BOXSHADE..." );

print "Creating result file: $aln_path" . "_boxshade.rtf\n";

$box_mech->follow_link( text_regex => qr/here is your output/i );
$box_mech->save_content( "$aln_path" . "_boxshade.rtf" );

print "Creating output file: $aln_path" . "_boxshade_frames.rtf\n";

open( $input, '<', $aln_path . "_boxshade.rtf" )
  || die "Can't open $aln_path" . "_boxshade.rtf: $!\n";
my @boxshade_in = <$input>;
close($input);

my @boxshade_out = boxshade_frames( \@boxshade_in, \@box_frames );

open( $output, '>', $aln_path . "_boxshade_frames.rtf" )
  || die "Can't open $aln_path" . "_boxshade_frames.rtf: $!\n";
print $output @boxshade_out;
close($output);

my $clustaln_file = $aln_path . '_frames.txt';
my $boxshade_file = $aln_path . '_boxshade.rtf';
my $boxframe_file = $aln_path . '_boxshade_frames.rtf';

print <<"END";
===== BOXSHADE DONE =====

Done!

The following are some of the useful files that were produced:

1. The numbered Clustal alignment is found in $aln_path.aln-clustal_num.clustal_num
2. The version with frame markers is found in $clustaln_file
3. A numberless version without markers is found in $aln_path.aln-clustal.clustal
4. The BOXSHADE alignment is found in $boxshade_file
5. The BOXSHADE alignment with frame markers is found in $boxframe_file

END

#-------------------- Pod Documentation --------------------

=head1 NAME

Usage of frames.pl:

=head1 SYNOPSIS

perl frames.pl [ -d DIR ] file1 file2 [ file3 ... ]

=head2 Options:

    --help              Will show this help page.

    -d          DIR     The value of DIR will be the directory where all
                        files produced will be located, and will be the
                        prefix in all of their names.
                        Default: ALN

=head1 DESCRIPTION

frames.pl accepts any number (greater than 0) of input files, each of which should contain the sequences of all exons for a given gene. Each exon should begin with a FASTA identification line, with no extra lines in the file. In other words, for each gene, you want a file that looks something like this, with only FASTA lines and exon sequence data:

    >exon1
    ATGGACTG...................................................
    >exon2
    TACAT......................................................
    >exon3
    TATTACGACGTA...............................................

It doesn't matter what the > id line says, as it will be stripped off, but each exon should have one to separate them from one another. frames.pl ignores lower-case letters in the sequences, so keep that in mind when designing the input files. frames.pl will then ask the user for a name for each gene. These names will be the labels on the left side of the alignment files generated.

The DIR value supplied to the -d option will be the name of the folder where all the files are located, and the prefix for all the file names. The default value if you don't supply one is "aln".

The markers generated by `frames.pl` are `<`, `>`, and `^`, which indicate how the codon for the corresponding amino acid relates to an exon border at that point:

    <   1 nucleotide left of border, 2 nt right of border
    >   2 nt left of border, 1 nt right of border
    ^   0 left, 3 right (i.e. it's in-frame; the codon begins where the exon does)

=head3 Clustal

frames.pl will send the genes you supplied to Clustal Omega at:

    https://www.ebi.ac.uk/Tools/msa/clustalo/

This will produce a number of files, including (assuming the -d option was not given and all files are prefixed with "aln"):

    aln.aln-clustal.clustal         <- the alignment with no line numbers
    aln.aln-clustal_num.clustal_num <- the alignment with line numbers
    aln_frames.txt                  <- with numbers and frame markers
    aln_box_frames.txt              <- just frame markers, no alignment

=head3 BOXSHADE

Next, frames.pl will send aln.aln-clustal.clustal to BOXSHADE at:

    http://www.ch.embnet.org/software/BOX_form.html

This will produce two new files:

    aln_boxshade.rtf                <- the alignment boxshaded
    aln_boxxshade_frames.rtf        <- the BOXSHADE version with frame markers

There will be several additional files produced by Clustal Omega and frames.pl which may be of some interest, but are not the focus of this documentation. Feel free to browse them.

=head1 DEPENDENCIES

    SOAP::Lite
    WWW::Mechanize

=head1 SEE ALSO

This is a command line client for the (hopefully) included ExonFrames
module. It should be located, relative to this file, in the ./lib
directory, and contains further documentation on the functions used.

=head1 AUTHOR

Andrew Trivette - adt dot pseudologic at gmail dot com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Andrew Trivette - adt dot pseudologic at gmail dot com

This module is free software; you can redistribute it and/or modify it under
the terms of the Artistic License 2.0. For details, see the full text of the
license in the file LICENSE.md.

This program is distributed in the hope that it will be
useful, but it is provided “as is” and without any express
or implied warranties. For details, see the full text of
the license in the file LICENSE.md.

=cut

