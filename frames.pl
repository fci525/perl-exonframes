#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use WWW::Mechanize;

use lib dirname(__FILE__) . '/lib';
use ExonFrames ':all';

my $dir = dirname(__FILE__);

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

unless ( -e uc($file_prefix) || mkdir( uc($file_prefix) ) ) {
    die "Could not create directory: " . uc($file_prefix);
}

my $aln_path = uc($file_prefix) . "/$file_prefix";

my @aa_seqs;
my @framelines;
my @gene_names;

foreach my $file (@files) {
    print "\nEnter a gene name for file $file: ";
    my $gene_name = <STDIN>;
    $gene_name =~ s/\s//g;

    open( my $infile, '<', $file ) || die "Can't open $file: $!\n";
    my @lines = <$infile>;
    close($infile);

    @lines = exon_split(@lines);
    my $dna_seq = join( '', @lines );
    my $aa_seq = translate($dna_seq);

    my $frameline = frame_line( $dna_seq, $aa_seq, @lines );

    push( @aa_seqs,    $aa_seq );
    push( @framelines, $frameline );
    push( @gene_names, $gene_name );
}

my @fasta_output = fasta_format( \@aa_seqs, \@gene_names );

open( my $output, '>', "$aln_path.sequence.txt" )
  || die "Can't open $aln_path.sequence.txt: $!\n";
print $output @fasta_output;
close($output);

print "\nEnter an email address for Clustal Omega (you will receive your "
  . "results immediately, but Clustal requires an email address be given): ";
my $email = <STDIN>;
$email =~ s/\s//g;

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

print "===== BOXSHADE DONE =====\n";

print "\nDone!\n\n",
  "The following are some of the useful files that were produced:\n\n",
  "1. The numbered Clustal alignment is found in "
  . "$aln_path.aln-clustal_num.clustal_num\n",
  "2. The version with frame markers is found in "
  . "$aln_path"
  . "_frames.txt\n",
  "3. A numberless version without markers is found in "
  . "$aln_path.aln-clustal.clustal\n",
  "4. The BOXSHADE alignment is found in $aln_path" . "_boxshade.rtf\n",
  "5. The BOXSHADE alignment with frame markers is foun in "
  . "$aln_path"
  . "_boxshade_frames.rtf\n\n";

#-------------------- Pod Documentation --------------------

=head1 NAME

Usage of frames.pl:

=head1 SYNOPSIS

perl frames.pl [-d DIR] file1 file2 [file3 ...]

=head2 Options:

    --help              Will show this help page.

    -d          DIR     The value of DIR will be the directory where all
                        files produced will be located, and will be the
                        prefix in all of their names.
                        Default: ALN

=head1 DESCRIPTION

frames.pl accepts any number (greater than 0) of input files, each of which should contain the sequences of all exons for a given gene (inputs for the --boxshade option are different; read on). Each exon should begin with a fasta identification line, with no extra lines in the file. In other words, for each gene, you want a file that looks something like this, with only fasta lines and exon sequence data:

    >exon1
    ATGGACTG...................................................
    >exon2
    TACAT......................................................
    >exon3
    TATTACGACGTA...............................................

It doesn't matter what the > id line says, as it will be stripped off, but each exon should have one to separate them from one another. By default, frames.pl ignores lower-case letters in the sequences. You can override this behavior using the --ignorecase option. frames.pl will then ask the user for a name for each gene. When using the --noalign option, these names will be used to name the output files for each gene. Otherwise, they will be the labels on the left side of the alignment files generated.

The DIR value supplied to the -d option will be the name of the folder where all the files are located, and the prefix for all the file names. The default value if you don't supply one is "ALN".

=head3 Clustal

frames.pl will send the genes you supplied to Clustal Omega at:

    https://www.ebi.ac.uk/Tools/msa/clustalo/

This will produce a number of files, including (assuming the -d option was not given and all files are prefixed with "ALN"):

    ALN.aln-clustal.clustal         <- the alignment with no line numbers
    ALN.aln-clustal_num.clustal_num <- the alignnment with line numbers
    ALN_frames.txt                  <- with numbers and frame markers
    ALN_box_frames.txt              <- just frame markers, no alignment

=head3 BOXSHADE

Next, frames.pl will send aln.aln-clustal.clustal to BOXSHADE at:

    http://www.ch.embnet.org/software/BOX_form.html

This will produce two new files:

    ALN_boxshade.rtf                <- the alignment boxshaded
    ALN_boxxshade_frames.rtf        <- the BOXSHADE version with frame markers

There will be several additional files produced by Clustal Omega and frames.pl which may be of some interest, but are not the focus of this documentation. Feel free to browse them.

=head1 DEPENDENCIES

    SOAP::Lite
    WWW::Mechanize

=head1 SEE ALSO

This is a command line client for the (hopefully) included ExonFrames
module. It should be located, relative to this file, in the ./lib
directory, and contains further documentation on the functions used.

=head1 AUTHOR

Andrew Trivette, adt.pseudologic@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Andrew Trivette

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.18.2 or,
at your option, any later version of Perl 5 you may have available.

=cut

