# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl ExonFrames.t'

#########################

use strict;
use warnings;

use Test::More tests => 22;
use File::Basename;
use lib dirname(__FILE__) . '/../lib';

BEGIN { use_ok('ExonFrames') };

#########################

use ExonFrames ':all';

##### Test that exons are being extracted properly

my @input = exon_split(

    ( ">this is line 1\n",
      "ACTAGCCATTACGAT\n",
      ">this is line 2\n",
      "TACCAT\n",
      "TACCAGGATTAAAC\n",
      ">this is line 3\n",
      "AAAAAAAAAAAAAAACCCCCTA\r",
      ">this is line 4\n",
      "abcdACT   AGaggATC566\n",
      "\n\n\n")

);

my @output = ( "ACTAGCCATTACGAT",
               "TACCATTACCAGGATTAAAC",
               "AAAAAAAAAAAAAAACCCCCTA",
               "ACTAGATC" );

is_deeply( \@input, \@output, 'exon_split() various' );

##### Test that translation is correct

my @dna_seqs = ( 'ATGAACTTGAGATAGCCA',
                 'TTAACGAGAGTACCAGTA',
                 'TACTACCATTAGTACAAGTAC',
                 'ATACATGAGGACCTAGA',
                 'ACTACCTCACGGTATTGGGT',
                 'CGTCGCTCCTCAACCTTAT');

my @aa_seqs = ( 'MNLR',
                'LTRVPV',
                'YYH',
                'IHEDL',
                'TTSRYW',
                'RRSSTL');

for ( my $i = 0; $i < scalar(@dna_seqs); ++$i ) {
    is( translate($dna_seqs[$i]), $aa_seqs[$i], "translate() $aa_seqs[$i]" );
}

##### Test frame line generation

my @exons = (

    [ "ACTAGCCATTACGA",
      "TTACCATTACCAGGATTATA",
      "CGAAAAAAAAAAAAAAACCCCCT",
      "ACTAGATCA",
    ],

    [ "ACTAGCCATTACGA",
      "TTACCATTACCAGGATTATA",
      "CGAAAAAAAAAAAAAAACCCCCT",
      "ACTAGATCAT",
    ],
    [ "ACTAGCCATTACGA",
      "TTACCATTACCAGGATTATA",
      "CGAAAAAAAAAAAAAAACCCCCT",
      "ACTAGATCATC",
    ],
    [ "ATGACCTATCAGACCAGTC" ],
);

my @framelines = ( "^   >      <       ^  ",
                   "^   >      <       ^  ",
                   "^   >      <       ^  ",
                   "^     " );

for ( my $i = 0; $i < scalar(@exons); ++$i ) {
    my @exon_array = @{ $exons[$i] };
    my $dna_seq = join( '', @exon_array );
    my $aa_seq = translate($dna_seq);

    is( frame_line($dna_seq, $aa_seq, @exon_array), $framelines[$i],
        "frame_line() $framelines[$i]" );
}

##### Test the fasta wrapping

my $aa_seq1 = "RSNTVWSMFESFS";
my $aa_seq2 = "MPVVNHEDSEFHLSHTEEDKLNEFQVITNFPPEDLPDVVRLLRNHGWQLEPALSRYFDGEWKGEPDQMGEPTQTSTPMAETLVPPALGPRPLLFTASLPVVRPLPANFRNDFRTIGLNGRSNTVWSMFESFSYDGNPFLFILLLIPRIINRLSATIFTFFCTLLSLHSISGGGNSGKPKISKVPKAPTRETHIPLAEILGDTKDKDAFCELKSFKPDISFNEALRIAKEEFKFMLLILVGDTYDTDTDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLGVRNTPEIFLIANVANKASHSETLPSQRLSILGKLKVNSLNRFLQSLTNVVEKYTPELVVNKTEMHELRMSREIKKLQEDAYKKSLEMDRIKAIEKEKSLKHAQDLKLNSTARQLKWLKACIDEIQPFETTGKQATLQFRTSSGKRFVKKFPSMTTLYQIYQSIGCHIYLAVYSSDPAEWSNALQDKIRQLSADDDMLCFKEGQLETATATTIEELGHIINNELTSFDLERGKLEFDFELVSPFPKYTVHPNEHMSVDQVPQLWPNGSLLVEALDEEDEEDEENEEQ";

my @fasta1 = line_wrap(   60, $aa_seq1 );
my @fasta2 = line_wrap(   60, $aa_seq2 );
my @fasta3 = line_wrap(    4, $aa_seq1 );
my @fasta4 = line_wrap(    1, $aa_seq1 );
my @fasta5 = line_wrap(    0, $aa_seq2 );
my @fasta6 = line_wrap( 1000, $aa_seq1 . $aa_seq2 );

my @wrap1 = ( "RSNTVWSMFESFS\n", );
my @wrap2 = ( "MPVVNHEDSEFHLSHTEEDKLNEFQVITNFPPEDLPDVVRLLRNHGWQLEPALSRYFDGE\n",
              "WKGEPDQMGEPTQTSTPMAETLVPPALGPRPLLFTASLPVVRPLPANFRNDFRTIGLNGR\n",
              "SNTVWSMFESFSYDGNPFLFILLLIPRIINRLSATIFTFFCTLLSLHSISGGGNSGKPKI\n",
              "SKVPKAPTRETHIPLAEILGDTKDKDAFCELKSFKPDISFNEALRIAKEEFKFMLLILVG\n",
              "DTYDTDTDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLGVR\n",
              "NTPEIFLIANVANKASHSETLPSQRLSILGKLKVNSLNRFLQSLTNVVEKYTPELVVNKT\n",
              "EMHELRMSREIKKLQEDAYKKSLEMDRIKAIEKEKSLKHAQDLKLNSTARQLKWLKACID\n",
              "EIQPFETTGKQATLQFRTSSGKRFVKKFPSMTTLYQIYQSIGCHIYLAVYSSDPAEWSNA\n",
              "LQDKIRQLSADDDMLCFKEGQLETATATTIEELGHIINNELTSFDLERGKLEFDFELVSP\n",
              "FPKYTVHPNEHMSVDQVPQLWPNGSLLVEALDEEDEEDEENEEQ\n" );
my @wrap3 = ( "RSNT\n",
              "VWSM\n",
              "FESF\n",
              "S\n" );
my @wrap4 = ( "R\n", "S\n", "N\n", "T\n", "V\n", "W\n", "S\n", "M\n", "F\n",
              "E\n", "S\n", "F\n", "S\n" );
my @wrap5 = ( "" );
my @wrap6 = ( "RSNTVWSMFESFSMPVVNHEDSEFHLSHTEEDKLNEFQVITNFPPEDLPDVVRLLRNHGWQLEPALSRYFDGEWKGEPDQMGEPTQTSTPMAETLVPPALGPRPLLFTASLPVVRPLPANFRNDFRTIGLNGRSNTVWSMFESFSYDGNPFLFILLLIPRIINRLSATIFTFFCTLLSLHSISGGGNSGKPKISKVPKAPTRETHIPLAEILGDTKDKDAFCELKSFKPDISFNEALRIAKEEFKFMLLILVGDTYDTDTDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLGVRNTPEIFLIANVANKASHSETLPSQRLSILGKLKVNSLNRFLQSLTNVVEKYTPELVVNKTEMHELRMSREIKKLQEDAYKKSLEMDRIKAIEKEKSLKHAQDLKLNSTARQLKWLKACIDEIQPFETTGKQATLQFRTSSGKRFVKKFPSMTTLYQIYQSIGCHIYLAVYSSDPAEWSNALQDKIRQLSADDDMLCFKEGQLETATATTIEELGHIINNELTSFDLERGKLEFDFELVSPFPKYTVHPNEHMSVDQVPQLWPNGSLLVEALDEEDEEDEENEEQ\n" );

is_deeply( \@fasta1, \@wrap1, 'fasta_wrap() item < wrap' );
is_deeply( \@fasta2, \@wrap2, 'fasta_wrap() standard' );
is_deeply( \@fasta3, \@wrap3, 'fasta_wrap() short wrap' );
is_deeply( \@fasta4, \@wrap4, 'fasta_wrap() 1 wrap' );
is_deeply( \@fasta5, \@wrap5, 'fasta_wrap() 0 wrap' );
is_deeply( \@fasta6, \@wrap6, 'fasta_wrap() 1000 wrap' );

##### Test fasta formatting of clustal input file

my @fasta_aa = (
    'MRTNADIYFL',

    'ALKYLYDETDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLG'
    . 'VRNTPEIFLIANVANKASHSETLPSQRLSILGKLK',

    'VGDTYDTDTDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLGGARGGARG'
    . 'VRNTPEIFLIANVANKASHSETLPSQRLSILGKLKVNSLNRFLQSLTNVVEKYTPELV'
    . 'VNKTEMHELRMSREIKKLQEDAYKKSLEMDRIKAIEKEKSLKHAQDLKLNSTARQ'
);

my @fasta_gn = ( 'gene1', 'gene2', 'gene3' );

my @fasta_out = (
    ">gene1\n",
    "MRTNADIYFL\n",
    ">gene2\n",
    "ALKYLYDETDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLG\n",
    "VRNTPEIFLIANVANKASHSETLPSQRLSILGKLK\n",
    ">gene3\n",
    "VGDTYDTDTDTVDVNSKLLLEKILLNKKTLQYLRKIDNDLIIYLKCVHELEPWLVARQLG\n",
    "GARGGARGVRNTPEIFLIANVANKASHSETLPSQRLSILGKLKVNSLNRFLQSLTNVVEK\n",
    "YTPELVVNKTEMHELRMSREIKKLQEDAYKKSLEMDRIKAIEKEKSLKHAQDLKLNSTAR\n",
    "Q\n"
);

my @fasta_test = fasta_format( \@fasta_aa, \@fasta_gn );

is_deeply( \@fasta_test, \@fasta_out, 'fasta_format() short test' );

##### Test that frames are put into clustal alignment properly
##### and that the boxframes are being made correctly

my $path = '/home/andrew/school/bioinformatics/ExonFrames/t/';

open( my $clustal, '<', "$path" . 'FULL/full.aln-clustal_num.clustal_num' ) ||
    die "ERROR: Can't open clustal_num file: $!\n";
my @clustal_lines = <$clustal>;
close($clustal);

open( my $frames, '<', "$path" . 'FULL/full_orig_frames.txt' ) ||
    die "ERROR: Can't open orig: $!\n";
my @frame_lines = <$frames>;
chomp @frame_lines;
close($frames);

my ( $clust, $boxframelines ) = clustal_frames( \@clustal_lines, \@frame_lines );

my @test_lines = @{$clust};
my @box_framelines = @{$boxframelines};

open( my $clustal_wframes, '<', "$path" . "FULL/full_frames.txt" ) ||
    die "ERROR: Can't open full_frames.txt: $!\n";
@clustal_lines = <$clustal_wframes>;
close($clustal_wframes);

open( my $box_test, '<', "$path" . 'FULL/full_box_frames.txt' ) ||
    die "ERROR: Can't open full_box_frames.txt: $!\n";
my @box_lines = <$box_test>;
chomp @box_lines;
close($box_test);

is_deeply( \@test_lines, \@clustal_lines, 'clustal_frames() vs full' );
is_deeply( \@box_framelines, \@box_lines, 'clustal_frames() boxframes' );

##### Test that frames are put into boxshade output properly

open( my $input, '<', "$path" . 'FULL/full_boxshade.rtf' )
  || die "Can't open $path" . "FULL/full_boxshade.rtf: $!\n";
my @boxshade_lines = <$input>;
close($input);

my @boxframes = boxshade_frames( \@boxshade_lines, \@box_framelines );

open( $input, '<', "$path" . 'FULL/full_boxshade_frames.rtf' )
  || die "Can't open $path" . "FULL/full_boxshade_frames.rtf: $!\n";
my @fullframes = <$input>;
close($input);

is_deeply( \@boxframes, \@fullframes, 'boxshade_frames() vs full' );
