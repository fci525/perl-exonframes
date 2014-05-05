#!/usr/bin/env perl

package ExonFrames;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
    'all' => [
        qw(
          exon_split
          translate
          frame_line
          line_wrap
          fasta_format
          clustal_frames
          boxshade_frames
          )
    ],
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
);

our $VERSION = '0.8';

#<1>========== Subs and Things ==========

my %amino_acids = (
    A => [ "GCT", "GCC", "GCA", "GCG" ],
    C => [ "TGT", "TGC" ],
    D => [ "GAT", "GAC" ],
    E => [ "GAA", "GAG" ],
    F => [ "TTT", "TTC" ],
    G => [ "GGT", "GGC", "GGA", "GGG" ],
    H => [ "CAT", "CAC" ],
    I => [ "ATT", "ATC", "ATA" ],
    K => [ "AAA", "AAG" ],
    L => [ "CTT", "CTC", "CTA", "CTG", "TTA", "TTG" ],
    M => ["ATG"],
    N => [ "AAT", "AAC" ],
    P => [ "CCT", "CCC", "CCA", "CCG" ],
    Q => [ "CAA", "CAG" ],
    R => [ "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" ],
    S => [ "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" ],
    T => [ "ACT", "ACC", "ACA", "ACG" ],
    U => ["TGA"],
    V => [ "GTT", "GTC", "GTA", "GTG" ],
    W => ["TGG"],
    Y   => [ "TAT", "TAC" ],
    "*" => [ "TAA", "TAG", "TGA" ],
);

# takes lines from a FASTA file with exon seqs for a gene, returns an array
# of the sequences joined and stripped of the FASTA lines
sub exon_split {
    my @lines = @_;
    my @exons;

    my $in_line;
    my $cur_exon;
    my $invalid;

    foreach my $line (@lines) {
        if ( $line =~ /^>/ ) {
            if ($cur_exon) {
                push( @exons, $cur_exon );
            }
            $cur_exon = '';
        }
        elsif ( $line =~ /[A-Z]/ ) {
            $in_line = $line;

            $in_line =~ s/\s{1,2}$//g;
            $line =~ s/[^ACTG]//g;

            unless ( $in_line eq $line ) {
                $invalid = 1;
            }

            $cur_exon .= $line;
        }
    }

    push( @exons, $cur_exon );

    if ($invalid) {
        print "WARNING: invalid characters in exon sequence were removed\n";
    }

    return @exons;
}

# takes a DNA seq, returns its AA translation
sub translate {
    my $dna_seq = shift;

    my @codons;
    my $count  = 0;
    my $aa_seq = '';

  STOP: while ( $count <= length($dna_seq) - 3 ) {
        my $codon = substr( $dna_seq, $count, length($dna_seq) > 3 ? 3 : 0 );
        foreach my $stop ( @{ $amino_acids{"*"} } ) {
            if ( $codon eq $stop ) {
                last STOP;
            }
        }

        push( @codons, $codon );
        $count += 3;
    }

    foreach my $codon (@codons) {
      MATCH: foreach my $acid ( keys(%amino_acids) ) {
            foreach my $alt ( @{ $amino_acids{$acid} } ) {
                if ( $codon eq $alt ) {
                    $aa_seq .= $acid;
                    last MATCH;
                }
            }
        }
    }

    return $aa_seq;
}

# given a DNA sequence, its translation, and the sequences of its exons,
# this returns a line with markers (<, >, or ^) indicating how the codons
# line up with the exon borders.
sub frame_line {
    my ( $dna_seq, $aa_seq, @exons ) = @_;

    my %frames = (
        0 => "^",
        1 => "<",
        2 => ">"
    );

    my $gene_frame_line = " " x length($aa_seq);
    my @gene_frame_array = split( '', $gene_frame_line );

    my $dna_start;
    my $aa_start;
    my $frame;
    my $offset;

    foreach my $exon (@exons) {
        if ( $dna_seq =~ /$exon/ ) {
            $dna_start = $-[0];
            $offset    = $dna_start % 3;

            $frame = $frames{$offset};
            $dna_start -= $offset;
        }
        else {
            die "ERROR: Exon not found in sequence\n";
        }

        $aa_start = $dna_start / 3;
        $gene_frame_array[$aa_start] = $frame;
    }

    $gene_frame_line = join( '', @gene_frame_array );

    return $gene_frame_line;
}

# takes a number and a string, and returns an array with the string wrapped
# to that number of characters. designed for sequences, so it makes no attempt
# at word identification. every line gets a "\n" character at the end.
sub line_wrap {
    my ( $line_wrap, $dna_seq ) = @_;

    return ('') if $line_wrap <= 0;

    my @seq_list;

    while ( length($dna_seq) > $line_wrap ) {

        my $new_length = length($dna_seq) - $line_wrap;
        my $seq_slice = substr( $dna_seq, 0, $line_wrap );

        push( @seq_list, $seq_slice . "\n" );
        $dna_seq = substr( $dna_seq, $line_wrap, $new_length );

    }

    push( @seq_list, $dna_seq . "\n" );

    return @seq_list;
}

# takes two array refs, one with the AA seqs for a set of proteins, the other
# with the names of those genes. returns an array with lines for a FASTA output
# file with all protein sequences for Clustal.
sub fasta_format {
    my ( $aa_seq_ref, $gene_name_ref ) = @_;

    my @aa_seqs    = @{$aa_seq_ref};
    my @gene_names = @{$gene_name_ref};

    my @aa_wrap;

    my @fasta_lines;

    for ( my $i = 0 ; $i < scalar(@aa_seqs) ; ++$i ) {
        push( @fasta_lines, ">$gene_names[$i]\n" );

        @aa_wrap = line_wrap( 60, $aa_seqs[$i] );

        foreach my $line (@aa_wrap) {
            push( @fasta_lines, $line );
        }
    }

    return @fasta_lines;
}

# takes as input two array references and a file prefix. the first array
# contains the lines of a Clustal Omega alignment output file, and the second
# contains the framelines generated by &frame_line(). the file prefix will tell
# ExonFrames where to put the new alignment with framelines added.
sub clustal_frames {
    my ( $clustal_line_ref, $frame_line_ref ) = @_;

    my @clustal_lines = @{$clustal_line_ref};               # deref
    my @frame_lines   = @{$frame_line_ref};                 # deref
    my $aln_re        = qr/(^[^ ]+ +)([A-Z-]+)\s+[0-9]+$/;  # regex for AA lines

    my $line_prefix;    # number of spaces before alignment lines
    my $line_wrap;      # number of characters in alignment lines

    my $gene_cnt = 0;   # keep track of what gene we're on
    my $frame;          # the frameline after exploding, adding spaces, and join
    my $frame_seg;      # the part of $frame we actually modified

    my @ex_line;        # the exploded alignment line
    my @ex_frame;       # the exploded frameline

    my @output;         # the clustal alignment + frame markers
    my @box_frames;     # just frame markers, for boxshade

    # initialize @box_frames with empty strings for each frame line
    foreach my $line (@frame_lines) {
        push( @box_frames, '' );
    }

    # the first 3 lines are always blank, so just take them
    for ( my $i = 0 ; $i < 3 ; ++$i ) {
        push( @output, shift(@clustal_lines) );
    }

    foreach my $line (@clustal_lines) {

        # if the line contains AAs or dashes
        if ( $line =~ /$aln_re/ ) {

            # get the $line_prefix and $line_wrap values ($line_wrap is
            # acquired every time because the last line is always shorter,
            # and causes loop problems otherwise
            $line_prefix = length($1);
            $line_wrap   = length($2);

            # explode the lines to loop over them
            @ex_line  = split( '', $2 );
            @ex_frame = split( '', $frame_lines[$gene_cnt] );

            # for every dash in the AA line, put a space in the frameline
            for ( my $i = 0 ; $i < scalar(@ex_line) ; ++$i ) {
                if ( $ex_line[$i] eq '-' ) {
                    splice( @ex_frame, $i, 0, ' ' );
                }
            }

            $frame = join( '', @ex_frame );

            # take the modified part of $frame, or the whole last line's-worth
            if ( length($frame) >= $line_wrap ) {
                $frame_seg = substr( $frame, 0, $line_wrap );
            }
            else {
                $frame_seg = $frame;
            }
            $box_frames[$gene_cnt] .= $frame_seg;

            # push the lines to the output array
            push( @output, $line );
            push( @output, ( ' ' x $line_prefix ) . $frame_seg . "\n" );

            # and remove the already done part from the frameline
            $frame_lines[$gene_cnt] = substr( $frame, $line_wrap );

            $gene_cnt += 1;
        }
        else {
            # just take the blank lines, and reset $gene_cnt
            push( @output, $line );
            $gene_cnt = 0;
        }
    }

    # and return the boxframes
    return ( \@output, \@box_frames );
}

sub boxshade_frames {
    my ( $boxshade_line_ref, $box_frame_ref ) = @_;

    my @boxshade_lines = @{$boxshade_line_ref};    # deref
    my @box_frames     = @{$box_frame_ref};        # deref

    my @output;    # will store lines for the output file

    my @boxframe_array;    # array of refs to "@frameline_array"s

    # regex vars for various repeated things
    my $aa_line   = qr/(\\highlight.\\cf. )([A-Z-]+)$/;
    my $pre_line  = qr/(\\highlight.\\cf. )([^ ]+ +[0-9]* )$/;
    my $lbrk_line = qr/\\highlight.\\cf. \\line$/;
    my $pbrk_line = qr/^\\page$/;

    my $cur_gene;          # used to pull out the frame array for a gene
    my @cur_gene_array;    # deref that array
    my $cur_line;          # and get the right line from it

    my $frame_prefix;      # prefix with coloring for framelines
    my $line_prefix;       # prefix with coloring for line breaks
    my $box_wrap;          # the line wrap that boxshade picked
    my $page_break;        # how many groups fit on a page (if any)
    my $tot_line_cnt;      # how many lines in each @frameline_array

    my $break_now = 0;     # indicates when to break the page (1) or not (0)
    my $prev_line = 0;     # marks what kind of data was on the last line
    my $line_cnt  = 0;     # which line of each frameline array we're on
    my $gene_cnt  = 0;     # which gene we're on

    my $grp_lines = ( scalar(@box_frames) * 2 ) + 2;

    # all used to calculate how many groups of genes fit on one page
    my $font_size;         # font size
    my $points;            # points (font spec) that a group would take up
    my $inches;            # points converted to inches
    my $margins;           # top and bottom margins added together
    my $space;             # space available on a page in inches

    # loop through every line of the boxshade file and pick out format stuff
    foreach my $line (@boxshade_lines) {
        if ( $line =~ /\\margt([0-9]+)\\margb([0-9]+)\\/ ) {
            $margins = $1 + $2;
        }
        if ( $line =~ /\\fs([0-9]+)$/ ) {
            $font_size = $1 / 2;
        }
        if ( $line =~ /$pre_line/ ) {
            $line_prefix = $1;
            $frame_prefix = $1 . ( ' ' x length($2) );
        }
        if ( $line =~ /$aa_line/ ) {
            $box_wrap += length($2);
        }
        if ( $line =~ /$lbrk_line/ && $box_wrap ) {
            last;
        }
    }

    # how much space will a group take up?
    $points = $font_size * $grp_lines;
    $inches = $points * 0.0138;

    # how much space is on a page?
    $margins = ( $margins * 0.0006944444 ) + 0.5;
    $space   = 11 - $margins;

    # unless 1 block of genes is too big for a page, we turn on page_break, which
    # says the number of blocks which can fully fit on one page.
    if ( not $inches > $space ) {
        $page_break = int( $space / $inches );
    }

    # line_wrap() all the framelines, and send that array to @boxframe_array
    foreach my $frameline (@box_frames) {
        my @frameline_array = line_wrap( $box_wrap, $frameline );
        push( @boxframe_array, \@frameline_array );
        $tot_line_cnt = scalar(@frameline_array);
    }

    # just take the first 9 lines; they're for formatting etc.
    for ( my $i = 0 ; $i < 9 ; ++$i ) {
        push( @output, shift(@boxshade_lines) );
    }

    foreach my $line (@boxshade_lines) {

        # if the line has a gene name (prefix) on it, just take it
        if ( $line =~ /$pre_line/ ) {
            push( @output, $line );
            $prev_line = 0;    # this was a prefix
        }

        # if it has amino acids on it, just take it
        elsif ( $line =~ /$aa_line/ ) {
            push( @output, $line );
            $prev_line = 1;    # this was AAs
        }

        # if it has a .rtf linebreak, we need to do stuff
        elsif ( $line =~ /$lbrk_line/ ) {

            # if the previous line was amino acids
            if ( $prev_line == 1 ) {

                # deref the proper gene frameline array and print the right line
                $cur_gene       = $boxframe_array[$gene_cnt];
                @cur_gene_array = @{$cur_gene};
                $cur_line       = $cur_gene_array[$line_cnt];

                push( @output, $line );
                push( @output, "$frame_prefix$cur_line" );
                push( @output, "$line_prefix\\line\n" );

                $gene_cnt += 1;    # next gene
                $prev_line = 2;    # so we know this was a linebreak
            }

            # if the previous line was another linebreak
            elsif ( $prev_line == 2 ) {

                # page break if necessary
                if ( $page_break && $break_now ) {
                    push( @output, "\\page\n" );
                    $break_now = 0;
                }

                # otherwise, newlines unless this is the end of the file
                elsif ( not( ( $line_cnt + 1 ) >= $tot_line_cnt ) ) {
                    push( @output, "$line_prefix\\line\n" );
                    push( @output, "$line_prefix\\line\n" );
                }

                $gene_cnt = 0;    # back to first gene
                $line_cnt += 1;   # next line for frameline arrays
                $prev_line = 3;   # this was a 'middle of nowhere' break
            }
        }

        # if it was a pagebreak instead
        elsif ( $line =~ /$pbrk_line/ ) {

            # let it break if we were due for one
            if ( $page_break && $break_now ) {
                $break_now = 0;
                push( @output, $line );
            }

            # otherwise, newlines unless this is the end of the file
            elsif ( not( ( $line_cnt + 1 ) >= $tot_line_cnt ) ) {
                push( @output, "$line_prefix\\line\n" );
                push( @output, "$line_prefix\\line\n" );
            }

            $gene_cnt = 0;    # back to first gene
            $line_cnt += 1;   # next line for frameline arrays
            $prev_line = 4;   # this was a pagebreak
        }

        # if it was just a \n and nothing else, take it
        else {
            push( @output, $line );
        }

        # turn on $break_now if no more blocks will fit on the page UNLESS
        # this is the end of the file, where we don't want a pagebreak
        if ( ( ( $line_cnt + 1 ) % $page_break ) == 0 ) {
            $break_now = ( $line_cnt + 1 ) >= $tot_line_cnt ? 0 : 1;
        }
        else {
            $break_now = 0;
        }
    }

    return @output;
}

__END__

#<3>========== Pod Documentation ==========

=head1 NAME

ExonFrames - a set of functions to notate the reading frames of codons in protein multiple alignments relative to exon borders

=head1 SYNOPSIS

Use frames.pl, included with this module. Run it with the --help option for usage information.

    perl frames.pl --help

=head1 DESCRIPTION

This module contains various functions to notate the arrangement of the reading frames of the exons for a set of genes, and display that information in a Clustal Omega alignment of the resulting proteins. The functions are as follows:

=head2 exon_split()

This takes as input an @array, meant to be the lines of an input file containing all exons for a given gene in FASTA format. It also returns an array, which contains just the sequences from the input, stripped of their FASTA lines and line wrapping. This funciton will ignore lower-case letters and any base which is not A, C, T, or G. Ambiguity code is not handled.

    my @exon_seqs = exon_split(@exon_file);

=head2 translate()

Translates a DNA sequence into an amino acid sequence. It does not try different reading frames or accept non-ACTG bases. It is assumed that the input sequence is just a join() of what came out of the exon_split() sub, as it is used in the frames.pl client.

    my $aa_seq = translate($dna_seq);

=head2 frame_line()

Takes as input a DNA sequence, a corresponding amino acid sequence (probably using the above translate()) and an array of exon sequences (generated with exon_split()) for the same gene. Outputs a string marking the reading frames for the exons in the gene. For example:

exon sequences:

    ATGCTTG ATTT CTCT AGACCT

joined into dna_seq with translation below:

    ATGCTTGATTTCTCTAGACCT
     M  L  D  F  S  R  P

But since the exons were not all divided into 3-base codons evenly, the '<', '>', and '^' symbols mark how the exon borders relate to the actual codons for the amino acids, e.g.:

    ATGCTTGATTTCTCTAGACCT
     M  L  D  F  S  R  P
     ^     <  >     ^

So the '^' means that it's in frame, the '<' means the codon had 1 base left of the border and 2 right, and '>' means 2 left 1 right. The actual translation and frameline are not spaced out to match codons, so it would be:

    my @exon_seqs = exon_split(@exon_file);
    my $dna_seq   = join( '', @exon_seqs );
    my $aa_seq    = translate($dna_seq);

    my $frameline = frame_line( $dna_seq, $aa_seq, @exon_seqs );

And you get:

    $aa_seq:    'MLDFSRP'
    $frameline: '^ <> ^ '

=head2 line_wrap()

Takes a string as input and simply wraps it to a given length, returning an @array of lines with "\n" characters at the end of each.

    my @lines = line_wrap( $wrap_length, $string );

=head2 fasta_format()

This is used to produce a FASTA file of several proteins to send to Clustal Omega for alignment. As input, it takes two \@array_refs. The first is to an array of the amino acid seqs for all genes to be aligned, and the second is an array of the names for those genes.The names must be acquired separately. frames.pl just asks the user for them.

As output, it returns an @array of lines to be stored to the output file.

    my @fasta_lines = fasta_format( \@aa_seqs, \@gene_names );

=head2 clustal_frames()

Once a Clustal alignment has been generated, this sub inserts the framelines into that alignment so that one can see not only the sequence homology, but how the exon borders for the aligned genes are similar or different.

As input, it takes two \@array_refs. The first is to an @array of the lines from the Clustal alignment, and the second @array contains the framelines for each gene that was aligned, in the same order. As output, you get two more \@array_refs. The first contains the lines to be output into the new alignment file, and the second contains framelines that have been spaced out to match the alignment. This second array is to be used with the next function, boxshade_frames().

    my ( $c_ref, $b_ref ) = clustal_frames( \@c_lines, \@framelines );

=head2 boxshade_frames()

This one does much the same as clustal_frames(), except that it's adding the framelines to a BOXSHADEd version of the Clustal alignment. The BOXSHADE version is a .rtf file, though, to allow for shading. 

As input, it takes two \@array_refs. The first contains the lines of the BOXSHADE output .rtf file, and the second contains the second array from clustal_frames() (@framelines in the example above).

As output, you get an @array of lines to be printed to a new output .rtf file, which has been reformatted so that the page breaks make sense with the added framelines.

    my @box_out_lines = boxshade_frames( \@box_lines, \@framelines );

=head1 DEPENDENCIES

This module itself does not require any non-core modules. The client included, however, does. frames.pl requires the following:

    SOAP::Lite
    WWW::Mechanize

As well as the (INCLUDED) Clustal Omega client provided by EMBL-EBI:

    http://www.ebi.ac.uk/Tools/webservices/download_clients/perl/soaplite/clustalo_soaplite.pl

=head1 AUTHOR

Andrew Trivette, adt.pseudologic@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Andrew Trivette

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself, either Perl version 5.18 or, at your option, any later version of Perl 5 you may have available.

=cut

