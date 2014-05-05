ExonFrames
==========

ExonFrames is a set of functions which can notate the reading frames
of codons relative to exon borders in a multiple alignment of some number
of proteins. Included is `frames.pl`, a small client which puts this into
practice, as well as some sample gene files and sample outputs.

The files in the `exons/` directory each contain all of the exons for
a human gene, named in the filename, in FASTA format. They can be
used as input files for `frames.pl`.

The files in the `TEST/` directory are an example of the output you
would get from `frames.pl` by running:

    perl frames.pl -d test exons/*

To display the same usage information listed below, use:

    perl frames.pl --help

You can also try running `perldoc` on the `ExonFrames.pm` file found
in the `lib/` directory for more information about the functions in
ExonFrames:

    perldoc ExonFrames.pm

### INSTALLATION

For now don't try to install this as a perl module. You can just download
the files using the `Download ZIP` button on the right and run them from
the unzipped folder.

### USAGE (of `frames.pl`)

    perl frames.pl [-d DIR] file1 file2 [file3 ...]

#### OPTIONS

    --help              Will show this help page.

    -d          DIR     The value of DIR will be the directory where all
                        files produced will be located, and will be the
                        prefix in all of their names.
                        Default: ALN

#### DESCRIPTION

`frames.pl` accepts any number (greater than 0) of input files, each of which should contain the sequences of all exons for a given gene. Each exon should begin with a FASTA identification line, with no extra lines in the file. In other words, for each gene, you want a file that looks something like this, with only FASTA lines and exon sequence data:

    >exon1
    ATGGACTG...................................................
    >exon2
    TACAT......................................................
    >exon3
    TATTACGACGTA...............................................

It doesn't matter what the `>exon whatever` id line says, as it will be stripped off, but each exon should have one to separate them from one another. `frames.pl` ignores lower-case letters in the sequences, so keep that in mind when creating the input files. `frames.pl` will then ask the user for a name for each gene. These names will be the labels on the left side of the alignment files generated.

The `DIR` value supplied to the `-d` option will be the name of the folder where all the files are located, and the prefix for all the file names. The default value if you don't supply one is "aln".

The markers generated by `frames.pl` are `<`, `>`, and `^`, which indicate how the codon for the corresponding amino acid relates to an exon border at that point:

    <   1 nucleotide left of border, 2 nt right of border
    >   2 nt left of border, 1 nt right of border
    ^   0 left, 3 right (i.e. it's in-frame; the codon begins where the exon does)

##### CLUSTAL

`frames.pl` will send the genes you supplied to Clustal Omega at https://www.ebi.ac.uk/Tools/msa/clustalo/.

This will produce a number of files, including (assuming the `-d` option was not given and all files are prefixed with "aln"):

    aln.aln-clustal.clustal         <- the alignment with no line numbers
    aln.aln-clustal_num.clustal_num <- the alignnment with line numbers
    aln_frames.txt                  <- with numbers and frame markers
    aln_box_frames.txt              <- just frame markers, no alignment

##### BOXSHADE

Next, `frames.pl` will send `aln.aln-clustal.clustal` to BOXSHADE at http://www.ch.embnet.org/software/BOX_form.html.

This will produce two new files:

    aln_boxshade.rtf                <- the alignment boxshaded
    aln_boxxshade_frames.rtf        <- the BOXSHADE version with frame markers

There will be several additional files produced by Clustal Omega and `frames.pl` which may be of some interest, but are not the focus of this documentation. Feel free to browse the examples in the `TEST/` directory.

### DEPENDENCIES

The ExonFrames module does not require any non-core modules, but the included
frames.pl client does. They are these:

    SOAP::Lite
    WWW::Mechanize

### COPYRIGHT AND LICENCE

Copyright (C) 2014 by Andrew Trivette - adt dot pseudologic at gmail dot com

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.18.2 or,
at your option, any later version of Perl 5 you may have available.

