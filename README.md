ExonFrames
==========

ExonFrames is a set of functions which can notate the reading frames
of exons in a multiple alignment of some number of genes. Included
is `frames.pl`, a small client which puts this into practice, as well
as some sample gene files and sample outputs.

The files in the exons/ directory each contain all of the exons for
a human gene, named in the filename, in FASTA format. They can be
used as input files for frames.pl.

The files in the TEST/ directory are an example of the output you
would get from frames.pl by running

    perl frames.pl -d test exons/*

For more usage information, try:

    perl frames.pl --help

You can also try running perldoc on the ExonFrames.pm file found
in the lib/ directory:

    perldoc ExonFrames.pm

### INSTALLATION (optional)

You can install this module by typing the following:

    perl Makefile.PL
    make
    make test
    make install

Or you can just download the files and run them from the unzipped
folder.

### DEPENDENCIES

This module does not require any non-core modules, but the included
frames.pl client does. They are these:

    SOAP::Lite
    WWW::Mechanize

COPYRIGHT AND LICENCE

Copyright (C) 2014 by Andrew Trivette, adt.pseudologic@gmail.com

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.18.2 or,
at your option, any later version of Perl 5 you may have available.

