use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'ExonFrames',
    VERSION_FROM      => 'lib/ExonFrames.pm', # finds $VERSION
    PREREQ_PM         => {SOAP::Lite,
                          WWW::Mechanize}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/ExonFrames.pm', # retrieve abstract from module
       AUTHOR         => 'Andrew <adt.pseudologic@gmail.com>') : ()),
);
