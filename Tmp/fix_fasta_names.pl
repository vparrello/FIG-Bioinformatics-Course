#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use gjoseqlib;

my %name_of = map { chomp; m/^(\d+\.\d+)\t(.*)$/ ? ($1 => $2) : () } `cat $ARGV[0]`;
#die Dumper(\%name_of);
    
my ($id, $comment, $seq);
while (($id, $comment, $seq) = gjoseqlib::next_fasta()) {
    #die Dumper($id, $comment, $seq);
    my $name;
    if (($id =~ m/^fig\|(\d+\.\d+)/) && defined($name = $name_of{$1})) {
	print STDOUT (qq(>$id $name\n), $seq, qq(\n));
    }
}
