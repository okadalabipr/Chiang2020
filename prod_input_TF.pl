#!/usr/bin/perl -w
  
use strict;
use warnings;

my $line;

my @paraName = qw( head
		k1
		km1
		k2
		omI2
		omE2
		k3
		km3 );

my $paraIdx = 0;
my %paraTerm = ();
open( FHtm, "para_dep/temp_param_TF" ) or die "Cannot open temp_param_TF: $!";
while ( <FHtm> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    if ( $paraName[$paraIdx] eq "head"
      || $paraName[$paraIdx] eq "k1"
      || $paraName[$paraIdx] eq "km1"
      || $paraName[$paraIdx] eq "k3"
      || $paraName[$paraIdx] eq "km3" )
    {
	$paraTerm{$paraName[$paraIdx]} = $_;
    } else {
	$_ =~ /^(\S+\s+\S+)\s+\S+\s+\S+\s+\S+/;
	$paraTerm{$paraName[$paraIdx]} = $1;
    }

    $paraIdx++;
}
close( FHtm );

my %tfTerm = ();
$tfTerm{'k2'} = $paraTerm{'k2'};
$tfTerm{'omI2'} = $paraTerm{'omI2'};
$tfTerm{'omE2'} = $paraTerm{'omE2'};

my @tf = ();
my %tf_para = ();
open( FHtf, "para_dep/para_TF" ) or die "Cannot open para_TF: $!";
$line = <FHtf>;
while ( <FHtf> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    my @temp = split( /\t/, $_ );

    push(@tf, $temp[0]);
    $tf_para{$temp[0]}{'ka'}[0] = 'P';
    $tf_para{$temp[0]}{'ka'}[1] = $temp[1];
    $tf_para{$temp[0]}{'ka'}[2] = $temp[2];
    $tf_para{$temp[0]}{'omI'}[0] = $temp[3];
    $tf_para{$temp[0]}{'omI'}[1] = $temp[4];
    $tf_para{$temp[0]}{'omI'}[2] = $temp[5];
    $tf_para{$temp[0]}{'omE'}[0] = $temp[6];
    $tf_para{$temp[0]}{'omE'}[1] = $temp[7];
    $tf_para{$temp[0]}{'omE'}[2] = $temp[8];
}
close( FHtf );

my $fileOutput;
for ( my $i = 0; $i <= $#tf; $i++  ) {
    $paraTerm{'k2'} = $tfTerm{'k2'} . "\t" . $tf_para{$tf[$i]}{'ka'}[0] . "\t" . $tf_para{$tf[$i]}{'ka'}[1] . "\t" . $tf_para{$tf[$i]}{'ka'}[2];
    $paraTerm{'omI2'} = $tfTerm{'omI2'} . "\t" . $tf_para{$tf[$i]}{'omI'}[0] . "\t" . $tf_para{$tf[$i]}{'omI'}[1] . "\t" . $tf_para{$tf[$i]}{'omI'}[2];
    $paraTerm{'omE2'} = $tfTerm{'omE2'} . "\t" . $tf_para{$tf[$i]}{'omE'}[0] . "\t" . $tf_para{$tf[$i]}{'omE'}[1] . "\t" . $tf_para{$tf[$i]}{'omE'}[2];

    my $nline = 0;
    $fileOutput = $tf[$i] . "_param";
    open( SH, ">input_TF/" . $fileOutput ) or die "Cannot open $fileOutput: $!";
    foreach ( @paraName ) {
	if ( $nline < 1 ) {
	    print SH "$paraTerm{$_}";
	} else {
	    print SH "\n$paraTerm{$_}";
	}
	$nline++;
    }
    close( SH );
}


