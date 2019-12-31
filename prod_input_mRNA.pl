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
		km2
		k3
		km3
		k4
		h4
		k5
		k6
		km6
		k7
		omI7
		omE7
		km7
		k8
		km8
		k9
		h9
		k10
		xi4
		xi9
		k11
		k12
		k19
		k20 );

my $paraIdx = 0;
my %paraTerm = ();
open( FHtm, "para_dep/temp_param_mRNA" ) or die "Cannot open temp_param_mRNA: $!";
while ( <FHtm> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    if ( $paraName[$paraIdx] eq "head"
      || $paraName[$paraIdx] eq "k1"
      || $paraName[$paraIdx] eq "km1"
      || $paraName[$paraIdx] eq "km2"
      || $paraName[$paraIdx] eq "k3"
      || $paraName[$paraIdx] eq "km3"
      || $paraName[$paraIdx] eq "h4"
      || $paraName[$paraIdx] eq "k6"
      || $paraName[$paraIdx] eq "km6"
      || $paraName[$paraIdx] eq "km7"
      || $paraName[$paraIdx] eq "k8"
      || $paraName[$paraIdx] eq "km8"
      || $paraName[$paraIdx] eq "h9" )
    {
	$paraTerm{$paraName[$paraIdx]} = $_;
    } else {
	$_ =~ /^(\S+\s+\S+\s+\S+)\s+\S+\s+\S+/;
	$paraTerm{$paraName[$paraIdx]} = $1;
    }

    $paraIdx++;
}
close( FHtm );

my %synrTerm = ();
$synrTerm{'k4'} = $paraTerm{'k4'};
$synrTerm{'k5'} = $paraTerm{'k5'};
$synrTerm{'k9'} = $paraTerm{'k9'};
$synrTerm{'k10'} = $paraTerm{'k10'};
$synrTerm{'xi4'} = $paraTerm{'xi4'};
$synrTerm{'xi9'} = $paraTerm{'xi9'};
$synrTerm{'k11'} = $paraTerm{'k11'};
$synrTerm{'k12'} = $paraTerm{'k12'};
my %tfTerm = ();
$tfTerm{'k2'} = $paraTerm{'k2'};
$tfTerm{'omI2'} = $paraTerm{'omI2'};
$tfTerm{'omE2'} = $paraTerm{'omE2'};
$tfTerm{'k7'} = $paraTerm{'k7'};
$tfTerm{'omI7'} = $paraTerm{'omI7'};
$tfTerm{'omE7'} = $paraTerm{'omE7'};
my %mrnaTerm = ();
$mrnaTerm{'k19'} = $paraTerm{'k19'};
$mrnaTerm{'k20'} = $paraTerm{'k20'};

my $high_set_ind = 1;
my %synr_para = ();
open( FHsp, "para_dep/para_synergism" ) or die "Cannot open para_synergism: $!";
$line = <FHsp>;
while ( <FHsp> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    my @temp = split( /\t/, $_ );

    if ( $high_set_ind > 0 ) {
	$synr_para{$temp[0]}{'k4'}[0] = $temp[1];
	$synr_para{$temp[0]}{'k4'}[1] = $temp[2];
	$synr_para{$temp[0]}{'k5'}[0] = $temp[3];
	$synr_para{$temp[0]}{'k5'}[1] = $temp[4];
	$synr_para{$temp[0]}{'k9'}[0] = $temp[5];
	$synr_para{$temp[0]}{'k9'}[1] = $temp[6];
	$synr_para{$temp[0]}{'k10'}[0] = $temp[7];
	$synr_para{$temp[0]}{'k10'}[1] = $temp[8];
	$synr_para{$temp[0]}{'xi4'}[0] = $temp[9];
	$synr_para{$temp[0]}{'xi4'}[1] = $temp[10];
	$synr_para{$temp[0]}{'xi9'}[0] = $temp[11];
	$synr_para{$temp[0]}{'xi9'}[1] = $temp[12];
	$synr_para{$temp[0]}{'k11'}[0] = $temp[13];
	$synr_para{$temp[0]}{'k11'}[1] = $temp[14];
	$synr_para{$temp[0]}{'k12'}[0] = $temp[15];
	$synr_para{$temp[0]}{'k12'}[1] = $temp[16];
	$high_set_ind = 0;
    } else {
	$synr_para{$temp[0]}{'k4'}[2] = $temp[1];
	$synr_para{$temp[0]}{'k4'}[3] = $temp[2];
	$synr_para{$temp[0]}{'k5'}[2] = $temp[3];
	$synr_para{$temp[0]}{'k5'}[3] = $temp[4];
	$synr_para{$temp[0]}{'k9'}[2] = $temp[5];
	$synr_para{$temp[0]}{'k9'}[3] = $temp[6];
	$synr_para{$temp[0]}{'k10'}[2] = $temp[7];
	$synr_para{$temp[0]}{'k10'}[3] = $temp[8];
	$synr_para{$temp[0]}{'xi4'}[2] = $temp[9];
	$synr_para{$temp[0]}{'xi4'}[3] = $temp[10];
	$synr_para{$temp[0]}{'xi9'}[2] = $temp[11];
	$synr_para{$temp[0]}{'xi9'}[3] = $temp[12];
	$synr_para{$temp[0]}{'k11'}[2] = $temp[13];
	$synr_para{$temp[0]}{'k11'}[3] = $temp[14];
	$synr_para{$temp[0]}{'k12'}[2] = $temp[15];
	$synr_para{$temp[0]}{'k12'}[3] = $temp[16];
	$high_set_ind = 1;
    }
}
close( FHsp );

my @tf_cat = ();
my %tf_para = ();
open( FHtf, "para_dep/para_cateTF" ) or die "Cannot open para_cateTF: $!";
$line = <FHtf>;
while ( <FHtf> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    my @temp = split( /\t/, $_ );

    push(@tf_cat, $temp[0]);
    $tf_para{$temp[0]}{'ka'}[0] = $temp[1];
    $tf_para{$temp[0]}{'ka'}[1] = $temp[2];
    $tf_para{$temp[0]}{'omI'}[0] = $temp[3];
    $tf_para{$temp[0]}{'omI'}[1] = $temp[4];
    $tf_para{$temp[0]}{'omE'}[0] = $temp[5];
    $tf_para{$temp[0]}{'omE'}[1] = $temp[6];
}
close( FHtf );

@tf_cat = reverse @tf_cat;

my %mrna_para = ();
open( FHmr, "para_dep/para_mRNA" ) or die "Cannot open para_mRNA: $!";
$line = <FHmr>;
while ( <FHmr> ) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    my @temp = split( /\t/, $_ );

    $mrna_para{$temp[0]}{'k19'}[0] = $temp[1];
    $mrna_para{$temp[0]}{'k19'}[1] = $temp[2];
    $mrna_para{$temp[0]}{'k20'}[0] = $temp[3];
    $mrna_para{$temp[0]}{'k20'}[1] = $temp[4];
}
close( FHmr );

my $fileOutput;
foreach my $synr ( qw( A C B ) ) {
    $paraTerm{'k4'} = $synrTerm{'k4'} . "\t" . $synr_para{$synr}{'k4'}[0] . "\t" . $synr_para{$synr}{'k4'}[1];
    $paraTerm{'k5'} = $synrTerm{'k5'} . "\t" . $synr_para{$synr}{'k5'}[0] . "\t" . $synr_para{$synr}{'k5'}[1];
    $paraTerm{'k9'} = $synrTerm{'k9'} . "\t" . $synr_para{$synr}{'k9'}[0] . "\t" . $synr_para{$synr}{'k9'}[1];
    $paraTerm{'k10'} = $synrTerm{'k10'} . "\t" . $synr_para{$synr}{'k10'}[0] . "\t" . $synr_para{$synr}{'k10'}[1];
    $paraTerm{'xi4'} = $synrTerm{'xi4'} . "\t" . $synr_para{$synr}{'xi4'}[0] . "\t" . $synr_para{$synr}{'xi4'}[1];
    $paraTerm{'xi9'} = $synrTerm{'xi9'} . "\t" . $synr_para{$synr}{'xi9'}[0] . "\t" . $synr_para{$synr}{'xi9'}[1];
    $paraTerm{'k11'} = $synrTerm{'k11'} . "\t" . $synr_para{$synr}{'k11'}[0] . "\t" . $synr_para{$synr}{'k11'}[1];
    $paraTerm{'k12'} = $synrTerm{'k12'} . "\t" . $synr_para{$synr}{'k12'}[0] . "\t" . $synr_para{$synr}{'k12'}[1];

    foreach my $mrna ( qw( M10 M15 M30 M60 M120 M180 ) ) {
	$paraTerm{'k19'} = $mrnaTerm{'k19'} . "\t" . $mrna_para{$mrna}{'k19'}[0] . "\t" . $mrna_para{$mrna}{'k19'}[1];
	$paraTerm{'k20'} = $mrnaTerm{'k20'} . "\t" . $mrna_para{$mrna}{'k20'}[0] . "\t" . $mrna_para{$mrna}{'k20'}[1];

	for ( my $i = 0; $i <= $#tf_cat; $i++  ) {
	    $paraTerm{'k2'} = $tfTerm{'k2'} . "\t" . $tf_para{$tf_cat[$i]}{'ka'}[0] . "\t" . $tf_para{$tf_cat[$i]}{'ka'}[1];
	    $paraTerm{'omI2'} = $tfTerm{'omI2'} . "\t" . $tf_para{$tf_cat[$i]}{'omI'}[0] . "\t" . $tf_para{$tf_cat[$i]}{'omI'}[1];
	    $paraTerm{'omE2'} = $tfTerm{'omE2'} . "\t" . $tf_para{$tf_cat[$i]}{'omE'}[0] . "\t" . $tf_para{$tf_cat[$i]}{'omE'}[1];
	    for ( my $j = $i; $j <= $#tf_cat; $j++ ) {
		$paraTerm{'k7'} = $tfTerm{'k7'} . "\t" . $tf_para{$tf_cat[$j]}{'ka'}[0] . "\t" . $tf_para{$tf_cat[$j]}{'ka'}[1];
		$paraTerm{'omI7'} = $tfTerm{'omI7'} . "\t" . $tf_para{$tf_cat[$j]}{'omI'}[0] . "\t" . $tf_para{$tf_cat[$j]}{'omI'}[1];
		$paraTerm{'omE7'} = $tfTerm{'omE7'} . "\t" . $tf_para{$tf_cat[$j]}{'omE'}[0] . "\t" . $tf_para{$tf_cat[$j]}{'omE'}[1];

		my $nline = 0;
		$fileOutput = $synr . "_" . $mrna . "_" . $tf_cat[$i] . $tf_cat[$j] . "_param";
		open( SH, ">input_mRNA/" . $fileOutput ) or die "Cannot open $fileOutput: $!";
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
	}
    }
}

