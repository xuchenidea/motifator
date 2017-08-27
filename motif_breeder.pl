#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Motifator;
use FAlite;
use Getopt::Std;
use vars qw($opt_v $opt_2 $opt_k $opt_t $opt_T $opt_g $opt_n $opt_p $opt_b $opt_s $opt_f $opt_o);

my $KMER = 9;
my $T1 = 0.50;
my $T2 = 0.75;
my $GEN  = 100;
my $POP  = 100;
my $BPOP = 50;
my $SNP  = 0.1;
my $OS   = 1;

# background instead of sequence 2?

getopts('v2k:t:T:g:p:b:n:s:f:o:');

die "
usage: motif_breeder.pl <sequence.fa>
options:
  -v         verbose
  -2         count both strands
  -b <file>  background sequence file [default is 1 sequence]
  -k <int>   motif size [$KMER]
  -t <float> starting threshold [$T1]
  -T <float> ending threshold [$T2]
  -g <int>   generations [$GEN]
  -p <int>   population size [$POP]
  -n <int>   number of parents [$BPOP]
  -s <float> mutation probability [$SNP]
  -f <file>  file of d25 motifs
  -o <int>   optimization strategy [$OS]
" unless @ARGV == 1;

my ($file1) = @ARGV;
$KMER = $opt_k if $opt_k;
$T1   = $opt_t if $opt_t;
$T2   = $opt_T if $opt_T;
$GEN  = $opt_g if $opt_g;
$POP  = $opt_p if $opt_p;
$BPOP = $opt_n if $opt_n;
$SNP  = $opt_s if $opt_s;
$OS   = $opt_o if $opt_o;
my $VERBOSE = $opt_v;
my $COUNT2 = $opt_2;
my $SEEDS = $opt_f;
my $BKGD = $opt_b;

# get the sequence(s) and set background frequencies
my $SEQ1 = read_sequences($file1);
my $LEN1; # total length of SEQ1 set
foreach my $seq (@$SEQ1) {$LEN1 += length($seq)}
my $SEQ2; # used when background file is present
my %BKGD; # no background file given
my @PNT;  # no background file given
if ($BKGD) {
	$SEQ2 = read_sequences($BKGD);
} else {
	my %count;
	my $total;
	foreach my $seq (@$SEQ1) {
		$count{A} += $seq =~ tr/A/A/;
		$count{C} += $seq =~ tr/C/C/;
		$count{G} += $seq =~ tr/G/G/;
		$count{T} += $seq =~ tr/T/T/;
		$total += length($seq);
	}
	foreach my $nt (keys %count) {
		$BKGD{$nt} = $count{$nt} / $total;
	}
	@PNT = ($BKGD{A}, $BKGD{C}, $BKGD{G}, $BKGD{T});
	#print "@PNT";
}
	
my ($RC1, $RC2);
if ($COUNT2) {
	foreach my $seq (@$SEQ1) {push @$RC1, rc($seq)}
	if ($BKGD) {
		foreach my $seq (@$SEQ2) {push @$RC2, rc($seq)}
	}
}

# initialize motif library
Motifator::d25_init(0.97, 0.70, 0.49, 0.40, 0.33);
my @alphabet = Motifator::d25_alphabet();

# create random initial population
my @pop;
for (my $i = 0; $i < $POP; $i++) {
	my $motif;
	for (my $i = 0; $i < $KMER; $i++) {$motif .= $alphabet[rand(@alphabet)]}
	my ($fit, $s1, $m1, $s2, $m2) = fitness($motif, $T1);
	push @pop, {
		motif => $motif,
		fitness => $fit,
		s1 => $s1, # sequence count 1
		m1 => $m1, # motif count 1
		s2 => $s2,
		m2 => $m2,
	};
}

# file of favorite motifs?
if ($SEEDS) {
	open(my $fh, $SEEDS) or die;
	my $n = 0;
	while (my $motif = <$fh>) {
		chomp $motif;
		if (length($motif) != $KMER) {die "motif length mismatch to k"}
		my ($fit, $s1, $m1, $s2, $m2) = fitness($motif, $T1);
		$pop[$n] = {
			motif => $motif,
			fitness => $fit,
			s1 => $s1,
			m1 => $m1,
			s2 => $s2,
			m2 => $m2,
		};
		$n++;
	};
	close $fh;
}

# main loop
for (my $g = 0; $g < $GEN; $g++) {
	@pop = sort {$b->{fitness} <=> $a->{fitness}} @pop;
	my $t = $T1 + (($T2 - $T1) / $GEN ) * $g;
		
	printf "g%d t=%.3f:  %s %d %d %d %d %.3f   %s %d %d %d %d %.3f\n",
		$g, $t,
		$pop[0]{motif},
		$pop[0]{s1}, $pop[0]{m1}, 
		$pop[0]{s2}, $pop[0]{m2}, 
		$pop[0]{fitness},
		
		$pop[$BPOP]{motif},
		$pop[$BPOP]{s1}, $pop[$BPOP]{m1}, 
		$pop[$BPOP]{s2}, $pop[$BPOP]{m2}, 
		$pop[$BPOP]{fitness};
	
	$t = $t ** $KMER;
	
	for (my $i = $BPOP; $i < @pop; $i++) {
		my $p1 = int rand($BPOP);
		my $p2 = int rand($BPOP);
		$pop[$i] = mate($pop[$p1], $pop[$p2], $t);
	}
	
}

###################################

sub mate {
	my ($p1, $p2, $t) = @_;
	
	# recombination
	my @motif;
	for (my $i = 0; $i < $KMER; $i++) {
		if (rand() < 0.5) {push @motif, substr($p1->{motif}, $i, 1)}
		else              {push @motif, substr($p2->{motif}, $i, 1)}
	}
	
	# substitution
	for (my $i = 0; $i < @motif; $i++) {
		if (rand() < $SNP) {
			$motif[$i] = $alphabet[rand(@alphabet)];
		}
	}
	
	# frame-shift
	if (rand() < $SNP) {
		if (rand() < 0.5) {
			push @motif, $alphabet[rand(@alphabet)];
			shift @motif;
		} else {
			unshift @motif, $alphabet[rand(@alphabet)];
			pop @motif;
		}
	}
	
	# fitness
	my $motif = join("", @motif);
	my ($fit, $s1, $m1, $s2, $m2) = fitness($motif, $t);

	# child
	return {
		motif => $motif,
		fitness => $fit,
		s1 => $s1,
		m1 => $m1,
		s2 => $s2,
		m2 => $m2,
	};
	
}

sub fitness {
	my ($motif, $t) = @_;
	
	my ($s1, $m1, $s2, $m2) = count_motifs($motif, $t);
	my $rs = $s1 / ($s2 + 1); # ratio of sequences
	my $rm = $m1 / ($m2 + 1); # ratio of motifs
	my $zs = ($s1 - $s2) / sqrt($s1 +1); # z-score diff in sequences
	my $zm = ($m1 - $m2) / sqrt($m1 +1); # z-score diff in motifs
	
	my $fit;	
	if    ($OS == 1) {$fit = $rs}
	elsif ($OS == 2) {$fit = $rm}
	elsif ($OS == 3) {$fit = $zs}
	elsif ($OS == 4) {$fit = $zm}
	elsif ($OS == 5) {$fit = $rs * $zs}
	elsif ($OS == 6) {$fit = $rm * $zm}
	elsif ($OS == 7) {$fit = $rs * $zm}
	elsif ($OS == 8) {$fit = $rm * $zs}
	else {die "unknown optimization strategy $OS"}
	
	return $fit, $s1, $m1, $s2, $m2;
}

sub count_motifs {
	my ($motif, $t) = @_;
	my ($s1, $m1, $s2, $m2) = (0, 0, 0, 0);
	
	# count group 1
	for (my $i = 0; $i < @$SEQ1; $i++) {
		my $count = Motifator::d25_count_motif($motif, $SEQ1->[$i], $t);
		if ($COUNT2) {
			$count += Motifator::d25_count_motif($motif, $RC1->[$i], $t);
		}
		$s1++ if $count;
		$m1 += $count;
	}
	
	# count group 2
	if ($BKGD) {
		for (my $i = 0; $i < @$SEQ2; $i++) {
			my $count = Motifator::d25_count_motif($motif, $SEQ2->[$i], $t);
			if ($COUNT2) {
				$count += Motifator::d25_count_motif($motif, $RC2->[$i], $t);
			}
			$s2++ if $count;
			$m2 += $count;
		}
	} else {
		my $p = Motifator::d25_seq_prob($motif, @PNT);
		$m2 = int($p * $LEN1);
		$s2 = $m2 < @$SEQ1 ? $m2 : @$SEQ1;
	}

	return $s1, $m1, $s2, $m2;
}

sub rc {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ACGT/TGCA/;
	return $seq;
}

sub read_sequences {
	my ($file) = @_;
	my @seq;
	open(my $fh, $file) or die "file not found: $file";
	my $fasta = new FAlite($fh);
	while (my $entry = $fasta->nextEntry) {
		push @seq, uc $entry->seq;
	}
	return \@seq;
}



