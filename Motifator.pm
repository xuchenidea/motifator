package Motifator;
use strict; 
use warnings;
use Inline C => 'DATA';

sub d25_alphabet {
	return split("", "AaCcGgTtRrYyMmKkWwSsBDHVN");
}


1;

__DATA__
__C__

static double ** Prob;
static int Initialized = 0;

void d25_init (double A, double a, double B, double b, double C) {
	int i;
	int j;
	
	Prob = malloc(128 * sizeof(double*));
	for (i = 0; i < 128; i++) {
		Prob[i] = malloc(128 * sizeof(double));
		for (j = 0; j < 128; j++) {
			Prob[i][j] = -1;
		}
	}
	
	Prob['A']['A'] = A;
	Prob['A']['C'] = (1 - A) / 3;
	Prob['A']['G'] = (1 - A) / 3;
	Prob['A']['T'] = (1 - A) / 3;
	
	Prob['a']['A'] = a;
	Prob['a']['C'] = (1 - a) / 3;
	Prob['a']['G'] = (1 - a) / 3;
	Prob['a']['T'] = (1 - a) / 3;
	
	Prob['C']['A'] = (1 - A) / 3;
	Prob['C']['C'] = A;
	Prob['C']['G'] = (1 - A) / 3;
	Prob['C']['T'] = (1 - A) / 3;
	
	Prob['c']['A'] = (1 - a) / 3;
	Prob['c']['C'] = a;
	Prob['c']['G'] = (1 - a) / 3;
	Prob['c']['T'] = (1 - a) / 3;
	
	Prob['G']['A'] = (1 - A) / 3;
	Prob['G']['C'] = (1 - A) / 3;
	Prob['G']['G'] = A;
	Prob['G']['T'] = (1 - A) / 3;
	
	Prob['g']['A'] = (1 - a) / 3;
	Prob['g']['C'] = (1 - a) / 3;
	Prob['g']['G'] = a;
	Prob['g']['T'] = (1 - a) / 3;
	
	Prob['T']['A'] = (1 - A) / 3;
	Prob['T']['C'] = (1 - A) / 3;
	Prob['T']['G'] = (1 - A) / 3;
	Prob['T']['T'] = A;
	
	Prob['t']['A'] = (1 - a) / 3;
	Prob['t']['C'] = (1 - a) / 3;
	Prob['t']['G'] = (1 - a) / 3;
	Prob['t']['T'] = a;
	
	Prob['R']['A'] = B;
	Prob['R']['C'] = 0.5 - B;
	Prob['R']['G'] = B;
	Prob['R']['T'] = 0.5 - B;
	
	Prob['r']['A'] = b;
	Prob['r']['C'] = 0.5 - b;
	Prob['r']['G'] = b;
	Prob['r']['T'] = 0.5 - b;
	
	Prob['Y']['A'] = 0.5 - B;
	Prob['Y']['C'] = B;
	Prob['Y']['G'] = 0.5 - B;
	Prob['Y']['T'] = B;
	
	Prob['y']['A'] = 0.5 - b;
	Prob['y']['C'] = b;
	Prob['y']['G'] = 0.5 - b;
	Prob['y']['T'] = b;
	
	Prob['M']['A'] = B;
	Prob['M']['C'] = B;
	Prob['M']['G'] = 0.5 - B;
	Prob['M']['T'] = 0.5 - B;
	
	Prob['m']['A'] = b;
	Prob['m']['C'] = b;
	Prob['m']['G'] = 0.5 - b;
	Prob['m']['T'] = 0.5 - b;
	
	Prob['K']['A'] = 0.5 - B;
	Prob['K']['C'] = 0.5 - B;
	Prob['K']['G'] = B;
	Prob['K']['T'] = B;
	
	Prob['k']['A'] = 0.5 - b;
	Prob['k']['C'] = 0.5 - b;
	Prob['k']['G'] = b;
	Prob['k']['T'] = b;
	
	Prob['W']['A'] = B;
	Prob['W']['C'] = 0.5 - B;
	Prob['W']['G'] = 0.5 - B;
	Prob['W']['T'] = B;
	
	Prob['w']['A'] = b;
	Prob['w']['C'] = 0.5 - b;
	Prob['w']['G'] = 0.5 - b;
	Prob['w']['T'] = b;
	
	Prob['S']['A'] = 0.5 - B;
	Prob['S']['C'] = B;
	Prob['S']['G'] = B;
	Prob['S']['T'] = 0.5 - B;
	
	Prob['s']['A'] = 0.5 - b;
	Prob['s']['C'] = b;
	Prob['s']['G'] = b;
	Prob['s']['T'] = 0.5 - b;
	
	Prob['B']['A'] = 1 - C*3;
	Prob['B']['C'] = C;
	Prob['B']['G'] = C;
	Prob['B']['T'] = C;
	
	Prob['D']['A'] = C;
	Prob['D']['C'] = 1 - C*3;
	Prob['D']['G'] = C;
	Prob['D']['T'] = C;
		
	Prob['H']['A'] = C;
	Prob['H']['C'] = C;
	Prob['H']['G'] = 1 - C*3;
	Prob['H']['T'] = C;
	
	Prob['V']['A'] = C;
	Prob['V']['C'] = C;
	Prob['V']['G'] = C;
	Prob['V']['T'] = 1 - C*3;
	
	Prob['N']['A'] = 0.25;
	Prob['N']['C'] = 0.25;
	Prob['N']['G'] = 0.25;
	Prob['N']['T'] = 0.25;
	
	Initialized = 1;
	
}

double d25_seq_prob (const char *motif, double a, double c, double g, double t) {
	int i;
	char s;
	double p, prob = 1;
	
	for (i = 0; i < strlen(motif); i++) {
		s = motif[i];
		p = Prob[s]['A'] * a + Prob[s]['C'] * c + Prob[s]['G'] * g + Prob[s]['T'] * t;
		/*switch (s) {
			case 'R': case 'Y': case 'r': case 'y': p *= 2; break;
			case 'M': case 'm': case 'K': case 'k': p *= 2; break;
			case 'S': case 's': case 'W': case 'w': p *= 2; break;
			case 'B': case 'D': case 'H': case 'V': p *= 3; break;
			case 'b': case 'd': case 'h': case 'v': p *= 3; break;
			case 'N': p *= 4; break;
		}*/
		prob *= p;
	}
	
	return prob;
}

void d25_print_motif (char *motif) {
	int i;
	double *f;
	
	if (!Initialized) {
		fprintf(stderr, "Motifator::init not called\n");
		exit(1);
	}
	
	for (i = 0; i < strlen(motif); i++) {
		switch (motif[i]) {
			case 'A': f = Prob['A']; break;
			case 'a': f = Prob['a']; break;
			case 'C': f = Prob['C']; break;
			case 'c': f = Prob['c']; break;
			case 'G': f = Prob['G']; break;
			case 'g': f = Prob['g']; break;
			case 'T': f = Prob['T']; break;
			case 't': f = Prob['t']; break;
			case 'R': f = Prob['R']; break;
			case 'r': f = Prob['r']; break;
			case 'Y': f = Prob['Y']; break;
			case 'y': f = Prob['y']; break;
			case 'M': f = Prob['M']; break;
			case 'm': f = Prob['m']; break;
			case 'K': f = Prob['K']; break;
			case 'k': f = Prob['k']; break;
			case 'W': f = Prob['W']; break;
			case 'w': f = Prob['w']; break;
			case 'S': f = Prob['S']; break;
			case 's': f = Prob['s']; break;
			case 'B': f = Prob['B']; break;
			case 'D': f = Prob['D']; break;
			case 'H': f = Prob['H']; break;
			case 'V': f = Prob['V']; break;
			case 'N': f = Prob['N']; break;
			default: ik_exit(1, "switch not possible");
		}
		printf("%c %.3f %.3f %.3f %.3f\n", motif[i],
			f['A'], f['C'], f['G'], f['T']);
	}
	
}

int d25_count_motif (char *motif, char *seq, double t) {
	int i, j, k, l, n, nt, mt;
	double p;
	
	if (!Initialized) {
		fprintf(stderr, "Motifator::init not called\n");
		exit(1);
	}
	
	k = strlen(motif);
	l = strlen(seq);
	n = 0; // number of motifs found over threshold
	
	for (i = 0; i < l -k + 1; i++) {
		p = 1;
		for (j = 0; j < k; j++) {
			nt = seq[i+j]; // the current nucleotide
			mt = motif[j]; // the current motif symbol
			p *= Prob[mt][nt];
		}
		if (p > t) n++;
	}
	
	return n;
}


