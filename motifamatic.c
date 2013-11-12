/*****************************************************************************\
 motifamatic.c

Copyright (C) Ian Korf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define MAX_LENGTH 64

struct motif {
	int     length;
	double  matrix[MAX_LENGTH][4];
};
typedef struct motif * Motif;

static Motif new_motif (int w) {
	Motif  m;
	int i, j;
	
	m = malloc(sizeof(struct motif));
	m->length = w;
	
	for (i = 0; i < w; i++) {
		for (j = 0; j < 4; j++) {
			m->matrix[i][j] = 0;
		}
	}
	
	return m;
}

static Motif read_motif (const char *filename) {
	char    line[256], name[64], length[64], value[64];
	double  count, total;
	int     i, j, len;
	FILE    *stream;
	Motif   m = NULL;
	
	stream = fopen(filename, "r");
	
	fgets(line, sizeof(line), stream);
	sscanf(line, "%s length=%s", name, length);
	len = atoi(length);
	
	m = new_motif(len);
	
	/* read counts */
	for (i = 0; i < 4; i++) {
		fscanf(stream, "%s", value);
		if (!(strcmp("A", value) == 0 || strcmp("C", value) == 0 ||
			strcmp("G", value) == 0 || strcmp("T", value) == 0 ))
			fprintf(stderr, "incorrect motif format\n");
		
		for (j = 0; j < len; j++) {
			fscanf(stream, "%s", value);
			count = atof(value) + 1; /* pseudocount of 1 */
			m->matrix[j][i] = count;
		}
	}
	
	/* transform to probability */
	for (i = 0; i < len; i++) {
		total = 0;
		for (j = 0; j < 4; j++) total += m->matrix[i][j];
		for (j = 0; j < 4; j++) m->matrix[i][j] /= total;
	}
	
	return m;
}

int main (int argc, char ** argv) {
	int     i, j, len, max_i;
	double  s, score, max_s, t, expect;
	Motif   motif;
	FILE   *stream;
	char    seq[MAX_LENGTH];
	
	if (argc != 5) {
		fprintf(stderr, "usage: %s <motif> <seq file> <seq len> <threshold>\n", argv[0]);
		exit(1);
	}
	
	motif = read_motif(argv[1]);
	stream = fopen(argv[2], "r");
	len = atoi(argv[3]);
	expect = pow(0.25, motif->length);
	t = atof(argv[4]);
	if (len >= MAX_LENGTH -1) {
		fprintf(stderr, "sequences too long, < %d\n", MAX_LENGTH);
	}
	
	/* main loop */
	i = 0; j = 0;
	while (fscanf(stream, "%s", seq) == 1) {
		max_i = -1;
		max_s = -1;
		for (i = 0; i < len - motif->length +1; i++) {
			score = 1;
			for (j = 0; j < motif->length; j++) {
				switch (seq[i+j]) {
					case 'A': s = motif->matrix[j][0]; break;
					case 'C': s = motif->matrix[j][1]; break;
					case 'G': s = motif->matrix[j][2]; break;
					case 'T': s = motif->matrix[j][3]; break;
					case 'N': s = 0.1;                 break;
					default:
						fprintf(stderr, "impossible %s : >%c< %d", seq, seq[i+j], i+j);
						exit(1);
				}
				score *= s;
			}
			if (score > max_s) {
				max_s = score;
				max_i = i;
			}
		}
		
		max_s = log(max_s / expect) / log(2);
		if (max_s > t) printf("%s %d %f\n", seq, max_i, max_s);
	}
	
	return 0;
}



