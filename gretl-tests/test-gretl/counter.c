#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int gcount = 0;
int fcount = 0;

int read_output (const char *fname) 
{
    FILE *fp;
    char line[512];
    int n, err = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't open %s\n", fname);
	return 1;
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (!strncmp(line, "Function evaluations:", 21)) {
	    if (sscanf(line + 22, "%d", &n) == 1) {
		fcount += n;
	    } else {
		fprintf(stderr, "Couldn't get fcount\n");
		err = 1;
	    }
	} else if (!strncmp(line, "Evaluations of gradient:", 24)) {
	    if (sscanf(line + 25, "%d", &n)) {
		gcount += n;
	    } else {
		fprintf(stderr, "Couldn't get gcount\n");
		err = 1;
	    }		
	} 
    }

    fclose(fp);

    return err;
}

int main (void) 
{
    FILE *fp;
    char line[512];
    char script[48];
    char outname[64];
    char *s;
    int gcount0, fcount0;
    int err = 0;

    fp = fopen("ps.list", "r");
    if (fp == NULL) {
	fputs("Couldn't open ps.list\n", stderr);
	exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line != '#' && !isspace(*line)) {
	    sscanf(line, "%47s", script);
	    sprintf(outname, "output/%s", script);
	    s = strstr(outname, ".inp");
	    if (s != NULL) {
		strcpy(s, ".out");
		err = read_output(outname);
	    } else {
		err = 1;
	    }
	}
    }

    gcount0 = gcount;
    fcount0 = fcount;
    gcount = fcount = 0;

    rewind(fp);

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line != '#' && !isspace(*line)) {
	    sscanf(line, "%47s", script);
	    sprintf(outname, "newout/%s", script);
	    s = strstr(outname, ".inp");
	    if (s != NULL) {
		strcpy(s, ".out");
		err = read_output(outname);
	    } else {
		err = 1;
	    }
	}
    }

    fclose(fp);

    printf("                     output newout\n"
	   "BFGS gradient evals: %6d %6d\n"
	   "BFGS function evals: %6d %6d\n",
	   gcount0, gcount, fcount0, fcount);    

    return 0;
}
