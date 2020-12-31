#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{
    FILE *fp;
    char *p, line[1024];
    int oldfcount = 0;
    int oldgcount = 0;
    int newfcount = 0;
    int newgcount = 0;
    int c;

    fp = fopen("diffs", "r");
    if (fp == NULL) {
	fprintf(stderr, "%s: couldn't open diffs file\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof line, fp)) {
	p = strstr(line, "Function evaluations:");
	if (p != NULL) {
	    sscanf(p + 21, "%d", &c);
	    if (*line == '<') {
		oldfcount += c;
	    } else {
		newfcount += c;
	    }
	} else {
	    p = strstr(line, "Evaluations of gradient:");
	    if (p != NULL) {
		sscanf(p + 24, "%d", &c);
		if (*line == '<') {
		    oldgcount += c;
		} else {
		    newgcount += c;
		}
	    }	
	}
    }

    printf("Total 'function evaluation' counts\n");
    printf("  baseline %d, new %d\n", oldfcount, newfcount);
    printf("Total 'gradient evaluation' counts\n");
    printf("  baseline %d, new %d\n", oldgcount, newgcount);

    fclose(fp);

    return 0;
}
