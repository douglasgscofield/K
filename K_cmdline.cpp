#include "K.h"

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Routines for handling command-line arguments                  */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////// */
int
cmdline_args(
	KConfig K,
	int argc,
	char *argv[])
{
	const char *thisfunction = "cmdline_args";
	int i;
	if (K->genotypes != 1) {
		not_implemented(thisfunction, "genotypes != 1");
	}
	i = 1;
	while (argv[i]) {
		if (!strcmp(argv[i], "--help")) {
			cmdline_usage(K);
			exit(0);
		} else if (!strcmp(argv[i], "-nested")) {
			/* request for nested, but this is the unnested cmdline_args! */
			char buf[300];
			sprintf(buf, "Unnested simulation, option invalid: <%s>\n",
					argv[i]);
			fatal(buf);
		} else if (!strcmp(argv[i], "-unnested")) {
			/* request for unnested, which is OK, and ignored here */
			i += 1;
		} else if (!strcmp(argv[i], "-U")) {
			/* single argument for genomic mutation rate */
			K->U = (double) atof(argv[i + 1]);
			i += 2;
		} else if (!strcmp(argv[i], "-s")) {
			/* single argument for selection coefficient */
			K->fit_s = (double) atof(argv[i + 1]);
			i += 2;
		} else if (!strcmp(argv[i], "-h")) {
			/* single argument for dominance coefficient */
			K->fit_h = (double) atof(argv[i + 1]);
			i += 2;
		} else if (!strcmp(argv[i], "-S")) {
			/* single argument for selfing rate */
			K->S[0] = (double) atof(argv[i + 1]);
			i += 2;
		} else if (!strcmp(argv[i], "-truncate")) {
			K->option_truncate = 1;
			i += 1;
		} else if (!strcmp(argv[i], "-nolethal")) {
			K->option_nolethal = 1;
			i += 1;
		} else if (!strcmp(argv[i], "-tableheadingonly")) {
			/* print the 'table' heading and exit */
			stats_print_table_heading(K);
			exit(0);
		} else if (!strcmp(argv[i], "-table")) {
			/* print stats as 'table' */
			K->option_table = 1;
			i += 1;
		} else if (!strcmp(argv[i], "-load_savefile")) {
			K->load_savefile = 1;
			i += 1;
		} else if (!strcmp(argv[i], "-load_savefile_name")) {
			K->load_savefile_name =
				(char *) calloc((size_t) sizeof(char),
								strlen(argv[i + 1]) + 1);
			strcpy(K->load_savefile_name, argv[i + 1]);
			i += 2;
		} else if (!strcmp(argv[i], "-save_savefile")) {
			K->save_savefile = 1;
			i += 1;
		} else if (!strcmp(argv[i], "-save_savefile_name")) {
			K->save_savefile_name =
				(char *) calloc((size_t) sizeof(char),
								strlen(argv[i + 1]) + 1);
			strcpy(K->save_savefile_name, argv[i + 1]);
			i += 2;
			/*
			   } else if (!strcmp( argv[i], "-inita1" )) {
			   // initial frequency of allele 1 
			   INITIAL_A1 = (double)atof( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-n" )) {
			   N = atol( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-fecm" )) {
			   FECM = atol( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-fecf" )) {
			   FECF = atol( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-gens" )) {
			   GENS = atol( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-selmod" )) {
			   selection_model = atoi( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-popmod" )) {
			   selection_model = atoi( argv[i+1] );
			   i += 2;
			   } else if (!strcmp( argv[i], "-selfsurv" )) {
			   self_survival = (double)atof( argv[i+1] );
			   i += 2;
			 */
		} else {
			char buf[300];
			sprintf(buf, "Unrecognized argument: <%s>\n", argv[i]);
			fatal(buf);
		}
	}
	if (i == 1) {
		/* there were no arguments */
		/* return 1; */
	}
	return 0;
}

void
cmdline_usage(
	KConfig K)
{
	const char *thisfunction = "cmdline_usage";
	fprintf(stdout, "K [args]\n\
\n\
  -FLAG <ARGUMENT>  EXPLANATION [CURRENT VALUE]\n\
  ---------------------------------------------\n\
  --help            This help message\n\
  -U <float>        Genomic mutation rate [%lg]\n\
  -s <float>        Selection coefficient [%lg]\n\
  -h <float>        Dominance coefficient [%lg]\n\
  -S <float>        Selfing rate [%lg]\n\
  -truncate         Turn on truncate [%d]\n\
  -nolethal         Do not shortcut computations with lethals [%d]\n\
  -table            Print stats output in single-line format\n\
  -tableheadingonly Print the heading for a table and exit\n\
  -load_savefile    Load a saved file of frequencies [%d]\n\
  -load_savefile_name <string>  File from which to load [%s]\n\
  -save_savefile    Save a file of frequencies [%d]\n\
  -save_savefile_name <string>  File to which to save [%s]\n\
", K->U, K->fit_s, K->fit_h, K->S[0], K->option_truncate, K->option_nolethal, K->load_savefile, K->load_savefile_name, K->save_savefile, K->save_savefile_name);
}
