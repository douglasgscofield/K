#include "K_n.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for handling command-line arguments                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
int         cmdline_args_n  (KConfig_n KN, int argc, char *argv[])
{
    const char* thisfunction = "cmdline_args_n";
    int i;
    i = 1;
    while (argv[i]) {
        if (!strcmp(argv[i], "--help")) {
            cmdline_usage_n(KN);
            exit(0);
        } else if (!strcmp( argv[i], "-nested" )) {
            /* request for nested, which is OK, and ignored here */
            i += 1;
        } else if (!strcmp( argv[i], "-unnested" )) {
            /* request for unnested, but this is the nested cmdline_args! */
            char buf[300];
            sprintf( buf, "Nested simulation, option invalid: <%s>\n",
                argv[i] );
            fatal( buf );
        } else if (!strcmp( argv[i], "-U0" )) {
            KN->U[0] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-U1" )) {
            KN->U[1] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-s0" )) {
            KN->fit_s[0] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-s1" )) {
            KN->fit_s[1] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-h0" )) {
            KN->fit_h[0] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-h1" )) {
            KN->fit_h[1] = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-S" )) {
            KN->S = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-A" )) {
            KN->A = (double)atof( argv[i+1] );
            i += 2;
        } else if (!strcmp( argv[i], "-truncate" )) {
            KN->option_truncate = 1;
            i += 1;
        } else if (!strcmp( argv[i], "-nolethal" )) {
            KN->option_nolethal = 1;
            i += 1;
        } else if (!strcmp( argv[i], "-tableheadingonly" )) {
            stats_print_table_heading_n(KN);
            exit(0);
        } else if (!strcmp( argv[i], "-table" )) {
            KN->option_table = 1;
            i += 1;
        } else if (!strcmp( argv[i], "-load_savefile" )) {
            KN->load_savefile = 1;
            i += 1;
        } else if (!strcmp( argv[i], "-load_savefile_name" )) {
            KN->load_savefile_name = (char*)calloc((size_t)sizeof(char),
                                                   strlen(argv[i+1])+1);
            strcpy(KN->load_savefile_name, argv[i+1]);
            i+= 2;
        } else if (!strcmp( argv[i], "-save_savefile" )) {
            KN->save_savefile = 1;
            i += 1;
        } else if (!strcmp( argv[i], "-save_savefile_name" )) {
            KN->save_savefile_name = (char*)calloc((size_t)sizeof(char),
                                                   strlen(argv[i+1])+1);
            strcpy(KN->save_savefile_name, argv[i+1]);
            i+= 2;
        } else {
            char buf[300];
            sprintf( buf, "Unrecognized argument: <%s>\n", argv[i] );
            fatal( buf );
        }
    }
    if (i == 1) {
        /* there were no arguments */
        /* return 1; */
    }
    return 0;
}

void        cmdline_usage_n     (KConfig_n KN)
{
    const char* thisfunction = "cmdline_usage_n";
    fprintf(stderr, "K [args]\n\
\n\
  -FLAG <ARGUMENT>  EXPLANATION [CURRENT VALUE]\n\
  ---------------------------------------------\n\
  --help            This help message\n\
  -U0 <float>       Genomic mutation rate mut class 0 [%lg]\n\
  -U1 <float>       Genomic mutation rate mut class 1 [%lg]\n\
  -s0 <float>       Selection coefficient mut class 0 [%lg]\n\
  -s1 <float>       Selection coefficient mut class 1 [%lg]\n\
  -h0 <float>       Dominance coefficient mut class 0 [%lg]\n\
  -h1 <float>       Dominance coefficient mut class 1 [%lg]\n\
  -S <float>        Selfing rate [%lg]\n\
  -A <float>        Apomixis rate [%lg]\n\
  -truncate         Turn on truncate [%d]\n\
  -nolethal         Do not shortcut computations with lethals [%d]\n\
  -table            Print stats out put in single-line format\n\
  -tableheadingonly Print the heading for a table and exit\n\
  -load_savefile    Load a saved file of frequencies [%d]\n\
  -load_savefile_name <string>  File from which to load [%s]\n\
  -save_savefile    Save a file of frequencies [%d]\n\
  -save_savefile_name <string>  File to which to save [%s]\n\
",
    KN->U[0],
    KN->U[1],
    KN->fit_s[0],
    KN->fit_s[1],
    KN->fit_h[0],
    KN->fit_h[1],
    KN->S,
    KN->A,
    KN->option_truncate,
    KN->option_nolethal,
    KN->load_savefile,
    KN->load_savefile_name,
    KN->save_savefile,
    KN->save_savefile_name
    );
}

