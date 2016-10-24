#include "K.h"
#include "SimpleOpt.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for handling command-line arguments                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
int         cmdline_args        (KConfig K, int argc, char *argv[])
{
    const char* thisfunction = "cmdline_args";
    int i;
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    enum {
        o_help,
        o_progress,
        o_generation_cutoff,
        o_nested,
        o_unnested,
        o_U_mutation,
        o_s_selection,
        o_h_dominance,
        o_S_selfing_rate,
        o_A_apomixis_rate,
        o_truncate,
        o_no_lethal,
        o_table_heading_only,
        o_table,
        o_load_savefile,
        o_load_savefile_name,
        o_save_savefile,
        o_save_savefile_name,
    };
    CSimpleOpt::SOption K_options[] = {
        { o_help, "-?", SO_NONE },
        { o_help, "--help", SO_NONE },
        { o_progress, "--progress", SO_REQ_SEP },
        { o_generation_cutoff, "--generation-cutoff", SO_REQ_SEP },
        { o_nested, "--nested", SO_NONE },
        { o_unnested, "--unnested", SO_NONE },
        { o_U_mutation, "-U", SO_REQ_SEP },
        { o_s_selection, "-s", SO_REQ_SEP },
        { o_h_dominance, "-h", SO_REQ_SEP },
        { o_S_selfing_rate, "-S", SO_REQ_SEP },
        { o_A_apomixis_rate, "-A", SO_REQ_SEP },
        { o_truncate, "--truncate", SO_NONE },
        { o_no_lethal, "--no-lethal-shortcut", SO_NONE },
        { o_table_heading_only, "--table-heading-only", SO_NONE },
        { o_table, "--table", SO_NONE },
        { o_load_savefile, "--load-savefile", SO_NONE },
        { o_load_savefile_name, "--load-savefile-name", SO_REQ_SEP },
        { o_save_savefile, "--save-savefile", SO_NONE },
        { o_save_savefile_name, "--save-savefile-name", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, K_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << "invalid argument '" << args.OptionText() << "'" << endl;
            cmdline_usage(K);
            exit(0);
        }
        switch (args.OptionId()) {
            case o_help:
                cmdline_usage(K); exit(0); break;
            case o_progress:
                K->progress = atol(args.OptionArg()); break;
            case o_generation_cutoff:
                GENERATION_CUTOFF = atoi(args.OptionArg()); break;
            case o_nested:
                cerr << "Unnested simulation, option invalid '" << args.OptionText() << "'" << endl; break;
                exit(1);
            case o_unnested:
                /* ignored here */ break;
            case o_U_mutation:
                K->U = atof(args.OptionArg()); break;
            case o_s_selection:
                K->fit_s = atof(args.OptionArg()); break;
            case o_h_dominance:
                K->fit_h = atof(args.OptionArg()); break;
            case o_S_selfing_rate:
                K->S[0] = atof(args.OptionArg()); break;
            case o_A_apomixis_rate:
                K->A[0] = atof(args.OptionArg()); break;
            case o_truncate:
                K->option_truncate = 1; break;
            case o_no_lethal:
                K->option_nolethal = 1; break;
            case o_table_heading_only:
                stats_print_table_heading(K); exit(0); break;
            case o_table:
                K->option_table = 1; break;
            case o_load_savefile:
                K->load_savefile = 1; break;
            case o_load_savefile_name:
                K->load_savefile_name = args.OptionArg(); break;
            case o_save_savefile:
                K->save_savefile = 1; break;
            case o_save_savefile_name:
                K->save_savefile_name = args.OptionArg(); break;
            default:
                cerr << "Invalid option '" << args.OptionText() << "'" << endl; exit(1); break;

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
        }
    }
    return 0;
}

void        cmdline_usage           (KConfig K)
{
    fprintf(stderr, "K [args]\n\
\n\
  -FLAG <ARGUMENT>  EXPLANATION [CURRENT VALUE]\n\
  ---------------------------------------------\n\
  --help            This help message\n\
  -U <float>        Genomic mutation rate [%lg]\n\
  -s <float>        Selection coefficient [%lg]\n\
  -h <float>        Dominance coefficient [%lg]\n\
  -S <float>        Selfing rate [%lg]\n\
  -A <float>        Apomixis rate [%lg]\n\
  --truncate         Turn on truncate [%d]\n\
  --no-lethal         Do not shortcut computations with lethals [%d]\n\
  --table            Print stats output in single-line format\n\
  --table-heading-only Print the heading for a table and exit\n\
  --load-savefile    Load a saved file of frequencies [%d]\n\
  --load-savefile-name <string>  File from which to load [%s]\n\
  --save-savefile    Save a file of frequencies [%d]\n\
  --save-savefile-name <string>  File to which to save [%s]\n\
",
    K->U,
    K->fit_s,
    K->fit_h,
    K->S[0],
    K->A[0],
    K->option_truncate,
    K->option_nolethal,
    K->load_savefile,
    K->load_savefile_name,
    K->save_savefile,
    K->save_savefile_name
    );
}

