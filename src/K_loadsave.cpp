#include "K.h"

////////////////////////////////////////////////////////////////
void        load_loadfile   (KConfig K, KArray& a)

// use savefile

{
    //const char* thisfunction = "load_loadfile";
    FILE* fp;
    KInt i, j, g;
    KScalar val;
    if ((fp = fopen(K->loadfile,"r")) != NULL) {
        fprintf(stderr, "loading from %s...", K->loadfile);
        while (fscanf(fp, "%d\t%d\t%d\t%lg\n",
                      &i, &j, &g, &val) != EOF) {
            if (i < 0 || i > K->MI) fatal("invalid input");
            if (j < 0 || j > K->MJ) fatal("invalid input");
            if (g < 0 || g > K->genotypes) fatal("invalid input");
            if (val < 0.0 || val > 1.0) fatal("invalid input");
            a[i][j][g] = val;
        }
        fclose(fp);
        fprintf(stderr, "done\n");
    } else {
        char buf[200];
        sprintf(buf, "could not find %s", K->loadfile);
        fatal(buf);
    }
}

////////////////////////////////////////////////////////////////
void        save_savefile   (KConfig K, KArray& a)

// create savefile

{
    //const char* thisfunction = "save_savefile";
    FILE* fp;
    if (K->savefile == (char*)0) {
        char buf[200];
        sprintf(buf, "no savefile named");
        fatal(buf);
    }
    fp = fopen(K->savefile, "w+");
    fprintf(stderr, "saving to %s...\n", K->savefile);
    dump_values_KArray(K, fp, a, K->MI, K->MJ, K->genotypes);
}

