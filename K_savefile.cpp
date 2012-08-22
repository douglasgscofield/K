#include "K.h"

/* ///////////////////////////////////////////////////////////// */
void
load_savefile(
	KConfig K,
	KArray & a)
/*
** use savefile
*/
{
	const char *thisfunction = "load_savefile";
	FILE *fp;
	KInt i, j, g;
	KScalar val;
	char *savefile;
	savefile = create_load_savefile_name(K);
	if ((fp = fopen(savefile, "r")) != NULL) {
		printf("loading from %s...", SAVEFILE);
		while (fscanf(fp, "%d\t%d\t%d\t%lg\n", &i, &j, &g, &val) != EOF) {
			if (i < 0 || i > K->MI)
				fatal("invalid input");
			if (j < 0 || j > K->MJ)
				fatal("invalid input");
			if (g < 0 || g > K->genotypes)
				fatal("invalid input");
			if (val < 0.0 || val > 1.0)
				fatal("invalid input");
			a[i][j][g] = val;
		}
		fclose(fp);
		printf("done\n");
	} else {
		char buf[200];
		sprintf(buf, "could not find %s", savefile);
		fatal(buf);
	}
}

/* ///////////////////////////////////////////////////////////// */
void
save_savefile(
	KConfig K,
	KArray & a)
/*
** create savefile
*/
{
	const char *thisfunction = "save_savefile";
	char *savefile;
	FILE *fp;
	savefile = create_save_savefile_name(K);
	fp = fopen(savefile, "w+");
	printf("saving to %s...\n", savefile);
	dump_values_KArray(K, fp, a, K->MI, K->MJ, K->genotypes);
}

/* ///////////////////////////////////////////////////////////// */
char *
create_load_savefile_name(
	KConfig K)
/*
** generate name for savefile from information in K
*/
{
	const char *thisfunction = "create_load_savefile_name";
	char *ans;
	char buf[1000];
	if (K->load_savefile_name != (char *) 0)
		return K->load_savefile_name;
	strcpy(buf, SAVEFILEBASE);
	strcat(buf, SAVEFILESUFFIX);
	ans = (char *) calloc((size_t) sizeof(char), strlen(buf) + 1);
	strcpy(ans, buf);
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
char *
create_save_savefile_name(
	KConfig K)
/*
** generate name for savefile from information in K
*/
{
	const char *thisfunction = "create_save_savefile_name";
	char *ans;
	char buf[1000];
	if (K->save_savefile_name != (char *) 0)
		return K->save_savefile_name;
	strcpy(buf, SAVEFILEBASE);
	strcat(buf, SAVEFILESUFFIX);
	ans = (char *) calloc((size_t) sizeof(char), strlen(buf) + 1);
	strcpy(ans, buf);
	return ans;
}
