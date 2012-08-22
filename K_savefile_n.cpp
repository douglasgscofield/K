#include "K.h"

/* ///////////////////////////////////////////////////////////// */
void
load_savefile_n(
	KConfig_n KN,
	KArray_n & a)
/*
** use savefile
*/
{
	const char *thisfunction = "load_savefile_n";
	FILE *fp;
	KInt i0, j0, i1, j1;
	KScalar val;
	char *savefile;
	savefile = create_load_savefile_name_n(KN);
	if ((fp = fopen(savefile, "r")) != NULL) {
		printf("loading from %s... ", savefile);
		while (fscanf(fp, "%d\t%d\t%d\t%d\t%lg\n", &i0, &j0, &i1, &j1, &val) !=
			   EOF) {
			if (i0 < 0 || i0 > KN->MI0)
				fatal("invalid input");
			if (j0 < 0 || j0 > KN->MJ0)
				fatal("invalid input");
			if (i1 < 0 || i1 > KN->MI1)
				fatal("invalid input");
			if (j1 < 0 || j1 > KN->MJ1)
				fatal("invalid input");
			if (val < 0.0 || val > 1.0)
				fatal("invalid input");
			a[i0][j0][i1][j1] = val;
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
save_savefile_n(
	KConfig_n KN,
	KArray_n & a)
/*
** create savefile
*/
{
	const char *thisfunction = "save_savefile_n";
	char *savefile;
	FILE *fp;
	savefile = create_save_savefile_name_n(KN);
	fp = fopen(savefile, "w+");
	printf("saving to %s... ", savefile);
	dump_values_KArray_n(KN, fp, a, KN->MI0, KN->MJ0, KN->MI1, KN->MJ1);
	printf("done\n");
}

/* ///////////////////////////////////////////////////////////// */
char *
create_load_savefile_name_n(
	KConfig_n KN)
/*
** generate name for savefile from information in KN
*/
{
	const char *thisfunction = "create_load_savefile_name_n";
	char *ans;
	char buf[1000];
	if (KN->load_savefile_name != (char *) NULL)
		return KN->load_savefile_name;
	strcpy(buf, SAVEFILEBASE_n);
	strcat(buf, SAVEFILESUFFIX_n);
	ans = (char *) calloc((size_t) sizeof(char), strlen(buf) + 1);
	strcpy(ans, buf);
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
char *
create_save_savefile_name_n(
	KConfig_n KN)
/*
** generate name for savefile from information in KN
*/
{
	const char *thisfunction = "create_save_savefile_name_n";
	char *ans;
	char buf[1000];
	if (KN->save_savefile_name != (char *) NULL)
		return KN->save_savefile_name;
	strcpy(buf, SAVEFILEBASE_n);
	strcat(buf, SAVEFILESUFFIX_n);
	ans = (char *) calloc((size_t) sizeof(char), strlen(buf) + 1);
	strcpy(ans, buf);
	return ans;
}
