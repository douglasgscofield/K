#define SAVEFILE_n       "savefile_n.txt"
#define SAVEFILEBASE_n   "savefile_n"
#define SAVEFILESUFFIX_n ".txt"

void        load_savefile_n (KConfig_n KN, KArray_n& a);
void        save_savefile_n (KConfig_n KN, KArray_n& a);
char*       create_load_savefile_name_n (KConfig_n KN);
char*       create_save_savefile_name_n (KConfig_n KN);

