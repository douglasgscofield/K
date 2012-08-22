#define SAVEFILE         "savefile.txt"
#define SAVEFILEBASE     "savefile"
#define SAVEFILESUFFIX   ".txt"

void        load_savefile   (KConfig K, KArray& a);
void        save_savefile   (KConfig K, KArray& a);
char*       create_load_savefile_name   (KConfig K);
char*       create_save_savefile_name   (KConfig K);

