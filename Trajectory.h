#ifndef TRAJECTORY_H
#define TRAJECTORY_H

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Class for tracking genotype trajectories through time
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

class Trajectory {
    private:
        bool           active_flag;
        bool           debug_flag;
        int            type;     // type of trajectory
#define  TRAJECTORY_ERROR        (-1)
#define  TRAJECTORY_END_OF_GEN    1
        static const struct struct_types {
            string       name;
            string       define_name;
            int               type;
        } types[];
        static const int num_types;
        int            freq;     // frequency of recording
        int            gen_started;  // gen trajectory started
        int            gen_to_end; // not implemented
        int            last_gen;  // gen last recorded
        int            next_gen;  // gen to record next
        int            num_recorded; // number of times
        ofstream  file;     // the file stream
        string    filename; // filename
        bool           dropzero_flag;
        bool           writenow_flag;

    public:
        Trajectory();
        ~Trajectory();

        void start(KConfig K, 
                   const string& file, 
                   const string& type, 
                   int freq);
        void check(KConfig K);
        void write(KConfig K);
        void debug(bool flag) { debug_flag = flag; }
        int  lastgen() { return(last_gen); }
        bool active() { return(active_flag); }
        void stop();
        void reset();

        static void                 dump_type  ();
        static int                  match_type (const string& type);
        static const string&   match_type (int type);
};

/////////////////////////////////////////////////////////////////
inline Trajectory::Trajectory()
    : active_flag(false), debug_flag(false), dropzero_flag(true), writenow_flag(false)
{ /* EMPTY */ }

/////////////////////////////////////////////////////////////////
inline Trajectory::~Trajectory()
{ /* EMPTY */ }

/////////////////////////////////////////////////////////////////
inline void Trajectory::check(KConfig K) {
    const char* thisfunction = "Trajectory::check";
    const int this_gen = K->generation - 1; // because *_nextgen increments it
    // will need to change with more than one trajectory
    // and with different types of trajectories
    if (debug_flag) {
        cerr << thisfunction << " in gen " << this_gen
            << " with next_gen == " << next_gen << endl;
    }
    if (! active_flag) {
        cerr << thisfunction << ": no active trajectory" << endl;
        fatal(NULL);
    }
    if (this_gen != next_gen && ! writenow_flag) {
        if (debug_flag) {
            cerr << thisfunction << " no write, returning" << endl;
        }
        return;
    }
    if (writenow_flag) writenow_flag = false;
    // write the trajectory information
    write(K);
    // set up for next time
    next_gen += freq;
    ++ num_recorded;
    if (debug_flag) {
        cerr << thisfunction << " in gen " << this_gen
            << " finished writing" << endl;
    }
}

/////////////////////////////////////////////////////////////////
inline void Trajectory::write(KConfig K) {
    //const char* thisfunction = "Trajectory::write";
    int this_gen = K->generation;
    if (type == TRAJECTORY_END_OF_GEN) {
        this_gen -= 1;  // see Trajectory::check
        KArray& a = K->x;
        file << endl << "Trajectory" 
            << "\tgen:" << this_gen
            << "\ttype:" << type 
            << "\tfreq:" << freq 
            << "\twritenow:" << writenow_flag
            << "\tdropzero:" << dropzero_flag
            << "\tatequilibrium:" << is_equilibrium(K)
            << "\theterozygotes:" << "0" << "," << K->MI
            << "\thomozygotes:" << "0" << "," << K->MJ
            << endl;
        file << "heterozygote\thomozygote\tgenotype\tfrequency" << endl;
        for (int i = 0; i <= K->MI; i++) {
            for (int j = 0; j <= K->MJ; j++) {
                for (int g = 0; g < K->genotypes; g++) {
                    if (dropzero_flag && a[i][j][g] == 0.0) {
                        continue;
                    }
                    file << i << "\t" << j << "\t" << g << "\t" 
                        << a[i][j][g] << endl;
                }
            }
        }
    }
}

#endif // TRAJECTORY_H

