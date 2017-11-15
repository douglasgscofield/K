#include "K.h"
#include "Trajectory.h"

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

const struct Trajectory::struct_types Trajectory::types[] = {
    {"*error*",     "TRAJECTORY_ERROR",       TRAJECTORY_ERROR},
    {"end_of_gen",  "TRAJECTORY_END_OF_GEN",  TRAJECTORY_END_OF_GEN},
};
const int Trajectory::num_types = sizeof(Trajectory::types) / sizeof(Trajectory::struct_types);


/////////////////////////////////////////////////////////////////
void Trajectory::start     (KConfig K, 
                            const string& fil, 
                            const string& typ,
                            int frq)
{
    const char* thisfunction = "Trajectory::start";
    if (debug_flag) {
        cerr << thisfunction << " gen=" << K->generation << " file='" << fil << "' type='" << typ << "' freq=" << frq << endl;
    }
    filename = fil;
    file.open(filename.c_str());
    if (file.fail()) {
        cerr << thisfunction << ": could not open trajectory file " << filename << endl;
        fatal(NULL);
    }
    if ((type = match_type(typ)) == TRAJECTORY_ERROR) {
        cerr << thisfunction << ": unrecognized trajectory type '" << type << "'" << endl;
        fatal(NULL);
    }
    if ((freq = frq) <= 0) {
        cerr << thisfunction << ": trajectory frequency (" << freq << ") must be > 0" << endl;
        fatal(NULL);
    }
    // write trajectory header
    active_flag = true;
    gen_started = K->generation;
    num_recorded = 0;
    writenow_flag = true;  // write a trajectory immediately
}


/////////////////////////////////////////////////////////////////
void Trajectory::stop ()
{
    const char* thisfunction = "Trajectory::stop";
    if (debug_flag) {
        cerr << thisfunction << endl;
    }
    if (! active_flag) {
        cerr << thisfunction << ": no active trajectory" << endl;
        fatal(NULL);
    }
    reset();
    active_flag = false;
}

/* should there be a stack of trajectories? */
/*
 * routine should include means of specifying genotypes to track,
 * as well as frequency of tracking.  specify trajectory file, or
 * to stdout.
 */

/////////////////////////////////////////////////////////////////
// static functions for managing the type of trajectory
void Trajectory::dump_type ()
{
    //const char* thisfunction = "Trajectory::dump_type";
    for (int i = 0; i < num_types; ++i) {
        cout << "type " << i 
            << " name=<" << types[i].name
            << "> define_name=<" << types[i].define_name
            << "> type=" << types[i].type << endl;
    }
}

/////////////////////////////////////////////////////////////////
const string& Trajectory::match_type (int type)
{
    //const char* thisfunction = "Trajectory::match_type";
    for (int i = 0; i < num_types; ++i) {
        if (type == types[i].type) {
            return(types[i].name);
        }
    }
    return( types[0].name );
}

/////////////////////////////////////////////////////////////////
int Trajectory::match_type (const string& type)
{
    //const char* thisfunction = "Trajectory::match_type";
    for (int i = 0; i < num_types; ++i) {
        if (type == types[i].name) {
            return(types[i].type);
        }
    }
    return( TRAJECTORY_ERROR );
}

/////////////////////////////////////////////////////////////////
void Trajectory::reset() 
{
    //const char* thisfunction = "Trajectory::reset";
    active_flag = false;
    filename = "";
    if (file)
        file.close();
    type = TRAJECTORY_ERROR;
    freq = 0;
    num_recorded = 0;
    next_gen = 0;
    gen_started = 0;
    gen_to_end = 0;
    dropzero_flag = true;
    writenow_flag = false;
    debug_flag = false;
}

