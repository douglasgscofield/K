
OBJ  = K.o \
	   K_cmdline.o \
	   K_debug.o \
	   K_equilibrium.o \
	   K_genotypes.o \
	   K_initiate.o \
	   K_math.o \
	   K_mutation.o \
	   K_nextgen.o \
	   K_reproduction.o \
	   K_savefile.o \
	   K_selection.o \
	   K_stats.o \
	   K_util.o \
	   main.o \
	   Trajectory.o

HEAD = K.h \
	   K_cmdline.h \
	   K_debug.h \
	   K_equilibrium.h \
	   K_genotypes.h \
	   K_initiate.h \
	   K_math.h \
	   K_mutation.h \
	   K_nextgen.h \
	   K_reproduction.h \
	   K_savefile.h \
	   K_selection.h \
	   K_stats.h \
	   K_util.h \
	   Trajectory.h

OBJN = K_n.o \
	   K_cmdline_n.o \
	   K_debug_n.o \
	   K_equilibrium_n.o \
	   K_initiate_n.o \
	   K_mutation_n.o \
	   K_nextgen_n.o \
	   K_prevgen_n.o \
	   K_reproduction_n.o \
	   K_savefile_n.o \
	   K_selection_n.o \
	   K_stats_n.o \
	   K_util_n.o

HEADN = K_n.h \
	   K_cmdline_n.h \
	   K_debug_n.h \
	   K_equilibrium_n.h \
	   K_genotypes.h \
	   K_initiate_n.h \
	   K_math.h \
	   K_mutation_n.h \
	   K_nextgen_n.h \
	   K_prevgen_n.h \
	   K_reproduction_n.h \
	   K_savefile_n.h \
	   K_selection_n.h \
	   K_stats_n.h \
	   K_util.h \
	   Trajectory.h

LIBS = -lm

BIN  = K Kn

CXXFLAGS = -g3

CFLAGS = -g3

all: K


clean:
	rm -f $(OBJ) $(OBJN) $(BIN)

K: $(OBJ)
	g++ $(OBJ) -o "K" $(LIBS)

Kn: $(OBJ) $(OBJN)
	g++ $(OBJ) $(OBJN) -o "Kn" $(LIBS)

$(OBJ):	$(HEAD)

.cpp.o:	$(HEAD) $(HEADN)

