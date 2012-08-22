
OBJ  = K.o K_cmdline.o K_debug.o K_equilibrium.o K_genotypes.o K_initiate.o K_math.o K_mutation.o K_nextgen.o K_reproduction.o K_savefile.o K_selection.o K_stats.o K_util.o main.o K_dropin.o Trajectory.o $(RES)
HEAD = K.h K_cmdline.h K_debug.h K_equilibrium.h K_genotypes.h K_initiate.h K_math.h K_mutation.h K_nextgen.h K_reproduction.h K_savefile.h K_selection.h K_stats.h K_util.h Trajectory.h
LIBS = -lm
BIN  = K
CXXFLAGS = -g3
CFLAGS = -g3

all: K


clean:
	rm -f $(OBJ) $(BIN)

$(BIN): $(OBJ)
	g++ $(OBJ) -o "K" $(LIBS)

.cpp.o:	$(HEAD)

