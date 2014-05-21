CXX		=	g++
O_FLAG	=	-O0
D_FLAG	=	-D_WITH_DEBUG -ggdb -g3 -fvar-tracking-assignments -fno-inline -fno-inline-small-functions -fno-eliminate-unused-debug-types
P_FLAG	=
CXXFLAGS=	-Wall -D_FILE_OFFSET_BITS=64 $(O_FLAG) $(D_FLAG) $(P_FLAG)
LIBS	=	-lm

BINDIR  =   ../bin

PROG	=	K

OBJS = K.o \
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

OBJS_N = K_n.o \
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

HEAD_N = K_n.h \
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

# Makefile uses recursive $(MAKE) to build separate versions.
#
#   all (default):   debug version unoptimised (O_FLAG=-O0)
#   opt:             debug version optimised (O_FLAG=-O2)
#   profile:         default with gprof profiling (P_FLAG=-pg)
#   release:         no debug optimised (D_FLAG= O_FLAG=-O2)
#   release-profile: no debug optimised (D_FLAG= O_FLAG=-O2 P_FLAG=-pg)


all: $(PROG)
	mkdir -p $(BINDIR)
	mv $? $(BINDIR)

K: $(OBJS) version.h
	$(CXX) $(LDFLAGS) $(OBJS) -o $@ $(LIBS)

Kn: $(OBJS) $(OBJS_N) version.h
	$(CXX) $(LDFLAGS) $(OBJS) $(OBJS_N) -o $@ $(LIBS) 

$(OBJS): $(HEAD)
$(OBJS_N): $(HEAD_N)

.cpp.o:	$(HEAD) $(HEAD_N)

K_cmdline.h:	SimpleOpt.h

K.h K_n.h: version.h

version.h: .FORCE
	./git-getversion.sh > version.h
	echo "#define CXX_VERSION \""`$(CXX) --version | head -n 1`"\"" >> version.h
	echo "#define CXXFLAGS \"$(CXXFLAGS)\"" >> version.h

.FORCE:

opt: .FORCE
	$(MAKE) clean
	$(MAKE) O_FLAG=-O2

K-opt:

profile: .FORCE
	$(MAKE) clean
	$(MAKE) P_FLAG=-pg


#---------------------------  Other targets


clean:
	rm -f gmon.out $(OBJS) $(OBJS_N) $(PROG) version.h


clean-all: clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


