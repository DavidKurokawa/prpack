IGRAPH_SUPPORT ?= 0

CXX = g++
CXXFLAGS = -Wall -g -fopenmp 
LDFLAGS =
OBJS = prpack_utils.o \
    prpack_base_graph.o \
    prpack_preprocessed_ge_graph.o \
    prpack_preprocessed_gs_graph.o \
    prpack_preprocessed_schur_graph.o \
    prpack_preprocessed_scc_graph.o \
    prpack_solver.o \
    prpack_solver_ge.o \
    prpack_solver_sccgs.o \
    prpack_solver_inout.o \
    prpack_result.o \
    prpack_driver.o \
    prpack_driver_benchmark.o
PROG = prpack_driver

ifeq ($(IGRAPH_SUPPORT),1)
	OBJS += prpack_igraph_graph.o
	CXXFLAGS += $(shell pkg-config igraph --cflags) -DPRPACK_IGRAPH_SUPPORT
	LDFLAGS += $(shell pkg-config igraph --libs)
endif

all: ${PROG}
	
${PROG}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LDFLAGS}
	
test: $(PROG)
	./prpack_driver data/jazz.smat --output=- 2>/dev/null | \
	  python test/checkprvec.py data/jazz.smat -
	./prpack_driver data/jazz.smat --output=- -m sgs 2>/dev/null | \
	  python test/checkprvec.py data/jazz.smat -
	./prpack_driver data/jazz.smat --output=- -m sccgs 2>/dev/null| \
	  python test/checkprvec.py data/jazz.smat - 
	./prpack_driver data/jazz.smat --output=- -m sccgs -w 2>/dev/null| \
	  python test/checkprvec.py data/jazz.smat - 
	./prpack_driver data/power.smat -w --output=- -m sccgs 2>/dev/null| \
	  python test/checkprvec.py data/power.smat - 
	./prpack_driver data/netscience.smat -w --output=- -m sccgs 2>/dev/null| \
	  python test/checkprvec.py data/netscience.smat - 
	./prpack_driver data/wb-cs.stanford.smat --output=- -m sccgs_uv \
		  -a 0.5 -v test/csstan-v.vec -u test/csstan-u.vec 2>/dev/null \
	  | python test/checkprvec.py data/wb-cs.stanford.smat - \
		  -a 0.5 -v test/csstan-v.vec -u test/csstan-u.vec
test_inout: $(PROG)
	./prpack_driver data/jazz.smat --output=- -m inout 2>/dev/null| \
	  python test/checkprvec.py data/jazz.smat - 
	./prpack_driver data/wb-cs.stanford.smat --output=- -m inout 2>/dev/null| \
	  python test/checkprvec.py data/wb-cs.stanford.smat - 
	./prpack_driver data/wb-cs.stanford.smat --output=- -m inout \
		  -a 0.5 -v test/csstan-v.vec  2>/dev/null \
	  | python test/checkprvec.py data/wb-cs.stanford.smat - \
		  -a 0.5 -v test/csstan-v.vec 	  
	./prpack_driver data/wb-cs.stanford.smat --output=- -m inout \
		  -a 0.5 -u test/csstan-u.vec 2>/dev/null \
	  | python test/checkprvec.py data/wb-cs.stanford.smat - \
		  -a 0.5 -u test/csstan-u.vec
	./prpack_driver data/wb-cs.stanford.smat --output=- -m inout \
		  -a 0.5 -v test/csstan-v.vec -u test/csstan-u.vec 2>/dev/null \
	  | python test/checkprvec.py data/wb-cs.stanford.smat - \
		  -a 0.5 -v test/csstan-v.vec -u test/csstan-u.vec
	
perf: $(PROG)
	./prpack_driver ?

matlab:	
	cd matlab; make

clean:
	$(RM) *.o ${PROG}
	

.PHONY: all clean test matlab
