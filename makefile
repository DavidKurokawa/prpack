CXX = g++
CXXFLAGS = -Wall -O3 -fopenmp
OBJS = prpack_utils.o prpack_adjacency_list.o prpack_preprocessed_gs_graph.o prpack_preprocessed_schur_graph.o prpack_preprocessed_scc_graph.o prpack_solver.o prpack_result.o prpack_driver.o
PROG = prpack_driver

all: ${PROG}
	
${PROG}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS}
	
test: $(PROG)
	./prpack_driver data/jazz.smat --output=- 2>/dev/null | \
	  python test/checkprvec.py data/jazz.smat -
	./prpack_driver data/jazz.smat --output=- -m sgs 2>/dev/null | \
	  python test/checkprvec.py data/jazz.smat -
	./prpack_driver data/jazz.smat --output=- -m sccgs 2>/dev/null| \
	  python test/checkprvec.py data/jazz.smat - 

clean:
	rm *.o ${PROG} -f

.PHONY: all clean test