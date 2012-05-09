CXXFLAGS = -Wall -O3 -fopenmp 

FILES = utils.cpp ../*.cpp ../*.h

all : use_pagerank_solver.mex* create_pagerank_solver.mex*

use_pagerank_solver.mex* : use_pagerank_solver.cpp $(FILES)
	mex -v -I../ use_pagerank_solver.cpp utils.cpp ../*.cpp CXXFLAGS="\$$CXXFLAGS $(CXXFLAGS)" LDFLAGS="\$$LDFLAGS $(CXXFLAGS)"

create_pagerank_solver.mex* : create_pagerank_solver.cpp $(FILES)
	mex -v -I../ create_pagerank_solver.cpp utils.cpp ../*.cpp CXXFLAGS="\$$CXXFLAGS $(CXXFLAGS)" LDFLAGS="\$$LDFLAGS $(CXXFLAGS)"
    
clean:
	$(RM) *.mex*
    
test: 
	matlab -nodesktop -nosplash -r test_main

.PHONY : all clean test
