#plik make używający skryptu komilujacego pakietu octave

CXX=g++
CXXFLAGS=-Wall -std=c++0x -c

OCT_CXX=mkoctfile
OCT_LDFLAGS=-lfftw3

NLS_solver.oct : NLS_solver.cpp solve.o input_test.cpp
	$(OCT_CXX) $(OCT_LDFLAGS) NLS_solver.cpp solve.o input_test.cpp

clean : 
	rm solve.o NLS_solver.o input_test.o NLS_solver.oct

solve.o: solve.cpp
	$(CXX) $(CXXFLAGS) solve.cpp

install :
	mv NLS_solver.oct ../octave_scripts
