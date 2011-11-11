#plik make używający skryptu komilujacego pakietu octave

NLS_sovler.oct : NLS_solver.cpp input_test.cpp solve.cpp
	mkoctfile NLS_solver.cpp input_test.cpp solve.cpp

clean : 
	rm solve.o input_test.o NLS_solver.o

install :
	mv NLS_solver.oct ../octave_scripts
