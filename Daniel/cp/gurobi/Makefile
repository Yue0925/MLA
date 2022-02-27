FLAGS = -g -Wall -Wextra 
GRBPATH= /opt/gurobi811/


all:CuttingPlanesEngine.o main.cpp 
	g++ $(FLAGS) CuttingPlanesEngine.o main.cpp -o main -I$(GRBPATH)/linux64/include -L$(GRBPATH)/linux64/lib -lgurobi_c++ -lgurobi81
CuttingPlanesEngine.o : CuttingPlanesEngine.cpp CuttingPlanesEngine.h
	g++ -c $(FLAGS) CuttingPlanesEngine.cpp -o CuttingPlanesEngine.o -I$(GRBPATH)/linux64/include -L$(GRBPATH)/linux64/lib -lgurobi_c++ -lgurobi81
howToExportLibPath: 
	write this in terminal "export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(GRBPATH)/linux64/lib/"
	write this in terminal "export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/opt/gurobi811/linux64/lib/"
	./main
clean: 
	rm -f main *.o *~
run:
	export LD_LIBRARY_PATH=/opt/gurobi811/linux64/lib/
	./main
