



FLAGS  = -O3 -Wall -Wno-sign-compare
INCS = -I/usr/include $(shell root-config --cflags)
LIBS = $(shell root-config --glibs) -lRooFit -lRooFitCore -lRooStats
OBJ= $(SRC:%.C=%.o)
SRC=PLRun3.C
EXE=runMass

all: $(EXE)

$(EXE): $(OBJ) $(SRC) 
	g++ $(OBJ) -o $@ $(FLAGS) $(INCS) $(LIBS)

$(OBJ): %.o: %.C Makefile
	g++ $(FLAGS) -c $(INCS) -o $@ $<

.PHONY : clean
clean: 
	rm -f $(OBJ) $(EXE) *.so *.d
