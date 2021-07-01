LIBS 	 := -L/usr/lib64/liblapack.a -llapack -L/usr/lib/atlas -L/<hypre_path>/src/lib/ -lHYPRE -atlas -lblas -lg2c -lrt 
CXXFLAGS := -DENV_MPI -O3 -Wall -DMPICH_IGNORE_CXX_SEEK -I/<hypre_path>/src/hypre/include/
CXX      := mpicxx
EXT      := cpp
TARGET   := ../bin/fspai-1.0
SOURCES  := $(wildcard *.$(EXT) sequential/*.$(EXT) parallel/*.$(EXT))
OBJECTS  := $(patsubst %.$(EXT), %.o, $(SOURCES))

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

$(OBJECTS): %.o: %.$(EXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf *.o
