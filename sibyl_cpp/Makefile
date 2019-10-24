ROOTLIBS := $(shell root-config --cflags --libs)
RATLIBS := -L$(RATROOT)/lib -lRATEvent -I$(RATROOT)/include

all:
	g++ -shared -Wl,-soname,libfastrat.so $(ROOTLIBS) $(RATLIBS) -o libfastrat.so -fPIC fastrat.cpp -lRATEvent

macos:
	g++ -shared  $(ROOTLIBS) $(RATLIBS) -o libfastrat.so -Wl,-rpath,$(shell root-config --libdir) -fPIC fastrat.cpp

clean:
	rm -f libfastrat.so
