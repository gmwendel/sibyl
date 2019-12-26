ROOTLIBS := $(shell root-config --cflags --libs)
RATLIBS := -L$(RATROOT)/lib -lRATEvent -I$(RATROOT)/include

all:
	g++ $(ROOTLIBS) $(RATLIBS) -o treescene generateTreeScene.cpp
