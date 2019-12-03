.DELETE_ON_ERROR:

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

INC_DIR   := -I.

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11 #-g -O0
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)

all: checkdirs slib/libBSA_cls.so run_BSA

checkdirs: dict slib

dict slib:
	@mkdir -p $@

%: %.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -Lslib -lBSA_cls

%.o: %.cxx 
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(INC_DIR)


dict/dictBSA.cxx: BSA_survey_cls.h 
	rootcling -f $@ -c $(ROOTCFLAGS) $(HIPOCFLAGS) -p $^


slib/libBSA_cls.so: BSA_survey_cls.cxx
	g++ -shared -fPIC -o $@ $(ROOTLDFLAGS) $(ROOTCFLAGS) $(HIPOCFLAGS) $(INC_DIR) $(ROOTCFLAGS) $^
#	cp dict/datadict_rdict.pcm slib/.

