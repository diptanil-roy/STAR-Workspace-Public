#include "Particle.h"
#include <vector>

#ifdef __ROOTCLING__
#pragma link C++ class vector<Particle>+;
#endif

//rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h
//g++ -shared -o libmydict.so mydict.cxx `root-config --cflags --libs`
