#include <vector>
#ifdef __ROOTCLING__
// #pragma link C++ class treeclass+;
#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;
#endif

//rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h
//g++ -shared -o libmydict.so mydict.cxx `root-config --cflags --libs`