R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/FASTJET/fastjet-install/lib/libfastjet.dylib);

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/PseudoJet.hh"

using namespace fastjet;
using namespace std;

void FastJetReader(){
    
    PseudoJet p(   2.0,  0.1,  0, 3.0);
    cout << "Constituent " << p.pt() << endl;

}