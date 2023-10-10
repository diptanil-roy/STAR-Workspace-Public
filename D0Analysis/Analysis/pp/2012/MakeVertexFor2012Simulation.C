#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TList.h"
#include "TRandom3.h"

void MakeVertexFor2012Simulation(){
	TString FileName = "runlistRun12pp_P12id";

	ifstream runnumberfile;

	runnumberfile.open(FileName.Data());

	string STRING;

	TRandom3 *r3 = new TRandom3();

	while(!runnumberfile.eof()) // To get you all the lines.
    {
        getline(runnumberfile,STRING); // Saves the line in STRING.

        TString outfilename = "VertexFiles/";
        outfilename += STRING;
        // outfilename += ".txt";

        ofstream vertexfile;
        vertexfile.open(outfilename.Data(), ios::out | ios::app | ios::binary); 

        for (int i = 1; i <= 1000; i++){
        	vertexfile << STRING << "\t" << i << "\t" << 0.1 << "\t" << -0.1 << "\t" << r3->Uniform(-6., 6.) << endl;
        }
    }


}