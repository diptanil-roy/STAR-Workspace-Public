#include "D0Filter.h"

#include "StarGenerator/EVENT/StarGenEvent.h"
#include "StarGenerator/EVENT/StarGenParticle.h"


//______________________________________________________________________________________________
D0Filter::D0Filter( const char* name ) : StarFilterMaker(name)
{

}
//______________________________________________________________________________________________
void D0Filter::AddTrigger( int id, double mnpt, double mxpt, double mneta, double mxeta, int idm )
{
  mTriggers.push_back( Trigger_t( id, mnpt, mxpt, mneta, mxeta, idm ) );
}
//______________________________________________________________________________________________
int D0Filter::Init() {

  return kStOK;

}

//______________________________________________________________________________________________
int D0Filter::Filter( StarGenEvent *_event ) 
{
  cout << "Filter Called" << endl;
  // Get a reference to the current event
  StarGenEvent& event = (_event)? *_event : *mEvent;

  int npart = event.GetNumberOfParticles();

  cout << "Number of Particles = " << npart << endl; 

  // Ensure a D0 is w/in our acceptance
  bool go = false;
  for ( int ipart=1 /*skip header*/; 
	ipart<npart; 
	ipart++ )
    {
      StarGenParticle *part = event[ipart];
      if ( part->GetStatus() != StarGenParticle::kFinal ) continue; // Must be a final state particle

      //      part->Print();

      TLorentzVector momentum = part->momentum();
      double pT  = momentum.Perp();
      double eta = momentum.Eta();

      // Make sure we're looking at something with no color
      if ( TMath::Abs(part->GetId()) != 421 ) continue;
      if ( pT < 1.00 || TMath::Abs(eta) > 1 ) continue;

      go = true;
      break;

    }

  //  if ( false == go ) LOG_INFO << "No candidate D0, reject event" << endm;
  

  if ( false == go ) return StarGenEvent::kReject;

  LOG_INFO << "Event Accepted." << endm;

  return StarGenEvent::kAccept;
}
