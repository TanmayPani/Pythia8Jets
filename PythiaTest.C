R__ADD_INCLUDE_PATH(/Users/tanmaypani/PYTHIA8/pythia8309/include)
//R__ADD_INCLUDE_PATH(/usr/local/include)
R__LOAD_LIBRARY(/Users/tanmaypani/PYTHIA8/pythia8309/lib/libpythia8.dylib)
//R__LOAD_LIBRARY(/usr/local/lib/libfastjet.dylib)

#include "Pythia8/Pythia.h" // access to Pythia objects.
//#include "fastjet/ClusterSequence.hh"
#include "FastJet_Input_Struct.h"

using namespace Pythia8; // allow simplified notation.
//using namespace fastjet;
void PythiaTest(int PtHatLow = 20, int PtHatHigh = 25, string decay_switch = "off", bool StoreFullEvent = false, string MergeWeightsFile = "testMergeOut.dat") {

	string EventStorage = "PartialEvents";
	if(StoreFullEvent)EventStorage="FullEvents";

	int maxEvents = 0;

	if(PtHatLow == 15)maxEvents = 1500000;
	if(PtHatLow == 20)maxEvents = 1500000;
	if(PtHatLow == 25)maxEvents = 1000000;
	if(PtHatLow == 35)maxEvents = 1000000;
	if(PtHatLow == 45)maxEvents = 500000;
	if(PtHatLow == 55)maxEvents = 500000;

	TFile *f = new TFile(Form("Pythia8Gluons_STARTune_pSet17_PtHat%i_%i_ALLWkDec%s_%s.root", PtHatLow, PtHatHigh, decay_switch.c_str(), EventStorage.c_str()), "RECREATE");
	TTree *t = new TTree("Pythia8", "Pythia8");

	FastJet_Input FJ_Input;

	if(!StoreFullEvent){
		t->Branch("EventID", &FJ_Input.EventId, "EventId/I");
		t->Branch("PtHat", &FJ_Input.PtHat, "PtHat/F");
		t->Branch("Weight", &FJ_Input.Weight, "Weight/F");
		t->Branch("NFinalState", &FJ_Input.NFinal, "NFinal/F");
		t->Branch("ParticlePID", &FJ_Input.ParticlePID);
		t->Branch("ParticleCharge", &FJ_Input.ParticleCharge);
		t->Branch("ParticleMass", &FJ_Input.ParticleMass);
		t->Branch("Particle4P", "vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&FJ_Input.Particle4Mom, 128000, 0);//Dont split this into components
	}

	Pythia *pythia = new Pythia(); // Define Pythia object.
	Event *MCEvent = static_cast<Event*>(&(pythia->event));

	if(StoreFullEvent)t->Branch("PythiaEvent", &MCEvent);

	pythia->readFile("detroit.cmnd");

	pythia->readString(Form("PhaseSpace:pTHatMin = %i", PtHatLow));           // low pthat cut
	pythia->readString(Form("PhaseSpace:pTHatMax = %i" ,PtHatHigh));           // high pthat cut

	//Switch off weak decay of the following particles if comparing with Embedding as they decay before they reach the detector 
 	pythia->readString(Form("111:mayDecay =  %s", decay_switch.c_str()));//pi0                                                                                                                                      
	pythia->readString(Form("211:mayDecay =  %s", decay_switch.c_str()));//pi+                                                                                                                                      
  	pythia->readString(Form("221:mayDecay =  %s", decay_switch.c_str()));//eta                                                                                                                                      
  	pythia->readString(Form("321:mayDecay =  %s", decay_switch.c_str()));//K+                                                                                                                                       
  	pythia->readString(Form("310:mayDecay =  %s", decay_switch.c_str()));//Kshort                                                                                                                                   
  	pythia->readString(Form("130:mayDecay =  %s", decay_switch.c_str()));//Klong                                                                                                                                    
  	pythia->readString(Form("3122:mayDecay = %s", decay_switch.c_str()));//Lambda0                                                                                                                                  
  	pythia->readString(Form("3212:mayDecay = %s", decay_switch.c_str()));//Sigma0                                                                                                                                   
  	pythia->readString(Form("3112:mayDecay = %s", decay_switch.c_str()));//Sigma-                                                                                                                                   
  	pythia->readString(Form("3222:mayDecay = %s", decay_switch.c_str()));//Sigma+                                                                                                                                   
  	pythia->readString(Form("3312:mayDecay = %s", decay_switch.c_str()));//Xi-                                                                                                                                      
  	pythia->readString(Form("3322:mayDecay = %s", decay_switch.c_str()));//Xi0                                                                                                                                      
  	pythia->readString(Form("3334:mayDecay = %s", decay_switch.c_str()));//Omega-

  	pythia->readString(Form("Main:numberOfEvents = %i", maxEvents));//Maximum events to be generated

	pythia->init(); // Initialize
	
	Pythia8::Info *GenInfo = const_cast<Pythia8::Info*>(&(pythia->info));

	int NEvents = pythia->mode("Main:numberOfEvents"); // The number of events to run.

	cout<<"Generating "<<NEvents<<" in pthat "<<PtHatLow<<" to "<<PtHatHigh<<endl;
	// --- The event loop ---
	for(int iEvent = 0; iEvent < NEvents; iEvent++){

		if(!pythia->next()) {
			continue;
		}if(StoreFullEvent){
			t->Fill();
			continue;
		}else{
		FJ_Input = {0};
		FJ_Input.EventId = iEvent;
		float PtHat =  GenInfo->pTHat();
		float W =  GenInfo->weight(0);
		float NFinal = MCEvent->nFinal();

		FJ_Input.PtHat = PtHat;
		FJ_Input.Weight = W;
		FJ_Input.NFinal = NFinal;

		for(int i=1; i<MCEvent->size(); i++){
			//delete particle;
			Particle *particle = static_cast<Particle*>(&(MCEvent->at(i)));

			if(!particle->isFinal())continue;
			if(fabs(particle->eta())>1.0)continue; //For fair comparision to detector level...
		//	if(particle->pT()>30.0)continue;

			int PID = particle->id();
			short Charge = particle->charge();
			float Mass = particle->m();
			float Px = particle->px();
			float Py = particle->py();
			float Pz = particle->pz();
			float E = particle->e();

		//	cout<<Px<<" "<<Py<<" "<<Pz<<" "<<E<<endl;

			FJ_Input.ParticlePID.push_back(PID);
			FJ_Input.ParticleCharge.push_back(Charge);
			FJ_Input.ParticleMass.push_back(Mass);
			
			PxPyPzEVector P(Px, Py, Pz, E);
	//		cout<<particle->eta()<<"  "<<P.Eta()<<endl;
			FJ_Input.Particle4Mom.push_back(P);
		}
//		cout<<MCEvent->nFinal()<<" "<<FJ_Input.Particle4Mom.size()<<endl;
		t->Fill();
		}
	} // End event loop.
	// --- Calculate final statistics ---
	pythia->stat();

//	f->cd();
	t->Write();
	f->Write();
	f->Close();

//	delete t;
//	delete f;

	ofstream out(MergeWeightsFile.c_str(), ios_base::app);
	out<<GenInfo->weightSum()<<" "<<GenInfo->sigmaGen()<<" "<<GenInfo->sigmaErr()<<endl;
}


