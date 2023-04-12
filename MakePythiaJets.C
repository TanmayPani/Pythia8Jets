#include "FastJetMaker.h"
#include "PythiaMaker.h"

#include "Pythia8Plugins/FastJet3.h" // to bridge between Pythia8 particle and FastJet3 PseudoJet classes

bool WeakDecaysOff = true;
float R = 0.4;
float minConstPt = 2.0;
JetInitPartons jet_partons = kAll;
Tune tune = kDetroit;
bool hadronize = true;

void MakePythiaJets(string name = "test_", int PtHatLow = 20, int PtHatHigh = 25, int NEvs = 1000) {

	FastJetMaker *FJ = new FastJetMaker(antikt_algorithm, R, BIpt2_scheme); //Define FastJet JetDefinition
	FJ->SetConstituentPtMin(minConstPt);
	FJ->SetConstituentAbsEtaMax(1.0);
	FJ->SetJetPtMin(0.0);
	FJ->SetJetAbsEtaMax(1.0-R);

	PythiaMaker *pyt = new PythiaMaker(tune, hadronize, NEvs); // Define Pythia object.
	pyt->SetPtHatRange(PtHatLow, PtHatHigh);
	pyt->TurnOnHardProcesses(jet_partons);
	if(WeakDecaysOff)pyt->SetWeakDecaysOff();
	pyt->Initialize(); // Initialize

	TreeMaker *tree = new TreeMaker(name+to_string(PtHatLow)+"_"+to_string(PtHatHigh)+hardQCD[jet_partons]);
	tree->AddBranch()
	
	int NEvents = pyt->GetNEvents(); // The number of events to run.

	for(int iEvent = 0; iEvent < NEvents; iEvent++){//Event Loop
		FJ->EmptyAllVectors();

		if(!pyt->NextEvent()) {
			continue;
		}
		for(int i=1; i<pyt->GetNParticles(); i++){//Particle loop
			//delete particle;
			Particle *particle = pyt->GetParticle(i);
			if(!particle->isFinal())continue;
			//cout<<"Input# "<<i<<" "<<particle->pT()<<" "<<particle->eta()<<" "<<particle->phi()<<endl;
			PseudoJet fj_input = *particle;
			fj_input.set_user_info(new UserInfo(i, particle->charge()));
			FJ->InputForClustering(fj_input);
			//cout<<Px<<" "<<Py<<" "<<Pz<<" "<<E<<endl;
		}

		vector<PseudoJet>* Jets = static_cast<vector<PseudoJet>*>(FJ->GetClusteredJets());
		for(int j = 0; j<Jets->size(); j++){//jet loop
			PseudoJet *jet = static_cast<PseudoJet*>(&(Jets->at(j)));
			//cout<<jet->pt()<<" "<<jet->eta()<<" "<<jet->phi()<<endl;

			vector<PseudoJet>* Constituents = FJ->GetConstituents(jet);
			//cout<<Constituents->size()<<endl;
			for(int c = 0; c<Constituents->size(); c++){//Constituent loop
				PseudoJet *constit = static_cast<PseudoJet*>(&(Constituents->at(c)));
				//cout<<"    "<<constit->user_info<Particle>().index()<<" "<<constit->pt()<<" "<<constit->eta()<<" "<<constit->phi()<<endl;
			}
		}
	} // End event loop.
	// --- Calculate final statistics ---
	pyt->ShowStats();


}


