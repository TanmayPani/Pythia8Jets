#define JetMaker_cxx
#include "JetMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void JetMaker::Loop()
{
	//JetDefinition jet_def(algo, R, reco_sch, Best);
	Selector eta_acceptance = SelectorAbsEtaMax(1.0-R);
//	Selector rapidity_acceptance = SelectorAbsRapMax(1.2);
	Selector is_ghost = SelectorIsPureGhost();
//	Selector jet_cut = (eta_acceptance && rapidity_acceptance) && !is_ghost;
	Selector jet_cut = eta_acceptance;
	//Selector jet_cut = eta_acceptance && rapidity_acceptance;

	DeclareHistograms();

	//ClusterSequenceArea *CS = nullptr;
	ClusterSequence *CS = nullptr;
	vector<PseudoJet> all_jets;
	vector<PseudoJet> jets;
	vector<PseudoJet>constits; 

   	vector<PseudoJet> all_particles;

   	bool doTest = false;

	//int NEventsTest[NPtHatBins] = {300000, 300000, 200000, 2000000, 100000, 100000};

	BookJetTree();

	for(int pth = 0; pth<NPtHatBins; pth++){
	//for(int pth = 3; pth<4; pth++){
		float Wt = W[pth]/W[NPtHatBins-1];
		if(doTest) Wt = xsecs[pth]/xsecs[NPtHatBins-1];

		cout<<"Doing jets in P_T_HAT "<< PtHats[pth]<<" to "<<PtHats[pth+1]<<endl;
		cout<<"PtHat weight = "<<Wt<<endl;
   
		TFile *fin = new TFile(Form("Pythia8Gluons_STARTune_pSet17_PtHat%i_%i_ALLWkDecoff_PartialEvents.root", PtHats[pth], PtHats[pth+1]));

		TTree *tree = (TTree*)fin->Get("Pythia8");

   		Init(tree);

  		Long64_t nevents = tree->GetEntriesFast();
   		 if(doTest)nevents = 10000;
	
		cout<<"Total events read: "<<nevents<<endl;

   		for (Long64_t event=0; event<nevents; event++) {
        		tree->GetEntry(event);  
             
			all_particles.clear();

			NParticles = Particle4P->size();
	
	//	cout<<"#Particles in event: "<<NParticles<<endl;
		//int Nacc = 0;
				Weight = Wt;

			for(int iPart = 0; iPart<NParticles; iPart++){
				PxPyPzEVector v = Particle4P->at(iPart);
				if(v.Pt() < 0.2)continue;
				short Ch = ParticleCharge->at(iPart);
				//Nacc++;
				//all_particles.push_back(PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()).set_user_info(new ParticleInfo(Ch)));

				     TrackE[iPart] = v.E();
				     TrackPt[iPart] = v.Pt();
					TrackEta[iPart] = v.Eta();
					if(v.Phi()>=0)TrackPhi[iPart] = v.Phi();
					else TrackPhi[iPart] = 2*pi+v.Phi();
				 TrackCharge[iPart] = Ch;
			TrackJetID[iPart] = 0;	
			if(v.Pt() < 2)continue;
                        //	cout<<iPart<<" "<<TrackJetID[iPart]<<" "<<TrackPt[iPart]<<" "<<TrackEta[iPart]<<" "<<TrackPhi[iPart]<<" "<<Ch<<endl;
				PseudoJet particle(v.Px(), v.Py(), v.Pz(), v.E());

				particle.set_user_info(new ParticleInfo(iPart, Ch));
				//cout<<v.Px()<<" "<<v.Py()<<" "<<v.Pz()<<" "<<v.E()<<endl;
				all_particles.push_back(move(particle));
				//cout<<v.Eta()<<" "<<particle.eta()<<endl;
				//if(ParticleCharge->at(iPart) != particle.user_info<ParticleInfo>().getCharge())cout<<"Someting's not OK with PseudoJet user info!! "<<endl;

			}
	
			if(event%10000 == 0)
				cout<<" Events read: "<<event<<endl;

			delete CS;
			//all_jets.clear();
			jets.clear();

		//	CS = new ClusterSequenceArea(all_particles, *jet_def, *area_def);
			CS = new ClusterSequence(all_particles, *jet_def);
			//all_jets = CS->inclusive_jets(10);
			jets = jet_cut(CS->inclusive_jets(10));

			NJets = jets.size();
	//		cout<<"# jets clustered = "<<NJets<<endl;
			for(int iJet = 0; iJet<NJets; iJet++){
				constits.clear();
				constits = jets[iJet].constituents();
				int NConst = constits.size();
				//if(NConst<2)continue;
				
				float jpt = jets[iJet].pt();
				float jeta = jets[iJet].eta();
				float jphi = jets[iJet].phi();
			//	float jarea = jets[iJet].area();

				float pt2sum = 0;
				float ptsum = 0;
				float ptdR = 0;
				float ptLead = 0;
				float ptSubLead = 0;
				
				//float LeSub = 0;
				int Nch = 0;
				for(int iConst = 0; iConst<NConst; iConst++){
					if(!constits[iConst].has_user_info())continue;//no ghosts

					short cch = constits[iConst].user_info<ParticleInfo>().getCharge();
					float cE = constits[iConst].e();
					float cpt = constits[iConst].pt();
					float ceta = constits[iConst].eta();
					float cphi = constits[iConst].phi();

					int uid = constits[iConst].user_info<ParticleInfo>().getIndex();

					TrackJetID[uid] = iJet+1;
	//					cout<<iJet<<" "<<TrackJetID[uid]<<" "<<uid<<" "<<cpt<<" "<<ceta<<" "<<cphi<<" "<<cch<<endl;
					if(cch==0)continue;
					Nch++;

					float dR = GetdR(jeta, jphi, ceta, cphi);

					pt2sum += cpt*cpt;
					ptsum += cpt;

					ptdR += cpt*dR;

					if(cpt>ptLead){
						ptSubLead = ptLead;
						ptLead = cpt;
					}else if(cpt>ptSubLead){
						ptSubLead = cpt;
					}
				}

				//if(Nch<2)continue;

				float ptd = sqrt(pt2sum)/ptsum;
				float girth = ptdR/jpt;
				//sort(pt_sorted.begin(), pt_sorted.end());
				float lesub = ptLead-ptSubLead;

				JetE[iJet] = jets[iJet].e();
				JetPt[iJet] = jpt;
				JetEta[iJet] = jeta;
				JetPhi[iJet] = jphi;
			//	JetArea[iJet] = jarea;
				JetNConst[iJet] = NConst;
				JetNCh[iJet] = Nch;
				JetPtD[iJet] = ptd;
				JetLeSub[iJet] = lesub;
				JetGirth[iJet] = girth;

				treeout->Fill();

				//if(jarea<0.15)continue;

				hPt->Fill(jpt, Wt);
				hEta->Fill(jeta, Wt);
				hPhi->Fill(jphi, Wt);
//				hArea->Fill(jarea, Wt);
				h2PtvPtD->Fill(jpt, ptd, Wt);
				h2PtvGirth->Fill(jpt, girth, Wt);
				h2PtvLeSub->Fill(jpt, lesub, Wt);
			}	

   		}
	}

   WriteOutput();
}
