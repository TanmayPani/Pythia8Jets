R__ADD_INCLUDE_PATH(/usr/local/include)
R__LOAD_LIBRARY(/usr/local/lib/libfastjet.dylib)

#ifndef JetMaker_h
#define JetMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "Math/Vector4D.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

using namespace fastjet;
using namespace ROOT::Math;

class ParticleInfo: public PseudoJet::UserInfoBase{//User Info class to include non-kinematic details of the particles into the pseudojet objects...
public:
	ParticleInfo(const int indx, const short & ch):_charge(ch), _index(indx){}
   short getIndex() const {return _index;}
	short getCharge() const {return _charge;}
protected:
   int _index;
	short _charge;
};

class JetMaker {
private :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           EventID;
   Float_t         PtHat;
   //Float_t         Weight;
   //Float_t         NFinalState;
   vector<int>     *ParticlePID;
   vector<short>   *ParticleCharge;
   vector<float>   *ParticleMass;
   vector<LorentzVector<PxPyPzE4D<double> > > *Particle4P;

   // List of branches
   TBranch        *b_EventId;   //!
   TBranch        *b_PtHat;   //!
   //TBranch        *b_Weight;   //!
   //TBranch        *b_NFinal;   //!
   TBranch        *b_ParticlePID;   //!
   TBranch        *b_ParticleCharge;   //!
   TBranch        *b_ParticleMass;   //!
   TBranch        *b_Particle4P;   //!

   JetDefinition *jet_def;
   //JetDefinition jet_def;
   AreaDefinition *area_def;
   GhostedAreaSpec *area_spec;
   JetAlgorithm	algo;
   RecombinationScheme reco_sch;
   AreaType area_type;

   static const int NPtHatBins = 6; 

   int PtHats[NPtHatBins+1] = {15, 20, 25, 35, 45, 55, -1};
   float xsecs[NPtHatBins] = {0};
   float xsecs_err[NPtHatBins] = {0};
   float NEvents[NPtHatBins] = {0};
   float W[NPtHatBins] = {0};
	
   float R;
   float MaxRap;
   float GhostArea; 
   float etaRange;
   float minJetPt;
   float maxJetPt;
   float minTrackPt;
   float minJetNCh;

   TFile *treeoutFile;
   TTree *treeout;

   float Weight    ;
   int NJets;
   int NParticles;

   float JetE[20]     ;
   float JetPt[20]     ;
   float JetEta[20]    ;
   float JetPhi[20]    ;
   float JetArea[20]    ;
   float JetPtD [20]   ;
   float JetLeSub[20]  ;
   float JetGirth[20]  ;
   int JetNConst[20] ;
   int JetNCh[20] ;

   float TrackE[5000]   ;
   float TrackPt[5000]   ;
   float TrackEta[5000]  ;
   float TrackPhi[5000]  ;
   short TrackCharge[5000];
   int TrackJetID[5000];

   TH1F *hPt;
   TH1F *hEta;
   TH1F *hPhi;
   TH1F *hArea;

   TH2F *h2PtvPtD;
   TH2F *h2PtvGirth;
   TH2F *h2PtvLeSub;

public:   
   JetMaker();
   virtual ~JetMaker();
   virtual void Loop();

protected:
   void Init(TTree *tree);
   void BookJetTree();
   float GetdR(float a, float b, float c, float d);
   void DeclareHistograms();
   void WriteOutput();
};

#endif

#ifdef JetMaker_cxx
JetMaker::JetMaker(){
	//Define jet clustering parameters
	algo = antikt_algorithm;
   //reco_sch = E_scheme;
	reco_sch = BIpt2_scheme;
	area_type = active_area_explicit_ghosts;
	R = 0.4;
	MaxRap = 1.2;
//	GhostArea = 0.01;

	//Setup jet clustering algorithms
	area_spec = new GhostedAreaSpec(MaxRap);
	area_def = new AreaDefinition(area_type, *area_spec);
	jet_def = new JetDefinition(algo, R, reco_sch, Best);

	cout<<jet_def->description()<<endl;

	ifstream weights("XSecGenWeightspSet17Gluons.dat", ios::in);
	for(int i = 0; i<NPtHatBins; i++){
		weights>>NEvents[i]>>xsecs[i]>>xsecs_err[i];
		
		W[i] = xsecs[i]/NEvents[i];
	}

}

JetMaker::~JetMaker()
{
                       
                                   
}

void JetMaker::Init(TTree *tree)
{
   // Set object pointer
   ParticlePID = 0;
   ParticleCharge = 0;
   ParticleMass = 0;
   Particle4P = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;

   tree->SetBranchAddress("EventID", &EventID, &b_EventId);
   tree->SetBranchAddress("PtHat", &PtHat, &b_PtHat);
   //tree->SetBranchAddress("Weight", &Weight, &b_Weight);
   //tree->SetBranchAddress("NFinalState", &NFinalState, &b_NFinal);
   tree->SetBranchAddress("ParticlePID", &ParticlePID, &b_ParticlePID);
   tree->SetBranchAddress("ParticleCharge", &ParticleCharge, &b_ParticleCharge);
   tree->SetBranchAddress("ParticleMass", &ParticleMass, &b_ParticleMass);
   tree->SetBranchAddress("Particle4P", &Particle4P, &b_Particle4P);
}

void JetMaker::BookJetTree(){
   treeoutFile = new TFile("Pythia8_R04GluonJets_with_Tracks_Full.root", "RECREATE");
   treeout = new TTree("Pythia8_Jets", "Pythia8_Jets");

   treeout->Branch("Weight",      &Weight    , "Weight/F"); 
   treeout->Branch("NJets",      &NJets    , "NJets/I"); 
   treeout->Branch("NParticles",      &NParticles    , "NParticles/I"); 

   treeout->Branch("JetE",       &JetE     , "JetE[NJets]/F"); 
   treeout->Branch("JetPt",       &JetPt     , "JetPt[NJets]/F"); 
   treeout->Branch("JetEta",      &JetEta    , "JetEta[NJets]/F"); 
   treeout->Branch("JetPhi",      &JetPhi    , "JetPhi[NJets]/F"); 
   treeout->Branch("JetArea",      &JetArea    , "JetArea[NJets]/F"); 
   treeout->Branch("JetPtD",      &JetPtD    , "JetPtD[NJets]/F"); 
   treeout->Branch("JetLeSub",    &JetLeSub  , "JetLeSub[NJets]/F"); 
   treeout->Branch("JetGirth",    &JetGirth  , "JetGirth[NJets]/F"); 
   treeout->Branch("JetNConst",   &JetNConst , "JetNConst[NJets]/I"); 
   treeout->Branch("JetNCh",   &JetNCh , "JetNCh[NJets]/I"); 
  
   treeout->Branch("TrackE",     &TrackE    ,   "TrackE[NParticles]/F");
   treeout->Branch("TrackPt",     &TrackPt    ,   "TrackPt[NParticles]/F");
   treeout->Branch("TrackEta",    &TrackEta   ,  "TrackEta[NParticles]/F");
   treeout->Branch("TrackPhi",    &TrackPhi   ,  "TrackPhi[NParticles]/F");
   treeout->Branch("TrackCharge", &TrackCharge, "TrackCharge[NParticles]/S");
   treeout->Branch("TrackJetID",  &TrackJetID, "TrackJetID[NParticles]/I");

}

float JetMaker::GetdR(float jeta, float jphi, float ceta, float cphi){
	float deta = fabs(jeta-ceta);
	float dphi = fabs(jphi-cphi);

	if(dphi>pi)dphi = twopi-dphi;

	return sqrt(deta*deta + dphi*dphi);
}

void JetMaker::DeclareHistograms(){
	static const int NPtDBins = 9;
	static const int NLeSubBins = 5;
	static const int NGirthBins = 14;
   static const int NJetPtBins = 32;

   double JetPtBins[NJetPtBins+1] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 45, 50, 60};
	double PtDBins[NPtDBins+1] = {0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 1.0};
	double LeSubBins[NLeSubBins+1] = {0, 5, 10, 15, 20, 30};
	double GirthBins[NGirthBins+1] = {0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.2};

	hPt = new TH1F("Jet_Pt", "N_{jets} vs p_{T}^{Jet}", 40, 10, 50);
	hEta = new TH1F("Jet_Eta", "N_{jets} vs #eta_{Jet}", 40, -1, 1);
	hPhi = new TH1F("Jet_Phi", "N_{jets} vs #phi_{Jet}", 40, 0, 2.0*pi);
	hArea = new TH1F("Jet_Area", "N_{jets} vs A_{Jet}", 40, 0, pi*0.09);

	h2PtvPtD = new TH2F("Jet_Pt_v_PtD", "p_{T}^{D} vs p_{T}^{Jet}", NJetPtBins, JetPtBins, NPtDBins, PtDBins);
        h2PtvGirth = new TH2F("Jet_Pt_v_Girth", "Girth vs p_{T}^{Jet}", NJetPtBins, JetPtBins, NGirthBins, GirthBins);
        h2PtvLeSub = new TH2F("Jet_Pt_v_LeSub", "LeSub vs p_{T}^{Jet}", NJetPtBins, JetPtBins, NLeSubBins, LeSubBins);
}

void JetMaker::WriteOutput(){
	TFile *fout = new TFile("Output_Pythia8_JetHistos.root", "RECREATE");

	hPt->Write();
	hEta->Write();
	hPhi->Write();
	hArea->Write();

	h2PtvPtD->Write();
	h2PtvLeSub->Write();
	h2PtvGirth->Write();
	
	fout->Write();
	fout->Close();

   treeoutFile->cd();
   treeout->Write();
   treeoutFile->Write();
   treeoutFile->Close();
}


#endif // #ifdef JetMaker_cxx
