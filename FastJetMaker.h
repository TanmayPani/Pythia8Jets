//Wrapper class for using the FastJet library...
//Uses ROOT6 classes and compiler setup for a few functionalities, but can be made ROOT independent easily

//Load fastjet libraries and main fastjet folders
R__ADD_INCLUDE_PATH(/usr/local/include)
R__LOAD_LIBRARY(/usr/local/lib/libfastjet.dylib)

#ifndef FastJetMaker_h
#define FastJetMaker_h

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

using namespace fastjet;

class FastJetMaker{
    public:
        FastJetMaker(); //default constructor, antikt, R=0.4, BIpt2 recombination
        FastJetMaker(JetAlgorithm algo, float R, RecombinationScheme reco_sch);
        virtual ~FastJetMaker();

        void FastJetAreaMaker(AreaType area_type, float MaxRap); //Will define an area for the jet cone, 
        //Calculating areas will increase computation time by a LOT, do this only for data and detector level
        //NO AREA NEEDED for particle level from MC simulations

        inline void InputForClustering(PseudoJet constit); 
       // inline void InputForClustering(int i, short ch, float px, float py, float pz, float E);
        inline void InputForClustering(int i, short ch, float pt, float eta, float phi, float M);

        vector<PseudoJet>* GetClusteredJets(); //returns a pointer to the vector<PseudoJet> jets defined below
        vector<PseudoJet>* GetConstituents(PseudoJet* jet); //returns a pointer to the vector<PseudoJet> for constituents for a jet

        //Kinematic cuts for constituents, can add more as needed...
        void SetConstituentPtMin(float pt);
        void SetConstituentPtMax(float pt);
        void SetConstituentEtaMin(float eta);
        void SetConstituentEtaMax(float eta);
        void SetConstituentAbsEtaMax(float eta);

        void ApplyConstituentKinematicCuts(); //Apply all constituent kinematic cuts that are called

        //Kinematic cuts for jets, can add more as needed...
        void SetJetPtMin(float pt);
        void SetJetPtMax(float pt);
        void SetJetEtaMin(float eta);
        void SetJetEtaMax(float eta);
        void SetJetAbsEtaMax(float eta);

        void ApplyJetKinematicCuts(); //Apply all the jet kinematic cuts that are called

        bool IsConstituentGhost(PseudoJet *constit);

        void EmptyAllVectors();
    
    private:
        JetDefinition *jet_def = nullptr; 
        //JetDefinition jet_def;
        bool doJetAreas = false;
        AreaDefinition *area_def = nullptr;
        GhostedAreaSpec *area_spec = nullptr;

        ClusterSequenceArea *CS_Area = nullptr;
	    ClusterSequence *CS = nullptr;

        Selector no_ghost = !SelectorIsPureGhost();
        vector<Selector> const_kin_cuts;
        vector<Selector> jet_kin_cuts;

        vector<PseudoJet> constituents_before_cuts;
        vector<PseudoJet> jet_constituents;
        vector<PseudoJet> all_jets;

    protected:
        vector<PseudoJet> jets; //Final jets of the events are contained here

};

class UserInfo: public PseudoJet::UserInfoBase{//User Info class to include non-kinematic details of the particles into the pseudojet objects...
public:
	UserInfo(const int indx, const short & ch):_charge(ch), _index(indx){}
    short getIndex() const {return _index;}
	short getCharge() const {return _charge;}
protected:
    int _index;
	short _charge;
};

#endif