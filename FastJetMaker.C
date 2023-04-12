#define FastJetMaker_cxx

#include "FastJetMaker.h"

FastJetMaker::FastJetMaker(){
    jet_def = new JetDefinition(antikt_algorithm, 0.4, BIpt2_scheme, Best);
    cout<<jet_def->description()<<endl;
}

FastJetMaker::FastJetMaker(JetAlgorithm algo, float R, RecombinationScheme reco_sch){
    jet_def = new JetDefinition(algo, R, reco_sch, Best);
    cout<<jet_def->description()<<endl;
}

FastJetMaker::~FastJetMaker(){

}

void FastJetMaker::FastJetAreaMaker(AreaType area_type, float MaxRap){
    area_spec = new GhostedAreaSpec(MaxRap);
	area_def = new AreaDefinition(area_type, *area_spec);
    doJetAreas = true;
}

void FastJetMaker::InputForClustering(PseudoJet constit){
    constituents_before_cuts.push_back(move(constit));
}

// void FastJetMaker::InputForClustering(int i, short ch, float px, float py, float pz, float E){
//     PseudoJet constit(px, py, pz, E);
//     constit.set_user_info(new UserInfo(i, ch));
//     constituents_before_cuts.push_back(move(constit));
// }

void FastJetMaker::InputForClustering(int i, short ch, float pt, float eta, float phi, float M){
    float px = pt * cos(phi);
    float py = pt * sin(phi);
    float pz = pt * sinh(eta);
    float E = sqrt(px*px + py*py + pz*pz + M);
    PseudoJet constit(px, py, pz, E);
    constit.set_user_info(new UserInfo(i, ch));
    constituents_before_cuts.push_back(move(constit));
}

vector<PseudoJet>* FastJetMaker::GetClusteredJets(){
    ApplyConstituentKinematicCuts();
    if(!doJetAreas){
        CS = new ClusterSequence(jet_constituents, *jet_def);
        all_jets = CS->inclusive_jets();
    }else{
        CS_Area = new ClusterSequenceArea(jet_constituents, *jet_def, *area_def);
        all_jets = no_ghost(CS_Area->inclusive_jets());
    }
    ApplyJetKinematicCuts();
    return static_cast<vector<PseudoJet>*>(&jets);
}

vector<PseudoJet>* FastJetMaker::GetConstituents(PseudoJet *jet){
    jet_constituents.clear();
    jet_constituents =  jet->constituents();
    return static_cast<vector<PseudoJet>*>(&jet_constituents);
}

bool FastJetMaker::IsConstituentGhost(PseudoJet *constit){
    if(!doJetAreas)return false;
    else return !constit->has_user_info();
}

void FastJetMaker::SetConstituentPtMin(float cpt){
    const_kin_cuts.push_back(SelectorPtMin(cpt));
}

void FastJetMaker::SetConstituentPtMax(float cpt){
    const_kin_cuts.push_back(SelectorPtMax(cpt));
}

void FastJetMaker::SetConstituentEtaMin(float ceta){
    const_kin_cuts.push_back(SelectorEtaMin(ceta));
}

void FastJetMaker::SetConstituentEtaMax(float ceta){
    const_kin_cuts.push_back(SelectorEtaMax(ceta));
}

void FastJetMaker::SetConstituentAbsEtaMax(float ceta){
    const_kin_cuts.push_back(SelectorAbsEtaMax(ceta));
}

void FastJetMaker::ApplyConstituentKinematicCuts(){
    Selector Selector_All = const_kin_cuts[0];
    for(int i = 1; i<const_kin_cuts.size(); i++){
        Selector_All = Selector_All*const_kin_cuts[i];
    }
    jet_constituents = Selector_All(constituents_before_cuts); 
}

void FastJetMaker::SetJetPtMin(float jpt){
    jet_kin_cuts.push_back(SelectorPtMin(jpt));
}

void FastJetMaker::SetJetPtMax(float jpt){
    jet_kin_cuts.push_back(SelectorPtMax(jpt));
}

void FastJetMaker::SetJetEtaMin(float jeta){
    jet_kin_cuts.push_back(SelectorEtaMin(jeta));
}

void FastJetMaker::SetJetEtaMax(float jeta){
    jet_kin_cuts.push_back(SelectorEtaMax(jeta));
}

void FastJetMaker::SetJetAbsEtaMax(float jeta){
    jet_kin_cuts.push_back(SelectorAbsEtaMax(jeta));
}

void FastJetMaker::ApplyJetKinematicCuts(){
    Selector Selector_All = jet_kin_cuts[0];
    for(int i = 1; i<jet_kin_cuts.size(); i++){
        Selector_All = Selector_All*jet_kin_cuts[i];
    }
    jets = Selector_All(all_jets); 
}

void FastJetMaker::EmptyAllVectors(){
    constituents_before_cuts.clear();
    jet_constituents.clear();
    all_jets.clear();
    jets.clear();
}

