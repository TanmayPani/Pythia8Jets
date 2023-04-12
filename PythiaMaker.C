#define PythiaMaker_cxx

#include "PythiaMaker.h"

PythiaMaker::PythiaMaker(){
    pythia = new Pythia();
    MCEvent = static_cast<Pythia8::Event*>(&(pythia->event));

    pythia->readString("Main:timesAllowErrors = 3");
    pythia->readString("Init:showChangedSettings = on");
    pythia->readString("Next:numberCount = 10000");
    pythia->readString("Stat:showPartonLevel = on");
    pythia->readString("Beams:idA = 2212");
    pythia->readString("Beams:idB = 2212");
    pythia->readString("Beams:eCM = 200");
    pythia->readString("HardQCD:all = on");

    GenInfo = const_cast<Pythia8::Info*>(&(pythia->info));
}

PythiaMaker::PythiaMaker(Tune tune, bool hadronize, int MaxEvents){
    pythia = new Pythia();
    MCEvent = static_cast<Pythia8::Event*>(&(pythia->event)); 
    
    if(tune == kDetroit) pythia->readFile("detroit.cmnd");

    //pythia->readString(Form("HardQCD:%s = on", hardQCD[jetps]));
    if(!hadronize)pythia->readString("HadronLevel:Hadronize = off");

    pythia->readString(Form("Main:numberOfEvents = %i", MaxEvents));

    GenInfo = const_cast<Pythia8::Info*>(&(pythia->info));
}

PythiaMaker::~PythiaMaker(){
    
}

void PythiaMaker::Initialize(){
    pythia->init();
}

void PythiaMaker::SetPtHatRange(int pthmin, int pthmax){
    pythia->readString(Form("PhaseSpace:pTHatMin = %i", pthmin)); // low pthat cut
	pythia->readString(Form("PhaseSpace:pTHatMax = %i" ,pthmax)); // high pthat cut
}

void PythiaMaker::SetWeakDecaysOff(){
    pythia->readString("111:mayDecay =  off");//pi0                                                                                                                                      
	pythia->readString("211:mayDecay =  off");//pi+                                                                                                                                      
  	pythia->readString("221:mayDecay =  off");//eta                                                                                                                                      
  	pythia->readString("321:mayDecay =  off");//K+                                                                                                                                       
  	pythia->readString("310:mayDecay =  off");//Kshort                                                                                                                                   
  	pythia->readString("130:mayDecay =  off");//Klong                                                                                                                                    
  	pythia->readString("3122:mayDecay = off");//Lambda0                                                                                                                                  
  	pythia->readString("3212:mayDecay = off");//Sigma0                                                                                                                                   
  	pythia->readString("3112:mayDecay = off");//Sigma-                                                                                                                                   
  	pythia->readString("3222:mayDecay = off");//Sigma+                                                                                                                                   
  	pythia->readString("3312:mayDecay = off");//Xi-                                                                                                                                      
  	pythia->readString("3322:mayDecay = off");//Xi0                                                                                                                                      
  	pythia->readString("3334:mayDecay = off");//Omega-
}

//A little redundancy here, but may be useful incase we have multiple hardprocesses for quark and gluon jets
//Right now, the below definitions are equivalent, as kQuark is the same as qq2qq and kGluon as gg2gg
void PythiaMaker::TurnOnHardProcesses(JetInitPartons p){
   pythia->readString(Form("HardQCD:%s = on", hardQCD[p].c_str())); 
}

void PythiaMaker::TurnOnHardProcesses(string process){
   pythia->readString(Form("HardQCD:%s = on", process.c_str())); 
}

void PythiaMaker::SetHadronization(bool h){
    if(h)pythia->readString("HadronLevel:Hadronize = on");
    if(!h)pythia->readString("HadronLevel:Hadronize = off");
}

int PythiaMaker::GetNEvents(){
    return pythia->mode("Main:numberOfEvents");
}

bool PythiaMaker::NextEvent(){
    bool nxt = pythia->next();
    return nxt;
}

float PythiaMaker::GetPtHat(){
    return GenInfo->pTHat();
}

float PythiaMaker::GetEventWeight(int i){
    return GenInfo->weight(i);
}

int PythiaMaker::GetNParticles(){
    return MCEvent->size();
}

int PythiaMaker::GetNFinalStateParticles(){
    return MCEvent->nFinal();
}

Particle* PythiaMaker::GetParticle(int i){
    return static_cast<Particle*>(&(MCEvent->at(i)));
}

double PythiaMaker::GetCrossSection(){
    return GenInfo->sigmaGen();
}

double PythiaMaker::GetWeightSum(){
    return GenInfo->weightSum();
}

void PythiaMaker::ShowStats(){
    pythia->stat();
}