//Wrapper class to use PYTHIA8, some of it uses ROOT6 compiler setup, but can be made ROOT independent easily
//Will make it even more streamlined as I use PYTHIA more :D
//Load PYTHIA8 libraries
R__ADD_INCLUDE_PATH(/Users/tanmaypani/PYTHIA8/pythia8309/include)
R__LOAD_LIBRARY(/Users/tanmaypani/PYTHIA8/pythia8309/lib/libpythia8.dylib)

#ifndef PythiaMaker_h
#define PythiaMaker_h

#include "Pythia8/Pythia.h" // access to Pythia objects.

using namespace Pythia8;

enum JetInitPartons {kAll, kQuarks, kGluons};

string hardQCD[3] = {"all", "qq2qq", "gg2gg"};

enum Tune {kMonash, kDetroit};

class PythiaMaker{
    public:
        PythiaMaker();
        PythiaMaker(Tune tune, bool hadronize, int MaxEvents);
        virtual ~PythiaMaker();

        void SetWeakDecaysOff();
        inline void TurnOnHardProcesses(JetInitPartons p);
        inline void TurnOnHardProcesses(string process);
        void SetHadronization(bool h);
        void SetPtHatRange(int pthmin, int pthmax);

        void Initialize();

        int GetNEvents();
        bool NextEvent();
        double GetCrossSection();
        float GetEventWeight(int i);
        float GetPtHat();
        double GetWeightSum();
        void ShowStats();

        int GetNParticles();
        int GetNFinalStateParticles();
        Particle* GetParticle(int i);

    private:
        Pythia *pythia;
        Pythia8::Event *MCEvent;
        Pythia8::Info *GenInfo;

};

#endif