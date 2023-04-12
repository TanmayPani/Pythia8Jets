#ifndef TreeMaker_h
#define TreeMaker_h

class TreeMaker{
    public:
        TreeMaker(string name);
        virtual ~TreeMaker();

        template <class T> void AddBranch(string name, T variable);
        template <class T> void AddVectorBranch(string name, T variable);

        void FillTree();
        void WriteTreeToFile();

    private:
        TFile *f;
        TTree *t;
};

#endif