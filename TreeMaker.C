#define TreeMaker_cxx

TreeMaker::TreeMaker(string name){
    f = new TFile(name.c_str(), "RECREATE");
    t = new TTree("TTreeMaker", "TTreeMaker");
}

template <class T> void TreeMaker::AddBranch(string name, T variable){
    t->Branch(name.c_str(), &variable);
}

void TreeMaker::FillTree(){
    t->Fill();
}

void TreeMaker::WriteTreeToFile(){
    f->cd();
    t->Write();
    f->Write();
    f->Close();
}

