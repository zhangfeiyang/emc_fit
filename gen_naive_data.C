{
	TFile *f = new TFile("new.root","read");
		
	TTree *evt = (TTree*)f->Get("evt");
	
	TFile *f1 = new TFile("new1.root","recreate");	
	TNtuple *t = new TNtuple("evt","","totalPE:edep");
	int N = evt->GetEntries();
	for(int i=0;i<N;i++){
		evt->GetEntry(i);
		float edep = evt->GetLeaf("edep")->GetValue(0);
		float totalPE = evt->GetLeaf("totalPE")->GetValue(0);
		//totalPE = gRandom->Gaus(edep*1252.63,1252.63*edep*sqrt(0.6617/edep)*3.7327429678e-2);
		float scale = 145.321*(edep-0.6617) + 1252.63;
		totalPE = gRandom->Gaus(edep*scale,scale*edep*sqrt(0.6617/edep)*3.7327429678e-2);
		t->Fill(totalPE,edep);
	}

	t->Write();
	f1->Close();

}
