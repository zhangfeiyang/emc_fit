{
	TFile *f = new TFile("new.root","read");
		
	TTree *evt = (TTree*)f->Get("evt");
		
	TH1F *h0 = new TH1F("h0","",200,500,900);
	TH1F *h1 = new TH1F("h1","",200,500,900);
		
	int N = evt->GetEntries();
	for(int i=0;i<N;i++){
		evt->GetEntry(i);
		float edep = evt->GetLeaf("edep")->GetValue(0);
		float totalPE = evt->GetLeaf("totalPE")->GetValue(0);
		if(edep<0.661){
			h0->Fill(totalPE);
			h1->Fill(gRandom->Gaus(edep*1252.63,1252.63*edep*0.03762*sqrt(0.6617/edep)));
		}
	}
	h0->SetLineColor(kRed);
	h0->Draw();
	h1->Draw("same");
}
