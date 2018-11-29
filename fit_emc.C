#include "rootheader.h"
double myfun(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
	double lambda = par[4];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
	double exp_c1 = (1-alpha)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
	double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);	

    return c*(gaus + exp_c1*exp1);
}

double myfun2(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
	double lambda = par[4];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);

    return c*(gaus);
}

double myfun3(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
	double lambda = par[4];
    double pi = 3.141592654;

	double exp_c1 = (1-alpha)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
	double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);	
    return c*(exp_c1*exp1);
}


double myfun4(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
    double beta = par[4];
    double lambda = par[5];
    double gamma = par[6];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);
    double constant = (gamma)/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));

    return c*(exp_c1*exp1)/(alpha+beta+gamma);
}


bool get_data(int argc,char **argv,double *pars,double *epars,double& chi2,double &ndf,double& idealmean,double& idealsigma,string dirname,int index)
{

	string source = argv[1];
	string R = argv[2];	
	string Z = argv[3];	

	TCanvas *c1 = new TCanvas();
    	gStyle->SetOptFit(1);
    	gStyle->SetStatX(0.4);
    	gStyle->SetStatY(0.9);
    	gStyle->SetStatW(0.15);
    	gStyle->SetStatH(0.15);
    	TFile *file;
    	TF1 *f = new TF1("f",myfun,0,20000,5);

	TChain *t = new TChain("evt");

	t->Add("new.root");
	
	int total_entries = t->GetEntries();	

	int max_entries;
	if(source=="Ge68" || source == "K40")
		max_entries = 400000;
	else
		max_entries = 40000;

	if(source=="Cs137" || source== "Mn54")
		max_entries = 5235226;

	if(max_entries*(index+1) > total_entries) return false;

    	float edep;
    	if(source == "Ge68") edep = 1.0219978;
    	if(source == "Cs137") edep = 0.661657;
    	if(source == "Mn54") edep = 0.834848;
    	if(source == "Co60") edep = 2.5057385;
    	if(source == "K40") edep = 1.4608;
	
	t->Draw("totalPE>>h1(200,0,0)",Form("totalPE>0 && TMath::Abs(edep-%f)>0.001",edep),"",max_entries,max_entries*index);
	TH1F *h1 = (TH1F*)gDirectory->Get("h1");	
	
	t->Draw("totalPE>>h(200,0,0)",Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries,max_entries*index);
	TH1F *h = (TH1F*)gDirectory->Get("h");	
	
	h->Fit("gaus");

    	TF1* fun = h->GetFunction("gaus");
    	double C = fun->GetParameter(0);
    	double mean = fun->GetParameter(1);
    	double emean = fun->GetParError(2);
    	double sigma = fun->GetParameter(2);
    	idealmean = mean;
    	idealsigma = sigma;

	cout << max_entries*(index+1) << "\t" << max_entries*index << "\n";
	t->Draw("totalPE>>h(200,0,0)","totalPE>0","",max_entries,max_entries*index);
	h = (TH1F*)gDirectory->Get("h");

	int maxbin = h->GetMaximumBin();
	double maxbincenter = h->GetBinCenter(maxbin);

	t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),"","",max_entries,max_entries*index);	
	h = (TH1F*)gDirectory->Get("h");
	h->Fit("gaus");

     fun = h->GetFunction("gaus");

     C = fun->GetParameter(0);
     mean = fun->GetParameter(1);
     emean = fun->GetParError(2);
     sigma = fun->GetParameter(2);
    

	pars[0] = C*sqrt(2*3.14159)*sigma;
	pars[1] = 0.9;
	pars[2] = mean;
	pars[3] = sigma;
	pars[4] = 0.1;
	f->SetParLimits(4,0.0001,0.001);
	f->SetParameters(pars);

	f->SetParNames("C","#alpha","#mu","#sigma","#lambda","#gamma");
	t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-10*sigma),int(mean-10*sigma)+int(5*sigma)*3),"","",max_entries,max_entries*index);	
	cout << "data range is " << int(5*sigma) << "\t" << int(mean-10*sigma) << "\t" << int(mean-10*sigma)+int(5*sigma)*3 << "\n";
	h = (TH1F*)gDirectory->Get("h");
	
	h->Fit(f,"M","");
	f->GetParameters(pars);

	if(source=="Ge68"){
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
	}
	else{
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"MQ","",int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","", int(mean-3*sigma),int(mean-5*sigma)+int(5*sigma)*2);
	}

	f->GetParameters(pars);
	double *tmps = f->GetParErrors();
	memcpy(epars,tmps,5*sizeof(double));

	ndf = f->GetNDF();
	chi2 = f->GetChisquare();

    	TF1 *f2 = new TF1("f",myfun2,0,20000,5);
    	TF1 *f3 = new TF1("f",myfun3,0,20000,5);
	f2->SetParameters(pars);
	f3->SetParameters(pars);

	f2->SetLineColor(kBlue);
	f3->SetLineColor(kBlack);
	f2->Draw("same");
	f3->Draw("same");

	h1->SetLineColor(kGreen);
	h1->SetLineWidth(2);
	h1->Draw("same");

	string filename = dirname+"result_emc.C";
	c1->SaveAs(&filename[0]);
	filename = dirname+"result_emc.png";
	c1->SaveAs(&filename[0]);
	
	for(int i=0;i<4;i++){
		cout << pars[i] <<"\t" << epars[i] <<"\t";
	}

	return true;

}
int main(int argc,char **argv){

	double *pars = new double[5];
	double *epars = new double[5];
	double chi2;
	double ndf;
	string source = argv[1];
	string R = argv[2];	
	string Z = argv[3];	

	string dirname = "./";
	string filename = dirname+"result_emc";
	ofstream fout(&filename[0]);
	int index = 0;


	ofstream fout2("test_data");
	double idealmean,idealsigma;
	while(get_data(argc,argv,pars,epars,chi2,ndf,idealmean,idealsigma,dirname,index) && index <1){
		
		fout<< chi2 <<"\t"<< ndf <<"\t";

		for(int i=0;i<5;i++){
			fout << pars[i] <<"\t" << epars[i] <<"\t";
			fout2 << pars[i] <<"\t" << epars[i] <<"\t";
		}

		fout << "\n";
		fout << "bias of mean is "<<100*(pars[2]/idealmean - 1) << "\n";
		fout << "bias of sigma is "<<100*(pars[3]/idealsigma - 1) << "\n";
		index++;
	}

	return 0;
}
