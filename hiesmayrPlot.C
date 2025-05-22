// skeleton code copied from caradonnaPlot.C to start constructing functions to plot Hiesmayr2019, eqn(15) etc.

#include "GetEnhancement.C"
#include "TMath.h"

using namespace TMath;

Double_t cos2phi(Double_t* x, Double_t* par)
{
	return par[0]*TMath::Cos(2*x[0]*TMath::DegToRad())+par[1];
}

void hiesmayrPlot(Float_t theta=10, Float_t theta_prime=82, Float_t phi=10)
{
	//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

// Hiesmayr and Moskal 2019, eqn. 15
	
	Float_t k_i=1;
//	Float_t theta=10;
	theta=theta/180*Pi();
//	Float_t theta_prime=82;
	theta_prime=theta_prime/180*Pi();
//	Float_t phi=10;
	phi=phi/180*Pi();
	Float_t r_0=2.82e-15;
	
	TF1* cosThetaTildeA=new TF1("cosThetaTildeA","Cos([0])*Cos([1])+Cos([2]-x)*Sin([0])*Sin([1])",0,2*Pi());
	TF1* cosThetaTildeB=new TF1("cosThetaTildeB","Cos(Pi()-[0])*Cos(Pi()-[1])+Cos([2]+Pi()-x)*Sin(Pi()-[0])*Sin(Pi()-[1])",0,2*Pi());
	cosThetaTildeA->SetTitle("cos#Theta^{#tilde}_{A}");
	cosThetaTildeB->SetTitle("cos#Theta^{#tilde}_{B}");
	cosThetaTildeA->SetParameters(theta_prime,theta,phi,k_i,r_0);
	cosThetaTildeB->SetParameters(theta_prime,theta,phi,k_i,r_0);
	
	TF1* k_primeA=new TF1("k_primeA","1/(1-(cosThetaTildeA(x))+1/[3])",0,2*Pi());
	TF1* k_primeB=new TF1("k_primeB","1/(1-(cosThetaTildeB(x))+1/[3])",0,2*Pi());
	k_primeA->SetTitle("k'_{A}");
	k_primeB->SetTitle("k'_{B}");
	k_primeA->SetParameters(theta_prime,theta,phi,k_i,r_0);
	k_primeB->SetParameters(theta_prime,theta,phi,k_i,r_0);

	TF1* gammaA=new TF1("gammaA","[3]/(k_primeA(x))+(k_primeA(x))/[3]",0,2*Pi());
	TF1* gammaB=new TF1("gammaB","[3]/(k_primeB(x))+(k_primeB(x))/[3]",0,2*Pi());
	gammaA->SetTitle("#gamma_{A}");
	gammaB->SetTitle("#gamma_{B}");
	gammaA->SetParameters(theta_prime,theta,phi,k_i,r_0);
	gammaB->SetParameters(theta_prime,theta,phi,k_i,r_0);
	
	TF1* V_A=new TF1("V_A","(1-(cosThetaTildeA(x))^2)/((gammaA(x))+(cosThetaTildeA(x))^2-1)",0,2*Pi());
	TF1* V_B=new TF1("V_B","(1-(cosThetaTildeB(x))^2)/((gammaB(x))+(cosThetaTildeB(x))^2-1)",0,2*Pi());
	V_A->SetTitle("V_{A}");
	V_B->SetTitle("V_{B}");
	V_A->SetParameters(theta_prime,theta,phi,k_i,r_0);
	V_B->SetParameters(theta_prime,theta,phi,k_i,r_0);

	TF1* F_A=new TF1("F_A","((k_primeA(x))/[3])^2*((gammaA(x))+(cosThetaTildeA(x))^2-1)",0,2*Pi());
	TF1* F_B=new TF1("F_B","((k_primeB(x))/[3])^2*((gammaB(x))+(cosThetaTildeB(x))^2-1)",0,2*Pi());
	F_A->SetTitle("F_{A}");
	F_A->SetTitle("F_{A}");
	F_A->SetParameters(theta_prime,theta,phi,k_i,r_0);
	F_B->SetParameters(theta_prime,theta,phi,k_i,r_0);

	TF2* hies_f=new TF2("hies_f",
		"[4]^2*(F_A(x))*(F_B(y))*(1./4)*(1-(V_A(x))*(V_B(y))*Cos(2*((x-[2])-(y-[2]))))"
	,0,2*Pi(),0,2*Pi());
	hies_f->SetParameters(theta_prime,theta,phi,k_i,r_0);
	hies_f->SetParName(0,"theta");
	hies_f->SetParName(1,"theta_prime");
	hies_f->SetParName(2,"phi");
	hies_f->SetParName(3,"k_i");
	hies_f->SetParName(4,"r_0");

	TCanvas* hies_c=new TCanvas("hies_c","Hiesmayr Plots");
	hies_c->Divide(3,2);

	hies_c->cd(1);
	cosThetaTildeA->Draw();
	cosThetaTildeA->SetLineColor(1);
	cosThetaTildeB->Draw("same");

	hies_c->cd(2);
	k_primeA->Draw();
	k_primeA->SetLineColor(1);
	k_primeB->Draw("same");
	
	hies_c->cd(3);
	gammaA->Draw();
	gammaA->SetLineColor(1);
	gammaB->Draw("same");
	
	hies_c->cd(4);
	V_A->Draw();
	V_A->SetLineColor(1);
	V_B->Draw("same");

	hies_c->cd(5);
	F_A->Draw();
	F_A->SetLineColor(1);
	F_B->Draw("same");

	hies_c->cd(6);
	hies_f->Draw("colz");
/*
// Caradonna et al
	TF3* cara_f=new TF3("cara_f","1./(16*pow(-2+Cos(x),3)*pow(-2+Cos(y),3))*	(9+9*pow(Cos(y),2)-3*pow(Cos(y),3)-3*pow(Cos(x),2)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))+pow(Cos(x),3)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))-4*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(-9+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))+	Cos(x)*(-9-9*pow(Cos(y),2)+3*pow(Cos(y),3)+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(9-Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))))",0,Pi(),0,Pi(),0,2*Pi());
//	TF3* cara_f=new TF3("cara_f","1./(16*pow(-2+Cos(x),3)*pow(-2+Cos(y),3))*	(9+9*pow(Cos(y),2)-3*pow(Cos(y),3)-3*pow(Cos(x),2)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))+pow(Cos(x),3)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))-4*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(-9+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))+	Cos(x)*(-9-9*pow(Cos(y),2)+3*pow(Cos(y),3)+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(9-Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))))",0,Pi(),0,Pi(),-Pi(),Pi());
////////////////////////////////////

// Pryce & Ward
	TF2* Ka_f=new TF2("Ka_f","(pow((1-TMath::Cos(x)),3)+2)*(pow((1-TMath::Cos(y)),3)+2)/(pow(2-TMath::Cos(x),3)*pow(2-TMath::Cos(y),3))",0,TMath::Pi(),0,TMath::Pi());
	TF2* Kb_f=new TF2("Kb_f","pow(TMath::Sin(x),2) * pow(TMath::Sin(y),2) / (pow(2-TMath::Cos(x),2) * pow(2-TMath::Cos(y),2))",0,TMath::Pi(),0,TMath::Pi());
	TF2* enhVsTh12_PW_f = new TF2("enhVsTh12_PW_f","(Ka_f+Kb_f)/(Ka_f-Kb_f)",0,TMath::Pi(),0,TMath::Pi());
	enhVsTh12_PW_f->GetXaxis()->SetTitle("#theta_{1}");
	enhVsTh12_PW_f->GetYaxis()->SetTitle("#theta_{2}");
	enhVsTh12_PW_f->SetTitle("Enhancement vs. #theta_{1} and #theta_{2} (Pryce & Ward)");

	TF1* cos2Dphi_f=new TF1("cos2Dphi_f","TMath::Cos(2*x)",-TMath::Pi(),TMath::Pi());
	
// Pryce & Ward 3D
	TF3* pryce_f = new TF3("pryce_f","1./16*(((pow((1-TMath::Cos(x)),3)+2)*(pow((1-TMath::Cos(y)),3)+2)/(pow(2-TMath::Cos(x),3)*pow(2-TMath::Cos(y),3)))-(pow(TMath::Sin(x),2) * pow(TMath::Sin(y),2) / (pow(2-TMath::Cos(x),2) * pow(2-TMath::Cos(y),2)))*TMath::Cos(2*z))",0,TMath::Pi(),0,TMath::Pi(),-TMath::Pi(),TMath::Pi());
	pryce_f->GetXaxis()->SetTitle("#theta_{1}");
	pryce_f->GetYaxis()->SetTitle("#theta_{2}");
////////////////////////////////////

// Independent polarised gammas (taken from Bohm & Aharonov)
	TF1* gamma_f=new TF1("gamma_f","2-TMath::Cos(x*TMath::DegToRad())+1/(2-TMath::Cos(x*TMath::DegToRad()))",0,180);
	TF1* polScatProbPara_f = new TF1("polScatProbPara_f","(2*(gamma_f)^2-4*gamma_f*(TMath::Sin(x*TMath::DegToRad()))^2+(TMath::Sin(x*TMath::DegToRad()))^4)",0,180);
	TF1* polScatProbPerp_f = new TF1("polScatProbPerp_f","(2*(gamma_f)^2-4*gamma_f*(TMath::Sin(x*TMath::DegToRad()))^2+3*(TMath::Sin(x*TMath::DegToRad()))^4)",0,180);
	TF1* polarised_f = new TF1("polarised_f",
	"(2*(gamma_f)^2-4*gamma_f*(TMath::Sin(x*TMath::DegToRad()))^2+3*(TMath::Sin(x*TMath::DegToRad()))^4) / (2*(gamma_f)^2-4*gamma_f*(TMath::Sin(x*TMath::DegToRad()))^2+(TMath::Sin(x*TMath::DegToRad()))^4)",0,180);
// Bohm & Aharonov Entanglement
	TF1* entScatProbPara_f = new TF1("entScatProbPara_f","(2*gamma_f*(gamma_f-2*TMath::Sin(x*TMath::DegToRad())^2))",0,180);
	TF1* entScatProbPerp_f = new TF1("entScatProbPerp_f","((gamma_f-2*TMath::Sin(x*TMath::DegToRad())^2)^2+gamma_f^2)",0,180);
	TF1* entangled_BA_f = new TF1("entangled_BA_f","entScatProbPerp_f/entScatProbPara_f",0,180);
////////////////////////////////////

// Snyder, Pasternack and Hornbostel
	TF1* g=new TF1("g","2-TMath::Cos(x)+1/(2-TMath::Cos(x))",0,TMath::Pi());
	TF2* snyder_f=new TF2("snyder_h",
	"(g(x)*g(y) - g(x)*pow(TMath::Sin(y),2) - g(y)*pow(TMath::Sin(x),2) +  2*pow(TMath::Sin(x),2)*pow(TMath::Sin(y),2)) / (g(x)*g(y)  - g(x)*pow(TMath::Sin(y),2) - g(y)*pow(TMath::Sin(x),2))"
	,0,TMath::Pi(),0,TMath::Pi());
////////////////////////////////////

	TCanvas* c1=new TCanvas("c1","Independent");
	c1->Divide(2,2);
	c1->cd(1);
	gamma_f->SetTitle("gamma");
	gamma_f->DrawCopy();
	c1->cd(2);
	polScatProbPara_f->SetTitle("Parallel Scattering Probability");
	polScatProbPara_f->DrawCopy("");
	c1->cd(3);
	polScatProbPerp_f->SetTitle("Perpendicular Scattering Probability");
	polScatProbPerp_f->DrawCopy("");
	c1->cd(4);
	polarised_f->SetTitle("Enhancement");
	polarised_f->DrawCopy();
	
	TCanvas* c2=new TCanvas("c2","Entangled (Bohm and Aharonov)");
	c2->Divide(2,2);
	c2->cd(1);
	gamma_f->DrawCopy();
	c2->cd(2);
	entScatProbPara_f->SetTitle("Parallel Scattering Probability");
	entScatProbPara_f->DrawCopy();
	c2->cd(3);
	entScatProbPerp_f->SetTitle("Perpendicular Scattering Probability");
	entScatProbPerp_f->DrawCopy();
	c2->cd(4);
	entangled_BA_f->SetTitle("Enhancement");
	entangled_BA_f->DrawCopy();
	
	TCanvas* c3=new TCanvas("c3","Comparison");
	entangled_BA_f->SetLineColor(4);
//	entangled_BA_f->SetTitle("Entangled");
	entangled_BA_f->GetHistogram()->SetTitle("");
	entangled_BA_f->GetXaxis()->SetTitle("#theta");
	entangled_BA_f->GetYaxis()->SetTitle("Enhancement factor");
	entangled_BA_f->GetXaxis()->SetTitleSize(0.04);
	entangled_BA_f->GetYaxis()->SetTitleSize(0.04);
	entangled_BA_f->Draw();
//	polarised_f->SetTitle("Independent");
	polarised_f->DrawCopy("same");
	TLegend* leg3 = new TLegend(0.65,0.7,0.9,0.9);
	leg3->AddEntry(entangled_BA_f,"Entangled");
	leg3->AddEntry(polarised_f,"Independent");
	leg3->Draw();
//	c3->BuildLegend(0.75,0.75,0.9,0.9);
	
	Double_t enh_ent=(entangled_BA_f->Integral(thetaMin,thetaMax))/(thetaMax-thetaMin);
	Double_t enh_pol=(polarised_f->Integral(thetaMin,thetaMax))/(thetaMax-thetaMin);
	
	Double_t enh_ent_fold=(entangled_BA_f->Integral(180-thetaMax,180-thetaMin))/(thetaMax-thetaMin);
	enh_ent_fold+=enh_ent;
	enh_ent_fold=enh_ent_fold/2;
	
//	TF1* cos2phi_f=new TF1("cos2phi_f","[0]*cos(2*x*TMath::DegToRad())+[1]",-180,180);
	TF1* cos2phi_ent=new TF1(Form("cos(2#phi) distributions for %.1f#circ<#theta<%.1f#circ",thetaMin,thetaMax),cos2phi,-180,180,2);
	TF1* cos2phi_pol=new TF1("cos(2#phi) distributions",cos2phi,-180,180,2);
	Double_t amp;
	Double_t off;
	
	TCanvas* c4=new TCanvas("c4","");
	c4->SetTitle("cos(2#phi) distributions");
	c4->cd();
	amp=(1-enh_ent)/(1+enh_ent);
//	off=1-amp;
	off=1;
	cos2phi_ent->SetParameter(0,amp);
	cos2phi_ent->SetParameter(1,off);
	cos2phi_ent->SetLineColor(4);
	cos2phi_ent->GetXaxis()->SetTitle("#phi (degrees)");
	cos2phi_ent->Draw("");

	
	amp=(1-enh_pol)/(1+enh_pol);
//	off=1-amp;
	cos2phi_pol->SetParameter(0,amp);
	cos2phi_pol->SetParameter(1,off);
	cos2phi_pol->SetLineColor(2);
	cos2phi_pol->Draw("same");
	
//	c4->BuildLegend(0.75,0.75,0.9,0.9);
	TLegend* leg4 = new TLegend(0.6,0.75,0.9,0.9);
	leg4->AddEntry(cos2phi_ent,Form("Entangled (enh.=%f)",enh_ent));
	leg4->AddEntry(cos2phi_pol,Form("Independent (enh.=%f)",enh_pol));
	leg4->Draw();

// Integrate enhancement function over various ranges of theta
	cout << endl << "d_theta	min.	max.	Enh." << endl;
	Int_t n=17;
	Double_t dTheta[n];
	Double_t enh[n];
	Double_t thetaLow, thetaHi;
	dTheta[0]=0;
	enh[0]=entangled_BA_f->Eval(82);
	for(int i=1;i<n;i++) {
		thetaLow=82-5*i;
		thetaHi=82+5*i;
		dTheta[i]=thetaHi-thetaLow;
		enh[i]=entangled_BA_f->Integral(thetaLow,thetaHi)/dTheta[i];
		cout << dTheta[i] << "	" << thetaLow << "	" << thetaHi << "	" << enh[i] << endl;
	}
	TGraph* enhVsDtheta_gr=new TGraph(n,dTheta,enh);
	enhVsDtheta_gr->SetTitle("Enhancement vs. #Delta#theta (symetric about 82#circ);#Delta#theta;Enhancement");
	enhVsDtheta_gr->SetMarkerStyle(21);
	TCanvas* enhVsDtheta_c=new TCanvas("enhVsDtheta_c","Enhancement vs. #Delta#theta");
	enhVsDtheta_gr->Draw("ALP");
	
	TCanvas* indTheta_c=new TCanvas("indTheta_c","Enh. for independent theta 1 and 2");
//	indTheta_c->Divide(2,1);
	indTheta_c->cd(1);
//	enhVsTh12_PW_f->SetTitle("Enhancement (Pryce & Ward);theta1;theta2");
	enhVsTh12_PW_f->Draw("surf1");
	
	Float_t bestTheta = 81.660564;

	cout << "Theta min = " << thetaMin*TMath::DegToRad() << ", Enh = " << enhVsTh12_PW_f->Eval(thetaMin*TMath::DegToRad(),thetaMin*TMath::DegToRad()) << endl;
	cout << "Theta max = " << thetaMax*TMath::DegToRad() << ", Enh = " << enhVsTh12_PW_f->Eval(thetaMax*TMath::DegToRad(),thetaMin*TMath::DegToRad()) << endl;
	
	Double_t enh_Pryce=(enhVsTh12_PW_f->Integral(thetaMin*TMath::DegToRad(),thetaMax*TMath::DegToRad(),thetaMin*TMath::DegToRad(),thetaMax*TMath::DegToRad()))/pow(thetaMax*TMath::DegToRad()-thetaMin*TMath::DegToRad(),2);
	cout << Form("Enh. from Pryce & Ward (with %.2f<theta<%.2f): %f",thetaMin,thetaMax,enh_Pryce) << endl;
	
	TCanvas* caraTest_c=new TCanvas("caraTest_c","Caradonna Plots");
	caraTest_c->Divide(2,3);
	
	Float_t thetaTest=95;
	
	Int_t bins=180;
	Float_t binWidth=360./bins;
	
	TH1F* pryceTest_h=new TH1F("pryceTest_h",Form("Pryce Ward #Delta#phi plot at %.2f",bestTheta),bins,-180,180);
	TH1F* pryceTest2_h=new TH1F("pryceTest2_h",Form("Pryce Ward #Delta#phi plot at %.2f",thetaTest),bins,-180,180);
	TH1F* pryceTest3_h=new TH1F("pryceTest3_h","Pryce Ward enhancement",bins,0,180);
//	TH1F* caraTest_h=new TH1F("caraTest_h",Form("Caradonna #Delta#phi plot at %.2f",bestTheta),bins,0,360);
//	TH1F* caraTest2_h=new TH1F("caraTest2_h",Form("Caradonna #Delta#phi plot at %.2f",thetaTest),bins,0,360);
	TH1F* caraTest_h=new TH1F("caraTest_h",Form("Caradonna #Delta#phi plot at %.2f",bestTheta),bins,-180,180);
	TH1F* caraTest2_h=new TH1F("caraTest2_h",Form("Caradonna #Delta#phi plot at %.2f",thetaTest),bins,-180,180);
	TH1F* caraTest3_h=new TH1F("caraTest3_h","Caradonna enhancement",bins,0,180);
	
	TH1F* bohmTest3_h=new TH1F("bohmTest3_h","Bohm and Aharonov enhancement",bins,0,180);
	
	TH1F* caraTests_h[8];
	for(int i=0;i<8;i++) caraTests_h[i]=new TH1F(Form("caraTests_%i_h",i),Form("Caradonna plot at #theta_1=#theta_2=%.0f",i*5.0+75.0),bins,-180,180);
	
	TH3F* cara3D_h=new TH3F("cara3D_h","Caradonna 3D",18,0,180,18,0,180,18,0,360);
	cara3D_h->GetXaxis()->SetTitle("theta1");
	cara3D_h->GetYaxis()->SetTitle("theta2");
	cara3D_h->GetZaxis()->SetTitle("phi");
	TH3F* cara3D_cut_h=new TH3F("cara3D_cut_h","Caradonna 3D",18,75,115,18,75,115,18,0,360);
	cara3D_cut_h->GetXaxis()->SetTitle("theta1");
	cara3D_cut_h->GetYaxis()->SetTitle("theta2");
	cara3D_cut_h->GetZaxis()->SetTitle("phi");
	TH3F* cara3D_cut2_h=new TH3F("cara3D_cut2_h","Caradonna 3D",18,70,140,18,70,140,18,110,270);
	cara3D_cut2_h->GetXaxis()->SetTitle("theta1");
	cara3D_cut2_h->GetYaxis()->SetTitle("theta2");
	cara3D_cut2_h->GetZaxis()->SetTitle("phi");
	
	TH3F* pryce3D_h=new TH3F("pryce3D_h","Pryce and Ward 3D",18,0,180,18,0,180,18,0,360);
	pryce3D_h->GetXaxis()->SetTitle("theta1");
	pryce3D_h->GetYaxis()->SetTitle("theta2");
	pryce3D_h->GetZaxis()->SetTitle("phi");
	TH3F* pryce3D_cut_h=new TH3F("pryce3D_cut_h","Pryce and Ward 3D",18,75,115,18,75,115,18,0,360);
	pryce3D_cut_h->GetXaxis()->SetTitle("theta1");
	pryce3D_cut_h->GetYaxis()->SetTitle("theta2");
	pryce3D_cut_h->GetZaxis()->SetTitle("phi");
	TH3F* pryce3D_cut2_h=new TH3F("pryce3D_cut2_h","Pryce and Ward 3D",18,70,140,18,70,140,18,110,270);
	pryce3D_cut2_h->GetXaxis()->SetTitle("theta1");
	pryce3D_cut2_h->GetYaxis()->SetTitle("theta2");
	pryce3D_cut2_h->GetZaxis()->SetTitle("phi");
	
	Float_t ninety, zero;
	for(int i=0; i<bins; i++) {
	
//		cout << i << " " << i+1 << " " << i*binWidth-180+binWidth/2 << " " << pryce_f->Eval(bestTheta*TMath::DegToRad(),bestTheta*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()) << endl;
	
		pryceTest_h->SetBinContent(i+1,pryce_f->Eval(bestTheta*TMath::DegToRad(),bestTheta*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()));
		pryceTest2_h->SetBinContent(i+1,pryce_f->Eval(thetaTest*TMath::DegToRad(),thetaTest*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()));
		caraTest_h->SetBinContent(i+1,cara_f->Eval(bestTheta*TMath::DegToRad(),bestTheta*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()));
		caraTest2_h->SetBinContent(i+1,cara_f->Eval(thetaTest*TMath::DegToRad(),thetaTest*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()));
		
		for(int i2=0;i2<8;i2++) caraTests_h[i2]->SetBinContent(i+1,cara_f->Eval((i2*5.0+75.0)*TMath::DegToRad(),(i2*5.0+75.0)*TMath::DegToRad(),(i*binWidth-180+binWidth/2)*TMath::DegToRad()));
		
		ninety = pryce_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),90*TMath::DegToRad());
		zero = pryce_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),0*TMath::DegToRad());
		pryceTest3_h->SetBinContent(i+1,ninety/zero);
		
		ninety = cara_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),90*TMath::DegToRad());
		zero = cara_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),0*TMath::DegToRad());
		caraTest3_h->SetBinContent(i+1,ninety/zero);
		
//		ninety = entangled_BA_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),90*TMath::DegToRad());
//		zero = entangled_BA_f->Eval((i+1)*binWidth/2*TMath::DegToRad(),(i+1)*binWidth/2*TMath::DegToRad(),0*TMath::DegToRad());
		
//		bohmTest2_h->SetBinContent(i+1,
		bohmTest3_h->SetBinContent(i+1,entangled_BA_f->Eval((i+1)*binWidth/2,(i+1)*binWidth/2));
		
//		for(int j=0;j<180;j++) for(int k=0;k<180;k++) {
//			if(i<90 && j<90 && k>0) continue;
//			else {		 cara3D_h->SetBinContent(i+1,j+1,k+1,cara_f->Eval((i*binWidth/2+binWidth/4)*TMath::DegToRad(),(j*binWidth/2+binWidth/4)*TMath::DegToRad(),(k*2-180+binWidth/2)*TMath::DegToRad()));
//			}
//		}
		
	}
	caraTest_c->cd(1);
//	caraTest_h->SetLineColor(2);
//	caraTest_h->SetLineStyle(3);
	caraTest_h->DrawCopy("");
	caraTest_c->cd(2);
	pryceTest_h->DrawCopy("hist");
	caraTest_c->cd(3);
	caraTest2_h->DrawCopy("");
	caraTest_c->cd(4);
	pryceTest2_h->DrawCopy("");
	caraTest_c->cd(5);
	caraTest3_h->DrawCopy("");
	caraTest_c->cd(6);
	pryceTest3_h->DrawCopy("");
	
	bins=18;
	binWidth=180/bins;
	
	Double_t th1, th2, ph;
	
	for(int i=0;i<bins;i++) for(int j=0;j<bins;j++) for(int k=0;k<bins;k++) {
		th1=(i*binWidth+binWidth/2);
		th2=(j*binWidth+binWidth/2);
		ph=(k*binWidth*2+binWidth);
//		if(th1>90 && th2>90) continue;
//		else {
		 cara3D_h->Fill(th1, th2, ph, cara_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
		 pryce3D_h->Fill(th1, th2, ph, pryce_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
//		}
	}
	
	Float_t binWidth2=(115.-75.)/bins;
	for(int i=0;i<bins;i++) for(int j=0;j<bins;j++) for(int k=0;k<bins;k++) {
		th1=(i*binWidth2+binWidth2/2+75);
		th2=(j*binWidth2+binWidth2/2+75);
		ph=(k*binWidth*2+binWidth);
		cara3D_cut_h->Fill(th1, th2, ph, 1/cara_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
		pryce3D_cut_h->Fill(th1, th2, ph, 1/pryce_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
	}
	
	binWidth2=(140.-70.)/bins;
	Float_t binWidth3=(80.-(-80.))/bins;
	for(int i=0;i<bins;i++) for(int j=0;j<bins;j++) for(int k=0;k<bins;k++) {
		th1=(i*binWidth2+binWidth2/2+75);
		th2=(j*binWidth2+binWidth2/2+75);
		ph=(k*binWidth3+110+binWidth3/2);
		cara3D_cut2_h->Fill(th1, th2, ph, cara_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
		pryce3D_cut2_h->Fill(th1, th2, ph, pryce_f->Eval( th1*TMath::DegToRad(), th2*TMath::DegToRad(), ph*TMath::DegToRad()));
	}
	
	TCanvas* cara3D_c=new TCanvas("cara3D_c","Caradonna 3D");
	cara3D_h->Draw("box2 z");
	TCanvas* cara3D_cut_c=new TCanvas("cara3D_cut_c","Caradonna 3D (cut)");
	cara3D_cut_h->Draw("box2 z");
//	cara3D_cut_h->Draw("iso fb");
	TCanvas* cara3D_cut2_c=new TCanvas("cara3D_cut2_c","Caradonna 3D (cut 2)");
	cara3D_cut2_h->Draw("iso FB");
	gPad->SetTheta(20);
	gPad->SetPhi(-135);
	
	TCanvas* pryce3D_c=new TCanvas("pryce3D_c","Pryce and Ward 3D");
	pryce3D_h->Draw("box2 z");
	TCanvas* pryce3D_cut_c=new TCanvas("pryce3D_cut_c","Pryce and Ward 3D (cut)");
	pryce3D_cut_h->Draw("box2 z");
//	pryce3D_cut_h->Draw("iso FB");
	TCanvas* pryce3D_cut2_c=new TCanvas("pryce3D_cut2_c","Pryce and Ward 3D (cut 2)");
	pryce3D_cut2_h->Draw("iso FB");
	gPad->SetTheta(20);
	gPad->SetPhi(-135);
	
//	cara3D_cut2_h->Draw("box2 z");
	//	cara3D_h->Draw("iso");
//	cara_f->SetClippingBoxOn(90,90,0);
//	cara_f->SetFillColor(30);
//	cara_f->SetLineColor(15);
//	cara_f->Draw("box");
	
	TCanvas* caraTest2_c=new TCanvas("caraTest2_c","Caradonna Plots");
	caraTest_h->SetStats(0);
	caraTest_h->Draw();
	for(int i=0;i<8;i++) {
		caraTests_h[i]->SetLineColor(i+1);
		caraTests_h[i]->SetStats(0);
		caraTests_h[i]->Draw("same");
	}
	caraTest2_c->BuildLegend(0.55,0.55,0.9,0.9);
	
	
	
	cout << "//////////////////////////////////////////////////////////" << endl;
	cout << "Enhancement at theta_1=theta_2=" << bestTheta << " using Caradonna (3D): "
		<< cara_f->Eval(bestTheta*TMath::DegToRad() , bestTheta*TMath::DegToRad() , 90.*TMath::DegToRad())/pryce_f->Eval(bestTheta*TMath::DegToRad() , bestTheta*TMath::DegToRad(),0.*TMath::DegToRad())
		<< endl;

	cout << "Enhancement at theta_1=theta_2=" << bestTheta << " using Pryce & Ward (3D): "
		<< pryce_f->Eval(bestTheta*TMath::DegToRad() , bestTheta*TMath::DegToRad() , 90.*TMath::DegToRad())/pryce_f->Eval(bestTheta*TMath::DegToRad() , bestTheta*TMath::DegToRad(),0.*TMath::DegToRad())
		<< endl;

	cout << "Enhancement at theta_1=theta_2=" << bestTheta << " using Pryce & Ward (2D): "
		<< enhVsTh12_PW_f->Eval(bestTheta/180*TMath::Pi(),bestTheta/180*TMath::Pi())
		<< endl;
	cout << "//////////////////////////////////////////////////////////" << endl;
	
	TCanvas* caraTest3_c = new TCanvas("caraTest3_c","Caradonna, Bohm & Aharonov Comparison");
	caraTest3_c->Divide(2,3);
	caraTest3_c->cd(1);
//	caraTest_h->SetLineColor(2);
//	caraTest_h->SetLineStyle(3);
	caraTest_h->DrawCopy("");
	caraTest3_c->cd(2);
//	pryceTest_h->DrawCopy("hist");
	caraTest3_c->cd(3);
	caraTest2_h->DrawCopy("");
	caraTest3_c->cd(4);
//	pryceTest2_h->DrawCopy("");
	caraTest3_c->cd(5);
	caraTest3_h->DrawCopy("");
	caraTest3_c->cd(6);
	bohmTest3_h->DrawCopy("");
	
	TCanvas* enh_comp_c = new TCanvas("enh_comp_c","Enhancement Formalism Comparison");
	enh_comp_c->Divide(2,2);
	enh_comp_c->cd(1);
	entangled_BA_f->SetTitle("Bohm and Aharonov enhancement plot");
	entangled_BA_f->DrawCopy();
	enh_comp_c->cd(2);
	pryceTest3_h->DrawCopy();
	enh_comp_c->cd(3);
	caraTest3_h->DrawCopy();
	enh_comp_c->cd(4);
//	entangled_BA_f->SetTitle("Comparison");
	entangled_BA_f->DrawCopy();
	pryceTest3_h->SetLineColor(2);
	pryceTest3_h->SetLineStyle(2);
	pryceTest3_h->DrawCopy("same");
	caraTest3_h->SetLineColor(3);
	caraTest3_h->SetLineStyle(3);
	caraTest3_h->DrawCopy("same");
	
	TCanvas* enh_comp2_c = new TCanvas("enh_comp2_c","Enhancement Formalism Comparison");
//	entangled_BA_f->DrawCopy();
//	bohmTest3_h->GetXaxis()->SetTitle("#theta");
//	bohmTest3_h->GetYaxis()->SetTitle("Enhancement factor");
//	bohmTest3_h->DrawCopy("C");
	
	pryceTest3_h->GetXaxis()->SetTitle("#theta");
	pryceTest3_h->GetYaxis()->SetTitle("Enhancement factor");
	pryceTest3_h->SetLineColor(4);
	pryceTest3_h->SetLineStyle(1);
	pryceTest3_h->DrawCopy("C");
	caraTest3_h->SetLineColor(2);
	caraTest3_h->SetMarkerStyle(5);
	caraTest3_h->SetMarkerColor(2);
	caraTest3_h->SetMarkerSize(0.5);
	caraTest3_h->DrawCopy("P same");
	enh_comp2_c->BuildLegend(0.6,0.75,0.9,0.9);
*/
}
