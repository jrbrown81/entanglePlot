#include "GetEnhancement.C"
#include "TMath.h"

using namespace TMath;

Double_t cos2phi(Double_t* x, Double_t* par)
{
	return par[0]*TMath::Cos(2*x[0]*TMath::DegToRad())+par[1];
}

//Double_t g(Double_t* x, Double_t* par)
//{
//	return 2-TMath::Cos(x*TMath::DegToRad())+1/(2-TMath::Cos(x*TMath::DegToRad()));
//}

void entanglePlot(Double_t thetaMin=67, Double_t thetaMax=97)
{

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

//	gROOT->ProcessLine(".L ~/work/root_macros/root_macros/GetEnhancement.C");

// Caradonna et al
	TF3* cara_f=new TF3("cara_f","1./(16*pow(-2+Cos(x),3)*pow(-2+Cos(y),3))*	(9+9*pow(Cos(y),2)-3*pow(Cos(y),3)-3*pow(Cos(x),2)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))+pow(Cos(x),3)*(-3+3*Cos(y)-3*pow(Cos(y),2)+pow(Cos(y),3))-4*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(-9+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))+	Cos(x)*(-9-9*pow(Cos(y),2)+3*pow(Cos(y),3)+2*Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2)+Cos(y)*(9-Cos(2*z)*pow(Sin(x),2)*pow(Sin(y),2))))",0,Pi(),0,Pi(),-Pi(),Pi());
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

	TFile *file = new TFile("entangledCZT_hist.root");
	TH1F* cztEnt_h=(TH1F*)file->Get("h4");
	cztEnt_h->SetTitle("CZT Block (entangled), 67#circ<#theta<97#circ");
//	cztEnt_h->SetTitle("QE-G4 Simulation");
	
	TFile *file2 = new TFile("unentangledCZT_hist.root");
	TH1F* cztUnent_h=(TH1F*)file2->Get("h4");
	cztUnent_h->SetTitle("CZT Block (unentangled), 67#circ<#theta<97#circ");
//	cztUnent_h->SetTitle("G4 Simulation");
	
	TCanvas* c5=new TCanvas("c5","CZT Block Simulation");
//	c5->Divide(2,2);
//	c5->cd(1);
//	cztEnt_h->Draw();
	Double_t norm=cztEnt_h->Integral(-180,180);
	Int_t nBins=cztEnt_h->GetNbinsX();
	norm=norm/nBins;
//	cout << norm << " " << nBins << endl;
//	c5->cd(2);
	TH1F* cztEntNorm_h=(TH1F*)cztEnt_h->Clone("cztEntNorm_h");
	cztEntNorm_h->Scale(1/norm);
	Double_t enh_cztEnt=GetEnhancement(cztEntNorm_h,"IL0Q");
//	cztEntNorm_h->SetLineColor(1);
	cztEntNorm_h->SetStats(0);
		
	norm=cztUnent_h->Integral(-180,180);
	nBins=cztUnent_h->GetNbinsX();
	norm=norm/nBins;
//	cout << norm << " " << nBins << endl;
//	c5->cd(2);
	TH1F* cztUnentNorm_h=(TH1F*)cztUnent_h->Clone("cztUnentNorm_h");
	cztUnentNorm_h->Scale(1/norm);
	Double_t enh_cztUnent=GetEnhancement(cztUnentNorm_h,"IL0Q");
	cztUnentNorm_h->SetLineColor(2);
	cztUnentNorm_h->SetMarkerColor(2);
	cztUnentNorm_h->SetStats(0);
	
	c5->cd();
	cztEntNorm_h->SetTitle("");
	cztEntNorm_h->SetMarkerStyle(47);
	cztEntNorm_h->SetMarkerColor(4);
	cztEntNorm_h->SetLineColor(4);
	cztEntNorm_h->SetFillColor(4);
	
	cztEntNorm_h->GetXaxis()->SetTitle("#Delta#phi (degrees)");
	cztEntNorm_h->GetYaxis()->SetTitle("Normalised Scattering Probability");
	
	cztEntNorm_h->GetXaxis()->SetLabelSize(0.05);
	cztEntNorm_h->GetXaxis()->SetTitleSize(0.05);
	cztEntNorm_h->GetXaxis()->SetTitleOffset(0.9);
	cztEntNorm_h->GetYaxis()->SetLabelSize(0.05);
	cztEntNorm_h->GetYaxis()->SetTitleSize(0.05);
	cztEntNorm_h->GetYaxis()->SetTitleOffset(0.9);
	
	cztEntNorm_h->Draw("e1");
	cztUnentNorm_h->SetMarkerStyle(23);
	cztUnentNorm_h->SetMarkerColor(2);
	cztUnentNorm_h->SetLineColor(2);
	cztUnentNorm_h->SetFillColor(2);
	cztUnentNorm_h->Draw("e1 same");
	cos2phi_ent->Draw("same");
	cos2phi_pol->Draw("same");
//	TLegend* leg5 = new TLegend(0.6,0.75,0.9,0.9);
//	TLegend* leg5 = new TLegend(0.56,0.1,0.84,0.25);
	TLegend* leg5 = new TLegend(0.4,0.87,0.99,0.99);
//	leg5->AddEntry(cztEntNorm_h,cztEnt_h->GetTitle());
//	leg5->AddEntry(cztEntNorm_h,Form("QE-G4 Simulation (enh.=%.3f)",enh_cztEnt));
//	leg5->AddEntry(cztUnentNorm_h,cztUnent_h->GetTitle());
//	leg5->AddEntry(cztUnentNorm_h,Form("G4 Simulation (enh.=%.3f)",enh_cztUnent));
//	leg5->AddEntry(cos2phi_ent,Form("Theoretical (Entangled) (enh.=%.3f)",enh_ent));
//	leg5->AddEntry(cos2phi_pol,Form("Theoretical (Independent) (enh.=%.3f)",enh_pol));
	leg5->SetNColumns(2);
	leg5->AddEntry(cztEntNorm_h,"QE-Geant4 Simulation","ep");
	leg5->AddEntry(cztUnentNorm_h,"Geant4 Simulation","ep");
	leg5->AddEntry(cos2phi_ent,"Theoretical (Entangled)","l");
	leg5->AddEntry(cos2phi_pol,"Theoretical (Independent)","l");
	leg5->Draw();
	
//	c5->cd(4);
//	cos2phi_ent->Draw("");
//	cos2phi_pol->Draw("same");

	TCanvas* c6=new TCanvas("c6","CZT Block Simulation, no errors");
	cztEntNorm_h->Draw("hist p");
	cztUnentNorm_h->Draw("hist p same");
	cos2phi_ent->Draw("same");
	cos2phi_pol->Draw("same");

	TF1* enhFunc=new TF1("enhFunc","(1-[0])/(1+[0])*TMath::Cos(2*TMath::DegToRad()*x)+1",-180,180);
	Double_t c_t[9];
	enhFunc->SetParName(0,"enh");
	enhFunc->SetParameter(0,enh_cztEnt);
	cout << "Enh. from fit to CZT block simulation: " << enh_cztEnt << endl;
	cout << "#Delta#phi	c_t" << endl;
	for(int i=0;i<9;i++) {
		c_t[i]=(enhFunc->Integral(i*20,(i+1)*20))/20;
		cout << i*20+10 << "	" << c_t[i] << endl;
	}
	
	cout << Form("Enh. from theory (with %.2f<theta<%.2f): %f",thetaMin,thetaMax,enh_ent) << endl;
	cout << Form("Enh. from theory (with %.2f<theta<%.2f and %.2f<theta<%.2f): %f",thetaMin,thetaMax,180-thetaMax, 180-thetaMin,enh_ent_fold) << endl;
	cout << "#Delta#phi	c_t" << endl;
	for(int i=0;i<9;i++) {
		c_t[i]=(cos2phi_ent->Integral(i*20,(i+1)*20))/20;
		cout << i*20+10 << "	" << c_t[i] << endl;
	}
	
	TCanvas* c7=new TCanvas("c7","CZT Block Simulation, error bands");
	c7->cd();
	
	TH1F* cztUnentNorm2_h=(TH1F*)cztUnent_h->Clone("cztUnentNorm2_h");
	cztUnentNorm2_h->Rebin(4);
	
//	cztEntNorm_h->DrawCopy("e4");
	cztUnentNorm2_h->SetFillColor(4);
	cztUnentNorm2_h->DrawCopy("e4");
//	cos2phi_ent->DrawCopy("same");
//	cos2phi_pol->DrawCopy("same");
	TLegend* leg7 = new TLegend(0.6,0.8,0.95,0.95);
//	leg7->AddEntry(cztEntNorm_h,cztEnt_h->GetTitle());
//	leg7->AddEntry(cztEntNorm_h,Form("QE-G4 Simulation (enh.=%.3f)",enh_cztEnt));
//	leg7->AddEntry(cztUnentNorm_h,cztUnent_h->GetTitle());
//	leg7->AddEntry(cztUnentNorm_h,Form("G4 Simulation (enh.=%.3f)",enh_cztUnent));
//	leg7->AddEntry(cos2phi_ent,Form("Theoretical (Entangled) (enh.=%.3f)",enh_ent));
//	leg7->AddEntry(cos2phi_pol,Form("Theoretical (Independent) (enh.=%.3f)",enh_pol));
	leg7->AddEntry(cztEntNorm_h,"QE-G4 Simulation","ep");
	leg7->AddEntry(cztUnentNorm_h,"G4 Simulation","ep");
	leg7->AddEntry(cos2phi_ent,"Theoretical (Entangled)","l");
	leg7->AddEntry(cos2phi_pol,"Theoretical (Independent)","l");
//	leg7->Draw();
	
	
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

	TCanvas* CZT_block_comp_c=new TCanvas("CZT_block_comp_c","CZT Block Simulation (Pryce & Ward)");
	CZT_block_comp_c->cd();
	cztEntNorm_h->Draw("e1");
	cztUnentNorm_h->Draw("e1 same");

	TF1* cos2phi_pryce=new TF1(Form("cos(2#phi) distributions for %.1f#circ<#theta<%.1f#circ",thetaMin,thetaMax),cos2phi,-180,180,2);
	amp=(1-enh_Pryce)/(1+enh_Pryce);
	off=1;
	cos2phi_pryce->SetParameter(0,amp);
	cos2phi_pryce->SetParameter(1,off);
	cos2phi_pryce->SetLineColor(4);
	cos2phi_pryce->GetXaxis()->SetTitle("#phi (degrees)");
	cos2phi_pryce->Draw("same");
	
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
	
	TH3F* cara3D_h=new TH3F("cara3D_h","Caradonna 3D",180,0,180,180,0,180,180,-180,180);
	
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
		
		for(int j=0;j<180;j++) for(int k=0;k<180;k++) cara3D_h->SetBinContent(i+1,j+1,k+1,cara_f->Eval((i*binWidth/2+binWidth/4)*TMath::DegToRad(),(j*binWidth/2+binWidth/4)*TMath::DegToRad(),(k*2-180+binWidth/2)*TMath::DegToRad()));
		
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
	
	TCanvas* cara3D_c=new TCanvas("cara3D","Caradonna 3D");
	cara3D_h->Draw("col");
	
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

}
