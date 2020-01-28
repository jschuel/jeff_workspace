#include<math.h>


void alpha_calibration(int first_run, int last_run) {


TChain *iiwi = new TChain("data");
TChain *honu = new TChain("data");
TChain *kohola = new TChain("data");
TChain *nene = new TChain("data");
TChain *tako = new TChain("data");
TChain *humu = new TChain("data");
TChain *palila = new TChain("data");
TChain *elepaio = new TChain("data");

 gStyle->SetTitleFontSize(0.10);
 gStyle->SetLabelSize(.10, "XY");
  
 cout << "Analyzing run " << first_run << " through " << last_run << " from all TPCs"<<endl;

  
  for(int iRun=first_run; iRun<=last_run; iRun++) {
    iiwi->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/iiwi/%d_iiwi.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "iiwi" << std::endl;
    
    honu->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/honu/%d_honu.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "honu" << std::endl;

    kohola->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/kohola/%d_kohola.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "kohola" << std::endl;

    nene->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/nene/%d_nene.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "nene" << std::endl;
    
    tako->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/tako/%d_tako.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "tako" << std::endl;

    humu->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/humu/%d_humu.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "humu" << std::endl;

    palila->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/palila/%d_palila.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "palila" << std::endl;

    
    elepaio->Add(Form("/Users/vahsengrouplaptop/data/post_phase2/elepaio/%d_elepaio.root",iRun));
    std::cout << "Added run " << iRun << " from module " << "elepaio" << std::endl;
  } 

  TCut alpha = "hitside==11 && theta<=92 && theta>=88 && phi <= 2 && phi >= -2";

  TCut palila_alpha = "hitside==111 && theta<=92 && theta>=88 && phi <= 2 && phi >= -2";

  TCanvas *c1 = new TCanvas("c1","Alpha Calibrations",800,600);
  gStyle->SetOptStat(0);
  c1->Divide(2,4);
  
  c1->cd(1);

  iiwi->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *iiwi_length = iiwi->GetV1();
  double *iiwi_energy = iiwi->GetV2();
  double *iiwi_uncalibrated_energy = iiwi->GetV3();
  TH1F *h_iiwi = new TH1F("h_iiwi", "After gain calibration", 100, 0, 2500);
  TH1F *h_iiwi_uncal = new TH1F("h_iiwi_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<iiwi->GetSelectedRows(); i++){
    h_iiwi->Fill(iiwi_energy[i]);
    h_iiwi_uncal->Fill(iiwi_uncalibrated_energy[i]);
  }
  THStack *hs_iiwi = new THStack("hs_iiwi","BWD 198");
  h_iiwi->SetLineWidth(3);
  h_iiwi->SetMaximum(100);
  h_iiwi_uncal->SetLineWidth(3);
  h_iiwi_uncal->SetLineColor(6);
  h_iiwi_uncal->SetLineStyle(7);
  h_iiwi_uncal->SetMaximum(100);
  hs_iiwi->Add(h_iiwi);
  hs_iiwi->Add(h_iiwi_uncal);
  hs_iiwi->Draw("nostack");

  c1->Update();
  
  TLine *l=new TLine (1430,0,1430,100);
  l->SetLineColor(2);
  l->SetLineStyle(7);
  l->Draw();

  
  c1->cd(2);

  honu->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *honu_length = honu->GetV1();
  double *honu_energy = honu->GetV2();
  double *honu_uncalibrated_energy = honu->GetV3();
  TH1F *h_honu = new TH1F("h_honu", "After gain calibration", 100, 0, 2500);
  TH1F *h_honu_uncal = new TH1F("h_honu_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<honu->GetSelectedRows(); i++){
    h_honu->Fill(honu_energy[i]);
    h_honu_uncal->Fill(honu_uncalibrated_energy[i]);
  }
  THStack *hs_honu = new THStack("hs_honu","BWD 270");
  h_honu->SetLineWidth(3);
  h_honu->SetMaximum(100);
  h_honu_uncal->SetLineWidth(3);
  h_honu_uncal->SetLineColor(6);
  h_honu_uncal->SetLineStyle(7);
  h_honu_uncal->SetMaximum(100);
  hs_honu->Add(h_honu);
  hs_honu->Add(h_honu_uncal);
  hs_honu->Draw("nostack");

  c1->Update();
  
  l->Draw();

  auto *leg = new TLegend(0.75, 0.75, 0.95, .95);
  leg->AddEntry(h_iiwi, "After gain calibration", "l");
  leg->AddEntry(h_iiwi_uncal, "Before gain calibration", "dl");
  leg->AddEntry(l, "Reference line for 1430 keV alpha", "dl");
  leg->Draw();

  c1->cd(3);
  
  kohola->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *kohola_length = kohola->GetV1();
  double *kohola_energy = kohola->GetV2();
  double *kohola_uncalibrated_energy = kohola->GetV3();
  TH1F *h_kohola = new TH1F("h_kohola", "After gain calibration", 100, 0, 2500);
  TH1F *h_kohola_uncal = new TH1F("h_kohola_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<kohola->GetSelectedRows(); i++){
    h_kohola->Fill(kohola_energy[i]);
    h_kohola_uncal->Fill(kohola_uncalibrated_energy[i]);
  }
  THStack *hs_kohola = new THStack("hs_kohola","BWD 18");
  h_kohola->SetLineWidth(3);
  h_kohola->SetMaximum(100);
  h_kohola_uncal->SetLineWidth(3);
  h_kohola_uncal->SetLineColor(6);
  h_kohola_uncal->SetLineStyle(7);
  h_kohola_uncal->SetMaximum(100);
  hs_kohola->Add(h_kohola);
  hs_kohola->Add(h_kohola_uncal);
  hs_kohola->Draw("nostack");

  c1->Update();

  l->Draw();

  

  c1->cd(4);

  nene->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *nene_length = nene->GetV1();
  double *nene_energy = nene->GetV2();
  double *nene_uncalibrated_energy = nene->GetV3();
  TH1F *h_nene = new TH1F("h_nene", "After gain calibration", 100, 0, 2500);
  TH1F *h_nene_uncal = new TH1F("h_nene_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<nene->GetSelectedRows(); i++){
    h_nene->Fill(nene_energy[i]);
    h_nene_uncal->Fill(nene_uncalibrated_energy[i]);
  }
  THStack *hs_nene = new THStack("hs_nene","BWD 90");
  h_nene->SetLineWidth(3);
  h_nene->SetMaximum(100);
  h_nene_uncal->SetLineWidth(3);
  h_nene_uncal->SetLineColor(6);
  h_nene_uncal->SetLineStyle(7);
  h_nene_uncal->SetMaximum(100);
  hs_nene->Add(h_nene);
  hs_nene->Add(h_nene_uncal);
  hs_nene->Draw("nostack");

  c1->Update();

  l->Draw();

  

  c1->cd(5);

  tako->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *tako_length = tako->GetV1();
  double *tako_energy = tako->GetV2();
  double *tako_uncalibrated_energy = tako->GetV3();
  TH1F *h_tako = new TH1F("h_tako", "After gain calibration", 100, 0, 2500);
  TH1F *h_tako_uncal = new TH1F("h_tako_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<tako->GetSelectedRows(); i++){
    h_tako->Fill(tako_energy[i]);
    h_tako_uncal->Fill(tako_uncalibrated_energy[i]);
  }
  THStack *hs_tako = new THStack("hs_tako","FWD 90");
  h_tako->SetLineWidth(3);
  h_tako->SetMaximum(100);
  h_tako_uncal->SetLineWidth(3);
  h_tako_uncal->SetLineColor(6);
  h_tako_uncal->SetLineStyle(7);
  h_tako_uncal->SetMaximum(100);
  hs_tako->Add(h_tako);
  hs_tako->Add(h_tako_uncal);
  hs_tako->Draw("nostack");

  c1->Update();

  l->Draw();

  c1->cd(6);

  humu->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *humu_length = humu->GetV1();
  double *humu_energy = humu->GetV2();
  double *humu_uncalibrated_energy = humu->GetV3();
  TH1F *h_humu = new TH1F("h_humu", "After gain calibration", 100, 0, 2500);
  TH1F *h_humu_uncal = new TH1F("h_humu_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<humu->GetSelectedRows(); i++){
    h_humu->Fill(humu_energy[i]);
    h_humu_uncal->Fill(humu_uncalibrated_energy[i]);
  }
  THStack *hs_humu = new THStack("hs_humu","FWD 22");
  h_humu->SetLineWidth(3);
  h_humu->SetMaximum(100);
  h_humu_uncal->SetLineWidth(3);
  h_humu_uncal->SetLineColor(6);
  h_humu_uncal->SetLineStyle(7);
  h_humu_uncal->SetMaximum(100);
  hs_humu->Add(h_humu);
  hs_humu->Add(h_humu_uncal);
  hs_humu->Draw("nostack");

  c1->Update();

  l->Draw();

  

  c1->cd(7);

  palila->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", palila_alpha, "goff"); //energies assuming gain of 2000
  double *palila_length = palila->GetV1();
  double *palila_energy = palila->GetV2();
  double *palila_uncalibrated_energy = palila->GetV3();
  TH1F *h_palila = new TH1F("h_palila", "After gain calibration", 100, 0, 2500);
  TH1F *h_palila_uncal = new TH1F("h_palila_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<palila->GetSelectedRows(); i++){
    h_palila->Fill(palila_energy[i]);
    h_palila_uncal->Fill(palila_uncalibrated_energy[i]);
  }
  THStack *hs_palila = new THStack("hs_palila","FWD 270");
  h_palila->SetLineWidth(3);
  h_palila->GetXaxis()->CenterTitle();
  h_palila->SetMaximum(100);
  h_palila_uncal->SetLineWidth(3);
  h_palila_uncal->SetLineColor(6);
  h_palila_uncal->SetLineStyle(7);
  h_palila_uncal->SetMaximum(100);
  hs_palila->Add(h_palila);
  hs_palila->Add(h_palila_uncal);
  hs_palila->Draw("nostack");

  c1->Update();

  l->Draw();

  c1->cd(8);

  elepaio->Draw("length:paper_calibrated_track_energy:uncalibrated_track_energy", alpha, "goff"); //energies assuming gain of 2000
  double *elepaio_length = elepaio->GetV1();
  double *elepaio_energy = elepaio->GetV2();
  double *elepaio_uncalibrated_energy = elepaio->GetV3();
  TH1F *h_elepaio = new TH1F("h_elepaio", "After gain calibration", 100, 0, 2500);
  TH1F *h_elepaio_uncal = new TH1F("h_elepaio_uncal", "Before gain calibration", 100, 0, 2500);
  for(int i=0; i<elepaio->GetSelectedRows(); i++){
    h_elepaio->Fill(elepaio_energy[i]);
    h_elepaio_uncal->Fill(elepaio_uncalibrated_energy[i]);
  }
  THStack *hs_elepaio = new THStack("hs_elepaio","FWD 202");
  h_elepaio->SetLineWidth(3);
  h_elepaio->SetMaximum(100);
  h_elepaio_uncal->SetLineWidth(3);
  h_elepaio_uncal->SetLineColor(6);
  h_elepaio_uncal->SetLineStyle(7);
  h_elepaio_uncal->SetMaximum(100);
  hs_elepaio->Add(h_elepaio);
  hs_elepaio->Add(h_elepaio_uncal);
  hs_elepaio->Draw("nostack");

  c1->Update();

  l->Draw();

  

}

 

