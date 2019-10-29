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

  iiwi->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *iiwi_length = iiwi->GetV1();
  double *iiwi_energy = iiwi->GetV2();
  TH1F *h_iiwi = new TH1F("h_iiwi", "Iiwi", 100, 0, 2500);
  for(int i=0; i<iiwi->GetSelectedRows(); i++){
    h_iiwi->Fill(iiwi_energy[i]);
  }
  h_iiwi->SetLineWidth(3);
  h_iiwi->GetXaxis()->SetTitle("E[keV]");
  h_iiwi->GetXaxis()->SetLabelSize(0.08);
  h_iiwi->GetYaxis()->SetLabelSize(0.08);
  h_iiwi->GetXaxis()->SetTitleSize(0.08);
  h_iiwi->GetXaxis()->SetTitleOffset(0.52);
  h_iiwi->GetXaxis()->CenterTitle();
  h_iiwi->SetMaximum(100);
  h_iiwi->Draw();

  c1->Update();
  TLine *l=new TLine(1430,0,1430,100);
  l->SetLineColor(2);
  l->SetLineStyle(9);
  l->Draw();
  
  c1->cd(3);
    
  honu->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *honu_length = honu->GetV1();
  double *honu_energy = honu->GetV2();
  TH1F *h_honu = new TH1F("h_honu", "Honu", 100, 0, 2500);
  for(int i=0; i<honu->GetSelectedRows(); i++){
    h_honu->Fill(honu_energy[i]);
  }
  h_honu->SetLineWidth(3);
  h_honu->GetXaxis()->SetTitle("E[keV]");
  h_honu->GetXaxis()->SetLabelSize(0.08);
  h_honu->GetYaxis()->SetLabelSize(0.08);
  h_honu->GetXaxis()->SetTitleSize(0.08);
  h_honu->GetXaxis()->SetTitleOffset(0.52);
  h_honu->GetXaxis()->CenterTitle();
  h_honu->SetMaximum(100);
  h_honu->Draw();
  c1->Update();
  l->Draw();

  c1->cd(5);
    
  kohola->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *kohola_length = kohola->GetV1();
  double *kohola_energy = kohola->GetV2();
  TH1F *h_kohola = new TH1F("h_kohola", "Kohola", 100, 0, 2500);
  for(int i=0; i<kohola->GetSelectedRows(); i++){
    h_kohola->Fill(kohola_energy[i]);
  }
  h_kohola->SetLineWidth(3);
  h_kohola->GetXaxis()->SetTitle("E[keV]");
  h_kohola->GetXaxis()->SetLabelSize(0.08);
  h_kohola->GetXaxis()->SetTitleSize(0.08);
  h_kohola->GetYaxis()->SetLabelSize(0.08);
  h_kohola->GetXaxis()->SetTitleOffset(0.52);
  h_kohola->GetXaxis()->CenterTitle();
  h_kohola->SetMaximum(100);
  h_kohola->Draw();
  c1->Update();
  l->Draw();

  c1->cd(7);

  nene->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *nene_length = nene->GetV1();
  double *nene_energy = nene->GetV2();
  TH1F *h_nene = new TH1F("h_nene", "Nene", 100, 0, 2500);
  for(int i=0; i<nene->GetSelectedRows(); i++){
    h_nene->Fill(nene_energy[i]);
  }
  h_nene->SetLineWidth(3);
  h_nene->GetXaxis()->SetTitle("E[keV]");
  h_nene->GetXaxis()->SetLabelSize(0.08);
  h_nene->GetYaxis()->SetLabelSize(0.08);
  h_nene->GetXaxis()->SetTitleSize(0.08);
  h_nene->GetXaxis()->SetTitleOffset(0.52);
  h_nene->GetXaxis()->CenterTitle();
  h_nene->SetMaximum(100);
  h_nene->Draw();
  c1->Update();
  l->Draw();

  c1->cd(2);

  tako->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *tako_length = tako->GetV1();
  double *tako_energy = tako->GetV2();
  TH1F *h_tako = new TH1F("h_tako", "Tako", 100, 0, 2500);
  for(int i=0; i<tako->GetSelectedRows(); i++){
    h_tako->Fill(tako_energy[i]);
  }
  h_tako->SetLineWidth(3);
  h_tako->GetXaxis()->SetTitle("E[keV]");
  h_tako->GetXaxis()->SetLabelSize(0.08);
  h_tako->GetYaxis()->SetLabelSize(0.08);
  h_tako->GetXaxis()->SetTitleSize(0.08);
  h_tako->GetXaxis()->SetTitleOffset(0.52);
  h_tako->GetXaxis()->CenterTitle();
  h_tako->SetMaximum(100);
  h_tako->Draw();
  c1->Update();
  l->Draw();

  c1->cd(4);

  humu->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *humu_length = humu->GetV1();
  double *humu_energy = humu->GetV2();
  TH1F *h_humu = new TH1F("h_humu", "Humu", 100, 0, 2500);
  for(int i=0; i<humu->GetSelectedRows(); i++){
    h_humu->Fill(humu_energy[i]);
  }
  h_humu->SetLineWidth(3);
  h_humu->GetXaxis()->SetTitle("E[keV]");
  h_humu->GetXaxis()->SetLabelSize(0.08);
  h_humu->GetYaxis()->SetLabelSize(0.08);
  h_humu->GetXaxis()->SetTitleSize(0.08);
  h_humu->GetXaxis()->SetTitleOffset(0.52);
  h_humu->GetXaxis()->CenterTitle();
  h_humu->SetMaximum(100);
  h_humu->Draw();
  c1->Update();
  l->Draw();

  c1->cd(6);

  palila->Draw("length:paper_calibrated_track_energy", palila_alpha); //energies assuming gain of 2000
  double *palila_length = palila->GetV1();
  double *palila_energy = palila->GetV2();
  TH1F *h_palila = new TH1F("h_palila", "Palila", 100, 0, 2500);
  for(int i=0; i<palila->GetSelectedRows(); i++){
    h_palila->Fill(palila_energy[i]);
  }
  h_palila->SetLineWidth(3);
  h_palila->GetXaxis()->SetTitle("E[keV]");
  h_palila->GetXaxis()->SetLabelSize(0.08);
  h_palila->GetYaxis()->SetLabelSize(0.08);
  h_palila->GetXaxis()->SetTitleSize(0.08);
  h_palila->GetXaxis()->SetTitleOffset(0.52);
  h_palila->GetXaxis()->CenterTitle();
  h_palila->SetMaximum(100);
  h_palila->Draw();
  c1->Update();
  l->Draw();

  c1->cd(8);

  elepaio->Draw("length:paper_calibrated_track_energy", alpha); //energies assuming gain of 2000
  double *elepaio_length = elepaio->GetV1();
  double *elepaio_energy = elepaio->GetV2();
  TH1F *h_elepaio = new TH1F("h_elepaio", "Elepaio", 100, 00, 2500);
  for(int i=0; i<elepaio->GetSelectedRows(); i++){
    h_elepaio->Fill(elepaio_energy[i]);
  }
  h_elepaio->SetLineWidth(3);
  h_elepaio->GetXaxis()->SetTitle("E[keV]");
  h_elepaio->GetXaxis()->SetLabelSize(0.08);
  h_elepaio->GetYaxis()->SetLabelSize(0.08);
  h_elepaio->GetXaxis()->SetTitleSize(0.08);
  h_elepaio->GetXaxis()->SetTitleOffset(0.52);
  h_elepaio->GetXaxis()->CenterTitle();
  h_elepaio->SetMaximum(100);
  h_elepaio->Draw();
  c1->Update();
  l->Draw();
}

 

