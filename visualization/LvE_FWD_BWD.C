//################ Data Chains for each TPC ##################
TChain *FWD = new TChain("data");
TChain *BWD = new TChain("data");
//############################################################

void LvE_FWD_BWD(int first_run, int last_run) {
  
  cout << "Generating summary ntuple using run " << first_run << " through " << last_run << " from all TPCs" <<endl;

//Populate Chains 
  for(int iRun=first_run; iRun<=last_run; iRun++) { 
    BWD->Add(Form("~/data/phase2/iiwi/%d_iiwi.root",iRun));
    cout << "Added run " << iRun << " from module iiwi" << endl;
 
    BWD->Add(Form("~/data/phase2/honu/%d_honu.root",iRun));
    cout << "Added run " << iRun << " from module honu" << endl;

    BWD->Add(Form("~/data/phase2/kohola/%d_kohola.root",iRun));
    cout << "Added run " << iRun << " from module kohola" << endl;

    BWD->Add(Form("~/data/phase2/nene/%d_nene.root",iRun));
    cout << "Added run " << iRun << " from module nene" << endl;

    FWD->Add(Form("~/data/phase2/tako/%d_tako.root",iRun));
    cout << "Added run " << iRun << " from module tako" << endl;
  
    FWD->Add(Form("~/data/phase2/humu/%d_humu.root",iRun));
    cout << "Added run " << iRun << " from module humu" << endl;

    FWD->Add(Form("~/data/phase2/palila/%d_palila.root",iRun));
    cout << "Added run " << iRun << " from module palila" << endl;

    FWD->Add(Form("~/data/phase2/elepaio/%d_elepaio.root",iRun));
    cout << "Added run " << iRun << " from module elepaio" << endl;
  }

TCut alpha = "hitside==11";

TCut palila_alpha = "hitside==111";

//TCut neutron = "hitside==0 && timestamp_start[0]>1000 && length > (800/(2.811*59.5)*calibrated_track_energy+350) && length < (18000/(2.811*210)*calibrated_track_energy + 1000) && length < 22000 && calibrated_track_energy > 2.811*35";

//TCut notneutron = "hitside==0 && timestamp_start[0]>1000 && (length < (800/(2.811*59.5)*calibrated_track_energy+350) || length > (18000/(2.811*210)*calibrated_track_energy + 1000)) || (length > 22000 || calibrated_track_energy < 2.811*35)"; 

//TCut neutron = "hitside==0 && length > (4*september_calibrated_track_energy+300) && length < ((17000/700)*september_calibrated_track_energy + 1571)  && september_calibrated_track_energy > 80";

TCut neutron = "hitside==0 && paper_calibrated_track_energy < (0.25*length-75) && paper_calibrated_track_energy > (0.0411764*length-64.688)  && paper_calibrated_track_energy > 80";
 
TCut xray = "hitside == 0 && length > 2500 && calibrated_track_energy < 2.811*22 && length > (18000/(2.811*210)*calibrated_track_energy + 1000)";

TCut other = "hitside!= 0 && hitside";

 TCut test = "hitside == 0";

  
TCanvas *c1 = new TCanvas("c1", "L vs. E for BWD TPCs" ,800,600);
  c1->cd(1);
  BWD->Draw("length:paper_calibrated_track_energy", "hitside == 0", "goff");
  //BWD->Draw("length:paper_calibrated_track_energy", neutron, "goff"); 
  double *l_BWD = BWD->GetV1();
  double *e_BWD = BWD->GetV2();

  double x_high[2] = {600, 40300};
  double y_high[2] = {75, 10000}; 
  double x_low[2] = {3513.857, 98714};
  double y_low[2] = {80, 4000};
  double x_flat[2] = {0, 20000};
  double y_flat[2] = {80, 80};
  
  TF1 *f_high = new TF1("f_high", "[0]+[1]*x",599,40301);
  TF1 *f_low = new TF1("f_low", "[0]+[1]*x",3512.857,98715);
  TF1 *f_flat = new TF1("f_flat", "[0]+[1]*x");
  
  TGraph *gr1 = new TGraph(BWD->GetSelectedRows(), l_BWD, e_BWD);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.05);
  gr1->SetLineWidth(0);
  gr1->SetTitle("BWD TPCs");

  FWD->Draw("length:paper_calibrated_track_energy", "hitside == 0", "goff");
  //FWD->Draw("length:paper_calibrated_track_energy", neutron, "goff");
  double *l_FWD = FWD->GetV1();
  double *e_FWD = FWD->GetV2();

  TGraph *gr2 = new TGraph(FWD->GetSelectedRows(), l_FWD, e_FWD);
  gr2->SetMarkerColor(6);
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.05);
  gr2->SetLineWidth(0);
  gr2->SetTitle("FWD TPCs");

  TGraph *gr_high = new TGraph(2, x_high, y_high);
  gr_high->SetMarkerSize(0);
  gr_high->Fit("f_high", "R");
  gr_high->SetMaximum(4000);
  gr_high->SetMinimum(80);
  gr_high->GetXaxis()->SetLimits(0,30000);
  gr_high->SetLineColor(2);
  gr_high->GetFunction("f_high")->SetLineWidth(4);
  //gr_high->GetFunction("f_high")->GetXaxis()->SetLimits(620,20000);
  
  TGraph *gr_low = new TGraph(2, x_low, y_low);
  gr_low->SetMarkerSize(0);
  gr_low->Fit("f_low", "R");
  gr_low->SetMaximum(4000);
  gr_low->SetMinimum(80);
  gr_low->GetXaxis()->SetLimits(0,30000);
  gr_low->GetFunction("f_low")->SetLineWidth(4);

  TGraph *gr_flat = new TGraph(2, x_flat, y_flat);
  gr_flat->SetMarkerSize(0);
  gr_flat->SetLineColor(1);
  gr_flat->SetLineWidth(7);
  gr_flat->Fit("f_flat", "S");
  gr_flat->SetMaximum(4000);
  gr_flat->SetMinimum(80);
  gr_flat->GetXaxis()->SetLimits(0,30000);
  gr_flat->GetFunction("f_flat")->SetLineColor(1);
  gr_flat->GetFunction("f_flat")->SetLineWidth(7);
  
  TGraph *gr_for_legend = new TGraph(2, x_low, y_low);
  gr_for_legend->SetMarkerStyle(20);
  gr_for_legend->SetMarkerColor(4);

  TGraph *gr_for_legend2 = new TGraph(2, x_low, y_low);
  gr_for_legend2->SetMarkerStyle(20);
  gr_for_legend2->SetMarkerColor(6);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr_high);
  mg->Add(gr_low);
  mg->Add(gr_flat);
  mg->SetTitle("Energy vs. Length; Length(um); Recoil Energy(keV)");
  mg->Draw("AP");
  gPad->Modified();
  mg->GetXaxis()->SetLimits(0,20000);
  mg->SetMaximum(2000);
  mg->SetMinimum(-50);

  TLegend *l = new TLegend(0.1,0.7,0.28,0.6);
  l->AddEntry(gr_for_legend, "BWD TPCs", "p");
  l->AddEntry(gr_for_legend2, "FWD TPCs", "p");
  l->AddEntry(gr_high, "dE/dx selection", "l");
  l->AddEntry(gr_flat, "80 keV Recoil Energy Threshold", "l");
  l->Draw();
  
}
