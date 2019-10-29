// USAGE:
//   root -l
//     .L event_viewer.C
//     event_viewer("module_id", run_num)
//     view_event(2)
//     view_event(30)
//     ch->Show("evt_num", cut), etc.

#include <string>

TCanvas *c2;
TChain *ch;
TGraph *gr;
TGraph *gr1;
TLine *line1;
TLine *line2;
TMultiGraph *mg;
TH2F *h1;
TH1F *h2;
TH1F *h3;
//std::string module_id;
//int run_num;

TCut alpha = "hitside==11";

//TCut neutron = "hitside==0 && length > (800/(2.811*59.5)*calibrated_track_energy+350) && length < (18000/(2.811*210)*calibrated_track_energy + 1000) && length < 22000 && calibrated_track_energy > 2.811*35"; //factors of 2.811 to scale energy and dE/dx cuts to give absolute energies

TCut neutron = "hitside==0 && paper_calibrated_track_energy < (0.25*length-75) && paper_calibrated_track_energy > (0.0411764*length-64.688)  && paper_calibrated_track_energy > 100";

TCut xray = "hitside == 0 && timestamp_start[0]>1000 && length > 2500 && calibrated_track_energy < 2.811*22 && length > (18000/(2.811*210)*calibrated_track_energy + 1000)";

TCut shower = "timestamp_start[0]>1000 && calibrated_track_energy>2.811*3500";

//TCut notneutron = "hitside==0 && timestamp_start[0]>1000 && (length < (800/(2.811*59.5)*calibrated_track_energy+ 350) || length > (18000/(2.811*210)*calibrated_track_energy + 1000)) && calibrated_track_energy > 2.811*11"; //factors of 2.811 for same reason as above

//TCut mesh = "hitside==0 && timestamp_start[0]>1000 && length > 22000 && calibrated_track_energy > 2.811*100"; 

void event_viewer(std::string module_id, int run_num) {
  gStyle->SetOptStat(0);

  ch = new TChain("data");

  ch->SetMarkerStyle(20);
  ch->SetMarkerSize(1.0);
  ch->Add(Form("~/data/phase2/%s/%d_%s.root", module_id.c_str(),run_num, module_id.c_str()));
  cout << "Added run " << run_num << " from module: " << module_id << endl;
}

void view_event(std::string module_id, int run_num, int entry) {
  std::stringstream ss_cut;
  ss_cut << "(Entry$==" << entry << ")";
  std::string entry_cut = ss_cut.str();
  ss_cut << "*(tot+1)";
  TCut evt_num = Form("event_number == %d+1", entry);

  TH2F *h1 = new TH2F("h1","2D event display", 80, 0, 80, 336, 0, 336);
  TH1F *h2 = new TH1F("h2","BCID distribution", 100, 0, 100);
  TH1F *h3 = new TH1F("h3","TOT distribution", 14,0,14);

  c2 = new TCanvas("c2","c2",1600,1000); 
  c2->Divide(2,3);

  c2->cd(1);
  ch->Draw("row:column>>h1", ss_cut.str().c_str(),"COLZ");

  c2->cd(2);
  ch->Draw("time_bin:row:column:(tot+1)", entry_cut.c_str());

  c2->cd(3);
  ch->Draw("time_bin>>h2", entry_cut.c_str());

  c2->cd(4);
  ch->Draw("(tot+1)>>h3", entry_cut.c_str());

  c2->cd(5);
  ch->Draw("length:paper_calibrated_track_energy", evt_num, "goff");
  double *single_length = ch->GetV1();
  double *single_energy = ch->GetV2();

  TGraph *gr = new TGraph(1, single_length, single_energy); //to be included in multigraph: plots event in red
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(2);

  c2->cd(6);
  ch->Draw("paper_calibrated_track_energy:length");
  double *length = ch->GetV2();
  double *energy = ch->GetV1();
  int total = ch->GetSelectedRows();
  
  TGraph *gr1 = new TGraph(total, length, energy); //to be included in multigraph
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.1);
  gr1->SetMarkerColor(4);

  //For selection line plotting
  double x_high[2] = {600, 40300};
  double y_high[2] = {75, 10000}; 
  double x_low[2] = {3513.857, 98714};
  double y_low[2] = {80, 4000};
  double x_flat[2] = {0, 30000};
  double y_flat[2] = {100, 100};
  
  TF1 *f_high = new TF1("f_high", "[0]+[1]*x",599,40301);
  TF1 *f_low = new TF1("f_low", "[0]+[1]*x",3512.857,98715);
  TF1 *f_flat = new TF1("f_flat", "[0]+[1]*x");

  TGraph *gr_high = new TGraph(2, x_high, y_high);
  gr_high->SetMarkerSize(0);
  gr_high->Fit("f_high", "R");
  gr_high->SetMaximum(4000);
  gr_high->SetMinimum(80);
  gr_high->GetXaxis()->SetLimits(0,30000);
  gr_high->SetLineColor(2);
  gr_high->GetFunction("f_high")->SetLineWidth(1);
  //gr_high->GetFunction("f_high")->GetXaxis()->SetLimits(620,20000);
  
  TGraph *gr_low = new TGraph(2, x_low, y_low);
  gr_low->SetMarkerSize(0);
  gr_low->Fit("f_low", "R");
  gr_low->SetMaximum(4000);
  gr_low->SetMinimum(80);
  gr_low->GetXaxis()->SetLimits(0,30000);
  gr_low->GetFunction("f_low")->SetLineWidth(1);

  TGraph *gr_flat = new TGraph(2, x_flat, y_flat);
  gr_flat->SetMarkerSize(0);
  gr_flat->SetLineColor(1);
  gr_flat->SetLineWidth(7);
  gr_flat->Fit("f_flat", "S");
  gr_flat->SetMaximum(4000);
  gr_flat->SetMinimum(80);
  gr_flat->GetXaxis()->SetLimits(0,30000);
  gr_flat->GetFunction("f_flat")->SetLineColor(1);
  gr_flat->GetFunction("f_flat")->SetLineWidth(1);


  //Plotting everything
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr);
  mg->Add(gr_high);
  mg->Add(gr_low);
  mg->Add(gr_flat);
  mg->SetTitle("Energy vs. Length; Length(um); Recoil Energy(keV)");
  mg->GetXaxis()->SetLimits(0,30000);
  mg->SetMaximum(2000);
  mg->SetMinimum(-50);
  mg->Draw("AP");
  gPad->Modified();

  c2->Update();
  //c2->SaveAs(Form("/Users/vahsengrouplaptop/Desktop/root_scans/above_section_boundary/%s/%s_run_%d_event_%d.png", module_id.c_str(), module_id.c_str(), run_num, entry));


}

