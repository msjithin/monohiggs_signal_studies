/////etau_analyzer.C
// For use with Ntuples made from ggNtuplizer
// Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
// To compile using rootcom to an executable named 'analyze':
//$ ./rootcom etau_analyzer analyze
//
// To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jmadhusu/LatestNtuples/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/etau/output.root -1 10000
//./analyze /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_ZZZ/180603_185329/0000/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/analyzer/output.root -1 10000
// Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
// and storing the resulting histograms in the file output.root.
//
// To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
// root[0] TFile *f = new TFile("output.root");
// root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
// root[2] efficiency->Draw("AP")
// root[3] efficiency->SetTitle("Single photon trigger efficiency")
// root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
// root[5] efficiency->Draw("AP")
//
#define etau_analyzer_cxx
#include "etau_analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>
#include "makeHisto.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include <numeric>
#include "commonFunctions.h"
#include "fractions.C"
//#include "ApplyFF.h"
#include "object_functions.h"
#include "fill_histograms.h"
#include "selections.h"
#include "get_met_systenatics.h"
//#include "myclass.cc"

using namespace std;
using std::vector;
int main(int argc, const char *argv[])
{
  TStopwatch sw;
  sw.Start();

  myMap1 = new map<string, TH1F *>();
  myMap2 = new map<string, TH2F*>();
  std::string SampleName = argv[7];
  std::string isMC = argv[6];
  std::string year = argv[5];
  std::string outputfile = argv[2];
  Long64_t maxEvents = atof(argv[3]);
  string sp = "0";
  cout << "argc = " << argc << endl;
  if (argc > 8)
    sp = string(argv[8]);
  cout << "sp = " << sp << endl;
  const char *signalpara = sp.c_str();

  if (maxEvents < -1LL)
  {
    std::cout << "Please enter a valid value for maxEvents (parameter 3)." << std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout << "Please enter a valid value for reportEvery (parameter 4) " << std::endl;
    return 1;
  }
  // std::string SampleName = argv[5];

  etau_analyzer t(argv[1], argv[2], isMC, SampleName, year, signalpara);
  t.Loop(maxEvents, reportEvery, SampleName);
  // delete myMap1;
  cout << " Outpt written to " << outputfile << endl;
  sw.Stop();
  sw.Print();
  return 0;
}

void etau_analyzer::BookHistos(const char *file1, const char *file2, string year)
{
  TFile *file_in = new TFile(file1, "READ");
  fileName = new TFile(file2, "RECREATE");

  // makeOutputTree(tree);
  fileName->cd();
  h_nEvents = (TH1F *)((TH1F *)file_in->Get("nEvents"))->Clone(TString("nEvents"));
  TH1F *h_luminosity = new TH1F("luminosity", "luminosity", 2, 0, 2);
  if (year == "2016")
    h_luminosity->SetBinContent(1, 39500.0);
  else if (year == "2017")
    h_luminosity->SetBinContent(1, 41520.0);
  else if (year == "2018")
    h_luminosity->SetBinContent(1, 59700.0);

  file_in->Close();
  Float_t Pt_Bins[36] = {0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 600, 800, 1000};

  Float_t MetBins[15] = {0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180., 200, 300., 400., 600.0, 800.0};
  Float_t TrMassBins[24] = {0.0, 20, 40, 60, 80, 100, 120, 140, 160, 180., 200, 220, 240, 260, 280, 300., 320, 340, 360, 380, 400., 600.0, 800.0, 1000.0};
}

// Fill the sequential histos at a particular spot in the sequence

void etau_analyzer::fillHistos(int histoNumber, double event_weight, int higgs_Index)
{

  //  h_HiggsPt[histoNumber]->Fill(mcPt->at(higgs_Index),event_weight);
}

void etau_analyzer::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{

  int nTotal;
  nTotal = 0;
  int report_ = 0;
  int report_test = 0;
  double numberOfEvents = 0;
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;
  nInspected_genWeighted = 0.0;
  debug = false;
  if (debug)
    cout << "******** debugging is on ******************" << endl;
  if (fChain == 0)
    return;

  TString sample = TString(SampleName);

  // TH1F* h_cutflow_Htt=new TH1F("cutflow_Htt", "cutflow_Htt", 11, 0, 11); h_cutflow_Htt->Sumw2();

  Long64_t nentries = fChain->GetEntries();
  if (is_MC == true)
    std::cout << ".... MC file ..... " << std::endl;
  else
    std::cout << ".... DATA file ..... " << std::endl;

  std::cout << "Coming in: " << std::endl;
  std::cout << "nentries:" << nentries << std::endl;
  // Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout << "Running over " << nTotal << " events." << std::endl;
  //// // cout<<__LINE__<<endl;
  for (Long64_t jentry = 0; jentry < nentriesToCheck; jentry++)
  {
    current_event_number = jentry;
    eleCand.clear();
    tauCand.clear();
    aisrtauCand.clear();
    eleCand_nom.clear();
    tauCand_nom.clear();
    aisrtauCand_nom.clear();
    jetCand.clear();
    jetCand_nom.clear();

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    double inspected_event_weight = 1.0;
    if (is_MC)
      fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight / fabs(genWeight) : inspected_event_weight = 0.0;
    nInspected_genWeighted += inspected_event_weight;
    nInspected += 1;
    double event_weight = 1.0;
    double weight = 1.0;
    applySf = 1.0;
    //// cout<<__LINE__<<endl;
    double pileup_sf = 1.0;
    passSingleTriggerPaths = false;
    passCrossTrigger = false;
    int report_i = 0;
    Ztt_selector = false;

    numberOfEvents += weight;
    if (is_MC)
      weight = inspected_event_weight;
    else
      weight = 1.0;

    //// cout<<__LINE__<<endl;
    // make_met_plots("a1");
    plotFill("prefiringweight_incoming", prefiringweight, 44, 0, 1.1, 1);
    plotFill("prefiringweight_incoming_up", prefiringweightup, 44, 0, 1.1, 1);
    plotFill("prefiringweight_incoming_dn", prefiringweightdown, 44, 0, 1.1, 1);
    // prefiringweight = 1.0;
    if (is_MC)
      pileup_sf = h_pileup->GetBinContent(h_pileup->GetXaxis()->FindBin(puTrue->at(0)));
    weight = weight * pileup_sf;
    if (is_MC)
      weight = weight * prefiringweight;

    if (isGoodVtx == false)
      continue;
    t_index = get_t_Cand();
    tbar_index = get_tbar_Cand();

    if (((HLTEleMuX >> 3 & 1 == 1)     // HLT_Ele27_WPTight_Gsf_v (HLTEleMuX>>3&1 == 1 )
         || (HLTEleMuX >> 61 & 1 == 1) // HLT_Ele32_WPTight_Gsf_v
         || (HLTEleMuX >> 5 & 1 == 1)  // HLT_Ele35_WPTight_Gsf_v
         ))
      passSingleTriggerPaths = true; //

    if ((HLTTau >> 1 & 1 == 1)) // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
      passCrossTrigger = true;

    /////
    if (debug)
      cout << "entry # : " << jentry << endl;
    eventNumber = jentry;

    if (!is_MC)
      event_weight = 1.0;
    else
      event_weight = weight;

    if (metFilters == 0 && (passSingleTriggerPaths || passCrossTrigger))
    {
      event_weight = 1.0;
      /// fill gen taus
      string hNumber = "0";
      int gen_tau_count = 0;
      // cout<<__LINE__<<endl;
      int count = 0;
      bool found_tau = false;
      bool found_etau = false;
      int gen_tau_parent1 = -99;
      int gen_tau_parent2 = -99;
      int gen_tau_e_idx = -99;
      int gen_tau_h_idx = -99;
      TLorentzVector reco_lep1, reco_lep2, gen_lep1, gen_lep2;
      for (int imc = 0; imc < nMC; imc++)
      {
        // cout<<__LINE__<<endl;
        if (mcPt->at(imc) > 15 && abs(mcEta->at(imc)) < 2.3 && abs(mcPID->at(imc)) == 15 && mcStatus->at(imc) == 2)
          gen_tau_count++;

        if (mcPt->at(imc) > 15 && abs(mcEta->at(imc)) < 2.3 && abs(mcMotherPID->at(imc)) == 15)
        {
          plotFill("gen_tau_children" + hNumber, abs(mcPID->at(imc)), 300, 0, 300, event_weight);
        }
      }
      for (int imc = 0; imc < nMC; imc++) //// loop to get gen hadronic tau
      {
        int mom_pid = 0;
        if (mcPt->at(imc) > 15 && abs(mcEta->at(imc)) < 2.3 && abs(mcPID->at(imc)) == 15 && mcStatus->at(imc) == 2) // found a tau
        {
          mom_pid = mcPID->at(imc);
          for (int jmc = imc; jmc < nMC; jmc++)
          {
            if (mcPt->at(jmc) > 15 && abs(mcEta->at(jmc)) < 2.3 && mcMotherPID->at(jmc) == mom_pid && abs(mcPID->at(jmc)) > 18 && found_tau != true)
            {
              found_tau = true;
              gen_tau_parent2 = imc;
              gen_tau_h_idx = imc;
            }
            if (gen_tau_h_idx >= 0)
              break;
          }
        }
        if (gen_tau_h_idx >= 0)
          break;
      }

      for (int imc = 0; imc < nMC; imc++) //// loop to get gen tau-e
      {
        int mom_pid = 0;
        if (mcPt->at(imc) > 15 && abs(mcEta->at(imc)) < 2.3 && abs(mcPID->at(imc)) == 15 && mcStatus->at(imc) == 2) // found a tau
        {
          mom_pid = mcPID->at(imc);
          for (int jmc = imc; jmc < nMC; jmc++)
          {
            if (mcPt->at(jmc) > 15 && abs(mcEta->at(jmc)) < 2.3 && mcMotherPID->at(jmc) == mom_pid && abs(mcPID->at(jmc)) == 11 && found_etau != true)
            {
              found_etau = true;
              gen_tau_e_idx = jmc;
              gen_tau_parent1 = imc;
            }
            if (gen_tau_e_idx >= 0)
              break;
          }
        }
        if (gen_tau_e_idx >= 0)
          break;
      }

      /// filling for reco
      nSingleTrgPassed += event_weight;
      tauCand.clear();
      eleCand.clear();
      tauCand = getTauCand_med(15.0, 2.3, 0);
      eleCand = getEleCand(15.0, 2.3, 0);
      int reco_tau_e_idx = -99;
      int reco_tau_h_idx = -99;
      make_plots("0", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
      // tauCand = getTauCand(30.0, 2.3, 0 ); // from analysis
      // cout<< "n reco taus "<< tauCand.size() <<endl;
      // cout<<__LINE__<<endl;

      if (eleCand.size() > 0 && tauCand.size() > 0)
      {
        // cout<<__LINE__<<endl;
        for (int i = 0; i < eleCand.size(); i++)
        {
          for (int j = 0; j < tauCand.size(); j++)
          {
            if (eleCharge->at(eleCand[i]) * tau_Charge->at(tauCand[j]) < 0)
            {
              reco_tau_e_idx = eleCand[i];
              reco_tau_h_idx = tauCand[j];
              break;
            }
          }
          // cout<<__LINE__<<endl;
          if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
            break;
        }
        // cout<<__LINE__<<endl;
        if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
        {
          make_plots("med_2", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
          // cout<<__LINE__<<endl;
          reco_lep1.SetPtEtaPhiE(elePt->at(reco_tau_e_idx), eleEta->at(reco_tau_e_idx), elePhi->at(reco_tau_e_idx), eleE->at(reco_tau_e_idx));
          reco_lep2.SetPtEtaPhiE(tau_Pt->at(reco_tau_h_idx), tau_Eta->at(reco_tau_h_idx), tau_Phi->at(reco_tau_h_idx), tau_Energy->at(reco_tau_h_idx));
          // cout<<__LINE__<<endl;
          if (thirdLeptonVeto(reco_tau_e_idx, reco_tau_h_idx))
          {
            // cout<<__LINE__<<endl;
            make_plots("med_3", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
            // cout<<__LINE__<<endl;
            if ((bJet_medium(reco_tau_e_idx, reco_tau_h_idx).size() == 0) && (bJet_loose(reco_tau_e_idx, reco_tau_h_idx).size() < 2))
            {
              // cout<<__LINE__<<endl;
              make_plots("med_4", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              // cout<<__LINE__<<endl;
              double deltaR = reco_lep1.DeltaR(reco_lep2);
              bool pass_dr = false;
              // cout<<__LINE__<<endl;
              if (reco_lep2.Pt() >= 100 && reco_lep1.Pt() >= 100)
                pass_dr = true;
              else if (deltaR > 0.5)
                pass_dr = true;
              make_plots("med_5", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              // cout<<__LINE__<<endl;
            }
          }
        }
      }
      // cout<<__LINE__<<endl;
      tauCand.clear();
      eleCand.clear();
      tauCand = getTauCand_vl(15.0, 2.3, 0);
      eleCand = getEleCand(15.0, 2.3, 0);
      reco_tau_e_idx = -99;
      reco_tau_h_idx = -99;
      if (eleCand.size() > 0 && tauCand.size() > 0)
      {
        // cout<<__LINE__<<endl;
        for (int i = 0; i < eleCand.size(); i++)
        {
          for (int j = 0; j < tauCand.size(); j++)
          {
            if (eleCharge->at(eleCand[i]) * tau_Charge->at(tauCand[j]) < 0)
            {
              reco_tau_e_idx = eleCand[i];
              reco_tau_h_idx = tauCand[j];
              break;
            }
          }
          // cout<<__LINE__<<endl;
          if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
            break;
        }
        // cout<<__LINE__<<endl;
        if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
        {
          // cout<<__LINE__<<endl;
          make_plots("vvvl_2", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
          reco_lep1.SetPtEtaPhiE(elePt->at(reco_tau_e_idx), eleEta->at(reco_tau_e_idx), elePhi->at(reco_tau_e_idx), eleE->at(reco_tau_e_idx));
          reco_lep2.SetPtEtaPhiE(tau_Pt->at(reco_tau_h_idx), tau_Eta->at(reco_tau_h_idx), tau_Phi->at(reco_tau_h_idx), tau_Energy->at(reco_tau_h_idx));
          // cout<<__LINE__<<endl;
          if (thirdLeptonVeto(reco_tau_e_idx, reco_tau_h_idx))
          {
            // cout<<__LINE__<<endl;
            make_plots("vvvl_3", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
            if ((bJet_medium(reco_tau_e_idx, reco_tau_h_idx).size() == 0) && (bJet_loose(reco_tau_e_idx, reco_tau_h_idx).size() < 2))
            {
              // cout<<__LINE__<<endl;
              make_plots("vvvl_4", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              double deltaR = reco_lep1.DeltaR(reco_lep2);
              bool pass_dr = false;
              // cout<<__LINE__<<endl;
              if (reco_lep2.Pt() >= 100 && reco_lep1.Pt() >= 100)
                pass_dr = true;
              else if (deltaR > 0.5)
                pass_dr = true;
              make_plots("vvvl_5", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              // cout<<__LINE__<<endl;
            }
          }
        }
      }
// cout<<__LINE__<<endl;
      tauCand.clear();
      eleCand.clear();
      tauCand = getTauCand_superloose(15.0, 2.3, 0);
      eleCand = getEleCand(15.0, 2.3, 0);
      reco_tau_e_idx = -99;
      reco_tau_h_idx = -99;
      if (eleCand.size() > 0 && tauCand.size() > 0)
      {
        // cout<<__LINE__<<endl;
        for (int i = 0; i < eleCand.size(); i++)
        {
          for (int j = 0; j < tauCand.size(); j++)
          {
            if (eleCharge->at(eleCand[i]) * tau_Charge->at(tauCand[j]) < 0)
            {
              reco_tau_e_idx = eleCand[i];
              reco_tau_h_idx = tauCand[j];
              break;
            }
          }
          // cout<<__LINE__<<endl;
          if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
            break;
        }
        if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
        {
          // cout<<__LINE__<<endl;
          make_plots("sl_2", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
          reco_lep1.SetPtEtaPhiE(elePt->at(reco_tau_e_idx), eleEta->at(reco_tau_e_idx), elePhi->at(reco_tau_e_idx), eleE->at(reco_tau_e_idx));
          reco_lep2.SetPtEtaPhiE(tau_Pt->at(reco_tau_h_idx), tau_Eta->at(reco_tau_h_idx), tau_Phi->at(reco_tau_h_idx), tau_Energy->at(reco_tau_h_idx));
          if (thirdLeptonVeto(reco_tau_e_idx, reco_tau_h_idx))
          {
            // cout<<__LINE__<<endl;
            make_plots("sl_3", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
            if ((bJet_medium(reco_tau_e_idx, reco_tau_h_idx).size() == 0) && (bJet_loose(reco_tau_e_idx, reco_tau_h_idx).size() < 2))
            {
              // cout<<__LINE__<<endl;
              make_plots("sl_4", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              double deltaR = reco_lep1.DeltaR(reco_lep2);
              bool pass_dr = false;
              // cout<<__LINE__<<endl;
              if (reco_lep2.Pt() >= 100 && reco_lep1.Pt() >= 100)
                pass_dr = true;
              else if (deltaR > 0.5)
                pass_dr = true;
              make_plots("sl_5", reco_tau_e_idx, reco_tau_h_idx, gen_tau_h_idx, gen_tau_e_idx);
              // cout<<__LINE__<<endl;
            }
          }
        }
      }
    }
    // cout<<__LINE__<<endl;
    // // cout<<__LINE__<<endl;
    previous_event_number = current_event_number;
    report_test = nentriesToCheck / 20;
    while (report_test > 10)
    {
      report_test = report_test / 10;
      report_i++;
    }
    if (nentriesToCheck > 20)
      reportEvery = report_test * pow(10, report_i);
    else
      reportEvery = 1;
    if (jentry % reportEvery == 0)
    {
      std::cout << "Finished entry " << jentry << "/" << (nentriesToCheck - 1) << std::endl;
    }
  }

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  if ((nentriesToCheck - 1) % reportEvery != 0)
    std::cout << "Finished entry " << (nentriesToCheck - 1) << "/" << (nentriesToCheck - 1) << std::endl;
  // sw.Stop();
  std::cout << "All events checked." << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << "******************Jithin's original*************************" << std::endl;
  std::cout << std::setw(20) << std::right << "Initial entries " << numberOfEvents << std::endl;
  std::cout << std::setw(20) << std::right << "Passing smikking " << nPassedSkimmed << std::endl;

  std::cout << std::setw(20) << std::right << "GoodMuonPassed " << nGoodMuonPassed << std::setw(10) << std::right << "   % change= " << (nSingleTrgPassed - nGoodMuonPassed) * 100 / nSingleTrgPassed << std::endl;
  std::cout << std::setw(20) << std::right << "GoodTauPassed " << nGoodTauPassed << std::setw(10) << std::right << "   % change= " << (nGoodMuonPassed - nGoodTauPassed) * 100 / nGoodMuonPassed << std::endl;
  //   std::cout<<std::setw(20) <<std::right <<"TauIsoPassed "<<nTauIsoPassed<<std::setw(10) <<std::right << "   % change= "<<(nGoodTauPassed-nTauIsoPassed)*100/nGoodTauPassed<<std::endl;
  // std::cout<<std::setw(20) <<std::right <<"TauDecayModePassed "<<nTauDecayModePassed<<std::setw(10) <<std::right << "   % change= "<<(nTauIsoPassed-nTauDecayModePassed)*100/nTauIsoPassed<<std::endl;

  std::cout << std::setw(20) << std::right << "opp charge " << nGoodMuTauPassed << std::setw(10) << std::right << "   % change= " << (nGoodTauPassed - nGoodMuTauPassed) * 100 / nGoodTauPassed << std::endl;

  std::cout << std::setw(20) << std::right << "PassedThirdLepVeto " << nPassedThirdLepVeto << std::setw(10) << std::right << "   % change= " << (nGoodMuTauPassed - nPassedThirdLepVeto) * 100 / nGoodMuTauPassed << std::endl;
  std::cout << std::setw(20) << std::right << "PassedBjetVeto " << nPassedBjetVeto << std::setw(10) << std::right << "   % change= " << (nPassedThirdLepVeto - nPassedBjetVeto) * 100 / nPassedThirdLepVeto << std::endl;
  std::cout << std::setw(20) << std::right << "DeltaRPassed " << nDeltaRPassed << std::setw(10) << std::right << "   % change= " << (nPassedBjetVeto - nDeltaRPassed) * 100 / nPassedBjetVeto << std::endl;

  std::cout << std::setw(20) << std::right << "Total change :" << (numberOfEvents - nDeltaRPassed) * 100 / numberOfEvents << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << "*******************************************" << std::endl;
  std::cout << std::setw(20) << std::right << "Number of events inspected: " << nInspected << std::endl;
  std::cout << std::setw(20) << std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl;
  cout << " \n \n " << endl;
  cout << "gen_ntaus = " << gen_ntaus << endl;
  cout << "gentau plus  = " << gen_ntaup << endl;
  cout << "gentau minus = " << gen_ntaum << endl;
}

void etau_analyzer::printP4values(string when = "")
{
  if (check_unc == false)
    return;
  printTabSeparated("entry # : ", eventNumber, "before/after = ", when, "\n",
                    "Shift ", unc_shift, "\n",
                    "electron pt", my_eleP4.Pt(),
                    "electron eta", my_eleP4.Eta(), "\n",
                    "tau pt", my_tauP4.Pt(),
                    "tau eta", my_tauP4.Eta(), "\n",
                    "MET", my_metP4.Pt(), "\n");
}
float etau_analyzer::average(std::vector<double> const &v)
{
  if (v.empty())
  {
    return 0;
  }
  auto const count = static_cast<double>(v.size());
  return std::accumulate(v.begin(), v.end(), 0.0) / count;
}

void etau_analyzer::make_plots(string category, int reco_tau_e_idx, int reco_tau_h_idx, int gen_tau_h_idx, int gen_tau_e_idx)
{
  // cout<<__LINE__<<" cat="<< category <<endl;
  float event_weight = 1.0;
  TLorentzVector reco_lep1, reco_lep2, gen_lep1, gen_lep2;
  if (gen_tau_h_idx >= 0 && gen_tau_e_idx >= 0)
  {
    float event_weight = 1.0;
    
    gen_tau1Pt = mcPt->at(gen_tau_e_idx);
    gen_tau1Eta = mcEta->at(gen_tau_e_idx);
    gen_tau1Phi = mcPhi->at(gen_tau_e_idx);
    // gen_tau1dm = mcTauDecayMode->at(gen_tau_e_idx);
    gen_tau1E = mcE->at(gen_tau_e_idx);
    plotFill("gen_tau1Pt_" + category, gen_tau1Pt, 80, 0, 400, event_weight);
    plotFill("gen_tau1Eta_" + category, gen_tau1Eta, 25, -2.5, 2.5, event_weight);
    plotFill("gen_tau1Phi_" + category, gen_tau1Phi, 30, -3.14, 3.14, event_weight);
    // plotFill("gen_tau1DecayMode_" + category, gen_tau1dm, 12, 0, 12, event_weight);
    gen_lep1.SetPtEtaPhiE(gen_tau1Pt, gen_tau1Eta, gen_tau1Phi, gen_tau1E);

    gen_tau2Pt = mcPt->at(gen_tau_h_idx);
    gen_tau2Eta = mcEta->at(gen_tau_h_idx);
    gen_tau2Phi = mcPhi->at(gen_tau_h_idx);
    gen_tau2dm = mcTauDecayMode->at(gen_tau_h_idx);
    gen_tau2E = mcE->at(gen_tau_h_idx);
    plotFill("gen_tau2Pt_" + category, gen_tau2Pt, 80, 0, 400, event_weight);
    plotFill("gen_tau2Eta_" + category, gen_tau2Eta, 25, -2.5, 2.5, event_weight);
    plotFill("gen_tau2Phi_" + category, gen_tau2Phi, 30, -3.14, 3.14, event_weight);
    plotFill("gen_tau2DecayMode_" + category, gen_tau2dm, 12, 0, 12, event_weight);
    gen_lep2.SetPtEtaPhiE(gen_tau2Pt, gen_tau2Eta, gen_tau2Phi, gen_tau2E);
  }
  // // cout<<__LINE__<<endl;
  if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0)
  {
    // // cout<<__LINE__<<endl;
    int tau1Index = reco_tau_e_idx;
    int tau2Index = reco_tau_h_idx;

    reco_tau1Pt = elePt->at(tau1Index);
    reco_tau1Eta = eleEta->at(tau1Index);
    reco_tau1Phi = elePhi->at(tau1Index);
    // reco_tau1dm  = tau_DecayMode->at(tau1Index);
    reco_tau1E = eleE->at(tau1Index);
    plotFill("reco_tau1Pt_" + category, reco_tau1Pt, 80, 0, 400, event_weight);
    plotFill("reco_tau1Eta_" + category, reco_tau1Eta, 25, -2.5, 2.5, event_weight);
    plotFill("reco_tau1Phi_" + category, reco_tau1Phi, 30, -3.14, 3.14, event_weight);
    // plotFill("reco_tau1DecayMode_"+category, reco_tau1dm , 12, 0, 12,  event_weight);
    plotFill("reco_tau1E_" + category, reco_tau1E, 80, 0, 400, event_weight);
    reco_lep1.SetPtEtaPhiE(reco_tau1Pt, reco_tau1Eta, reco_tau1Phi, reco_tau1E);

    // // cout<<__LINE__<<endl;
    reco_tau2Pt = tau_Pt->at(tau2Index);
    reco_tau2Eta = tau_Eta->at(tau2Index);
    reco_tau2Phi = tau_Phi->at(tau2Index);
    reco_tau2dm = tau_DecayMode->at(tau2Index);
    reco_tau2E = tau_Energy->at(tau2Index);
    plotFill("reco_tau2Pt_" + category, reco_tau2Pt, 80, 0, 400, event_weight);
    plotFill("reco_tau2Eta_" + category, reco_tau2Eta, 25, -2.5, 2.5, event_weight);
    plotFill("reco_tau2Phi_" + category, reco_tau2Phi, 30, -3.14, 3.14, event_weight);
    plotFill("reco_tau2DecayMode_" + category, reco_tau2dm, 12, 0, 12, event_weight);
    reco_lep2.SetPtEtaPhiE(reco_tau2Pt, reco_tau2Eta, reco_tau2Phi, reco_tau2E);
  }
  if (reco_tau_e_idx >= 0 && reco_tau_h_idx >= 0 && gen_tau_h_idx >= 0 && gen_tau_e_idx >= 0)
  {
    /// calcualte and fill deltaR(gen_tau, reco_tau)
    float dR_ee = gen_lep1.DeltaR(reco_lep1);
    float dR_tt = gen_lep2.DeltaR(reco_lep2);
    plotFill("deltaR_ee_" + category, dR_ee, 15, 0, 3, event_weight);
    plotFill("deltaR_tt_" + category, dR_tt, 15, 0, 3, event_weight);

    /// calculate and fill deltaPt (gen_tau , reco_tau)
    float delta_pt = (gen_lep2.Pt() - reco_lep2.Pt());
    float delta_pt_ratio = delta_pt / gen_lep2.Pt();
    plotFill("deltaPt_" + category, delta_pt, 100, -50, 50, event_weight);
    plotFill("deltaPtRatio_" + category, delta_pt_ratio, 40, -2, 2, event_weight);

    plotFill_2D("gen_reco_tau2Pt_"+category, gen_lep2.Pt(), reco_lep2.Pt(),
                  80, 0, 400, 80, 0, 400 );
    plotFill_2D("gen_reco_tau1Pt_"+category, gen_lep1.Pt(), reco_lep1.Pt(),
                  80, 0, 400, 80, 0, 400 );
  }
}
