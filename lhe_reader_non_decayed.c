//&>/dev/null;x="${0%.*}";[ ! "$x" -ot "$0" ]||(rm -f "$x";g++ -o "$x" "$0" -I`root-config --incdir` `root-config --libs`);exit

// Build: g++ lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
//#include </home/muhammad/root6/build/include/TLorentzVector.h>
#include <TLorentzVector.h>
using namespace std;

// Pour une description du format leshouches
// hep-ph/0609017
//
// pour chaque evenement
// une ligne générale : NbPart idprocess poids scale alpha_em alpha_s
// pour chaque particule : id status mere1 mere2 couleur1 couleur2 px py pz E m lifetime spin  
int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;

  TH1::SetDefaultSumw2(true);
// --- insert here histogram definition
  TH1F* hdphi_gj_gl  = new TH1F("dphi_gj_gl","dphi_gj_gl",50,-3.5,3.5) ;
  TH1F* hDPhi_tbarLep_topbJet = new TH1F("dphi_tbarlep_topbJet","dphi tbarlep topbJet",50,-3.5,3.5) ;
  TH1F* hDPhi_topLep_tbarbJet = new TH1F("dphi_topLep_tbarbJet","dphi topLep tbarbJet",50,-3.5,3.5) ;
  TH1F* hDPhi_Jets_Lep  = new TH1F("dPhi_Jets_Lep","dPhi Jets Lep",50,-3.5,3.5) ;
  TH1F* hDPhi_bJets_lJets  = new TH1F("DPhi_bJets_lJets","DPhi bJets lJets",50,-3.5,3.5) ;
  TH1F* hDPhi_topLep_tbarJet  = new TH1F("dPhi_topLep_tbarJet","dPhi topLep tbarJet",50,-3.5,3.5) ;
  TH1F* hDPhi_bJets_bbar  = new TH1F("hDPhi_bJets_bbar","hDPhi bJets bbar",50,-3.5,3.5) ;
  TH1F* hdR_bjet_bbar  = new TH1F("dR_bjet_bbar","dR_bjet_bbar",50,-2.0,10.0) ;
  TH1F* hdR_ljets_lep  = new TH1F("dR_gen_ljets_lep","dR gen ljets lep",50,-2.0,10.0) ;
  TH1F* hdR_b_lep   = new TH1F("dR_gen_b_lep","dR_gen b lep",50,-2.0,10.0) ;

  TH1F* hDPhi_nu_jets = new TH1F("dphi_nu_jets","dphi nu jets",50,-3.5,3.5) ;
  TH1F* hgen_met_pt = new TH1F("MET","MET",50,0,200) ;
  TH1F* hMtt = new TH1F("mtt","tt invariant mass",50,250,900) ;
  TH1F* hPt_ttbar = new TH1F("pt_ttbar","Pt ttbar",50,0,500) ;
  TH1F* hdphi_ttbar = new TH1F("dphi_ttbar","dphi ttbar",50,0,10) ;
  TH1F* hBoost_ttbar = new TH1F("boost_ttbar","Boost ttbar",50,-0.2,1.2);
  TH1F* hPt_top = new TH1F("top_pt","top pt",50,0.0,400) ;
  TH1F* hEta_top = new TH1F("top_eta","top eta",50,-6,6) ;
  TH1F* hPhi_top = new TH1F("top_phi","top phi",50,-3.16,3.16) ;
  TH1F* hPx_top = new TH1F("top_px","top px",50,-200,200) ;
  TH1F* hPy_top = new TH1F("top_py","top py",50,-200,200) ;
  TH2F* hPx_top_tbar = new TH2F("top_tbar_px","top_tbar px",50,-50,50,50,-50,50) ;
  TH2F* hPy_top_tbar = new TH2F("top_tbar_py","top_tbar py",50,-50,50,50,-50,50) ;
  TH1F* hPt_tbar = new TH1F("tbar_pt","tbar pt",50,-10,450) ;
  TH1F* hEta_tbar = new TH1F("tbar_eta","tbar eta",50,-6,6) ;
  TH1F* hPhi_tbar = new TH1F("tbar_phi","tbar phi",50,-3.16,3.16) ;
  TH1F* hPx_tbar = new TH1F("tbar_px","tbar px",50,-200,200) ;
  TH1F* hPy_tbar = new TH1F("tbar_py","tbar py",50,-200,200) ;
  TH1F* hNextrajets = new TH1F("njets","n extra jets",50,0,4);
  TH1F* hextrajet_flavour = new TH1F("jet_flavour","extra jets flavour",50,-5,26);
  TH1F* hPt_extrajet = new TH1F("extrajet_pt","extra jets pt",50,0,200) ;
  TH1F* hEta_extrajet = new TH1F("extrajet_eta","extra jets eta",50,-5,5) ;
  TH1F* hPhi_extrajet = new TH1F("extrajet_phi","extra jets phi",50,-3.15,3.15) ;
  TH1F* hPt_Zp = new TH1F("zprime_pt","Z' pt",50,-10,200) ;
  TH1F* hPt_lmax = new TH1F("lmax_pt","lmax pt",50,-10,300) ;
  TH1F* hEta_lmax = new TH1F("lmax_eta","lmax eta",50,-6,6) ;
  TH1F* hPhi_lmax = new TH1F("lmax_phi","lmax phi",50,-3.16,3.16) ;
  TH1F* hPt_lmin = new TH1F("lmin_pt","lmin pt",50,0,200) ;
  TH1F* hEta_lmin = new TH1F("lmin_eta","lmin eta",50,-6,6) ;
  TH1F* hPhi_lmin = new TH1F("lmin_phi","lmin phi",50,-3.16,3.16) ;
  TH1F* hPt_met = new TH1F("met_pt","met pt",50,0,300) ;
  TH1F* hEta_met = new TH1F("met_eta","met eta",50,-6,6) ;
  TH1F* hPhi_met = new TH1F("met_phi","met phi",50,-3.16,3.16) ;
  TH1F* hPz_met = new TH1F("met_pz","met pz",50,-200,200) ;
  TH1F* hPt_bmax = new TH1F("bmax_pt","bmax pt",50,0,350) ;
  TH1F* hEta_bmax = new TH1F("bmax_eta","bmax eta",50,-6,6) ;
  TH1F* hPhi_bmax = new TH1F("bmax_phi","bmax phi",50,-3.16,3.16) ;
  TH1F* hPt_bmin = new TH1F("bmin_pt","bmin pt",50,0,250) ;
  TH1F* hEta_bmin = new TH1F("bmin_eta","bmin eta",50,-6,6) ;
  TH1F* hPhi_bmin = new TH1F("bmin_phi","bmin phi",50,-3.16,3.16) ;

  TH1F* hPt_met_muon = new TH1F("met_muon_pt","met muon pt",50,0,200) ;
  TH1F* hM_met_muon = new TH1F("met_muon_mass","met muon mass",50,0,200) ;
  TH1F* hM_reco_wpz10 = new TH1F("reco_w_masspz10","reco w mass pz10",200,0,400) ;
  TH1F* hM_reco_wpz20 = new TH1F("reco_w_masspz20","reco w mass pz20",200,0,400) ;
  TH1F* hM_reco_wpz30 = new TH1F("reco_w_masspz30","reco w mass pz30",200,0,400) ;
  TH1F* hM_reco_wpz40 = new TH1F("reco_w_masspz40","reco w mass pz40",200,0,400) ;

  TH1F** hPt_q = new TH1F*[4];
  TH1F** hEta_q = new TH1F*[4];
  TH1F** hPhi_q = new TH1F*[4];
cout<<"I am here: "<<endl;
  for (int iquark=0 ; iquark<4 ; iquark++) {
    ostringstream oss;
    oss << "q" << iquark ;
    string ptstr = oss.str()+"_pt";
    hPt_q[iquark] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_q[iquark] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_q[iquark] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.16,3.16) ;
  }
  TH1F* hCoM_eta_top = new TH1F("top_eta_CoM","top_eta_CoM",50, -5,5);
  TH1F* hCoM_pt_top = new TH1F("top_pt_CoM","top_pt_CoM",50, 0, 500);
  TH1F* hCoM_theta_top = new TH1F("top_theta_CoM","top_theta_CoM",50, 0,3.16);
//  TH1F* hCoM_theta_top = new TH1F("top_theta_CoM","top_theta_CoM",50, -M_PI,M_PI);
// --- end histogram definition

  int nlept=0, nsemi=0, nhadr=0;
  ifstream ff(lhefname.c_str(),ios::in); //ouverture du fichier .lhe
  //ifstream ff("test.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/s1/madevent/Events/zp4000_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/QCD/madevent/Events/qcd_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  int negativeWeight = 0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;

    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;
      event++;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }
      /*weight = 1.;*/

      if (event%1==0) cout << "reading event "<< event << endl;
      int lmin=-1, lmax=-1, bmin=-1, met1=-1, met2=-1, bmax=-1, muon=-1, met=-1, topLep,topbJet, tbarLep, tbarbJet, top_jet, tbar_jet;
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
int n_lep=0, n_alep=0, n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bbar=0, n_tbar_nu=0, n_top_nu=0;
int  bj[2]={-1,-1};
int lp[2]={-1,-1};
      int q[4]={-1,-1,-1,-1};
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_top_lep, v_top_bJet, v_tbar_lep, v_tbar_bJet, v_top_jet, v_tbar_jet, v_tbar_nu, v_top_nu;
      // in lhe first line is number 1, so fill unused array [0] with a crazy value;
      Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      spin[0]= -99999;
     for (int i=1 ; i<npart+1 ; i++) { //start at one
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        if (Status[i]==-1) continue; // status -1 = initial quark ==> skip
        if (Id[i]==6)  top=i;
        if (Id[i]==-6) topbar=i;
        if (Id[i]>6000000) zprime=i;

//--------met for muon----------------------
        int id = abs(Id[i]);  
        Part_absId = abs(Id[i]);
        Moth_absId =abs( Id[Mother1[i]]);
        Part_Id = Id[i];
        Moth_Id = Id[Mother1[i]];
        if ( id==14 ) { met=i;}
//cout<<"this is met"<<Id[i]<<endl;
////----------------------------------------
        if (abs(Id[Mother1[i]])==24 || abs(Id[Mother1[i]])==6) { // mother = W
          int id = abs(Id[i]);
          if ( id==11 || id==13 || id==15 ) { // charged leptons
          if (lmax == -1){ lmax=i;}
          else if (v[i]->Pt() > v[lmax]->Pt()) {
            lmin=lmax; lmax=i;
            }
            else
            {lmin=i;}
//-------------------------only muons and met --------------------------//
          if ( id==13 ) { muon=i;}
//-------------------------------------------------------           
  if(lp[0]==-1) lp[0]=i; 
  else if(lp[1]==-1) lp[1]=i;
  else cout<<"Error: more than two leptons in the event"<<endl;
          }
          if ( id==12 || id==14 || id==16 ) { // neutrinos = MET 
            if (met1 == -1) met1=i;
            else if (met2==-1) met2=i;
            else cout << "ERROR : more than 2 neutrinos" << endl;
          }
          if (id<5) { // light quarks
            //cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " << endl;
            if (q[0] == -1) q[0]=i;
            else if (q[1]==-1) q[1]=i;
            else if (q[2]==-1) q[2]=i;
            else if (q[3]==-1) q[3]=i;
            else cout << "ERROR : more than 4 light quarks" << endl;
          }
        }

          if ( Status[i]==1 && Mother1[i]==1 && Mother2[i]==2 && (abs(Id[i])<6 || Id[i]==21) ) {
            // extra jet
            hextrajet_flavour->Fill(Id[i]+0.5); // in order to center the bin (a quark is in the next bin)
            hPt_extrajet->Fill(v[i]->Pt());
            hEta_extrajet->Fill(v[i]->Eta());
            hPhi_extrajet->Fill(v[i]->Phi());
          }
//---------------------------------------------------------
          else if (abs(Id[i])==5) { // bjets
            if (bmax == -1) bmax=i;
            else if (v[i]->Pt() > v[bmax]->Pt()) {
              bmin=bmax; bmax=i;
            }
            else
              bmin=i;
           if(bj[0]==-1) bj[0]=i;
           else if(bj[1]==-1) bj[1]=i;
           else cout<<"Error: more than two bjets in the event"<<endl;
           }

//------------------top and tbar lep, bJets-----------------------------------
          if (Id[Mother1[i]]==6 || Id[Mother1[i]]==24) { // mother = W+ coming from top
          if ( Id[i]==-11 || Id[i]==-13 || Id[i]==-15 ) { // charged anti leptons
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_top_lep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_alep++;
          }
          if ( Id[i]==12 || Id[i]==14 || Id[i]==16 ) {//neutrino
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_top_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_top_nu++;
          }
          }
          if ( Id[Mother1[i]]==6 || Id[Mother1[i]==24]) { // mother = top
          if (Id[i]==5) { // bjets
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_top_bJet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_topbjets++;
          }
          if ( Id[i]==2 || Id[i]==4 || Id[i]==-1 || Id[i]==-3) {  // jets coming from top
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_top_jet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_topjets++;
          }}
          if ( Id[Mother1[i]]==-6 || Id[Mother1[i]==-24]) {
          if ( Id[i]==11 || Id[i]==13 || Id[i]==15 ) { // charged leptons
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_tbar_lep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_lep++;
          }
          if ( Id[i]==-12 || Id[i]==-14 || Id[i]==-16 ) {// anit neutrino
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_tbar_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_tbar_nu++;
          }
}
          if ( Id[Mother1[i]]==-6 || Id[Mother1[i]]==-24) { // mother = tbar
          if (Id[i]==-5) { // bbarjets
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_tbar_bJet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_bbar++;
          }
          if (Id[i]==1 || Id[i]==3 || Id[i]==-2 || Id[i]==-4) {  // jets coming from tbar
//if (v[i]->Pt() <= 30 || abs(v[i]->Eta()) > 2.5) continue;
          v_tbar_jet.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
          n_tbarjets++;
          }}

} //loop of i

//if((n_topbjets+n_bbar)!=2)continue;
//if((n_topbjets+n_bbar)!=2 || (n_tbarjets+n_topjets)!=2 || (n_lep+n_alep)!=1)continue;//semileptonic conditions
//if((n_topbjets+n_bbar)!=2 || (n_tbarjets+n_topjets)!=4 || (n_lep+n_alep)!=0)continue; // hadronic conditions
//if((n_topbjets+n_bbar)!=2 || (n_tbarjets+n_topjets)!=0 || (n_lep+n_alep)!=2)continue; // hadronic conditions
//cout<<"n_topbjets:  "<<n_topbjets<<"  n_bbar:  "<<n_bbar<<"  n_tbarjets:  "<<n_tbarjets<<"  n_topjets:  "<<n_topjets<<"   n_lep:  "<<n_lep<<"  n_alep: "<<n_alep<<endl;
//if ( n_lep!=1)continue; if(n_jets!=2)continue;
//cout<<"this is no of nu :"<<n_tbar_nu<<"   "<<n_lep<<endl;
int n_bjets = 0;
if(bj[1]!=-1) n_bjets = 2;

int n_leptons = 0;
if(lp[0]!=-1) n_leptons = 1;
if(lp[1]!=-1) n_leptons = 2;
//--------------------------------------
      // sort the quarks
      int nquarks=0;
      if (q[1] != -1) nquarks = 2;
      if (q[3] != -1) nquarks = 4;
      for (int j=0; j<nquarks; j++) {
        double ptmax=-1.; int imax=-1;
        for (int k=j; k<nquarks; k++) {
          if (v[q[k]]->Pt()>ptmax) {
            ptmax = v[q[k]]->Pt();
            imax = k;
          }
        }
        int tmp=q[j];
        q[j]=q[imax];
        q[imax]=tmp;
        }
////-----------------  Jets cuts
          int n_jets_passing_cuts = 0;
          TLorentzVector v_selected_jets;

          for (int nq=0; nq<nquarks ; nq++) {
  //        if (v[q[nq]]->Pt() <= 30 || abs(v[q[nq]]->Eta()) > 2.5) continue;
          n_jets_passing_cuts++;
          v_selected_jets.SetPxPyPzE(v[q[nq]]->Px(),v[q[nq]]->Py(),v[q[nq]]->Pz(),v[q[nq]]->E());
          hPt_q[nq]->Fill( v_selected_jets.Pt(), weight );
          hEta_q[nq]->Fill( v_selected_jets.Eta(), weight );
          hPhi_q[nq]->Fill( v_selected_jets.Phi(), weight );
           }
//if(n_jets_passing_cuts==0)continue;
//--------------------------bjets_max cuts--------------------
      int n_bjmax_passing_cuts = 0;
      TLorentzVector v_selected_bjmax;
      for (int nbj=0; nbj < n_bjets; nbj++ ){
//      if(v[bmax]->Pt() <= 30 || abs(v[bmax]->Eta()) >= 2.5 ) continue;
      n_bjmax_passing_cuts++;
      v_selected_bjmax.SetPxPyPzE(v[bmax]->Px(),v[bmax]->Py(),v[bmax]->Pz(),v[bmax]->E());
      }
//      if(n_bjmax_passing_cuts==0)continue;
//------------------------bjet_min cuts----------------------
      int n_bjmin_passing_cuts = 0;
      TLorentzVector v_selected_bjmin;
      for (int nbj=0; nbj < n_bjets; nbj++ ){
//      if(v[bmin]->Pt() <= 30 || abs(v[bmin]->Eta()) >= 2.5 ) continue;
      n_bjmin_passing_cuts++;
      v_selected_bjmin.SetPxPyPzE(v[bmin]->Px(),v[bmin]->Py(),v[bmin]->Pz(),v[bmin]->E());
      }
//      if(n_bjmin_passing_cuts==0 )continue;
//------------------------Leptons cuts----------------------
      int n_lepmin_passing_cuts = 0;
      TLorentzVector v_selected_lepmin;
      for (int nlp=0; nlp < n_leptons; nlp++ ){
      if(lmin!=-1){
//      if(v[lmin]->Pt() <= 30 || abs(v[lmin]->Eta()) >= 2.5 ) continue;
      v_selected_lepmin.SetPxPyPzE(v[lmin]->Px(),v[lmin]->Py(),v[lmin]->Pz(),v[lmin]->E());
        hPt_lmin->Fill( v_selected_lepmin.Pt(), weight );
        hEta_lmin->Fill( v_selected_lepmin.Eta(), weight );
        hPhi_lmin->Fill( v_selected_lepmin.Phi(), weight );
      }
      n_lepmin_passing_cuts++;
      }
//      if (n_lepmin_passing_cuts == 0 ) continue;
//----------------
      int n_lepmax_passing_cuts = 0;
      TLorentzVector v_selected_lepmax;
      for (int nlp=0; nlp < n_leptons; nlp++ ){
      if(lmax!=-1){
//      if(v[lp[nlp]]->Pt() <= 30 || abs(v[lp[nlp]]->Eta()) >= 2.5 ) continue;
      v_selected_lepmax.SetPxPyPzE(v[lp[nlp]]->Px(),v[lp[nlp]]->Py(),v[lp[nlp]]->Pz(),v[lp[nlp]]->E());
        hPt_lmax->Fill( v_selected_lepmax.Pt(), weight );
        hEta_lmax->Fill( v_selected_lepmax.Eta(), weight );
        hPhi_lmax->Fill( v_selected_lepmax.Phi(), weight );
      }
      n_lepmax_passing_cuts++;
      }
//cout<<"n lepmax  "<<n_lepmax_passing_cuts+ n_lepmin_passing_cuts<<endl;
//if (n_lepmax_passing_cuts == 0 || n_lepmin_passing_cuts==0 || n_jets_passing_cuts==0 ) continue;
//    -------------Top and Tbar seperated--------------------------------------
// cout<<"No of jets:  "<<n_topjets+n_tbarjets<<endl;
//if(n_topbjets!=1 || n_bbar!=1 || (n_tbarjets!=2 && n_topjets!=2) || (n_lep!=1 && n_alep!=1) || ( n_lep + n_tbar_nu)!=2 && (n_alep+n_top_nu)!=2)continue;
//cout<<"n_topbjets:  "<<n_topbjets<<"  n_bbar:  "<<n_bbar<<"  n_tbarjets:  "<<n_tbarjets<<"  n_topjets:  "<<n_topjets<<"   n_lep:  "<<n_lep<<"+"<< n_tbar_nu<<"  n_alep: "<<n_alep<<"+"<<n_top_nu<<endl;

      float top_jets_phi, top_bjets_phi, top_lep_phi, bjets_eta;
      top_jets_phi = v_top_jet.Phi();
      top_bjets_phi = v_top_bJet.Phi();
      top_lep_phi = v_top_lep.Phi();
      bjets_eta = v_top_bJet.Eta();
//      float top_jets_eta=v_top_jet.Eta();
      float top_nu_phi=v_top_nu.Phi();

      float tbar_jets_phi, tbar_bjets_phi, tbar_lep_phi, bbar_eta;
      tbar_jets_phi = v_tbar_jet.Phi();
      tbar_bjets_phi = v_tbar_bJet.Phi();
      tbar_lep_phi = v_tbar_lep.Phi();
      bbar_eta = v_tbar_bJet.Eta();
      float tbar_nu_phi=v_tbar_nu.Phi();
//      float tbar_lep_eta = v_tbar_lep.Eta();
//-------neutrino------------
//if(weight==-1)continue;
float tbar_nu_pt=v_tbar_nu.Pt(), top_nu_pt=v_top_nu.Pt();

      float DPhi_tbarLep_topbJet,DPhi_topLep_tbarbJet, DPhi_Jets_Lep,DPhi_topLep_tbarJet, DPhi_tbarbJet_topbJet, DPhi_topljets_tbarbjets, DPhi_topbjets_tbarljets;
      float DPhi_tbar_nu_top_jets, DPhi_top_nu_tbar_jets, dphi_nu_jets;
      float dphi_jets_lep=0, dphi_bjets_lep=0, dphi_bjets_ljets=0, nu_pt=0;
      DPhi_tbarLep_topbJet=DeltaPhi(tbar_lep_phi,top_bjets_phi);
      DPhi_topLep_tbarbJet=DeltaPhi(top_lep_phi,tbar_bjets_phi);
      DPhi_Jets_Lep=DeltaPhi(top_jets_phi,tbar_lep_phi);
      DPhi_topLep_tbarJet=DeltaPhi(top_lep_phi,tbar_jets_phi);
      DPhi_tbarbJet_topbJet=DeltaPhi(top_bjets_phi,tbar_bjets_phi);
      DPhi_topljets_tbarbjets=DeltaPhi(top_jets_phi,tbar_bjets_phi);
      DPhi_topbjets_tbarljets=DeltaPhi(top_bjets_phi, tbar_jets_phi);

      DPhi_tbar_nu_top_jets=DeltaPhi(tbar_nu_phi,top_jets_phi);
      DPhi_top_nu_tbar_jets=DeltaPhi(top_nu_phi,tbar_jets_phi);
      if(n_lep==0 && n_topjets==0 && n_alep==1 && n_tbarjets==2)
        {dphi_jets_lep=DPhi_topLep_tbarJet;dphi_bjets_lep=DPhi_topLep_tbarbJet;
        dphi_bjets_ljets=DPhi_topbjets_tbarljets; dphi_nu_jets=DPhi_top_nu_tbar_jets;nu_pt=top_nu_pt;}
      if(n_lep==1 && n_topjets==2 && n_alep==0 && n_tbarjets==0)
        {dphi_jets_lep=DPhi_Jets_Lep;dphi_bjets_lep=DPhi_tbarLep_topbJet; dphi_bjets_ljets=DPhi_topljets_tbarbjets;
        dphi_nu_jets=DPhi_tbar_nu_top_jets;nu_pt=tbar_nu_pt;}


      hDPhi_tbarLep_topbJet->Fill(dphi_bjets_lep, weight);
      hDPhi_topLep_tbarbJet->Fill(dphi_bjets_lep, weight);
      hDPhi_Jets_Lep->Fill(dphi_jets_lep, weight);
      hdphi_gj_gl->Fill(DPhi_Jets_Lep, weight);//just for test
      hDPhi_topLep_tbarJet->Fill(dphi_jets_lep, weight);
      hDPhi_bJets_bbar->Fill(DPhi_tbarbJet_topbJet,weight);
      hDPhi_bJets_lJets->Fill(dphi_bjets_ljets,weight);

      hDPhi_nu_jets->Fill(dphi_nu_jets,weight);
      hgen_met_pt->Fill(nu_pt,weight);
/*      float dR_bjet_bbar, dR_ljets_lep, dR_b_lep;
      dR_bjet_bbar= DeltaR( bjets_eta, bbar_eta, top_bjets_phi, tbar_bjets_phi);
      dR_ljets_lep= DeltaR(top_jets_eta,tbar_lep_eta, top_jets_phi, tbar_lep_phi);
      dR_b_lep= DeltaR(bjets_eta,tbar_lep_eta,top_bjets_phi,tbar_lep_phi);
      hdR_bjet_bbar->Fill(dR_bjet_bbar,weight);
      hdR_ljets_lep->Fill(dR_ljets_lep,weight);
      hdR_b_lep->Fill(dR_b_lep,weight);*/
//------------------------------------------------------
       // bmax kinematics
        hPt_bmax->Fill(v_selected_bjmax.Pt(), weight );
        hEta_bmax->Fill(v_selected_bjmax.Eta(), weight );
        hPhi_bmax->Fill(v_selected_bjmax.Phi(), weight );
      // bmin kinematics
        hPt_bmin->Fill(v_selected_bjmin.Pt(), weight );
        hEta_bmin->Fill(v_selected_bjmin.Eta(), weight );
        hPhi_bmin->Fill(v_selected_bjmin.Phi(), weight );
//------------------------------------------------------
      if (top == -1 || topbar == -1) {
        std::cout << "Warning: no tt~ in this event (line " << line << ")" << std::endl;
        std::cout << buffer.str() << std::endl;
        continue;
      }

// --- insert here the code to fill the histograms
      // top kinematics
      float top_mass = v[top]->M();
      hPt_top->Fill( v[top]->Pt(), weight);
      hEta_top->Fill( v[top]->Eta(), weight);
      hPhi_top->Fill( v[top]->Phi(), weight);
      hPx_top->Fill( px[top], weight);
      hPy_top->Fill( py[top], weight);
      hPx_top_tbar->Fill( px[top],px[topbar], weight);
      hPy_top_tbar->Fill( py[top],py[topbar], weight);
      // tbar kinematics
      float tbar_mass = v[topbar]->M();
      hPt_tbar->Fill( v[topbar]->Pt(), weight);
      hEta_tbar->Fill( v[topbar]->Eta(), weight);
      hPhi_tbar->Fill( v[topbar]->Phi(), weight );
      hPx_tbar->Fill( px[topbar], weight);
      hPy_tbar->Fill( py[topbar], weight);
//std::cout<<"top and tbar masses :"<<top_mass<<"  :"<<tbar_mass<<endl;
      // Z' kinematics
      if (zprime != -1)
        hPt_Zp->Fill( v[zprime]->Pt(), weight );
      //hEta_Zp->Fill( v[zprime]->Eta() );
      //hPhi_Zp->Fill( v[zprime]->Phi() );
      // mtt
      double mtt = sqrt( m[top]*m[top] + m[topbar]*m[topbar]
              + 2*(E[top]*E[topbar] - px[top]*px[topbar] - py[top]*py[topbar] - pz[top]*pz[topbar]) );

      /*if (mtt > 500)*/
        /*weight = -1;*/
      /*else*/
        /*weight = 1;*/

      hMtt->Fill(mtt, weight);
      // boost ttbar system
      hBoost_ttbar->Fill((*(v[top])+*(v[topbar])).Beta(), weight);

      // go back to CoM
      TLorentzVector ttbar_CoM =  *v[top] + *v[topbar];
      //cout << ttbar_CoM.E()<< endl;
      TVector3 b(0,0,-1*ttbar_CoM.Pz()/ttbar_CoM.E());
      ttbar_CoM.Boost(b);
      //cout << ttbar_CoM.Pz() << endl;
      TLorentzVector vtop_com = *v[top];
      vtop_com.Boost(b);
      hCoM_eta_top->Fill(vtop_com.Eta(), weight);
      hCoM_pt_top->Fill(vtop_com.Pt(), weight);
      hCoM_theta_top->Fill(vtop_com.Theta(), weight);

      // ttbar system
//      double pt_Ttbar = sqrt(abs(px[top]*px[topbar]+py[top]*py[topbar]));
double pt_Ttbar = sqrt(pow(px[top]+px[topbar],2)+pow(py[top]+py[topbar],2));
      hPt_ttbar->Fill(pt_Ttbar, weight);
      float phi_top = v[top]->Phi(), phi_tbar = v[topbar]->Phi();

      float DPhi_ttbar;
      DPhi_ttbar=DeltaPhi(phi_top, phi_tbar);
      hdphi_ttbar->Fill(DPhi_ttbar);
      // n extra jets
      hNextrajets->Fill( (npart - 13)+0.5 ) ; //for decayed ttbar*/
      //hNextrajets->Fill( (npart - 5)+0.5 )  ; // for undecayed ttbar

      // lmax kinematics
      if (lmax != -1) {
        hPt_lmax->Fill( v[lmax]->Pt(), weight );
        hEta_lmax->Fill( v[lmax]->Eta(), weight );
        hPhi_lmax->Fill( v[lmax]->Phi(), weight );
      }
      // lmin kinematics
      if (lmin != -1) {
        hPt_lmin->Fill( v[lmin]->Pt() );
        hEta_lmin->Fill( v[lmin]->Eta() );
        hPhi_lmin->Fill( v[lmin]->Phi() );
      }
      // met kinematics
TLorentzVector  v_reco_wpz10,  v_reco_wpz20,  v_reco_wpz30,  v_reco_wpz40;
      float met_px = 0, met_py = 0,met_pz = 0, met_E = 0;
      if (met1 !=-1) {
        if (met2 !=-1 ) {
          hPt_met->Fill( (*v[met1]+*v[met2]).Pt() );
          hEta_met->Fill( (*v[met1]+*v[met2]).Eta() );
          hPhi_met->Fill( (*v[met1]+*v[met2]).Phi() );
          hPz_met->Fill( (*v[met1]+*v[met2]).Pz() );
       met_pz = (*v[met1]+*v[met2]).Pz();
        }
        else {
          hPt_met->Fill( v[met1]->Pt() );
          hEta_met->Fill( v[met1]->Eta() );
          hPhi_met->Fill( v[met1]->Phi() );
          hPz_met->Fill(  v[met1]->Pz() );
       met_pz = v[met1]->Pz();
        }
      }
//v_selected_jets.SetPxPyPzE(v[q[nq]]->Px(),v[q[nq]]->Py(),v[q[nq]]->Pz(),v[q[nq]]->E());
//----------------------------Reco W boson -----------------------------//

     if (met !=-1 && muon !=-1) {
         met_px = v[met]->Px(); met_py = v[met]->Py(); met_pz = v[met]->Pz(); met_E = v[met]->E();
 v_reco_wpz10.SetPxPyPzE(met_px,met_py,10.0,met_E);
 v_reco_wpz20.SetPxPyPzE(met_px,met_py,20.0,met_E);
 v_reco_wpz30.SetPxPyPzE(met_px,met_py,30.0,met_E);
 v_reco_wpz40.SetPxPyPzE(met_px,met_py,40.0,met_E);

         hPt_met_muon->Fill( (*v[met]+*v[muon]).Pt() );
         hM_met_muon->Fill( (*v[met]+*v[muon]).M() );
  hM_reco_wpz10->Fill( ( v_reco_wpz10+*v[muon]).M() );
  hM_reco_wpz20->Fill( ( v_reco_wpz20+*v[muon]).M() );
  hM_reco_wpz30->Fill( ( v_reco_wpz30+*v[muon]).M() );
  hM_reco_wpz40->Fill( ( v_reco_wpz40+*v[muon]).M() );

         }
//  cout<<"this is met pz"<<met_px<<endl;
/*      // bmax kinematics
      if (bmax != -1) {
        hPt_bmax->Fill( v[bmax]->Pt() );
        hEta_bmax->Fill( v[bmax]->Eta() );
        hPhi_bmax->Fill( v[bmax]->Phi() );
      }*/
      // bmin kinematics
      if (bmin != -1) {
//        hPt_bmin->Fill( v[bmin]->Pt() );
//        hEta_bmin->Fill( v[bmin]->Eta() );
//        hPhi_bmin->Fill( v[bmin]->Phi() );
      }
      // light quarks kinematics
      for (int nq=0; nq<nquarks ; nq++) {
        hPt_q[nq]->Fill( v[q[nq]]->Pt() );
        hEta_q[nq]->Fill( v[q[nq]]->Eta() );
        hPhi_q[nq]->Fill( v[q[nq]]->Phi() );
      }
      if (lmin !=-1) nlept++;
      else if (nquarks==4) nhadr++;
      else nsemi++;
// --- end filling the histograms
      ff>>tt;
      line++;
      //if (event==100)  break;
      delete Id;
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);



    }

}
  cout << " Total number of events --> " << event << endl;
  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");
  hDPhi_nu_jets->Write();
  hgen_met_pt->Write();

  hdphi_gj_gl->Write();
  hdR_bjet_bbar->Write();
  hdR_ljets_lep->Write();
  hdR_b_lep->Write();
  hDPhi_tbarLep_topbJet->Write();
  hDPhi_Jets_Lep->Write();
  hDPhi_bJets_lJets->Write();
  hDPhi_bJets_bbar->Write();
  hDPhi_topLep_tbarbJet->Write();
  hDPhi_topLep_tbarJet->Write();
  hMtt->Write();
  hPt_ttbar->Write();
  hdphi_ttbar->Write();
  hPt_top->Write();
  hEta_top->Write();
  hPhi_top->Write();
  hPx_top->Write();
  hPy_top->Write();
  hPx_top_tbar->Write();
  hPy_top_tbar->Write();
  hPt_tbar->Write();
  hEta_tbar->Write();
  hPhi_tbar->Write();
  hPx_tbar->Write();
  hPy_tbar->Write();
  hPt_Zp->Write();
  //hEta_Zp->Write();
  //hPhi_Zp->Write();
  hNextrajets->Write();
  hPt_lmax->Write();
  hEta_lmax->Write();
  hPhi_lmax->Write();
  hPt_lmin->Write();
  hEta_lmin->Write();
  hPhi_lmin->Write();
  hPt_bmax->Write();
  hEta_bmax->Write();
  hPhi_bmax->Write();
  hPt_bmin->Write();
  hEta_bmin->Write();
  hPhi_bmin->Write();
  hPt_met->Write();
  hEta_met->Write();
  hPhi_met->Write();
  hPz_met->Write();
  hPt_met_muon->Write();
  hM_met_muon->Write();
  hM_reco_wpz10->Write();
  hM_reco_wpz20->Write();
  hM_reco_wpz30->Write();
  hM_reco_wpz40->Write();

  for (int i=0;i<4;i++) {
    hPt_q[i]->Write();
    hEta_q[i]->Write();
    hPhi_q[i]->Write();
  }
  hextrajet_flavour->Write();
  hPt_extrajet->Write();
  hEta_extrajet->Write();
  hPhi_extrajet->Write();
  hBoost_ttbar->Write();
  hCoM_eta_top->Write();
  hCoM_theta_top->Write();
  hCoM_pt_top->Write();
  rootfile->Close();

  cout << "Events with negative weight: " << negativeWeight << endl;
  cout << "lept decay = " << nlept*1.0/event << endl;
  cout << "hadr decay = " << nhadr*1.0/event << endl;
  cout << "semi decay = " << nsemi*1.0/event << endl;
  exit(0);




}
float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}



