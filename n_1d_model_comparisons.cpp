/*
 * A root macro to plot various parameters for SBND-like events
 * The events are generated using the MiniBooNE flux 
 *
 * The plots will be normalised 
 * The normalisation will account for the fact that 
 * SBND will be closer than MiniBooNE to the target hall
 *
 * The normalisation is calculated using:
 *
 *
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : February 2017
 *
 * For now, this is a simple macro for neutrino-mode only
 * I will amend this entire system to be a class structure with 
 * constructors and destructors for neutrino and antineutrino-modes
 *
*/

#include "model_comparisons.h"

using namespace std; 

int n_1d_model_comparisons() {

    // -------------------------------------------------------------------------
    // Normalisation:
    // For each model configuration, the distributions need to be normalised
    // This is done using:
    //     N = ( sum_i(flux_i*x-sec_i) / sum_i(flux_i) ) * (targetdist_MB / targetdist_SBND)^2  
    //
    // The flux_i will be the same for all
    // The cross sections will vary by model configuration
    // The ratio 
    //      (targetdist_MB / targetdist_SBND)^2 = (450 / 110 )^2 = 16.7 
    //      this accounts for the different distances of MiniBooNE and SBND from the target
    //      when using the MiniBooNE flux
    //
    // Read in the files
    // File 1:
    //      Model     : Default
    //      Nu/NuBar  : Nu 
    //      Experiment: SBND
    //
    // File 2:
    //      Model     : Default+MEC
    //      Nu/NuBar  : Nu 
    //      Experiment: SBND
    //
    // File 3:
    //      Model     : G16_02b
    //      Nu/NuBar  : Nu 
    //      Experiment: SBND
    // -------------------------------------------------------------------------   
    // -------------------------------------------------------------------------
    //                              Open flux file
    // -------------------------------------------------------------------------
    
    TFile f_flux("/hepstore/rjones/Software/Genie/comparisons/builds/SL6.x86_64/vtrunk/data/fluxes/miniboone/miniboone_april07_baseline_rgen610.6_flux_pospolarity_fluxes.root");
    if (f_flux.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " MiniBooNE flux file is open " << endl;
    }
    
    // -------------------------------------------------------------------------
    //                         Open the cross section files
    // -------------------------------------------------------------------------
    
    // Default
    TFile f1_xsec("/hepstore/rjones/Exercises/Flavours/Default/xsec_files/total_xsec.root");
    if (f1_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default xsec file is open " << endl;
    }
    
    // Default+MEC
    TFile f2_xsec("/hepstore/rjones/Exercises/Flavours/Default+MEC/xsec_files/total_xsec.root");
    if (f2_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default+MEC xsec file is open " << endl;
    }
    
    // G16_02b
    TFile f3_xsec("/hepstore/rjones/Exercises/Flavours/G16_02b/xsec_files/total_xsec.root");
    if (f3_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_02b xsec file is open " << endl;
    }
    
    // G16_01a
    TFile f4_xsec("/hepstore/rjones/Exercises/Flavours/G16_01a/xsec_files/total_xsec.root");
    if (f4_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01a xsec file is open " << endl;
    }
    
    // G16_01b
    TFile f5_xsec("/hepstore/rjones/Exercises/Flavours/G16_01b/xsec_files/total_xsec.root");
    if (f5_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01b xsec file is open " << endl;
    }
    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open Default
    TFile f1("/hepstore/rjones/Exercises/Flavours/Default/sbnd/1M/gntp.10000.gst.root");
    if (f1.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default event file 1 is open " << endl;
    }

    // Open Default+MEC
    TFile f2("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if (f2.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default+MEC event file 1 is open " << endl;
    }

    // Open G16_02b
    TFile f3("/hepstore/rjones/Exercises/Flavours/G16_02b/sbnd/1M/gntp.10000.gst.root");
    if (f3.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_02b event file 1 is open " << endl;
    }

    // Open G16_01a
    TFile f4("/hepstore/rjones/Exercises/Flavours/G16_01a/sbnd/1M/gntp.10000.gst.root");
    if (f4.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01a event file 1 is open " << endl;
    }

    // Open G16_01b
    TFile f5("/hepstore/rjones/Exercises/Flavours/G16_01b/sbnd/1M/gntp.10000.gst.root");
    if (f5.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01b event file 1 is open " << endl;
    }

    // -------------------------------------------------------------------------
    //                     Calculate the reco energy difference
    // -------------------------------------------------------------------------
    
    // Get the trees we want from the root files
    TTree *gst1 = (TTree*) f1.Get("gst");
    TTree *gst2 = (TTree*) f2.Get("gst");
    TTree *gst3 = (TTree*) f3.Get("gst");
    TTree *gst4 = (TTree*) f4.Get("gst");
    TTree *gst5 = (TTree*) f5.Get("gst");
   
    // -------------------------------------------------------------------------
    //                          Get the normalisations
    // -------------------------------------------------------------------------
    
    vector< double > norms;
    norms.push_back( Norm(1000000, f1_xsec, f_flux) );
    norms.push_back( Norm(1000000, f2_xsec, f_flux) );
    norms.push_back( Norm(1000000, f3_xsec, f_flux) );
    norms.push_back( Norm(1000000, f4_xsec, f_flux) );
    norms.push_back( Norm(1000000, f5_xsec, f_flux) );
    
    // -------------------------------------------------------------------------
    //                          Canvases and histograms
    // -------------------------------------------------------------------------
    // Legend entries will be the same for all stacked histograms
    vector< string > leg_entries;
    
    leg_entries.push_back(" Default " );
    leg_entries.push_back(" Default + MEC " );
    leg_entries.push_back(" G17_02b " );
    leg_entries.push_back(" G17_01a " );
    leg_entries.push_back(" G17_01b " );

    // -------------------------------------------------------------------------
    //                     CC and NC 0pi numerical distributions
    //                              CC0pi NN 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > NN_abs_cc0pi_stack;

    // The NN histograms
    TH1D *h1_NN_abs_cc0pi = new TH1D("h1_NN_abs_cc0pi","N_{N} , CC0#pi",25,0,25);
    TH1D *h2_NN_abs_cc0pi = new TH1D("h2_NN_abs_cc0pi","N_{N} , CC0#pi",25,0,25);
    TH1D *h3_NN_abs_cc0pi = new TH1D("h3_NN_abs_cc0pi","N_{N} , CC0#pi",25,0,25);
    TH1D *h4_NN_abs_cc0pi = new TH1D("h4_NN_abs_cc0pi","N_{N} , CC0#pi",25,0,25);
    TH1D *h5_NN_abs_cc0pi = new TH1D("h5_NN_abs_cc0pi","N_{N} , CC0#pi",25,0,25);

    gst1->Draw("(nfp+nfn)>>h1_NN_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("(nfp+nfn)>>h2_NN_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("(nfp+nfn)>>h3_NN_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("(nfp+nfn)>>h4_NN_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("(nfp+nfn)>>h5_NN_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    NN_abs_cc0pi_stack.push_back(h1_NN_abs_cc0pi);
    NN_abs_cc0pi_stack.push_back(h2_NN_abs_cc0pi);
    NN_abs_cc0pi_stack.push_back(h3_NN_abs_cc0pi);
    NN_abs_cc0pi_stack.push_back(h4_NN_abs_cc0pi);
    NN_abs_cc0pi_stack.push_back(h5_NN_abs_cc0pi);

    HistStacker(NN_abs_cc0pi_stack, leg_entries, norms, "FS Nucleons, CC0#pi", "NN_CC0Pi_Plot.png", "N_{N}", "Number of SBND events");

    delete h1_NN_abs_cc0pi;
    delete h2_NN_abs_cc0pi;
    delete h3_NN_abs_cc0pi;
    delete h4_NN_abs_cc0pi;
    delete h5_NN_abs_cc0pi;

    // -------------------------------------------------------------------------
    //                              NC0pi NN 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > NN_abs_nc0pi_stack;

    // The NN histograms
    TH1D *h1_NN_abs_nc0pi = new TH1D("h1_NN_abs_nc0pi","N_{N} , NC0#pi",25,0,25);
    TH1D *h2_NN_abs_nc0pi = new TH1D("h2_NN_abs_nc0pi","N_{N} , NC0#pi",25,0,25);
    TH1D *h3_NN_abs_nc0pi = new TH1D("h3_NN_abs_nc0pi","N_{N} , NC0#pi",25,0,25);
    TH1D *h4_NN_abs_nc0pi = new TH1D("h4_NN_abs_nc0pi","N_{N} , NC0#pi",25,0,25);
    TH1D *h5_NN_abs_nc0pi = new TH1D("h5_NN_abs_nc0pi","N_{N} , NC0#pi",25,0,25);

    gst1->Draw("(nfp+nfn)>>h1_NN_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("(nfp+nfn)>>h2_NN_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("(nfp+nfn)>>h3_NN_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("(nfp+nfn)>>h4_NN_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("(nfp+nfn)>>h5_NN_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    NN_abs_nc0pi_stack.push_back(h1_NN_abs_nc0pi);
    NN_abs_nc0pi_stack.push_back(h2_NN_abs_nc0pi);
    NN_abs_nc0pi_stack.push_back(h3_NN_abs_nc0pi);
    NN_abs_nc0pi_stack.push_back(h4_NN_abs_nc0pi);
    NN_abs_nc0pi_stack.push_back(h5_NN_abs_nc0pi);

    HistStacker(NN_abs_nc0pi_stack, leg_entries, norms, "FS Nucleons, NC0#pi", "NN_NC0Pi_Plot.png", "N_{N}", "Number of SBND events");

    delete h1_NN_abs_nc0pi;
    delete h2_NN_abs_nc0pi;
    delete h3_NN_abs_nc0pi;
    delete h4_NN_abs_nc0pi;
    delete h5_NN_abs_nc0pi;
    
    // -------------------------------------------------------------------------
    //                              CC0pi Np 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Np_abs_cc0pi_stack;

    // The Np histograms
    TH1D *h1_Np_abs_cc0pi = new TH1D("h1_Np_abs_cc0pi","N_{p} , CC0#pi",15,0,15);
    TH1D *h2_Np_abs_cc0pi = new TH1D("h2_Np_abs_cc0pi","N_{p} , CC0#pi",15,0,15);
    TH1D *h3_Np_abs_cc0pi = new TH1D("h3_Np_abs_cc0pi","N_{p} , CC0#pi",15,0,15);
    TH1D *h4_Np_abs_cc0pi = new TH1D("h4_Np_abs_cc0pi","N_{p} , CC0#pi",15,0,15);
    TH1D *h5_Np_abs_cc0pi = new TH1D("h5_Np_abs_cc0pi","N_{p} , CC0#pi",15,0,15);

    gst1->Draw("nfp>>h1_Np_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("nfp>>h2_Np_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("nfp>>h3_Np_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("nfp>>h4_Np_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("nfp>>h5_Np_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Np_abs_cc0pi_stack.push_back(h1_Np_abs_cc0pi);
    Np_abs_cc0pi_stack.push_back(h2_Np_abs_cc0pi);
    Np_abs_cc0pi_stack.push_back(h3_Np_abs_cc0pi);
    Np_abs_cc0pi_stack.push_back(h4_Np_abs_cc0pi);
    Np_abs_cc0pi_stack.push_back(h5_Np_abs_cc0pi);

    HistStacker(Np_abs_cc0pi_stack, leg_entries, norms, "FS Protons, CC0#pi", "Np_CC0Pi_Plot.png", "N_{p}", "Number of SBND events");

    delete h1_Np_abs_cc0pi;
    delete h2_Np_abs_cc0pi;
    delete h3_Np_abs_cc0pi;
    delete h4_Np_abs_cc0pi;
    delete h5_Np_abs_cc0pi;

    // -------------------------------------------------------------------------
    //                              NC0pi Np 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Np_abs_nc0pi_stack;

    // The Np histograms
    TH1D *h1_Np_abs_nc0pi = new TH1D("h1_Np_abs_nc0pi","N_{p} , NC0#pi",15,0,15);
    TH1D *h2_Np_abs_nc0pi = new TH1D("h2_Np_abs_nc0pi","N_{p} , NC0#pi",15,0,15);
    TH1D *h3_Np_abs_nc0pi = new TH1D("h3_Np_abs_nc0pi","N_{p} , NC0#pi",15,0,15);
    TH1D *h4_Np_abs_nc0pi = new TH1D("h4_Np_abs_nc0pi","N_{p} , NC0#pi",15,0,15);
    TH1D *h5_Np_abs_nc0pi = new TH1D("h5_Np_abs_nc0pi","N_{p} , NC0#pi",15,0,15);

    gst1->Draw("nfp>>h1_Np_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("nfp>>h2_Np_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("nfp>>h3_Np_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("nfp>>h4_Np_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("nfp>>h5_Np_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Np_abs_nc0pi_stack.push_back(h1_Np_abs_nc0pi);
    Np_abs_nc0pi_stack.push_back(h2_Np_abs_nc0pi);
    Np_abs_nc0pi_stack.push_back(h3_Np_abs_nc0pi);
    Np_abs_nc0pi_stack.push_back(h4_Np_abs_nc0pi);
    Np_abs_nc0pi_stack.push_back(h5_Np_abs_nc0pi);

    HistStacker(Np_abs_nc0pi_stack, leg_entries, norms, "FS Protons, NC0#pi", "Np_NC0Pi_Plot.png", "N_{p}", "Number of SBND events");

    delete h1_Np_abs_nc0pi;
    delete h2_Np_abs_nc0pi;
    delete h3_Np_abs_nc0pi;
    delete h4_Np_abs_nc0pi;
    delete h5_Np_abs_nc0pi;

    // -------------------------------------------------------------------------
    //                              CC0pi Nn 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Nn_abs_cc0pi_stack;

    // The Nn histograms
    TH1D *h1_Nn_abs_cc0pi = new TH1D("h1_Nn_abs_cc0pi","N_{n} , CC0#pi",15,0,15);
    TH1D *h2_Nn_abs_cc0pi = new TH1D("h2_Nn_abs_cc0pi","N_{n} , CC0#pi",15,0,15);
    TH1D *h3_Nn_abs_cc0pi = new TH1D("h3_Nn_abs_cc0pi","N_{n} , CC0#pi",15,0,15);
    TH1D *h4_Nn_abs_cc0pi = new TH1D("h4_Nn_abs_cc0pi","N_{n} , CC0#pi",15,0,15);
    TH1D *h5_Nn_abs_cc0pi = new TH1D("h5_Nn_abs_cc0pi","N_{n} , CC0#pi",15,0,15);

    gst1->Draw("nfn>>h1_Nn_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("nfn>>h2_Nn_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("nfn>>h3_Nn_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("nfn>>h4_Nn_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("nfn>>h5_Nn_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Nn_abs_cc0pi_stack.push_back(h1_Nn_abs_cc0pi);
    Nn_abs_cc0pi_stack.push_back(h2_Nn_abs_cc0pi);
    Nn_abs_cc0pi_stack.push_back(h3_Nn_abs_cc0pi);
    Nn_abs_cc0pi_stack.push_back(h4_Nn_abs_cc0pi);
    Nn_abs_cc0pi_stack.push_back(h5_Nn_abs_cc0pi);

    HistStacker(Nn_abs_cc0pi_stack, leg_entries, norms, "FS Neutrons, CC0#pi", "Nn_CC0Pi_Plot.png", "N_{n}", "Number of SBND events");

    delete h1_Nn_abs_cc0pi;
    delete h2_Nn_abs_cc0pi;
    delete h3_Nn_abs_cc0pi;
    delete h4_Nn_abs_cc0pi;
    delete h5_Nn_abs_cc0pi;

    // -------------------------------------------------------------------------
    //                              NC0pi Nn 
    // -------------------------------------------------------------------------

    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Nn_abs_nc0pi_stack;

    // The Nn histograms
    TH1D *h1_Nn_abs_nc0pi = new TH1D("h1_Nn_abs_nc0pi","N_{n} , NC0#pi",15,0,15);
    TH1D *h2_Nn_abs_nc0pi = new TH1D("h2_Nn_abs_nc0pi","N_{n} , NC0#pi",15,0,15);
    TH1D *h3_Nn_abs_nc0pi = new TH1D("h3_Nn_abs_nc0pi","N_{n} , NC0#pi",15,0,15);
    TH1D *h4_Nn_abs_nc0pi = new TH1D("h4_Nn_abs_nc0pi","N_{n} , NC0#pi",15,0,15);
    TH1D *h5_Nn_abs_nc0pi = new TH1D("h5_Nn_abs_nc0pi","N_{n} , NC0#pi",15,0,15);

    gst1->Draw("nfn>>h1_Nn_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("nfn>>h2_Nn_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("nfn>>h3_Nn_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("nfn>>h4_Nn_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("nfn>>h5_Nn_abs_nc0pi", "nc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Nn_abs_nc0pi_stack.push_back(h1_Nn_abs_nc0pi);
    Nn_abs_nc0pi_stack.push_back(h2_Nn_abs_nc0pi);
    Nn_abs_nc0pi_stack.push_back(h3_Nn_abs_nc0pi);
    Nn_abs_nc0pi_stack.push_back(h4_Nn_abs_nc0pi);
    Nn_abs_nc0pi_stack.push_back(h5_Nn_abs_nc0pi);

    HistStacker(Nn_abs_nc0pi_stack, leg_entries, norms, "FS Neutrons, NC0#pi", "Nn_NC0Pi_Plot.png", "N_{n}", "Number of SBND events");

    delete h1_Nn_abs_nc0pi;
    delete h2_Nn_abs_nc0pi;
    delete h3_Nn_abs_nc0pi;
    delete h4_Nn_abs_nc0pi;
    delete h5_Nn_abs_nc0pi;

    return 0;
}

// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------
//
//                                                DEFINING THE FUNCTIONS
//
// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------
//                      normalisation function
// -------------------------------------------------------------------------
double Norm(int n_events, TFile &xsec_file, TFile &flux_file){
    // The constant values used throughout
    double sbnd_scale = 16.74;       // The ratio of MB:SBND distances^2 
    double Na         = 6.022e23;    // Avogadro
    double A_Ar       = 0.04;        // Argon mass # kg/mol
    double rho_Ar     = 1390;        // kg/m^3
    double M_fid      = rho_Ar * 55; // 55 m^3 for nu_mu, should be ~ 77000
    double xsec_scale = 1e-38;       // cm^2
    double e_bins     = 0.05;        // GeV
    double POT_sbnd   = 6.6e20;      // POT
    double n_sbnd;
    double norm_sbnd;

    // Get the flux histogram
    TH1D *h_flux = (TH1D*) flux_file.Get("flux_pos_pol_numu");
      
    // Loop over the entries of the histogram, multiply all by 0.05 MeV for the
    // binning
    int n_bins = h_flux->GetNbinsX(); 
    double tot_flux = 0; // Initiate cumulative flux
    double e_x_flux;     // Bin width * total flux
    
    for (int i = 0; i < n_bins; ++i){
        // For the current bin, get the content and add it to the current total
        tot_flux += h_flux->GetBinContent(i);
    }
    
    e_x_flux = tot_flux * e_bins;
    
    // Get the x_sec graphs
    TGraph *cc = (TGraph*) xsec_file.Get("nu_mu_Ar40/tot_cc");
    TGraph *nc = (TGraph*) xsec_file.Get("nu_mu_Ar40/tot_nc");
        
    // Integrate over the TGraphs
    double tot_xsec_cc = cc->Integral(0,-1);
    double tot_xsec_nc = nc->Integral(0,-1);
    double tot_xsec    = tot_xsec_cc + tot_xsec_nc;
    
    // Calculate the normalisation and return the value
    n_sbnd = e_x_flux * 16.74 * tot_xsec * xsec_scale * POT_sbnd * M_fid * Na * (1 / A_Ar );
    norm_sbnd = (n_sbnd / n_events);
    return norm_sbnd;
}

// -------------------------------------------------------------------------
//                      hist stacking function
// -------------------------------------------------------------------------
void HistStacker ( vector< TH1D* >   &hists,
                   vector< string >  &leg_entries,
                   vector< double >  &norm, 
                   const char* title,
                   const char* file_name,
                   const char* x_axis,
                   const char* y_axis ){
    
    // The Canvas, legend and empty histogram to print the title and axes labels
    // Canvas
    TCanvas *c   = new TCanvas( "c", title, 800, 600 );
    TLegend *leg = new TLegend( 0.68, 0.68, 0.88, 0.88 );

    // Get the sizes of the vectors
    int n_hists, n_leg, n_norm;
    n_hists = hists.size();
    n_leg   = leg_entries.size();
    n_norm  = norm.size();

    // Check they are the same, or return an error
    if ( n_hists != n_leg || n_hists != n_norm || n_leg != n_norm ) {
        cerr << " The vectors should all have the same number of entries " << endl;
        exit(1);
    }
    
    // Loop over the histograms 
    for ( int i = 0; i < n_hists; ++i ) {
    // Set their line colours, styles and widths
         if ( i < 4 || i == 5 ) {
             hists[i]->SetLineColor( i + 1 );
             hists[i]->SetLineStyle( 1 );
             hists[i]->SetLineWidth( 1 );
             hists[i]->SetMarkerColor( i + 1 );
         }
         else if ( i == 4 || i == 6 || i == 7 ){
             hists[i]->SetLineColor( i + 3 );
             hists[i]->SetLineStyle( 1 );
             hists[i]->SetLineWidth( 1 );
             hists[i]->SetMarkerColor( i + 3 );
         }
 
        // Normalise the histogram with the associated normalisation
        hists[i]->Scale( norm[i] );
    }

    for ( int i = 0; i < n_hists; ++i ) {

        // Don't draw G17_01b for now 
        if( i != 4 ){
            // Add in the legend entry
            leg->AddEntry( hists[i], leg_entries[i].c_str(), "l" );
        }
    }

    // Variables for calculating the maximum y value
    int i_max;
    double max = -1000; 
    double max_y;

    // Loop over the histograms and find which has the highest y-val
    // Print this first
    for ( int i = 0; i < n_hists; ++i ) {
        if ( hists[i]->GetMaximum() > max){
            max = hists[i]->GetMaximum();
        }
    }
    max_y = max + 0.1*max;
   
    /*
    // ----------------------------------------------------------------------
    //                          Statistical errors
    // ----------------------------------------------------------------------
    // Loop over histograms
    for ( int i = 0; i < n_hists; ++i ){
        // Loop over bins, get value, calculate sqrt(N)
        // Vector for each hist
        int n_bins = 0;
        n_bins = hists[i]->GetNbinsX();

        for ( int j = 0; j < n_bins; ++j ){
            double n_events = 0;
            n_events = hists[i]->GetBinContent(j);   

            double err_val = 0;
            err_val = sqrt(n_events);
            hists[i]->SetBinError(j,err_val);
        }   

    }
    */
    gStyle->SetEndErrorSize(1);
    gStyle->SetErrorX(0);
    
    // Draw the histograms
    hists[0]->Draw("e1");
    hists[0]->Draw("csame");
    for ( int i = 1; i < n_hists; ++i ) {        
        if( i != 4 ){
            // Draw the histograms 
            hists[i]->Draw( "samee1" );
            hists[i]->Draw( "csame" );
        }

    }

    hists[0]->GetXaxis()->SetTitle(x_axis);
    hists[0]->GetYaxis()->SetTitle(y_axis);
    hists[0]->SetAxisRange(0,max_y, "Y");
    hists[0]->SetTitleOffset(1.5, "Y");    
    hists[0]->SetStats(kFALSE);

    leg->Draw();
    c->SaveAs(file_name);

    delete c;
    delete leg;
        
}

// -------------------------------------------------------------------------
//                    reconstructed energy calculation
// -------------------------------------------------------------------------
void RecoNuE( TTree *event_tree,
              vector< double > &reco_E_CC, 
              vector< double > &reco_E_NC,
              vector< double > &MC_reco_E_CC, 
              vector< double > &MC_reco_E_NC ){

    // Get the branches to calculate reconstructed energy and MC energy
    TBranch *b_mu_e  = event_tree->GetBranch("El");
    TBranch *b_nu_e  = event_tree->GetBranch("Ev");
    TBranch *b_mu_p  = event_tree->GetBranch("pl");
    TBranch *b_theta = event_tree->GetBranch("cthl");
    TBranch *b_nfpi0 = event_tree->GetBranch("nfpi0");
    TBranch *b_nfpip = event_tree->GetBranch("nfpip");
    TBranch *b_nfpim = event_tree->GetBranch("nfpim");
    TBranch *b_cc    = event_tree->GetBranch("cc");
    TBranch *b_nc    = event_tree->GetBranch("nc");
    TBranch *b_pdgf  = event_tree->GetBranch("pdgf");
    
    // The variables from the branches and get the leaves
    double m_n   = 0.93828;   // Nucleon mass, GeV
    double m_mu  = 0.10566;   // Muon mass, GeV

    int n_values = event_tree->GetEntries(); // Number of entries to loop over
    
    // Loop over the leaves and calculate the reconstructed energy
    for( int i = 0; i < n_values; ++i){
        
        event_tree->GetEntry(i);

        double reco, reco_mc, e, p, cth;
        
        // For CC0pi
        if( b_cc->GetLeaf("cc")->GetValue() != 0 
            && b_nfpip->GetLeaf("nfpip")->GetValue()
             + b_nfpim->GetLeaf("nfpim")->GetValue()
             + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0
            && b_pdgf->GetLeaf("pdgf")->GetValue() == 13 ){
         
                // Get the values needed
                e   = b_mu_e->GetLeaf("El")->GetValue();
                p   = b_mu_p->GetLeaf("pl")->GetValue();
                cth = b_theta->GetLeaf("cthl")->GetValue(); 
            
                reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( 1 - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 
        
                reco_mc = TMath::Abs( reco - double(b_nu_e->GetLeaf("Ev")->GetValue()) );

                // Make the vectors of reconstructed and reconstructed-MC energy for CC0pi
                MC_reco_E_CC.push_back(reco_mc);
                reco_E_CC.push_back(reco);

        }
        // For NC0pi
        else if( b_nc->GetLeaf("nc")->GetValue() != 0 
                 && b_nfpip->GetLeaf("nfpip")->GetValue()
                  + b_nfpim->GetLeaf("nfpim")->GetValue()
                  + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0
                 && b_pdgf->GetLeaf("pdgf")->GetValue() == 13 ){
              
                // Get the values
                e   = b_mu_e->GetLeaf("El")->GetValue();
                p   = b_mu_p->GetLeaf("pl")->GetValue();
                cth = b_theta->GetLeaf("cthl")->GetValue(); 
            
                reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( 1 - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 
        
                reco_mc = TMath::Abs( reco - double(b_nu_e->GetLeaf("Ev")->GetValue()) );

                // Make the vectors of reconstructed and reconstructed-MC energy for CC0pi
                MC_reco_E_NC.push_back(reco_mc);
                reco_E_NC.push_back(reco);

        }
    }
}

// -------------------------------------------------------------------------
//                    Make final state particles map
// -------------------------------------------------------------------------
void FSPNumbers( TTree *event_tree,
                 m_map &n_fsp ){
    
    // Firstly, get out the trees we want to look at
    // All correspond to number of particles AFTER FSI
    //  - fspl  == 13 (muon)
    //  - fspl  == 14 (numu)
    //  - nfp   == protons
    //  - nfn   == neutrons
    //  - nfpip == pi+ } Collate these to 
    //  - nfpim == pi- } nfcpi (#final charged pions
    //  - nfpi0 == pi0
    // Set the branch addresses for these leaves
    TBranch *b_fspl  = event_tree->GetBranch("fspl");
    TBranch *b_nfp   = event_tree->GetBranch("nfp");
    TBranch *b_nfn   = event_tree->GetBranch("nfn");
    TBranch *b_nfpip = event_tree->GetBranch("nfpip");
    TBranch *b_nfpim = event_tree->GetBranch("nfpim");
    TBranch *b_nfpi0 = event_tree->GetBranch("nfpi0");

    // Create the variables to use as counters
    int nfmu   = 0;
    int nfe    = 0;
    int nfnumu = 0;
    int nfp    = 0;
    int nfn    = 0;
    int nfcpi  = 0;
    int nfpi0  = 0;
   
    // Get the number of events which contain final state muons
    int n_values = event_tree->GetEntries(); // Number of entries to loop over
    
    // Loop over the leaves and calculate the reconstructed energy
    for( int i = 0; i < n_values; ++i){
        
        // Get the current entry
        event_tree->GetEntry(i);
     
        // Count #final state leptons
        if( b_fspl->GetLeaf("fspl")->GetValue() == 13 ){
            ++nfmu;
        }
        else if( b_fspl->GetLeaf("fspl")->GetValue() == 14 ){
            ++nfnumu;
        }
        else if( b_fspl->GetLeaf("fspl")->GetValue() == 11 ){
            // Print the ntuple number for this electron
            cout << " NTuple entry : " << event_tree->GetEntry(i) <<  endl;
            cout << " NTuple event : " << event_tree->GetEvent(i) <<  endl;
            cout << " ----------------------- " << endl;
            ++nfe;
        }

        // Count #final state nucleons
        if( b_nfp->GetLeaf("nfp")->GetValue() == 1 ){
            ++nfp;
        }
        if( b_nfn->GetLeaf("nfn")->GetValue() == 1 ){
            ++nfn;
        }

        // Count #final state pions
        if( b_nfpip->GetLeaf("nfpip")->GetValue() == 1 || b_nfpim->GetLeaf("nfpim")->GetValue() == 1){
            ++nfcpi;           
        }
        if( b_nfpi0->GetLeaf("nfpi0")->GetValue() == 1 ){
            ++nfpi0;           
        }
    }

    // Now fill the map with the values
    // Do this in such a way that the key can be printed straight into a table
    n_fsp.insert( pair< string, int > (" Protons ", nfp));
    n_fsp.insert( pair< string, int > (" Neutrons ", nfn));
    n_fsp.insert( pair< string, int > (" Muons ", nfmu));
    n_fsp.insert( pair< string, int > (" Muon Neutrinos ", nfnumu));
    n_fsp.insert( pair< string, int > (" Charged Pions ", nfcpi));
    n_fsp.insert( pair< string, int > (" Neutral Pions ", nfpi0));
    n_fsp.insert( pair< string, int > (" Electrons ", nfe));
}

// -------------------------------------------------------------------------
//                  Make final state interactions map
// -------------------------------------------------------------------------
void FSINumbers( TTree *event_tree,
                 ostream &file,
                 double norm,
                 m_map &n_cc_fsi,
                 m_map &n_nc_fsi,
                 m_map &n_cc_mc_fsi,
                 m_map &n_nc_mc_fsi ){
    
    // Firstly, get out the trees we want to look at
    // Need both cc and nc with varying number of outgoing pions
    //      - cc : charged current FSI
    //      - nc : neutral current FSI
    // The conditions:
    //      Charged Current
    //      - cc0pi
    //      - cc1pi0
    //      - cc1pi+
    //      - cc1pi-
    //      - ccpi+pi-
    //      - cc2pi+
    //      - cc2pi-
    //      - cc2pi0
    //      - ccpi+pi0
    //      - ccpi-pi0
    //      - cc3pi
    //      - cccoh
    //      Neutral Current
    //      - nc0pi
    //      - nc1pi0
    //      - nc1pi+
    //      - nc1pi-
    //      - ncpi+pi-
    //      - nc2pi+
    //      - nc2pi-
    //      - nc2pi0
    //      - ncpi+pi0
    //      - ncpi-pi0
    //      - nc3pi
    //      - nccoh
    //
    // Set the branch addresses for these leaves
    TBranch *b_pdg = event_tree->GetBranch("pdgf");
    TBranch *b_cc  = event_tree->GetBranch("cc");
    TBranch *b_nc  = event_tree->GetBranch("nc");
    TBranch *b_coh = event_tree->GetBranch("coh");
    TBranch *b_pi0 = event_tree->GetBranch("nfpi0");
    TBranch *b_pip = event_tree->GetBranch("nfpip");
    TBranch *b_pim = event_tree->GetBranch("nfpim");

    // Charged current counters
    int ncc0pi    = 0;
    int ncc1pip   = 0;
    int ncc1pim   = 0;
    int ncc1pi0   = 0;
    int ncc2pip   = 0;
    int ncc2pim   = 0;
    int ncc2pi0   = 0;
    int nccpippim = 0;
    int nccpippi0 = 0;
    int nccpimpi0 = 0;
    int ncc3pi    = 0;
    int ncccoh    = 0;

    // Neutral current counters
    int nnc0pi    = 0;
    int nnc1pip   = 0;
    int nnc1pim   = 0;
    int nnc1pi0   = 0;
    int nnc2pip   = 0;
    int nnc2pim   = 0;
    int nnc2pi0   = 0;
    int nncpippim = 0;
    int nncpippi0 = 0;
    int nncpimpi0 = 0;
    int nnc3pi    = 0;
    int nnccoh    = 0;

    // Get the number of events which contain final state muons
    int n_values = event_tree->GetEntries(); // Number of entries to loop over
    
    // Loop and count for various conditions
    for( int i = 0; i < n_values; ++i){
    
        // Get the current entry
        event_tree->GetEntry(i);
    
        // Charged current
        // CC0Pi
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc0pi;
        }

        // CC1Pi+
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc1pip;
        }
     
        // CC1Pi-
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++ncc1pim;
        }
        
        // CC1Pi0
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc1pi0;
        }
        
        // CC2Pi+
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 2 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc2pip;
        }
     
        // CC2Pi-
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 2 ){
            
            ++ncc2pim;
        }
        
        // CC2Pi0
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 2 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc2pi0;
        }
        
        // CCPi+Pi-
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++nccpippim;
        }
     
        // CCPi+Pi0
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nccpippi0;
        }
        
        // CCPi-Pi0
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++nccpimpi0;
        }
        
        // CC >3Pi
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && ( b_pi0->GetLeaf("nfpi0")->GetValue() 
                + b_pip->GetLeaf("nfpip")->GetValue()
                + b_pim->GetLeaf("nfpim")->GetValue() ) >= 3 ){
            
            ++ncc3pi;
        }
        
        // CCCOH
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_coh->GetLeaf("coh")->GetValue() !=0 ){
            
            ++ncccoh;
        }
       
        //----------------------------------------------------------
        // Neutral current
        // NC0Pi
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nnc0pi;
        }

        // nc1Pi+
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nnc1pip;
        }
     
        // nc1Pi-
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++nnc1pim;
        }
        
        // nc1Pi0
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nnc1pi0;
        }
        
        // nc2Pi+
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 2 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nnc2pip;
        }
     
        // nc2Pi-
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 2 ){
            
            ++nnc2pim;
        }
        
        // nc2Pi0
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 2 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nnc2pi0;
        }
        
        // ncPi+Pi-
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++nncpippim;
        }
     
        // ncPi+Pi0
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 1 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++nncpippi0;
        }
        
        // ncPi-Pi0
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 1 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 1 ){
            
            ++nncpimpi0;
        }
        
        // nc >3Pi
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && ( b_pi0->GetLeaf("nfpi0")->GetValue() 
                + b_pip->GetLeaf("nfpip")->GetValue()
                + b_pim->GetLeaf("nfpim")->GetValue() ) >= 3 ){
            
            ++nnc3pi;
        }
        
        // NCCOH
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_coh->GetLeaf("coh")->GetValue() !=0 ){
            
            ++nnccoh;
        }
    }


    // Now fill the map with the values
    // Do this in such a way that the key can be printed straight into a table
    n_cc_fsi.insert( pair< string, int > (" CC0 $ \\pi$ ", ncc0pi));
    n_cc_fsi.insert( pair< string, int > (" CC1 $ \\pi^{+} $ ", ncc1pip));
    n_cc_fsi.insert( pair< string, int > (" CC1 $ \\pi^{-} $ ", ncc1pim));
    n_cc_fsi.insert( pair< string, int > (" CC1 $ \\pi^{0} $ ", ncc1pi0));
    n_cc_fsi.insert( pair< string, int > (" CC2 $ \\pi^{+} $ ", ncc2pip));
    n_cc_fsi.insert( pair< string, int > (" CC2 $ \\pi^{-} $ ", ncc2pim));
    n_cc_fsi.insert( pair< string, int > (" CC2 $ \\pi^{0} $ ", ncc2pi0));
    n_cc_fsi.insert( pair< string, int > (" CC  $ \\pi^{+}\\pi^{-} $ ", nccpippim));
    n_cc_fsi.insert( pair< string, int > (" CC  $ \\pi^{+}\\pi^{0} $ ", nccpippi0));
    n_cc_fsi.insert( pair< string, int > (" CC  $ \\pi^{-}\\pi^{0} $ ", nccpimpi0));
    n_cc_fsi.insert( pair< string, int > (" CC  $ \\gt 3 \\pi$ ", ncc3pi));
    n_cc_fsi.insert( pair< string, int > (" CCCOH ", ncccoh));

    n_nc_fsi.insert( pair< string, int > (" NC0 $ \\pi$ ", nnc0pi));
    n_nc_fsi.insert( pair< string, int > (" NC1 $ \\pi^{+} $ ", nnc1pip));
    n_nc_fsi.insert( pair< string, int > (" NC1 $ \\pi^{-} $ ", nnc1pim));
    n_nc_fsi.insert( pair< string, int > (" NC1 $ \\pi^{0} $ ", nnc1pi0));
    n_nc_fsi.insert( pair< string, int > (" NC2 $ \\pi^{+} $ ", nnc2pip));
    n_nc_fsi.insert( pair< string, int > (" NC2 $ \\pi^{-} $ ", nnc2pim));
    n_nc_fsi.insert( pair< string, int > (" NC2 $ \\pi^{0} $ ", nnc2pi0));
    n_nc_fsi.insert( pair< string, int > (" NC  $ \\pi^{+}\\pi^{-} $ ", nncpippim));
    n_nc_fsi.insert( pair< string, int > (" NC  $ \\pi^{+}\\pi^{0} $ ", nncpippi0));
    n_nc_fsi.insert( pair< string, int > (" NC  $ \\pi^{-}\\pi^{0} $ ", nncpimpi0));
    n_nc_fsi.insert( pair< string, int > (" NC  $ \\gt 3 \\pi$ ", nnc3pi));
    n_nc_fsi.insert( pair< string, int > (" NCCOH ", nnccoh));

    /*
    cout << " CC0Pi    : " << ncc0pi << endl;
    cout << " CC1Pi+   : " << ncc1pip << endl;
    cout << " CC1Pi-   : " << ncc1pim << endl;
    cout << " CC1Pi0   : " << ncc1pi0 << endl;
    cout << " CC2Pi+   : " << ncc2pip << endl;
    cout << " CC2Pi-   : " << ncc2pim << endl;
    cout << " CC2Pi0   : " << ncc2pi0 << endl;
    cout << " CCPi+Pi- : " << nccpippim << endl;
    cout << " CCPi+Pi0 : " << nccpippi0 << endl;
    cout << " CCPi-Pi0 : " << nccpimpi0 << endl;
    cout << " CCCOH    : " << ncccoh << endl;
    cout << " ------------------------- " << endl;
    cout << " NC0Pi    : " << nnc0pi << endl;
    cout << " NC1Pi+   : " << nnc1pip << endl;
    cout << " NC1Pi-   : " << nnc1pim << endl;
    cout << " NC1Pi0   : " << nnc1pi0 << endl;
    cout << " NC2Pi+   : " << nnc2pip << endl;
    cout << " NC2Pi-   : " << nnc2pim << endl;
    cout << " NC2Pi0   : " << nnc2pi0 << endl;
    cout << " NCPi+Pi- : " << nncpippim << endl;
    cout << " NCPi+Pi0 : " << nncpippi0 << endl;
    cout << " NCPi-Pi0 : " << nncpimpi0 << endl;
    cout << " NCCOH    : " << nnccoh << endl;
    cout << " ------------------------- " << endl;
    cout << " ------------------------- " << endl;
    */
}

// -------------------------------------------------------------------------
//                      Make final state tables
// -------------------------------------------------------------------------
void MakeTable( const m_outer &n_cc_map,
                const m_outer &n_nc_map,
                ostream &file){
    /*
    // Iterators 
    typedef m_map::const_iterator map_it;
    typedef m_outer::const_iterator outer_it;

    // Make a table from an input map - of - maps to print nice things
    // Loop over the outer map and write the first line of the table 
    // Number of columns to be made
    int n_models; 
    n_models = n_cc_map.size();
   
    vector< string > m_names;

    for (outer_it it = n_cc_map.begin(); it != n_cc_map.end(); it++ ){
        m_names.push_back(it->first);
    }

    // Begin the tabular environment in the LaTeX output file for n_predictions columns
    file << "\\begin{tabular}{| c || * {" << n_models << "}{c | } }" << endl;
    file << "\\hline" << endl;
 
    // Fill the first line with "Configurations" for all the right-hand columns
    file << "  & \\multicolumn{ " << n_models << " }{c|}{ \\textbf{Configurations} } \\\\" << endl;
    file << "\\hline" << endl;
 
 
    // Fill a line of the table with the prediction names
    file << " \\rotatebox{90}{ \\textbf{ Final Hadronic State } } & " ;
    for( unsigned int i = 0; i < n_models - 1; ++i ){
        file << "\\rotatebox{90}{ \\textbf{ " << m_names[i] << " } } & ";
    }
    file << "\\rotatebox{90}{ \\textbf{ " << m_names[n_models - 1] << " } } \\\\" << endl;
    file << " \\hline " << endl;
 
     // -------------------------------------------------------------------
     //      Loop over the outer maps and fill the rest of the table                                   
     // -------------------------------------------------------------------

     for( map_it it1 = n_cc_map->second.begin(); it1 != n_cc_map->second.end(); it1++ ){


        for( outer_it it2 = n_cc_map.begin(); it3 = --n_cc_map.end(); it2 != it3; ++it2 ){
            
            file << it1->first << " & " << it1->
        }
     }
 
     // -------------------------------------------------------------------
     //                    End the table                                   
     // -------------------------------------------------------------------
 
     // End the tabular environment
     file << "\\end{tabular}" << endl;
    */

}

