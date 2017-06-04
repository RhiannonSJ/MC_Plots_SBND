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

int model_comparisons() {

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
    //           Charged and neutral current maps to put into table
    // -------------------------------------------------------------------------
    m_outer cc_mc_model_ints;
    m_outer nc_mc_model_ints;
    m_outer cc_reco_model_ints;
    m_outer nc_reco_model_ints;
    
    vector< double > n_cc_fsi_1;
    vector< double > n_cc_fsi_2;
    vector< double > n_cc_fsi_3;
    vector< double > n_cc_fsi_4;
    vector< double > n_cc_fsi_5;
    vector< double > n_nc_fsi_1;
    vector< double > n_nc_fsi_2;
    vector< double > n_nc_fsi_3;
    vector< double > n_nc_fsi_4;
    vector< double > n_nc_fsi_5;
    
    vector< double > n_cc_mc_fsi_1;
    vector< double > n_cc_mc_fsi_2;
    vector< double > n_cc_mc_fsi_3;
    vector< double > n_cc_mc_fsi_4;
    vector< double > n_cc_mc_fsi_5;
    vector< double > n_nc_mc_fsi_1;
    vector< double > n_nc_mc_fsi_2;
    vector< double > n_nc_mc_fsi_3;
    vector< double > n_nc_mc_fsi_4;
    vector< double > n_nc_mc_fsi_5;
   
    ofstream file_n;
    file_n.open("n_interactions.txt");

    file_n << " Default " << endl;
    FSINumbers( gst1, file_n, norms[0], n_cc_fsi_1, n_nc_fsi_1, n_cc_mc_fsi_1, n_nc_mc_fsi_1  );
    cc_reco_model_ints.insert( pair< string, vector< double > >( "Default", n_cc_fsi_1 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "Default", n_nc_fsi_1 ) );
    cc_mc_model_ints.insert( pair< string, vector< double > >( "Default", n_cc_mc_fsi_1 ) );
    nc_mc_model_ints.insert( pair< string, vector< double > >( "Default", n_nc_mc_fsi_1 ) );
    
    file_n << " Default + MEC " << endl;
    FSINumbers( gst2, file_n, norms[1], n_cc_fsi_2, n_nc_fsi_2, n_cc_mc_fsi_2, n_nc_mc_fsi_2  );
    cc_reco_model_ints.insert( pair< string, vector< double > >( "Default + MEC", n_cc_fsi_2 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "Default + MEC", n_nc_fsi_2 ) );
    cc_mc_model_ints.insert( pair< string, vector< double > >( "Default + MEC", n_cc_mc_fsi_2 ) );
    nc_mc_model_ints.insert( pair< string, vector< double > >( "Default + MEC", n_nc_mc_fsi_2 ) );
    
    
    file_n << " G16_02b " << endl;
    FSINumbers( gst3, file_n, norms[2], n_cc_fsi_3, n_nc_fsi_3, n_cc_mc_fsi_3, n_nc_mc_fsi_3 );
    cc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_02b", n_cc_fsi_3 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_02b", n_nc_fsi_3 ) );
    cc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_02b", n_cc_mc_fsi_3 ) );
    nc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_02b", n_nc_mc_fsi_3 ) );
    
    
    file_n << " G16_01a " << endl;
    FSINumbers( gst4, file_n, norms[3], n_cc_fsi_4, n_nc_fsi_4, n_cc_mc_fsi_4, n_nc_mc_fsi_4 );
    cc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01a", n_cc_fsi_4 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01a", n_nc_fsi_4 ) );
    cc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_01a", n_cc_mc_fsi_4 ) );
    nc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_01a", n_nc_mc_fsi_4 ) );
    
    
    file_n << " G16_01b " << endl;
    FSINumbers( gst5, file_n, norms[4], n_cc_fsi_5, n_nc_fsi_5, n_cc_mc_fsi_5, n_nc_mc_fsi_5 );
    cc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_cc_fsi_5 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_nc_fsi_5 ) );
    cc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_cc_mc_fsi_5 ) );
    nc_mc_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_nc_mc_fsi_5 ) );
    
    // Make vector of final hadronic states
    // Using LaTeX syntax
    vector< string > FHS;

    FHS.push_back( "0 \\( \\pi \\)" );
    FHS.push_back( "1 \\( \\pi^+ \\)" );
    FHS.push_back( "1 \\( \\pi^- \\)" );
    FHS.push_back( "1 \\( \\pi^0 \\)" );
    FHS.push_back( "2 \\( \\pi^+ \\)" );
    FHS.push_back( "2 \\( \\pi^- \\)" );
    FHS.push_back( "2 \\( \\pi^0 \\)" );
    FHS.push_back( "\\( \\pi^+ \\pi^- \\)" );
    FHS.push_back( "\\( \\pi^+ \\pi^0 \\)" );
    FHS.push_back( "\\( \\pi^- \\pi^0 \\)" );
    FHS.push_back( "\\( > 3 \\pi \\)" );
    FHS.push_back( "Coherent" );

    ofstream file_reco;
    file_reco.open( "FSI_Reco_Table.tex" );
    MakeTable( cc_reco_model_ints, nc_reco_model_ints, FHS, file_reco );
   
    ofstream file_mc;
    file_mc.open( "FSI_MC_Table.tex" );
    MakeTable( cc_mc_model_ints, nc_mc_model_ints, FHS, file_mc );
    
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
    //                          Q2 Absolute norm CC0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_abs_cc0pi_stack;

    // The Q2 histograms
    TH1D *h1_Q2_abs_cc0pi = new TH1D("h1_Q2_abs_cc0pi","Q^{2} Absolute Normalisation, CC0#pi",50,0,2.5);
    TH1D *h2_Q2_abs_cc0pi = new TH1D("h2_Q2_abs_cc0pi","Q^{2} Absolute Normalisation, CC0#pi",50,0,2.5);
    TH1D *h3_Q2_abs_cc0pi = new TH1D("h3_Q2_abs_cc0pi","Q^{2} Absolute Normalisation, CC0#pi",50,0,2.5);
    TH1D *h4_Q2_abs_cc0pi = new TH1D("h4_Q2_abs_cc0pi","Q^{2} Absolute Normalisation, CC0#pi",50,0,2.5);
    TH1D *h5_Q2_abs_cc0pi = new TH1D("h5_Q2_abs_cc0pi","Q^{2} Absolute Normalisation, CC0#pi",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("Q2>>h2_Q2_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("Q2>>h3_Q2_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("Q2>>h4_Q2_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("Q2>>h5_Q2_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Q2_abs_cc0pi_stack.push_back(h1_Q2_abs_cc0pi);
    Q2_abs_cc0pi_stack.push_back(h2_Q2_abs_cc0pi);
    Q2_abs_cc0pi_stack.push_back(h3_Q2_abs_cc0pi);
    Q2_abs_cc0pi_stack.push_back(h4_Q2_abs_cc0pi);
    Q2_abs_cc0pi_stack.push_back(h5_Q2_abs_cc0pi);

    HistStacker(Q2_abs_cc0pi_stack, leg_entries, norms, "Q2 Absolute Normalisation, CC0#pi", "Q2_CC0Pi_Plot.png", "Q^{2}", "Number of SBND events");

    delete h1_Q2_abs_cc0pi;
    delete h2_Q2_abs_cc0pi;
    delete h3_Q2_abs_cc0pi;
    delete h4_Q2_abs_cc0pi;
    delete h5_Q2_abs_cc0pi;
    
    // -------------------------------------------------------------------------
    //                          Q2 Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_abs_cc1pip_stack;

    // The Q2 histograms
    TH1D *h1_Q2_abs_cc1pip = new TH1D("h1_Q2_abs_cc1pip","Q^{2} Absolute Normalisation, CC1#pi^{+}",50,0,2.5);
    TH1D *h2_Q2_abs_cc1pip = new TH1D("h2_Q2_abs_cc1pip","Q^{2} Absolute Normalisation, CC1#pi^{+}",50,0,2.5);
    TH1D *h3_Q2_abs_cc1pip = new TH1D("h3_Q2_abs_cc1pip","Q^{2} Absolute Normalisation, CC1#pi^{+}",50,0,2.5);
    TH1D *h4_Q2_abs_cc1pip = new TH1D("h4_Q2_abs_cc1pip","Q^{2} Absolute Normalisation, CC1#pi^{+}",50,0,2.5);
    TH1D *h5_Q2_abs_cc1pip = new TH1D("h5_Q2_abs_cc1pip","Q^{2} Absolute Normalisation, CC1#pi^{+}",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    
    // Push them onto the vector
    Q2_abs_cc1pip_stack.push_back(h1_Q2_abs_cc1pip);
    Q2_abs_cc1pip_stack.push_back(h2_Q2_abs_cc1pip);
    Q2_abs_cc1pip_stack.push_back(h3_Q2_abs_cc1pip);
    Q2_abs_cc1pip_stack.push_back(h4_Q2_abs_cc1pip);
    Q2_abs_cc1pip_stack.push_back(h5_Q2_abs_cc1pip);

    HistStacker(Q2_abs_cc1pip_stack, leg_entries, norms, "Q2 Absolute Normalisation, CC1#pi^{+}", "Q2_CC1PiP_Plot.png", "Q^{2}", "Number of SBND events");

    delete h1_Q2_abs_cc1pip;
    delete h2_Q2_abs_cc1pip;
    delete h3_Q2_abs_cc1pip;
    delete h4_Q2_abs_cc1pip;
    delete h5_Q2_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          Q2 Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_abs_cc1pi0_stack;

    // The Q2 histograms
    TH1D *h1_Q2_abs_cc1pi0 = new TH1D("h1_Q2_abs_cc1pi0","Q^{2} Absolute Normalisation, CC1#pi^{0}",50,0,2.5);
    TH1D *h2_Q2_abs_cc1pi0 = new TH1D("h2_Q2_abs_cc1pi0","Q^{2} Absolute Normalisation, CC1#pi^{0}",50,0,2.5);
    TH1D *h3_Q2_abs_cc1pi0 = new TH1D("h3_Q2_abs_cc1pi0","Q^{2} Absolute Normalisation, CC1#pi^{0}",50,0,2.5);
    TH1D *h4_Q2_abs_cc1pi0 = new TH1D("h4_Q2_abs_cc1pi0","Q^{2} Absolute Normalisation, CC1#pi^{0}",50,0,2.5);
    TH1D *h5_Q2_abs_cc1pi0 = new TH1D("h5_Q2_abs_cc1pi0","Q^{2} Absolute Normalisation, CC1#pi^{0}",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    
    // Push them onto the vector
    Q2_abs_cc1pi0_stack.push_back(h1_Q2_abs_cc1pi0);
    Q2_abs_cc1pi0_stack.push_back(h2_Q2_abs_cc1pi0);
    Q2_abs_cc1pi0_stack.push_back(h3_Q2_abs_cc1pi0);
    Q2_abs_cc1pi0_stack.push_back(h4_Q2_abs_cc1pi0);
    Q2_abs_cc1pi0_stack.push_back(h5_Q2_abs_cc1pi0);

    HistStacker(Q2_abs_cc1pi0_stack, leg_entries, norms, "Q2 Absolute Normalisation, CC1#pi^{0}", "Q2_CC1Pi0_Plot.png", "Q^{2}", "Number of SBND events");

    delete h1_Q2_abs_cc1pi0;
    delete h2_Q2_abs_cc1pi0;
    delete h3_Q2_abs_cc1pi0;
    delete h4_Q2_abs_cc1pi0;
    delete h5_Q2_abs_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          Q2 Absolute norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_abs_nc1pi0_stack;

    // The Q2 histograms
    TH1D *h1_Q2_abs_nc1pi0 = new TH1D("h1_Q2_abs_nc1pi0","Q^{2} Absolute Normalisation, NC1#pi^{0}",50,0,2.5);
    TH1D *h2_Q2_abs_nc1pi0 = new TH1D("h2_Q2_abs_nc1pi0","Q^{2} Absolute Normalisation, NC1#pi^{0}",50,0,2.5);
    TH1D *h3_Q2_abs_nc1pi0 = new TH1D("h3_Q2_abs_nc1pi0","Q^{2} Absolute Normalisation, NC1#pi^{0}",50,0,2.5);
    TH1D *h4_Q2_abs_nc1pi0 = new TH1D("h4_Q2_abs_nc1pi0","Q^{2} Absolute Normalisation, NC1#pi^{0}",50,0,2.5);
    TH1D *h5_Q2_abs_nc1pi0 = new TH1D("h5_Q2_abs_nc1pi0","Q^{2} Absolute Normalisation, NC1#pi^{0}",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    
    // Push them onto the vector
    Q2_abs_nc1pi0_stack.push_back(h1_Q2_abs_nc1pi0);
    Q2_abs_nc1pi0_stack.push_back(h2_Q2_abs_nc1pi0);
    Q2_abs_nc1pi0_stack.push_back(h3_Q2_abs_nc1pi0);
    Q2_abs_nc1pi0_stack.push_back(h4_Q2_abs_nc1pi0);
    Q2_abs_nc1pi0_stack.push_back(h5_Q2_abs_nc1pi0);

    HistStacker(Q2_abs_nc1pi0_stack, leg_entries, norms, "Q2 Absolute Normalisation, NC1#pi^{0}", "Q2_NC1Pi0_Plot.png", "Q^{2}", "Number of SBND events");

    delete h1_Q2_abs_nc1pi0;
    delete h2_Q2_abs_nc1pi0;
    delete h3_Q2_abs_nc1pi0;
    delete h4_Q2_abs_nc1pi0;
    delete h5_Q2_abs_nc1pi0;

    // -------------------------------------------------------------------------
    //                          Q2 Area norm CC0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_area_cc0pi_stack;
    vector< double > int_Q2_cc0pi;

    // The Q2 histograms
    TH1D *h1_Q2_area_cc0pi = new TH1D("h1_Q2_area_cc0pi","Q^{2} Area Normalised, CC0#pi",50,0,2.5);
    TH1D *h2_Q2_area_cc0pi = new TH1D("h2_Q2_area_cc0pi","Q^{2} Area Normalised, CC0#pi",50,0,2.5);
    TH1D *h3_Q2_area_cc0pi = new TH1D("h3_Q2_area_cc0pi","Q^{2} Area Normalised, CC0#pi",50,0,2.5);
    TH1D *h4_Q2_area_cc0pi = new TH1D("h4_Q2_area_cc0pi","Q^{2} Area Normalised, CC0#pi",50,0,2.5);
    TH1D *h5_Q2_area_cc0pi = new TH1D("h5_Q2_area_cc0pi","Q^{2} Area Normalised, CC0#pi",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("Q2>>h2_Q2_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("Q2>>h3_Q2_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("Q2>>h4_Q2_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("Q2>>h5_Q2_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push back the scaling
    double int1_Q2_cc0pi, int2_Q2_cc0pi, int3_Q2_cc0pi, int4_Q2_cc0pi, int5_Q2_cc0pi;

    int1_Q2_cc0pi = h1_Q2_area_cc0pi->Integral();
    int2_Q2_cc0pi = h2_Q2_area_cc0pi->Integral();
    int3_Q2_cc0pi = h3_Q2_area_cc0pi->Integral();
    int4_Q2_cc0pi = h4_Q2_area_cc0pi->Integral();
    int5_Q2_cc0pi = h5_Q2_area_cc0pi->Integral();
    
    int_Q2_cc0pi.push_back( 1 / int1_Q2_cc0pi );
    int_Q2_cc0pi.push_back( 1 / int2_Q2_cc0pi );
    int_Q2_cc0pi.push_back( 1 / int3_Q2_cc0pi );
    int_Q2_cc0pi.push_back( 1 / int4_Q2_cc0pi );
    int_Q2_cc0pi.push_back( 1 / int5_Q2_cc0pi );
    
    // Push them onto the vector
    Q2_area_cc0pi_stack.push_back(h1_Q2_area_cc0pi);
    Q2_area_cc0pi_stack.push_back(h2_Q2_area_cc0pi);
    Q2_area_cc0pi_stack.push_back(h3_Q2_area_cc0pi);
    Q2_area_cc0pi_stack.push_back(h4_Q2_area_cc0pi);
    Q2_area_cc0pi_stack.push_back(h5_Q2_area_cc0pi);

    HistStacker(Q2_area_cc0pi_stack, leg_entries, int_Q2_cc0pi, "Q2 Area Normalised, CC0#pi", "Q2_CC0Pi_Area_Plot.png", "Q^{2}", "SBND events, area normalised");

    delete h1_Q2_area_cc0pi;
    delete h2_Q2_area_cc0pi;
    delete h3_Q2_area_cc0pi;
    delete h4_Q2_area_cc0pi;
    delete h5_Q2_area_cc0pi;
    
    // -------------------------------------------------------------------------
    //                          Q2 Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_area_cc1pip_stack;
    vector< double > int_Q2_cc1pip;

    // The Q2 histograms
    TH1D *h1_Q2_area_cc1pip = new TH1D("h1_Q2_area_cc1pip","Q^{2} Area Normalised, CC1#pi^{+}",50,0,2.5);
    TH1D *h2_Q2_area_cc1pip = new TH1D("h2_Q2_area_cc1pip","Q^{2} Area Normalised, CC1#pi^{+}",50,0,2.5);
    TH1D *h3_Q2_area_cc1pip = new TH1D("h3_Q2_area_cc1pip","Q^{2} Area Normalised, CC1#pi^{+}",50,0,2.5);
    TH1D *h4_Q2_area_cc1pip = new TH1D("h4_Q2_area_cc1pip","Q^{2} Area Normalised, CC1#pi^{+}",50,0,2.5);
    TH1D *h5_Q2_area_cc1pip = new TH1D("h5_Q2_area_cc1pip","Q^{2} Area Normalised, CC1#pi^{+}",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
     
    // Push back the scaling
    double int1_Q2_cc1pip, int2_Q2_cc1pip, int3_Q2_cc1pip, int4_Q2_cc1pip, int5_Q2_cc1pip;

    int1_Q2_cc1pip = h1_Q2_area_cc1pip->Integral();
    int2_Q2_cc1pip = h2_Q2_area_cc1pip->Integral();
    int3_Q2_cc1pip = h3_Q2_area_cc1pip->Integral();
    int4_Q2_cc1pip = h4_Q2_area_cc1pip->Integral();
    int5_Q2_cc1pip = h5_Q2_area_cc1pip->Integral();
    
    int_Q2_cc1pip.push_back( 1 / int1_Q2_cc1pip );
    int_Q2_cc1pip.push_back( 1 / int2_Q2_cc1pip );
    int_Q2_cc1pip.push_back( 1 / int3_Q2_cc1pip );
    int_Q2_cc1pip.push_back( 1 / int4_Q2_cc1pip );
    int_Q2_cc1pip.push_back( 1 / int5_Q2_cc1pip );
    
    // Push them onto the vector
    Q2_area_cc1pip_stack.push_back(h1_Q2_area_cc1pip);
    Q2_area_cc1pip_stack.push_back(h2_Q2_area_cc1pip);
    Q2_area_cc1pip_stack.push_back(h3_Q2_area_cc1pip);
    Q2_area_cc1pip_stack.push_back(h4_Q2_area_cc1pip);
    Q2_area_cc1pip_stack.push_back(h5_Q2_area_cc1pip);

    HistStacker(Q2_area_cc1pip_stack, leg_entries, int_Q2_cc1pip, "Q2 Area Normalised, CC1#pi^{+}", "Q2_CC1PiP_Area_Plot.png", "Q^{2}", "SBND events, area normalised");

    delete h1_Q2_area_cc1pip;
    delete h2_Q2_area_cc1pip;
    delete h3_Q2_area_cc1pip;
    delete h4_Q2_area_cc1pip;
    delete h5_Q2_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          Q2 Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_area_cc1pi0_stack;
    vector< double > int_Q2_cc1pi0;

    // The Q2 histograms
    TH1D *h1_Q2_area_cc1pi0 = new TH1D("h1_Q2_area_cc1pi0","Q^{2} Area Normalised, CC1#pi^{0}",50,0,2.5);
    TH1D *h2_Q2_area_cc1pi0 = new TH1D("h2_Q2_area_cc1pi0","Q^{2} Area Normalised, CC1#pi^{0}",50,0,2.5);
    TH1D *h3_Q2_area_cc1pi0 = new TH1D("h3_Q2_area_cc1pi0","Q^{2} Area Normalised, CC1#pi^{0}",50,0,2.5);
    TH1D *h4_Q2_area_cc1pi0 = new TH1D("h4_Q2_area_cc1pi0","Q^{2} Area Normalised, CC1#pi^{0}",50,0,2.5);
    TH1D *h5_Q2_area_cc1pi0 = new TH1D("h5_Q2_area_cc1pi0","Q^{2} Area Normalised, CC1#pi^{0}",50,0,2.5);
    
    gst1->Draw("Q2>>h1_Q2_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
     
    // Push back the scaling
    double int1_Q2_cc1pi0, int2_Q2_cc1pi0, int3_Q2_cc1pi0, int4_Q2_cc1pi0, int5_Q2_cc1pi0;

    int1_Q2_cc1pi0 = h1_Q2_area_cc1pi0->Integral();
    int2_Q2_cc1pi0 = h2_Q2_area_cc1pi0->Integral();
    int3_Q2_cc1pi0 = h3_Q2_area_cc1pi0->Integral();
    int4_Q2_cc1pi0 = h4_Q2_area_cc1pi0->Integral();
    int5_Q2_cc1pi0 = h5_Q2_area_cc1pi0->Integral();
    
    int_Q2_cc1pi0.push_back( 1 / int1_Q2_cc1pi0 );
    int_Q2_cc1pi0.push_back( 1 / int2_Q2_cc1pi0 );
    int_Q2_cc1pi0.push_back( 1 / int3_Q2_cc1pi0 );
    int_Q2_cc1pi0.push_back( 1 / int4_Q2_cc1pi0 );
    int_Q2_cc1pi0.push_back( 1 / int5_Q2_cc1pi0 );
    
    // Push them onto the vector
    Q2_area_cc1pi0_stack.push_back(h1_Q2_area_cc1pi0);
    Q2_area_cc1pi0_stack.push_back(h2_Q2_area_cc1pi0);
    Q2_area_cc1pi0_stack.push_back(h3_Q2_area_cc1pi0);
    Q2_area_cc1pi0_stack.push_back(h4_Q2_area_cc1pi0);
    Q2_area_cc1pi0_stack.push_back(h5_Q2_area_cc1pi0);

    HistStacker(Q2_area_cc1pi0_stack, leg_entries, int_Q2_cc1pi0, "Q2 Area Normalised, CC1#pi^{0}", "Q2_CC1Pi0_Area_Plot.png", "Q^{2}", "SBND events, area normalised");

    delete h1_Q2_area_cc1pi0;
    delete h2_Q2_area_cc1pi0;
    delete h3_Q2_area_cc1pi0;
    delete h4_Q2_area_cc1pi0;
    delete h5_Q2_area_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          Q2 Area norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Q2_area_nc1pi0_stack;
    vector< double > int_Q2_nc1pi0;

    // The Q2 histograms
    TH1D *h1_Q2_area_nc1pi0 = new TH1D("h1_Q2_area_nc1pi0","Q^{2} Area Normalised, NC1#pi^{0}",50,0,2.5);
    TH1D *h2_Q2_area_nc1pi0 = new TH1D("h2_Q2_area_nc1pi0","Q^{2} Area Normalised, NC1#pi^{0}",50,0,2.5);
    TH1D *h3_Q2_area_nc1pi0 = new TH1D("h3_Q2_area_nc1pi0","Q^{2} Area Normalised, NC1#pi^{0}",50,0,2.5);
    TH1D *h4_Q2_area_nc1pi0 = new TH1D("h4_Q2_area_nc1pi0","Q^{2} Area Normalised, NC1#pi^{0}",50,0,2.5);
    TH1D *h5_Q2_area_nc1pi0 = new TH1D("h5_Q2_area_nc1pi0","Q^{2} Area Normalised, NC1#pi^{0}",50,0,2.5);

    gst1->Draw("Q2>>h1_Q2_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Q2>>h2_Q2_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Q2>>h3_Q2_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Q2>>h4_Q2_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Q2>>h5_Q2_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    
    // Push back the scaling
    double int1_Q2_nc1pi0, int2_Q2_nc1pi0, int3_Q2_nc1pi0, int4_Q2_nc1pi0, int5_Q2_nc1pi0;

    int1_Q2_nc1pi0 = h1_Q2_area_nc1pi0->Integral();
    int2_Q2_nc1pi0 = h2_Q2_area_nc1pi0->Integral();
    int3_Q2_nc1pi0 = h3_Q2_area_nc1pi0->Integral();
    int4_Q2_nc1pi0 = h4_Q2_area_nc1pi0->Integral();
    int5_Q2_nc1pi0 = h5_Q2_area_nc1pi0->Integral();
    
    int_Q2_nc1pi0.push_back( 1 / int1_Q2_nc1pi0 );
    int_Q2_nc1pi0.push_back( 1 / int2_Q2_nc1pi0 );
    int_Q2_nc1pi0.push_back( 1 / int3_Q2_nc1pi0 );
    int_Q2_nc1pi0.push_back( 1 / int4_Q2_nc1pi0 );
    int_Q2_nc1pi0.push_back( 1 / int5_Q2_nc1pi0 );
    
    // Push them onto the vector
    Q2_area_nc1pi0_stack.push_back(h1_Q2_area_nc1pi0);
    Q2_area_nc1pi0_stack.push_back(h2_Q2_area_nc1pi0);
    Q2_area_nc1pi0_stack.push_back(h3_Q2_area_nc1pi0);
    Q2_area_nc1pi0_stack.push_back(h4_Q2_area_nc1pi0);
    Q2_area_nc1pi0_stack.push_back(h5_Q2_area_nc1pi0);

    HistStacker(Q2_area_nc1pi0_stack, leg_entries, int_Q2_nc1pi0, "Q2 Area Normalised, NC1#pi^{0}", "Q2_NC1Pi0_Area_Plot.png", "Q^{2}", "SBND events, area normalised");

    delete h1_Q2_area_nc1pi0;
    delete h2_Q2_area_nc1pi0;
    delete h3_Q2_area_nc1pi0;
    delete h4_Q2_area_nc1pi0;
    delete h5_Q2_area_nc1pi0;

    // -------------------------------------------------------------------------
    //                          Ev Absolute norm CC0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_abs_cc0pi_stack;

    // The Ev histograms
    TH1D *h1_Ev_abs_cc0pi = new TH1D("h1_Ev_abs_cc0pi","E_{#nu} Absolute Normalisation, CC0#pi",50,0,3);
    TH1D *h2_Ev_abs_cc0pi = new TH1D("h2_Ev_abs_cc0pi","E_{#nu} Absolute Normalisation, CC0#pi",50,0,3);
    TH1D *h3_Ev_abs_cc0pi = new TH1D("h3_Ev_abs_cc0pi","E_{#nu} Absolute Normalisation, CC0#pi",50,0,3);
    TH1D *h4_Ev_abs_cc0pi = new TH1D("h4_Ev_abs_cc0pi","E_{#nu} Absolute Normalisation, CC0#pi",50,0,3);
    TH1D *h5_Ev_abs_cc0pi = new TH1D("h5_Ev_abs_cc0pi","E_{#nu} Absolute Normalisation, CC0#pi",50,0,3);

    gst1->Draw("Ev>>h1_Ev_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("Ev>>h2_Ev_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("Ev>>h3_Ev_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("Ev>>h4_Ev_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("Ev>>h5_Ev_abs_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push them onto the vector
    Ev_abs_cc0pi_stack.push_back(h1_Ev_abs_cc0pi);
    Ev_abs_cc0pi_stack.push_back(h2_Ev_abs_cc0pi);
    Ev_abs_cc0pi_stack.push_back(h3_Ev_abs_cc0pi);
    Ev_abs_cc0pi_stack.push_back(h4_Ev_abs_cc0pi);
    Ev_abs_cc0pi_stack.push_back(h5_Ev_abs_cc0pi);

    HistStacker(Ev_abs_cc0pi_stack, leg_entries, norms, "Ev Absolute Normalisation, CC0#pi", "Ev_CC0Pi_Plot.png", "E_{#nu}", "Number of SBND events");

    delete h1_Ev_abs_cc0pi;
    delete h2_Ev_abs_cc0pi;
    delete h3_Ev_abs_cc0pi;
    delete h4_Ev_abs_cc0pi;
    delete h5_Ev_abs_cc0pi;
    
    // -------------------------------------------------------------------------
    //                          Ev Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_abs_cc1pip_stack;

    // The Ev histograms
    TH1D *h1_Ev_abs_cc1pip = new TH1D("h1_Ev_abs_cc1pip","E_{#nu} Absolute Normalisation, CC1#pi^{+}",50,0,3);
    TH1D *h2_Ev_abs_cc1pip = new TH1D("h2_Ev_abs_cc1pip","E_{#nu} Absolute Normalisation, CC1#pi^{+}",50,0,3);
    TH1D *h3_Ev_abs_cc1pip = new TH1D("h3_Ev_abs_cc1pip","E_{#nu} Absolute Normalisation, CC1#pi^{+}",50,0,3);
    TH1D *h4_Ev_abs_cc1pip = new TH1D("h4_Ev_abs_cc1pip","E_{#nu} Absolute Normalisation, CC1#pi^{+}",50,0,3);
    TH1D *h5_Ev_abs_cc1pip = new TH1D("h5_Ev_abs_cc1pip","E_{#nu} Absolute Normalisation, CC1#pi^{+}",50,0,3);

    gst1->Draw("Ev>>h1_Ev_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    
    // Push them onto the vector
    Ev_abs_cc1pip_stack.push_back(h1_Ev_abs_cc1pip);
    Ev_abs_cc1pip_stack.push_back(h2_Ev_abs_cc1pip);
    Ev_abs_cc1pip_stack.push_back(h3_Ev_abs_cc1pip);
    Ev_abs_cc1pip_stack.push_back(h4_Ev_abs_cc1pip);
    Ev_abs_cc1pip_stack.push_back(h5_Ev_abs_cc1pip);

    HistStacker(Ev_abs_cc1pip_stack, leg_entries, norms, "Ev Absolute Normalisation, CC1#pi^{+}", "Ev_CC1PiP_Plot.png", "E_{#nu}", "Number of SBND events");

    delete h1_Ev_abs_cc1pip;
    delete h2_Ev_abs_cc1pip;
    delete h3_Ev_abs_cc1pip;
    delete h4_Ev_abs_cc1pip;
    delete h5_Ev_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          Ev Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_abs_cc1pi0_stack;

    // The Ev histograms
    TH1D *h1_Ev_abs_cc1pi0 = new TH1D("h1_Ev_abs_cc1pi0","E_{#nu} Absolute Normalisation, CC1#pi^{0}",50,0,3);
    TH1D *h2_Ev_abs_cc1pi0 = new TH1D("h2_Ev_abs_cc1pi0","E_{#nu} Absolute Normalisation, CC1#pi^{0}",50,0,3);
    TH1D *h3_Ev_abs_cc1pi0 = new TH1D("h3_Ev_abs_cc1pi0","E_{#nu} Absolute Normalisation, CC1#pi^{0}",50,0,3);
    TH1D *h4_Ev_abs_cc1pi0 = new TH1D("h4_Ev_abs_cc1pi0","E_{#nu} Absolute Normalisation, CC1#pi^{0}",50,0,3);
    TH1D *h5_Ev_abs_cc1pi0 = new TH1D("h5_Ev_abs_cc1pi0","E_{#nu} Absolute Normalisation, CC1#pi^{0}",50,0,3);

    gst1->Draw("Ev>>h1_Ev_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");

    // Push them onto the vector
    Ev_abs_cc1pi0_stack.push_back(h1_Ev_abs_cc1pi0);
    Ev_abs_cc1pi0_stack.push_back(h2_Ev_abs_cc1pi0);
    Ev_abs_cc1pi0_stack.push_back(h3_Ev_abs_cc1pi0);
    Ev_abs_cc1pi0_stack.push_back(h4_Ev_abs_cc1pi0);
    Ev_abs_cc1pi0_stack.push_back(h5_Ev_abs_cc1pi0);

    HistStacker(Ev_abs_cc1pi0_stack, leg_entries, norms, "Ev Absolute Normalisation, CC1#pi^{0}", "Ev_CC1Pi0_Plot.png", "E_{#nu}", "Number of SBND events");

    delete h1_Ev_abs_cc1pi0;
    delete h2_Ev_abs_cc1pi0;
    delete h3_Ev_abs_cc1pi0;
    delete h4_Ev_abs_cc1pi0;
    delete h5_Ev_abs_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          Ev Absolute norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_abs_nc1pi0_stack;

    // The Ev histograms
    TH1D *h1_Ev_abs_nc1pi0 = new TH1D("h1_Ev_abs_nc1pi0","E_{#nu} Absolute Normalisation, NC1#pi^{0}",50,0,3);
    TH1D *h2_Ev_abs_nc1pi0 = new TH1D("h2_Ev_abs_nc1pi0","E_{#nu} Absolute Normalisation, NC1#pi^{0}",50,0,3);
    TH1D *h3_Ev_abs_nc1pi0 = new TH1D("h3_Ev_abs_nc1pi0","E_{#nu} Absolute Normalisation, NC1#pi^{0}",50,0,3);
    TH1D *h4_Ev_abs_nc1pi0 = new TH1D("h4_Ev_abs_nc1pi0","E_{#nu} Absolute Normalisation, NC1#pi^{0}",50,0,3);
    TH1D *h5_Ev_abs_nc1pi0 = new TH1D("h5_Ev_abs_nc1pi0","E_{#nu} Absolute Normalisation, NC1#pi^{0}",50,0,3);

    gst1->Draw("Ev>>h1_Ev_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    
    // Push them onto the vector
    Ev_abs_nc1pi0_stack.push_back(h1_Ev_abs_nc1pi0);
    Ev_abs_nc1pi0_stack.push_back(h2_Ev_abs_nc1pi0);
    Ev_abs_nc1pi0_stack.push_back(h3_Ev_abs_nc1pi0);
    Ev_abs_nc1pi0_stack.push_back(h4_Ev_abs_nc1pi0);
    Ev_abs_nc1pi0_stack.push_back(h5_Ev_abs_nc1pi0);

    HistStacker(Ev_abs_nc1pi0_stack, leg_entries, norms, "Ev Absolute Normalisation, NC1#pi^{0}", "Ev_NC1Pi0_Plot.png", "E_{#nu}", "Number of SBND events");

    delete h1_Ev_abs_nc1pi0;
    delete h2_Ev_abs_nc1pi0;
    delete h3_Ev_abs_nc1pi0;
    delete h4_Ev_abs_nc1pi0;
    delete h5_Ev_abs_nc1pi0;
    
    // -------------------------------------------------------------------------
    //                          Ev Area norm CC0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_area_cc0pi_stack;
    vector< double > int_Ev_cc0pi;

    // The Ev histograms
    TH1D *h1_Ev_area_cc0pi = new TH1D("h1_Ev_area_cc0pi","E_{#nu} Area Normalised, CC0#pi",50,0,3);
    TH1D *h2_Ev_area_cc0pi = new TH1D("h2_Ev_area_cc0pi","E_{#nu} Area Normalised, CC0#pi",50,0,3);
    TH1D *h3_Ev_area_cc0pi = new TH1D("h3_Ev_area_cc0pi","E_{#nu} Area Normalised, CC0#pi",50,0,3);
    TH1D *h4_Ev_area_cc0pi = new TH1D("h4_Ev_area_cc0pi","E_{#nu} Area Normalised, CC0#pi",50,0,3);
    TH1D *h5_Ev_area_cc0pi = new TH1D("h5_Ev_area_cc0pi","E_{#nu} Area Normalised, CC0#pi",50,0,3);

    gst1->Draw("Ev>>h1_Ev_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst2->Draw("Ev>>h2_Ev_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst3->Draw("Ev>>h3_Ev_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst4->Draw("Ev>>h4_Ev_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    gst5->Draw("Ev>>h5_Ev_area_cc0pi", "cc && (nfpip + nfpim + nfpi0 == 0) ");
    
    // Push back the scaling
    double int1_Ev_cc0pi, int2_Ev_cc0pi, int3_Ev_cc0pi, int4_Ev_cc0pi, int5_Ev_cc0pi;

    int1_Ev_cc0pi = h1_Ev_area_cc0pi->Integral();
    int2_Ev_cc0pi = h2_Ev_area_cc0pi->Integral();
    int3_Ev_cc0pi = h3_Ev_area_cc0pi->Integral();
    int4_Ev_cc0pi = h4_Ev_area_cc0pi->Integral();
    int5_Ev_cc0pi = h5_Ev_area_cc0pi->Integral();
    
    int_Ev_cc0pi.push_back( 1 / int1_Ev_cc0pi );
    int_Ev_cc0pi.push_back( 1 / int2_Ev_cc0pi );
    int_Ev_cc0pi.push_back( 1 / int3_Ev_cc0pi );
    int_Ev_cc0pi.push_back( 1 / int4_Ev_cc0pi );
    int_Ev_cc0pi.push_back( 1 / int5_Ev_cc0pi );
    
    // Push them onto the vector
    Ev_area_cc0pi_stack.push_back(h1_Ev_area_cc0pi);
    Ev_area_cc0pi_stack.push_back(h2_Ev_area_cc0pi);
    Ev_area_cc0pi_stack.push_back(h3_Ev_area_cc0pi);
    Ev_area_cc0pi_stack.push_back(h4_Ev_area_cc0pi);
    Ev_area_cc0pi_stack.push_back(h5_Ev_area_cc0pi);

    HistStacker(Ev_area_cc0pi_stack, leg_entries, int_Ev_cc0pi, "Ev Area Normalised, CC0#pi", "Ev_CC0Pi_Area_Plot.png", "E_{#nu}", "SBND events, area normalised");

    delete h1_Ev_area_cc0pi;
    delete h2_Ev_area_cc0pi;
    delete h3_Ev_area_cc0pi;
    delete h4_Ev_area_cc0pi;
    delete h5_Ev_area_cc0pi;
    
    // -------------------------------------------------------------------------
    //                          Ev Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_area_cc1pip_stack;
    vector< double > int_Ev_cc1pip;

    // The Ev histograms
    TH1D *h1_Ev_area_cc1pip = new TH1D("h1_Ev_area_cc1pip","E_{#nu} Area Normalised, CC1#pi^{+}",50,0,3);
    TH1D *h2_Ev_area_cc1pip = new TH1D("h2_Ev_area_cc1pip","E_{#nu} Area Normalised, CC1#pi^{+}",50,0,3);
    TH1D *h3_Ev_area_cc1pip = new TH1D("h3_Ev_area_cc1pip","E_{#nu} Area Normalised, CC1#pi^{+}",50,0,3);
    TH1D *h4_Ev_area_cc1pip = new TH1D("h4_Ev_area_cc1pip","E_{#nu} Area Normalised, CC1#pi^{+}",50,0,3);
    TH1D *h5_Ev_area_cc1pip = new TH1D("h5_Ev_area_cc1pip","E_{#nu} Area Normalised, CC1#pi^{+}",50,0,3);

    gst1->Draw("Ev>>h1_Ev_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) ");
     
    // Push back the scaling
    double int1_Ev_cc1pip, int2_Ev_cc1pip, int3_Ev_cc1pip, int4_Ev_cc1pip, int5_Ev_cc1pip;

    int1_Ev_cc1pip = h1_Ev_area_cc1pip->Integral();
    int2_Ev_cc1pip = h2_Ev_area_cc1pip->Integral();
    int3_Ev_cc1pip = h3_Ev_area_cc1pip->Integral();
    int4_Ev_cc1pip = h4_Ev_area_cc1pip->Integral();
    int5_Ev_cc1pip = h5_Ev_area_cc1pip->Integral();
    
    int_Ev_cc1pip.push_back( 1 / int1_Ev_cc1pip );
    int_Ev_cc1pip.push_back( 1 / int2_Ev_cc1pip );
    int_Ev_cc1pip.push_back( 1 / int3_Ev_cc1pip );
    int_Ev_cc1pip.push_back( 1 / int4_Ev_cc1pip );
    int_Ev_cc1pip.push_back( 1 / int5_Ev_cc1pip );
    
    // Push them onto the vector
    Ev_area_cc1pip_stack.push_back(h1_Ev_area_cc1pip);
    Ev_area_cc1pip_stack.push_back(h2_Ev_area_cc1pip);
    Ev_area_cc1pip_stack.push_back(h3_Ev_area_cc1pip);
    Ev_area_cc1pip_stack.push_back(h4_Ev_area_cc1pip);
    Ev_area_cc1pip_stack.push_back(h5_Ev_area_cc1pip);

    HistStacker(Ev_area_cc1pip_stack, leg_entries, int_Ev_cc1pip, "Ev Area Normalised, CC1#pi^{+}", "Ev_CC1PiP_Area_Plot.png", "E_{#nu}", "SBND events, area normalised");

    delete h1_Ev_area_cc1pip;
    delete h2_Ev_area_cc1pip;
    delete h3_Ev_area_cc1pip;
    delete h4_Ev_area_cc1pip;
    delete h5_Ev_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          Ev Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_area_cc1pi0_stack;
    vector< double > int_Ev_cc1pi0;

    // The Ev histograms
    TH1D *h1_Ev_area_cc1pi0 = new TH1D("h1_Ev_area_cc1pi0","E_{#nu} Area Normalised, CC1#pi^{0}",50,0,3);
    TH1D *h2_Ev_area_cc1pi0 = new TH1D("h2_Ev_area_cc1pi0","E_{#nu} Area Normalised, CC1#pi^{0}",50,0,3);
    TH1D *h3_Ev_area_cc1pi0 = new TH1D("h3_Ev_area_cc1pi0","E_{#nu} Area Normalised, CC1#pi^{0}",50,0,3);
    TH1D *h4_Ev_area_cc1pi0 = new TH1D("h4_Ev_area_cc1pi0","E_{#nu} Area Normalised, CC1#pi^{0}",50,0,3);
    TH1D *h5_Ev_area_cc1pi0 = new TH1D("h5_Ev_area_cc1pi0","E_{#nu} Area Normalised, CC1#pi^{0}",50,0,3);
    
    gst1->Draw("Ev>>h1_Ev_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
     
    // Push back the scaling
    double int1_Ev_cc1pi0, int2_Ev_cc1pi0, int3_Ev_cc1pi0, int4_Ev_cc1pi0, int5_Ev_cc1pi0;

    int1_Ev_cc1pi0 = h1_Ev_area_cc1pi0->Integral();
    int2_Ev_cc1pi0 = h2_Ev_area_cc1pi0->Integral();
    int3_Ev_cc1pi0 = h3_Ev_area_cc1pi0->Integral();
    int4_Ev_cc1pi0 = h4_Ev_area_cc1pi0->Integral();
    int5_Ev_cc1pi0 = h5_Ev_area_cc1pi0->Integral();
    
    int_Ev_cc1pi0.push_back( 1 / int1_Ev_cc1pi0 );
    int_Ev_cc1pi0.push_back( 1 / int2_Ev_cc1pi0 );
    int_Ev_cc1pi0.push_back( 1 / int3_Ev_cc1pi0 );
    int_Ev_cc1pi0.push_back( 1 / int4_Ev_cc1pi0 );
    int_Ev_cc1pi0.push_back( 1 / int5_Ev_cc1pi0 );
    
    // Push them onto the vector
    Ev_area_cc1pi0_stack.push_back(h1_Ev_area_cc1pi0);
    Ev_area_cc1pi0_stack.push_back(h2_Ev_area_cc1pi0);
    Ev_area_cc1pi0_stack.push_back(h3_Ev_area_cc1pi0);
    Ev_area_cc1pi0_stack.push_back(h4_Ev_area_cc1pi0);
    Ev_area_cc1pi0_stack.push_back(h5_Ev_area_cc1pi0);

    HistStacker(Ev_area_cc1pi0_stack, leg_entries, int_Ev_cc1pi0, "Ev Area Normalised, CC1#pi^{0}", "Ev_CC1Pi0_Area_Plot.png", "E_{#nu}", "SBND events, area normalised");

    delete h1_Ev_area_cc1pi0;
    delete h2_Ev_area_cc1pi0;
    delete h3_Ev_area_cc1pi0;
    delete h4_Ev_area_cc1pi0;
    delete h5_Ev_area_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          Ev Area norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > Ev_area_nc1pi0_stack;
    vector< double > int_Ev_nc1pi0;

    // The Ev histograms
    TH1D *h1_Ev_area_nc1pi0 = new TH1D("h1_Ev_area_nc1pi0","E_{#nu} Area Normalised, NC1#pi^{0}",50,0,3);
    TH1D *h2_Ev_area_nc1pi0 = new TH1D("h2_Ev_area_nc1pi0","E_{#nu} Area Normalised, NC1#pi^{0}",50,0,3);
    TH1D *h3_Ev_area_nc1pi0 = new TH1D("h3_Ev_area_nc1pi0","E_{#nu} Area Normalised, NC1#pi^{0}",50,0,3);
    TH1D *h4_Ev_area_nc1pi0 = new TH1D("h4_Ev_area_nc1pi0","E_{#nu} Area Normalised, NC1#pi^{0}",50,0,3);
    TH1D *h5_Ev_area_nc1pi0 = new TH1D("h5_Ev_area_nc1pi0","E_{#nu} Area Normalised, NC1#pi^{0}",50,0,3);

    gst1->Draw("Ev>>h1_Ev_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst2->Draw("Ev>>h2_Ev_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst3->Draw("Ev>>h3_Ev_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst4->Draw("Ev>>h4_Ev_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    gst5->Draw("Ev>>h5_Ev_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) ");
    
    // Push back the scaling
    double int1_Ev_nc1pi0, int2_Ev_nc1pi0, int3_Ev_nc1pi0, int4_Ev_nc1pi0, int5_Ev_nc1pi0;

    int1_Ev_nc1pi0 = h1_Ev_area_nc1pi0->Integral();
    int2_Ev_nc1pi0 = h2_Ev_area_nc1pi0->Integral();
    int3_Ev_nc1pi0 = h3_Ev_area_nc1pi0->Integral();
    int4_Ev_nc1pi0 = h4_Ev_area_nc1pi0->Integral();
    int5_Ev_nc1pi0 = h5_Ev_area_nc1pi0->Integral();
    
    int_Ev_nc1pi0.push_back( 1 / int1_Ev_nc1pi0 );
    int_Ev_nc1pi0.push_back( 1 / int2_Ev_nc1pi0 );
    int_Ev_nc1pi0.push_back( 1 / int3_Ev_nc1pi0 );
    int_Ev_nc1pi0.push_back( 1 / int4_Ev_nc1pi0 );
    int_Ev_nc1pi0.push_back( 1 / int5_Ev_nc1pi0 );
    
    // Push them onto the vector
    Ev_area_nc1pi0_stack.push_back(h1_Ev_area_nc1pi0);
    Ev_area_nc1pi0_stack.push_back(h2_Ev_area_nc1pi0);
    Ev_area_nc1pi0_stack.push_back(h3_Ev_area_nc1pi0);
    Ev_area_nc1pi0_stack.push_back(h4_Ev_area_nc1pi0);
    Ev_area_nc1pi0_stack.push_back(h5_Ev_area_nc1pi0);

    HistStacker(Ev_area_nc1pi0_stack, leg_entries, int_Ev_nc1pi0, "Ev Area Normalised, NC1#pi^{0}", "Ev_NC1Pi0_Area_Plot.png", "E_{#nu}", "SBND events, area normalised");

    delete h1_Ev_area_nc1pi0;
    delete h2_Ev_area_nc1pi0;
    delete h3_Ev_area_nc1pi0;
    delete h4_Ev_area_nc1pi0;
    delete h5_Ev_area_nc1pi0;
    // -------------------------------------------------------------------------
    //                          CosThetaPi Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_abs_cc1pip_stack;

    // The Ev histograms
    TH1D *h1_CosThPi_abs_cc1pip = new TH1D("h1_CosThPi_abs_cc1pip","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h2_CosThPi_abs_cc1pip = new TH1D("h2_CosThPi_abs_cc1pip","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h3_CosThPi_abs_cc1pip = new TH1D("h3_CosThPi_abs_cc1pip","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h4_CosThPi_abs_cc1pip = new TH1D("h4_CosThPi_abs_cc1pip","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h5_CosThPi_abs_cc1pip = new TH1D("h5_CosThPi_abs_cc1pip","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);

    gst1->Draw("cthf>>h1_CosThPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst2->Draw("cthf>>h2_CosThPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst3->Draw("cthf>>h3_CosThPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst4->Draw("cthf>>h4_CosThPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst5->Draw("cthf>>h5_CosThPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    
    // Push them onto the vector
    CosThPi_abs_cc1pip_stack.push_back(h1_CosThPi_abs_cc1pip);
    CosThPi_abs_cc1pip_stack.push_back(h2_CosThPi_abs_cc1pip);
    CosThPi_abs_cc1pip_stack.push_back(h3_CosThPi_abs_cc1pip);
    CosThPi_abs_cc1pip_stack.push_back(h4_CosThPi_abs_cc1pip);
    CosThPi_abs_cc1pip_stack.push_back(h5_CosThPi_abs_cc1pip);

    HistStacker(CosThPi_abs_cc1pip_stack, leg_entries, norms, "CosThetaPi Absolute Normalisation, CC1#pi^{+}", "CosThPi_CC1PiP_Plot.png", "cos(#theta_{#pi})", "Number of SBND events");

    delete h1_CosThPi_abs_cc1pip;
    delete h2_CosThPi_abs_cc1pip;
    delete h3_CosThPi_abs_cc1pip;
    delete h4_CosThPi_abs_cc1pip;
    delete h5_CosThPi_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          CosThetaPi Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_abs_cc1pi0_stack;

    // The Ev histograms
    TH1D *h1_CosThPi_abs_cc1pi0 = new TH1D("h1_CosThPi_abs_cc1pi0","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThPi_abs_cc1pi0 = new TH1D("h2_CosThPi_abs_cc1pi0","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThPi_abs_cc1pi0 = new TH1D("h3_CosThPi_abs_cc1pi0","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThPi_abs_cc1pi0 = new TH1D("h4_CosThPi_abs_cc1pi0","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThPi_abs_cc1pi0 = new TH1D("h5_CosThPi_abs_cc1pi0","cos(#theta_{#pi}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);

    gst1->Draw("cthf>>h1_CosThPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("cthf>>h2_CosThPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("cthf>>h3_CosThPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("cthf>>h4_CosThPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("cthf>>h5_CosThPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    
    // Push them onto the vector
    CosThPi_abs_cc1pi0_stack.push_back(h1_CosThPi_abs_cc1pi0);
    CosThPi_abs_cc1pi0_stack.push_back(h2_CosThPi_abs_cc1pi0);
    CosThPi_abs_cc1pi0_stack.push_back(h3_CosThPi_abs_cc1pi0);
    CosThPi_abs_cc1pi0_stack.push_back(h4_CosThPi_abs_cc1pi0);
    CosThPi_abs_cc1pi0_stack.push_back(h5_CosThPi_abs_cc1pi0);

    HistStacker(CosThPi_abs_cc1pi0_stack, leg_entries, norms, "CosThetaPi Absolute Normalisation, CC1#pi^{0}", "CosThPi_CC1Pi0_Plot.png", "cos(#theta_{#pi})", "Number of SBND events");

    delete h1_CosThPi_abs_cc1pi0;
    delete h2_CosThPi_abs_cc1pi0;
    delete h3_CosThPi_abs_cc1pi0;
    delete h4_CosThPi_abs_cc1pi0;
    delete h5_CosThPi_abs_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          CosThetaPi Absolute norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_abs_nc1pi0_stack;

    // The Ev histograms
    TH1D *h1_CosThPi_abs_nc1pi0 = new TH1D("h1_CosThPi_abs_nc1pi0","cos(#theta_{#pi}) Absolute Normalisation, NC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThPi_abs_nc1pi0 = new TH1D("h2_CosThPi_abs_nc1pi0","cos(#theta_{#pi}) Absolute Normalisation, NC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThPi_abs_nc1pi0 = new TH1D("h3_CosThPi_abs_nc1pi0","cos(#theta_{#pi}) Absolute Normalisation, NC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThPi_abs_nc1pi0 = new TH1D("h4_CosThPi_abs_nc1pi0","cos(#theta_{#pi}) Absolute Normalisation, NC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThPi_abs_nc1pi0 = new TH1D("h5_CosThPi_abs_nc1pi0","cos(#theta_{#pi}) Absolute Normalisation, NC1#pi^{0}",40,-1,1);

    gst1->Draw("cthf>>h1_CosThPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("cthf>>h2_CosThPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("cthf>>h3_CosThPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("cthf>>h4_CosThPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("cthf>>h5_CosThPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    
    // Push them onto the vector
    CosThPi_abs_nc1pi0_stack.push_back(h1_CosThPi_abs_nc1pi0);
    CosThPi_abs_nc1pi0_stack.push_back(h2_CosThPi_abs_nc1pi0);
    CosThPi_abs_nc1pi0_stack.push_back(h3_CosThPi_abs_nc1pi0);
    CosThPi_abs_nc1pi0_stack.push_back(h4_CosThPi_abs_nc1pi0);
    CosThPi_abs_nc1pi0_stack.push_back(h5_CosThPi_abs_nc1pi0);

    HistStacker(CosThPi_abs_nc1pi0_stack, leg_entries, norms, "CosThetaPi Absolute Normalisation, NC1#pi^{0}", "CosThPi_NC1Pi0_Plot.png", "cos(#theta_{#pi})", "Number of SBND events");

    delete h1_CosThPi_abs_nc1pi0;
    delete h2_CosThPi_abs_nc1pi0;
    delete h3_CosThPi_abs_nc1pi0;
    delete h4_CosThPi_abs_nc1pi0;
    delete h5_CosThPi_abs_nc1pi0;
    
    // -------------------------------------------------------------------------
    //                      CosThetaPi Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_area_cc1pip_stack;
    vector< double > int_cth_pi_cc1pip;

    // The Ev histograms
    TH1D *h1_CosThPi_area_cc1pip = new TH1D("h1_CosThPi_area_cc1pip","cos(#theta_{#pi}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h2_CosThPi_area_cc1pip = new TH1D("h2_CosThPi_area_cc1pip","cos(#theta_{#pi}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h3_CosThPi_area_cc1pip = new TH1D("h3_CosThPi_area_cc1pip","cos(#theta_{#pi}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h4_CosThPi_area_cc1pip = new TH1D("h4_CosThPi_area_cc1pip","cos(#theta_{#pi}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h5_CosThPi_area_cc1pip = new TH1D("h5_CosThPi_area_cc1pip","cos(#theta_{#pi}) Area Normalised, CC1#pi^{+}",40,-1,1);

    gst1->Draw("cthf>>h1_CosThPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211 ");
    gst2->Draw("cthf>>h2_CosThPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst3->Draw("cthf>>h3_CosThPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst4->Draw("cthf>>h4_CosThPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst5->Draw("cthf>>h5_CosThPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
     
    // Push back the scaling
    double int1_cth_pi_cc1pip, int2_cth_pi_cc1pip, int3_cth_pi_cc1pip, int4_cth_pi_cc1pip, int5_cth_pi_cc1pip;

    int1_cth_pi_cc1pip = h1_CosThPi_area_cc1pip->Integral();
    int2_cth_pi_cc1pip = h2_CosThPi_area_cc1pip->Integral();
    int3_cth_pi_cc1pip = h3_CosThPi_area_cc1pip->Integral();
    int4_cth_pi_cc1pip = h4_CosThPi_area_cc1pip->Integral();
    int5_cth_pi_cc1pip = h5_CosThPi_area_cc1pip->Integral();
    
    int_cth_pi_cc1pip.push_back( 1 / int1_cth_pi_cc1pip );
    int_cth_pi_cc1pip.push_back( 1 / int2_cth_pi_cc1pip );
    int_cth_pi_cc1pip.push_back( 1 / int3_cth_pi_cc1pip );
    int_cth_pi_cc1pip.push_back( 1 / int4_cth_pi_cc1pip );
    int_cth_pi_cc1pip.push_back( 1 / int5_cth_pi_cc1pip );
    
    // Push them onto the vector
    CosThPi_area_cc1pip_stack.push_back(h1_CosThPi_area_cc1pip);
    CosThPi_area_cc1pip_stack.push_back(h2_CosThPi_area_cc1pip);
    CosThPi_area_cc1pip_stack.push_back(h3_CosThPi_area_cc1pip);
    CosThPi_area_cc1pip_stack.push_back(h4_CosThPi_area_cc1pip);
    CosThPi_area_cc1pip_stack.push_back(h5_CosThPi_area_cc1pip);

    HistStacker(CosThPi_area_cc1pip_stack, leg_entries, int_cth_pi_cc1pip, "CosThetaPi Area Normalised, CC1#pi^{+}", "CosThPi_CC1PiP_Area_Plot.png", "cos(#theta_{#pi})", "SBND events, area normalised");

    delete h1_CosThPi_area_cc1pip;
    delete h2_CosThPi_area_cc1pip;
    delete h3_CosThPi_area_cc1pip;
    delete h4_CosThPi_area_cc1pip;
    delete h5_CosThPi_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          CosThetaPi Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_area_cc1pi0_stack;
    vector< double > int_cth_pi_cc1pi0;

    // The Ev histograms
    TH1D *h1_CosThPi_area_cc1pi0 = new TH1D("h1_CosThPi_area_cc1pi0","cos(#theta_{#pi}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThPi_area_cc1pi0 = new TH1D("h2_CosThPi_area_cc1pi0","cos(#theta_{#pi}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThPi_area_cc1pi0 = new TH1D("h3_CosThPi_area_cc1pi0","cos(#theta_{#pi}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThPi_area_cc1pi0 = new TH1D("h4_CosThPi_area_cc1pi0","cos(#theta_{#pi}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThPi_area_cc1pi0 = new TH1D("h5_CosThPi_area_cc1pi0","cos(#theta_{#pi}) Area Normalised, CC1#pi^{0}",40,-1,1);
    
    gst1->Draw("cthf>>h1_CosThPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("cthf>>h2_CosThPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("cthf>>h3_CosThPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("cthf>>h4_CosThPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("cthf>>h5_CosThPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
     
    // Push back the scaling
    double int1_cth_pi_cc1pi0, int2_cth_pi_cc1pi0, int3_cth_pi_cc1pi0, int4_cth_pi_cc1pi0, int5_cth_pi_cc1pi0;

    int1_cth_pi_cc1pi0 = h1_CosThPi_area_cc1pi0->Integral();
    int2_cth_pi_cc1pi0 = h2_CosThPi_area_cc1pi0->Integral();
    int3_cth_pi_cc1pi0 = h3_CosThPi_area_cc1pi0->Integral();
    int4_cth_pi_cc1pi0 = h4_CosThPi_area_cc1pi0->Integral();
    int5_cth_pi_cc1pi0 = h5_CosThPi_area_cc1pi0->Integral();
    
    int_cth_pi_cc1pi0.push_back( 1 / int1_cth_pi_cc1pi0 );
    int_cth_pi_cc1pi0.push_back( 1 / int2_cth_pi_cc1pi0 );
    int_cth_pi_cc1pi0.push_back( 1 / int3_cth_pi_cc1pi0 );
    int_cth_pi_cc1pi0.push_back( 1 / int4_cth_pi_cc1pi0 );
    int_cth_pi_cc1pi0.push_back( 1 / int5_cth_pi_cc1pi0 );
    
    // Push them onto the vector
    CosThPi_area_cc1pi0_stack.push_back(h1_CosThPi_area_cc1pi0);
    CosThPi_area_cc1pi0_stack.push_back(h2_CosThPi_area_cc1pi0);
    CosThPi_area_cc1pi0_stack.push_back(h3_CosThPi_area_cc1pi0);
    CosThPi_area_cc1pi0_stack.push_back(h4_CosThPi_area_cc1pi0);
    CosThPi_area_cc1pi0_stack.push_back(h5_CosThPi_area_cc1pi0);

    HistStacker(CosThPi_area_cc1pi0_stack, leg_entries, int_cth_pi_cc1pi0, "CosThetaPi Area Normalised, CC1#pi^{0}", "CosThPi_CC1Pi0_Area_Plot.png", "cos(#theta_{#pi})", "SBND events, area normalised");

    delete h1_CosThPi_area_cc1pi0;
    delete h2_CosThPi_area_cc1pi0;
    delete h3_CosThPi_area_cc1pi0;
    delete h4_CosThPi_area_cc1pi0;
    delete h5_CosThPi_area_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          CosTheta Area norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThPi_area_nc1pi0_stack;
    vector< double > int_cth_pi_nc1pi0;

    // The Ev histograms
    TH1D *h1_CosThPi_area_nc1pi0 = new TH1D("h1_CosThPi_area_nc1pi0","cos(#theta_{#pi}) Area Normalised, NC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThPi_area_nc1pi0 = new TH1D("h2_CosThPi_area_nc1pi0","cos(#theta_{#pi}) Area Normalised, NC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThPi_area_nc1pi0 = new TH1D("h3_CosThPi_area_nc1pi0","cos(#theta_{#pi}) Area Normalised, NC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThPi_area_nc1pi0 = new TH1D("h4_CosThPi_area_nc1pi0","cos(#theta_{#pi}) Area Normalised, NC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThPi_area_nc1pi0 = new TH1D("h5_CosThPi_area_nc1pi0","cos(#theta_{#pi}) Area Normalised, NC1#pi^{0}",40,-1,1);

    gst1->Draw("cthf>>h1_CosThPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst2->Draw("cthf>>h2_CosThPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst3->Draw("cthf>>h3_CosThPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst4->Draw("cthf>>h4_CosThPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst5->Draw("cthf>>h5_CosThPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    
    // Push back the scaling
    double int1_cth_pi_nc1pi0, int2_cth_pi_nc1pi0, int3_cth_pi_nc1pi0, int4_cth_pi_nc1pi0, int5_cth_pi_nc1pi0;

    int1_cth_pi_nc1pi0 = h1_CosThPi_area_nc1pi0->Integral();
    int2_cth_pi_nc1pi0 = h2_CosThPi_area_nc1pi0->Integral();
    int3_cth_pi_nc1pi0 = h3_CosThPi_area_nc1pi0->Integral();
    int4_cth_pi_nc1pi0 = h4_CosThPi_area_nc1pi0->Integral();
    int5_cth_pi_nc1pi0 = h5_CosThPi_area_nc1pi0->Integral();
    
    int_cth_pi_nc1pi0.push_back( 1 / int1_cth_pi_nc1pi0 );
    int_cth_pi_nc1pi0.push_back( 1 / int2_cth_pi_nc1pi0 );
    int_cth_pi_nc1pi0.push_back( 1 / int3_cth_pi_nc1pi0 );
    int_cth_pi_nc1pi0.push_back( 1 / int4_cth_pi_nc1pi0 );
    int_cth_pi_nc1pi0.push_back( 1 / int5_cth_pi_nc1pi0 );
    
    // Push them onto the vector
    CosThPi_area_nc1pi0_stack.push_back(h1_CosThPi_area_nc1pi0);
    CosThPi_area_nc1pi0_stack.push_back(h2_CosThPi_area_nc1pi0);
    CosThPi_area_nc1pi0_stack.push_back(h3_CosThPi_area_nc1pi0);
    CosThPi_area_nc1pi0_stack.push_back(h4_CosThPi_area_nc1pi0);
    CosThPi_area_nc1pi0_stack.push_back(h5_CosThPi_area_nc1pi0);

    HistStacker(CosThPi_area_nc1pi0_stack, leg_entries, int_cth_pi_nc1pi0, "CosThetaPi Area Normalised, NC1#pi^{0}", "CosThPi_NC1Pi0_Area_Plot.png", "cos(#theta_{#pi})", "SBND events, area normalised");

    delete h1_CosThPi_area_nc1pi0;
    delete h2_CosThPi_area_nc1pi0;
    delete h3_CosThPi_area_nc1pi0;
    delete h4_CosThPi_area_nc1pi0;
    delete h5_CosThPi_area_nc1pi0;
    
    // -------------------------------------------------------------------------
    //                          PPi Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_abs_cc1pip_stack;

    // The Ev histograms
    TH1D *h1_PPi_abs_cc1pip = new TH1D("h1_PPi_abs_cc1pip","P_{#pi} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h2_PPi_abs_cc1pip = new TH1D("h2_PPi_abs_cc1pip","P_{#pi} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h3_PPi_abs_cc1pip = new TH1D("h3_PPi_abs_cc1pip","P_{#pi} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h4_PPi_abs_cc1pip = new TH1D("h4_PPi_abs_cc1pip","P_{#pi} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h5_PPi_abs_cc1pip = new TH1D("h5_PPi_abs_cc1pip","P_{#pi} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);

    gst1->Draw("pf>>h1_PPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst2->Draw("pf>>h2_PPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst3->Draw("pf>>h3_PPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst4->Draw("pf>>h4_PPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst5->Draw("pf>>h5_PPi_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    
    // Push them onto the vector
    PPi_abs_cc1pip_stack.push_back(h1_PPi_abs_cc1pip);
    PPi_abs_cc1pip_stack.push_back(h2_PPi_abs_cc1pip);
    PPi_abs_cc1pip_stack.push_back(h3_PPi_abs_cc1pip);
    PPi_abs_cc1pip_stack.push_back(h4_PPi_abs_cc1pip);
    PPi_abs_cc1pip_stack.push_back(h5_PPi_abs_cc1pip);

    HistStacker(PPi_abs_cc1pip_stack, leg_entries, norms, "PPi Absolute Normalisation, CC1#pi^{+}", "PPi_CC1PiP_Plot.png", "P_{#pi}", "Number of SBND events");

    delete h1_PPi_abs_cc1pip;
    delete h2_PPi_abs_cc1pip;
    delete h3_PPi_abs_cc1pip;
    delete h4_PPi_abs_cc1pip;
    delete h5_PPi_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          PPi Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_abs_cc1pi0_stack;

    // The Ev histograms
    TH1D *h1_PPi_abs_cc1pi0 = new TH1D("h1_PPi_abs_cc1pi0","P_{#pi} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h2_PPi_abs_cc1pi0 = new TH1D("h2_PPi_abs_cc1pi0","P_{#pi} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h3_PPi_abs_cc1pi0 = new TH1D("h3_PPi_abs_cc1pi0","P_{#pi} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h4_PPi_abs_cc1pi0 = new TH1D("h4_PPi_abs_cc1pi0","P_{#pi} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h5_PPi_abs_cc1pi0 = new TH1D("h5_PPi_abs_cc1pi0","P_{#pi} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);

    gst1->Draw("pf>>h1_PPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("pf>>h2_PPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("pf>>h3_PPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("pf>>h4_PPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("pf>>h5_PPi_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    
    // Push them onto the vector
    PPi_abs_cc1pi0_stack.push_back(h1_PPi_abs_cc1pi0);
    PPi_abs_cc1pi0_stack.push_back(h2_PPi_abs_cc1pi0);
    PPi_abs_cc1pi0_stack.push_back(h3_PPi_abs_cc1pi0);
    PPi_abs_cc1pi0_stack.push_back(h4_PPi_abs_cc1pi0);
    PPi_abs_cc1pi0_stack.push_back(h5_PPi_abs_cc1pi0);

    HistStacker(PPi_abs_cc1pi0_stack, leg_entries, norms, "PPi Absolute Normalisation, CC1#pi^{0}", "PPi_CC1Pi0_Plot.png", "P_{#pi}", "Number of SBND events");

    delete h1_PPi_abs_cc1pi0;
    delete h2_PPi_abs_cc1pi0;
    delete h3_PPi_abs_cc1pi0;
    delete h4_PPi_abs_cc1pi0;
    delete h5_PPi_abs_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          PPi Absolute norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_abs_nc1pi0_stack;

    // The Ev histograms
    TH1D *h1_PPi_abs_nc1pi0 = new TH1D("h1_PPi_abs_nc1pi0","P_{#pi} Absolute Normalisation, NC1#pi^{0}",40,0,1.5);
    TH1D *h2_PPi_abs_nc1pi0 = new TH1D("h2_PPi_abs_nc1pi0","P_{#pi} Absolute Normalisation, NC1#pi^{0}",40,0,1.5);
    TH1D *h3_PPi_abs_nc1pi0 = new TH1D("h3_PPi_abs_nc1pi0","P_{#pi} Absolute Normalisation, NC1#pi^{0}",40,0,1.5);
    TH1D *h4_PPi_abs_nc1pi0 = new TH1D("h4_PPi_abs_nc1pi0","P_{#pi} Absolute Normalisation, NC1#pi^{0}",40,0,1.5);
    TH1D *h5_PPi_abs_nc1pi0 = new TH1D("h5_PPi_abs_nc1pi0","P_{#pi} Absolute Normalisation, NC1#pi^{0}",40,0,1.5);

    gst1->Draw("pf>>h1_PPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("pf>>h2_PPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("pf>>h3_PPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("pf>>h4_PPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("pf>>h5_PPi_abs_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    
    // Push them onto the vector
    PPi_abs_nc1pi0_stack.push_back(h1_PPi_abs_nc1pi0);
    PPi_abs_nc1pi0_stack.push_back(h2_PPi_abs_nc1pi0);
    PPi_abs_nc1pi0_stack.push_back(h3_PPi_abs_nc1pi0);
    PPi_abs_nc1pi0_stack.push_back(h4_PPi_abs_nc1pi0);
    PPi_abs_nc1pi0_stack.push_back(h5_PPi_abs_nc1pi0);

    HistStacker(PPi_abs_nc1pi0_stack, leg_entries, norms, "PPi Absolute Normalisation, NC1#pi^{0}", "PPi_NC1Pi0_Plot.png", "P_{#pi}", "Number of SBND events");

    delete h1_PPi_abs_nc1pi0;
    delete h2_PPi_abs_nc1pi0;
    delete h3_PPi_abs_nc1pi0;
    delete h4_PPi_abs_nc1pi0;
    delete h5_PPi_abs_nc1pi0;
    
    // -------------------------------------------------------------------------
    //                      PPi Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_area_cc1pip_stack;
    vector< double > int_p_pi_cc1pip;

    // The Ev histograms
    TH1D *h1_PPi_area_cc1pip = new TH1D("h1_PPi_area_cc1pip","P_{#pi} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h2_PPi_area_cc1pip = new TH1D("h2_PPi_area_cc1pip","P_{#pi} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h3_PPi_area_cc1pip = new TH1D("h3_PPi_area_cc1pip","P_{#pi} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h4_PPi_area_cc1pip = new TH1D("h4_PPi_area_cc1pip","P_{#pi} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h5_PPi_area_cc1pip = new TH1D("h5_PPi_area_cc1pip","P_{#pi} Area Normalised, CC1#pi^{+}",40,0,1.5);

    gst1->Draw("pf>>h1_PPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211 ");
    gst2->Draw("pf>>h2_PPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst3->Draw("pf>>h3_PPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst4->Draw("pf>>h4_PPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
    gst5->Draw("pf>>h5_PPi_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 ) && pdgf == 211  ");
     
    // Push back the scaling
    double int1_p_pi_cc1pip, int2_p_pi_cc1pip, int3_p_pi_cc1pip, int4_p_pi_cc1pip, int5_p_pi_cc1pip;

    int1_p_pi_cc1pip = h1_PPi_area_cc1pip->Integral();
    int2_p_pi_cc1pip = h2_PPi_area_cc1pip->Integral();
    int3_p_pi_cc1pip = h3_PPi_area_cc1pip->Integral();
    int4_p_pi_cc1pip = h4_PPi_area_cc1pip->Integral();
    int5_p_pi_cc1pip = h5_PPi_area_cc1pip->Integral();
    
    int_p_pi_cc1pip.push_back( 1 / int1_p_pi_cc1pip );
    int_p_pi_cc1pip.push_back( 1 / int2_p_pi_cc1pip );
    int_p_pi_cc1pip.push_back( 1 / int3_p_pi_cc1pip );
    int_p_pi_cc1pip.push_back( 1 / int4_p_pi_cc1pip );
    int_p_pi_cc1pip.push_back( 1 / int5_p_pi_cc1pip );
    
    // Push them onto the vector
    PPi_area_cc1pip_stack.push_back(h1_PPi_area_cc1pip);
    PPi_area_cc1pip_stack.push_back(h2_PPi_area_cc1pip);
    PPi_area_cc1pip_stack.push_back(h3_PPi_area_cc1pip);
    PPi_area_cc1pip_stack.push_back(h4_PPi_area_cc1pip);
    PPi_area_cc1pip_stack.push_back(h5_PPi_area_cc1pip);

    HistStacker(PPi_area_cc1pip_stack, leg_entries, int_p_pi_cc1pip, "PPi Area Normalised, CC1#pi^{+}", "PPi_CC1PiP_Area_Plot.png", "P_{#pi}", "SBND events, area normalised");

    delete h1_PPi_area_cc1pip;
    delete h2_PPi_area_cc1pip;
    delete h3_PPi_area_cc1pip;
    delete h4_PPi_area_cc1pip;
    delete h5_PPi_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          PPi Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_area_cc1pi0_stack;
    vector< double > int_p_pi_cc1pi0;

    // The Ev histograms
    TH1D *h1_PPi_area_cc1pi0 = new TH1D("h1_PPi_area_cc1pi0","P_{#pi} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h2_PPi_area_cc1pi0 = new TH1D("h2_PPi_area_cc1pi0","P_{#pi} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h3_PPi_area_cc1pi0 = new TH1D("h3_PPi_area_cc1pi0","P_{#pi} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h4_PPi_area_cc1pi0 = new TH1D("h4_PPi_area_cc1pi0","P_{#pi} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h5_PPi_area_cc1pi0 = new TH1D("h5_PPi_area_cc1pi0","P_{#pi} Area Normalised, CC1#pi^{0}",40,0,1.5);
    
    gst1->Draw("pf>>h1_PPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst2->Draw("pf>>h2_PPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst3->Draw("pf>>h3_PPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst4->Draw("pf>>h4_PPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
    gst5->Draw("pf>>h5_PPi_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111 ");
     
    // Push back the scaling
    double int1_p_pi_cc1pi0, int2_p_pi_cc1pi0, int3_p_pi_cc1pi0, int4_p_pi_cc1pi0, int5_p_pi_cc1pi0;

    int1_p_pi_cc1pi0 = h1_PPi_area_cc1pi0->Integral();
    int2_p_pi_cc1pi0 = h2_PPi_area_cc1pi0->Integral();
    int3_p_pi_cc1pi0 = h3_PPi_area_cc1pi0->Integral();
    int4_p_pi_cc1pi0 = h4_PPi_area_cc1pi0->Integral();
    int5_p_pi_cc1pi0 = h5_PPi_area_cc1pi0->Integral();
    
    int_p_pi_cc1pi0.push_back( 1 / int1_p_pi_cc1pi0 );
    int_p_pi_cc1pi0.push_back( 1 / int2_p_pi_cc1pi0 );
    int_p_pi_cc1pi0.push_back( 1 / int3_p_pi_cc1pi0 );
    int_p_pi_cc1pi0.push_back( 1 / int4_p_pi_cc1pi0 );
    int_p_pi_cc1pi0.push_back( 1 / int5_p_pi_cc1pi0 );
    
    // Push them onto the vector
    PPi_area_cc1pi0_stack.push_back(h1_PPi_area_cc1pi0);
    PPi_area_cc1pi0_stack.push_back(h2_PPi_area_cc1pi0);
    PPi_area_cc1pi0_stack.push_back(h3_PPi_area_cc1pi0);
    PPi_area_cc1pi0_stack.push_back(h4_PPi_area_cc1pi0);
    PPi_area_cc1pi0_stack.push_back(h5_PPi_area_cc1pi0);

    HistStacker(PPi_area_cc1pi0_stack, leg_entries, int_p_pi_cc1pi0, "PPi Area Normalised, CC1#pi^{0}", "PPi_CC1Pi0_Area_Plot.png", "P_{#pi}", "SBND events, area normalised");

    delete h1_PPi_area_cc1pi0;
    delete h2_PPi_area_cc1pi0;
    delete h3_PPi_area_cc1pi0;
    delete h4_PPi_area_cc1pi0;
    delete h5_PPi_area_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          CosTheta Area norm NC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PPi_area_nc1pi0_stack;
    vector< double > int_p_pi_nc1pi0;

    // The Ev histograms
    TH1D *h1_PPi_area_nc1pi0 = new TH1D("h1_PPi_area_nc1pi0","P_{#pi} Area Normalised, NC1#pi^{0}",40,0,1.5);
    TH1D *h2_PPi_area_nc1pi0 = new TH1D("h2_PPi_area_nc1pi0","P_{#pi} Area Normalised, NC1#pi^{0}",40,0,1.5);
    TH1D *h3_PPi_area_nc1pi0 = new TH1D("h3_PPi_area_nc1pi0","P_{#pi} Area Normalised, NC1#pi^{0}",40,0,1.5);
    TH1D *h4_PPi_area_nc1pi0 = new TH1D("h4_PPi_area_nc1pi0","P_{#pi} Area Normalised, NC1#pi^{0}",40,0,1.5);
    TH1D *h5_PPi_area_nc1pi0 = new TH1D("h5_PPi_area_nc1pi0","P_{#pi} Area Normalised, NC1#pi^{0}",40,0,1.5);

    gst1->Draw("pf>>h1_PPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst2->Draw("pf>>h2_PPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst3->Draw("pf>>h3_PPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst4->Draw("pf>>h4_PPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    gst5->Draw("pf>>h5_PPi_area_nc1pi0", "nc && nfpi0 == 1 && ( nfpim + nfpip == 0 ) && pdgf == 111");
    
    // Push back the scaling
    double int1_p_pi_nc1pi0, int2_p_pi_nc1pi0, int3_p_pi_nc1pi0, int4_p_pi_nc1pi0, int5_p_pi_nc1pi0;

    int1_p_pi_nc1pi0 = h1_PPi_area_nc1pi0->Integral();
    int2_p_pi_nc1pi0 = h2_PPi_area_nc1pi0->Integral();
    int3_p_pi_nc1pi0 = h3_PPi_area_nc1pi0->Integral();
    int4_p_pi_nc1pi0 = h4_PPi_area_nc1pi0->Integral();
    int5_p_pi_nc1pi0 = h5_PPi_area_nc1pi0->Integral();
    
    int_p_pi_nc1pi0.push_back( 1 / int1_p_pi_nc1pi0 );
    int_p_pi_nc1pi0.push_back( 1 / int2_p_pi_nc1pi0 );
    int_p_pi_nc1pi0.push_back( 1 / int3_p_pi_nc1pi0 );
    int_p_pi_nc1pi0.push_back( 1 / int4_p_pi_nc1pi0 );
    int_p_pi_nc1pi0.push_back( 1 / int5_p_pi_nc1pi0 );
    
    // Push them onto the vector
    PPi_area_nc1pi0_stack.push_back(h1_PPi_area_nc1pi0);
    PPi_area_nc1pi0_stack.push_back(h2_PPi_area_nc1pi0);
    PPi_area_nc1pi0_stack.push_back(h3_PPi_area_nc1pi0);
    PPi_area_nc1pi0_stack.push_back(h4_PPi_area_nc1pi0);
    PPi_area_nc1pi0_stack.push_back(h5_PPi_area_nc1pi0);

    HistStacker(PPi_area_nc1pi0_stack, leg_entries, int_p_pi_nc1pi0, "PPi Area Normalised, NC1#pi^{0}", "PPi_NC1Pi0_Area_Plot.png", "P_{#pi}", "SBND events, area normalised");

    delete h1_PPi_area_nc1pi0;
    delete h2_PPi_area_nc1pi0;
    delete h3_PPi_area_nc1pi0;
    delete h4_PPi_area_nc1pi0;
    delete h5_PPi_area_nc1pi0;
    // -------------------------------------------------------------------------
    //                          CosThetaMu Absolute norm cc0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_abs_cc0pi_stack;

    // The Ev histograms
    TH1D *h1_CosThMu_abs_cc0pi = new TH1D("h1_CosThMu_abs_cc0pi","cos(#theta_{#mu}) Absolute Normalisation, CC0#pi",40,-1,1);
    TH1D *h2_CosThMu_abs_cc0pi = new TH1D("h2_CosThMu_abs_cc0pi","cos(#theta_{#mu}) Absolute Normalisation, CC0#pi",40,-1,1);
    TH1D *h3_CosThMu_abs_cc0pi = new TH1D("h3_CosThMu_abs_cc0pi","cos(#theta_{#mu}) Absolute Normalisation, CC0#pi",40,-1,1);
    TH1D *h4_CosThMu_abs_cc0pi = new TH1D("h4_CosThMu_abs_cc0pi","cos(#theta_{#mu}) Absolute Normalisation, CC0#pi",40,-1,1);
    TH1D *h5_CosThMu_abs_cc0pi = new TH1D("h5_CosThMu_abs_cc0pi","cos(#theta_{#mu}) Absolute Normalisation, CC0#pi",40,-1,1);

    gst1->Draw("cthl>>h1_CosThMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst2->Draw("cthl>>h2_CosThMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst3->Draw("cthl>>h3_CosThMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst4->Draw("cthl>>h4_CosThMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst5->Draw("cthl>>h5_CosThMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    
    // Push them onto the vector
    CosThMu_abs_cc0pi_stack.push_back(h1_CosThMu_abs_cc0pi);
    CosThMu_abs_cc0pi_stack.push_back(h2_CosThMu_abs_cc0pi);
    CosThMu_abs_cc0pi_stack.push_back(h3_CosThMu_abs_cc0pi);
    CosThMu_abs_cc0pi_stack.push_back(h4_CosThMu_abs_cc0pi);
    CosThMu_abs_cc0pi_stack.push_back(h5_CosThMu_abs_cc0pi);

    HistStacker(CosThMu_abs_cc0pi_stack, leg_entries, norms, "CosThetaMu Absolute Normalisation, CC0#pi", "CosThMu_cc0pi_Plot.png", "cos(#theta_{#mu})", "Number of SBND events");

    delete h1_CosThMu_abs_cc0pi;
    delete h2_CosThMu_abs_cc0pi;
    delete h3_CosThMu_abs_cc0pi;
    delete h4_CosThMu_abs_cc0pi;
    delete h5_CosThMu_abs_cc0pi;
    // -------------------------------------------------------------------------
    //                          CosThetaMu Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_abs_cc1pip_stack;

    // The Ev histograms
    TH1D *h1_CosThMu_abs_cc1pip = new TH1D("h1_CosThMu_abs_cc1pip","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h2_CosThMu_abs_cc1pip = new TH1D("h2_CosThMu_abs_cc1pip","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h3_CosThMu_abs_cc1pip = new TH1D("h3_CosThMu_abs_cc1pip","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h4_CosThMu_abs_cc1pip = new TH1D("h4_CosThMu_abs_cc1pip","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);
    TH1D *h5_CosThMu_abs_cc1pip = new TH1D("h5_CosThMu_abs_cc1pip","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{+}",40,-1,1);

    gst1->Draw("cthl>>h1_CosThMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst2->Draw("cthl>>h2_CosThMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst3->Draw("cthl>>h3_CosThMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst4->Draw("cthl>>h4_CosThMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst5->Draw("cthl>>h5_CosThMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    
    // Push them onto the vector
    CosThMu_abs_cc1pip_stack.push_back(h1_CosThMu_abs_cc1pip);
    CosThMu_abs_cc1pip_stack.push_back(h2_CosThMu_abs_cc1pip);
    CosThMu_abs_cc1pip_stack.push_back(h3_CosThMu_abs_cc1pip);
    CosThMu_abs_cc1pip_stack.push_back(h4_CosThMu_abs_cc1pip);
    CosThMu_abs_cc1pip_stack.push_back(h5_CosThMu_abs_cc1pip);

    HistStacker(CosThMu_abs_cc1pip_stack, leg_entries, norms, "CosThetaMu Absolute Normalisation, CC1#pi^{+}", "CosThMu_CC1PiP_Plot.png", "cos(#theta_{#mu})", "Number of SBND events");

    delete h1_CosThMu_abs_cc1pip;
    delete h2_CosThMu_abs_cc1pip;
    delete h3_CosThMu_abs_cc1pip;
    delete h4_CosThMu_abs_cc1pip;
    delete h5_CosThMu_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          CosThetaMu Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_abs_cc1pi0_stack;

    // The Ev histograms
    TH1D *h1_CosThMu_abs_cc1pi0 = new TH1D("h1_CosThMu_abs_cc1pi0","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThMu_abs_cc1pi0 = new TH1D("h2_CosThMu_abs_cc1pi0","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThMu_abs_cc1pi0 = new TH1D("h3_CosThMu_abs_cc1pi0","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThMu_abs_cc1pi0 = new TH1D("h4_CosThMu_abs_cc1pi0","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThMu_abs_cc1pi0 = new TH1D("h5_CosThMu_abs_cc1pi0","cos(#theta_{#mu}) Absolute Normalisation, CC1#pi^{0}",40,-1,1);

    gst1->Draw("cthl>>h1_CosThMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst2->Draw("cthl>>h2_CosThMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst3->Draw("cthl>>h3_CosThMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst4->Draw("cthl>>h4_CosThMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst5->Draw("cthl>>h5_CosThMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    
    // Push them onto the vector
    CosThMu_abs_cc1pi0_stack.push_back(h1_CosThMu_abs_cc1pi0);
    CosThMu_abs_cc1pi0_stack.push_back(h2_CosThMu_abs_cc1pi0);
    CosThMu_abs_cc1pi0_stack.push_back(h3_CosThMu_abs_cc1pi0);
    CosThMu_abs_cc1pi0_stack.push_back(h4_CosThMu_abs_cc1pi0);
    CosThMu_abs_cc1pi0_stack.push_back(h5_CosThMu_abs_cc1pi0);

    HistStacker(CosThMu_abs_cc1pi0_stack, leg_entries, norms, "CosThetaMu Absolute Normalisation, CC1#pi^{0}", "CosThMu_CC1Pi0_Plot.png", "cos(#theta_{#mu})", "Number of SBND events");

    delete h1_CosThMu_abs_cc1pi0;
    delete h2_CosThMu_abs_cc1pi0;
    delete h3_CosThMu_abs_cc1pi0;
    delete h4_CosThMu_abs_cc1pi0;
    delete h5_CosThMu_abs_cc1pi0;
    // -------------------------------------------------------------------------
    //                          CosTheta Area norm cc0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_area_cc0pi_stack;
    vector< double > int_cth_mu_cc0pi;

    // The Ev histograms
    TH1D *h1_CosThMu_area_cc0pi = new TH1D("h1_CosThMu_area_cc0pi","cos(#theta_{#mu}) Area Normalised, CC0#pi",40,-1,1);
    TH1D *h2_CosThMu_area_cc0pi = new TH1D("h2_CosThMu_area_cc0pi","cos(#theta_{#mu}) Area Normalised, CC0#pi",40,-1,1);
    TH1D *h3_CosThMu_area_cc0pi = new TH1D("h3_CosThMu_area_cc0pi","cos(#theta_{#mu}) Area Normalised, CC0#pi",40,-1,1);
    TH1D *h4_CosThMu_area_cc0pi = new TH1D("h4_CosThMu_area_cc0pi","cos(#theta_{#mu}) Area Normalised, CC0#pi",40,-1,1);
    TH1D *h5_CosThMu_area_cc0pi = new TH1D("h5_CosThMu_area_cc0pi","cos(#theta_{#mu}) Area Normalised, CC0#pi",40,-1,1);

    gst1->Draw("cthl>>h1_CosThMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst2->Draw("cthl>>h2_CosThMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst3->Draw("cthl>>h3_CosThMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst4->Draw("cthl>>h4_CosThMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst5->Draw("cthl>>h5_CosThMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    
    // Push back the scaling
    double int1_cth_mu_cc0pi, int2_cth_mu_cc0pi, int3_cth_mu_cc0pi, int4_cth_mu_cc0pi, int5_cth_mu_cc0pi;

    int1_cth_mu_cc0pi = h1_CosThMu_area_cc0pi->Integral();
    int2_cth_mu_cc0pi = h2_CosThMu_area_cc0pi->Integral();
    int3_cth_mu_cc0pi = h3_CosThMu_area_cc0pi->Integral();
    int4_cth_mu_cc0pi = h4_CosThMu_area_cc0pi->Integral();
    int5_cth_mu_cc0pi = h5_CosThMu_area_cc0pi->Integral();
    
    int_cth_mu_cc0pi.push_back( 1 / int1_cth_mu_cc0pi );
    int_cth_mu_cc0pi.push_back( 1 / int2_cth_mu_cc0pi );
    int_cth_mu_cc0pi.push_back( 1 / int3_cth_mu_cc0pi );
    int_cth_mu_cc0pi.push_back( 1 / int4_cth_mu_cc0pi );
    int_cth_mu_cc0pi.push_back( 1 / int5_cth_mu_cc0pi );
    
    // Push them onto the vector
    CosThMu_area_cc0pi_stack.push_back(h1_CosThMu_area_cc0pi);
    CosThMu_area_cc0pi_stack.push_back(h2_CosThMu_area_cc0pi);
    CosThMu_area_cc0pi_stack.push_back(h3_CosThMu_area_cc0pi);
    CosThMu_area_cc0pi_stack.push_back(h4_CosThMu_area_cc0pi);
    CosThMu_area_cc0pi_stack.push_back(h5_CosThMu_area_cc0pi);

    HistStacker(CosThMu_area_cc0pi_stack, leg_entries, int_cth_mu_cc0pi, "CosThetaMu Area Normalised, CC0#pi", "CosThMu_cc0pi_Area_Plot.png", "cos(#theta_{#mu})", "SBND events, area normalised");

    delete h1_CosThMu_area_cc0pi;
    delete h2_CosThMu_area_cc0pi;
    delete h3_CosThMu_area_cc0pi;
    delete h4_CosThMu_area_cc0pi;
    delete h5_CosThMu_area_cc0pi;
    // -------------------------------------------------------------------------
    //                      CosThetaMu Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_area_cc1pip_stack;
    vector< double > int_cth_mu_cc1pip;

    // The Ev histograms
    TH1D *h1_CosThMu_area_cc1pip = new TH1D("h1_CosThMu_area_cc1pip","cos(#theta_{#mu}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h2_CosThMu_area_cc1pip = new TH1D("h2_CosThMu_area_cc1pip","cos(#theta_{#mu}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h3_CosThMu_area_cc1pip = new TH1D("h3_CosThMu_area_cc1pip","cos(#theta_{#mu}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h4_CosThMu_area_cc1pip = new TH1D("h4_CosThMu_area_cc1pip","cos(#theta_{#mu}) Area Normalised, CC1#pi^{+}",40,-1,1);
    TH1D *h5_CosThMu_area_cc1pip = new TH1D("h5_CosThMu_area_cc1pip","cos(#theta_{#mu}) Area Normalised, CC1#pi^{+}",40,-1,1);

    gst1->Draw("cthl>>h1_CosThMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )  ");
    gst2->Draw("cthl>>h2_CosThMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst3->Draw("cthl>>h3_CosThMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst4->Draw("cthl>>h4_CosThMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst5->Draw("cthl>>h5_CosThMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
     
    // Push back the scaling
    double int1_cth_mu_cc1pip, int2_cth_mu_cc1pip, int3_cth_mu_cc1pip, int4_cth_mu_cc1pip, int5_cth_mu_cc1pip;

    int1_cth_mu_cc1pip = h1_CosThMu_area_cc1pip->Integral();
    int2_cth_mu_cc1pip = h2_CosThMu_area_cc1pip->Integral();
    int3_cth_mu_cc1pip = h3_CosThMu_area_cc1pip->Integral();
    int4_cth_mu_cc1pip = h4_CosThMu_area_cc1pip->Integral();
    int5_cth_mu_cc1pip = h5_CosThMu_area_cc1pip->Integral();
    
    int_cth_mu_cc1pip.push_back( 1 / int1_cth_mu_cc1pip );
    int_cth_mu_cc1pip.push_back( 1 / int2_cth_mu_cc1pip );
    int_cth_mu_cc1pip.push_back( 1 / int3_cth_mu_cc1pip );
    int_cth_mu_cc1pip.push_back( 1 / int4_cth_mu_cc1pip );
    int_cth_mu_cc1pip.push_back( 1 / int5_cth_mu_cc1pip );
    
    // Push them onto the vector
    CosThMu_area_cc1pip_stack.push_back(h1_CosThMu_area_cc1pip);
    CosThMu_area_cc1pip_stack.push_back(h2_CosThMu_area_cc1pip);
    CosThMu_area_cc1pip_stack.push_back(h3_CosThMu_area_cc1pip);
    CosThMu_area_cc1pip_stack.push_back(h4_CosThMu_area_cc1pip);
    CosThMu_area_cc1pip_stack.push_back(h5_CosThMu_area_cc1pip);

    HistStacker(CosThMu_area_cc1pip_stack, leg_entries, int_cth_mu_cc1pip, "CosThetaMu Area Normalised, CC1#pi^{+}", "CosThMu_CC1PiP_Area_Plot.png", "cos(#theta_{#mu})", "SBND events, area normalised");

    delete h1_CosThMu_area_cc1pip;
    delete h2_CosThMu_area_cc1pip;
    delete h3_CosThMu_area_cc1pip;
    delete h4_CosThMu_area_cc1pip;
    delete h5_CosThMu_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          CosThetaMu Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > CosThMu_area_cc1pi0_stack;
    vector< double > int_cth_mu_cc1pi0;

    // The Ev histograms
    TH1D *h1_CosThMu_area_cc1pi0 = new TH1D("h1_CosThMu_area_cc1pi0","cos(#theta_{#mu}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h2_CosThMu_area_cc1pi0 = new TH1D("h2_CosThMu_area_cc1pi0","cos(#theta_{#mu}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h3_CosThMu_area_cc1pi0 = new TH1D("h3_CosThMu_area_cc1pi0","cos(#theta_{#mu}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h4_CosThMu_area_cc1pi0 = new TH1D("h4_CosThMu_area_cc1pi0","cos(#theta_{#mu}) Area Normalised, CC1#pi^{0}",40,-1,1);
    TH1D *h5_CosThMu_area_cc1pi0 = new TH1D("h5_CosThMu_area_cc1pi0","cos(#theta_{#mu}) Area Normalised, CC1#pi^{0}",40,-1,1);
    
    gst1->Draw("cthl>>h1_CosThMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst2->Draw("cthl>>h2_CosThMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst3->Draw("cthl>>h3_CosThMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst4->Draw("cthl>>h4_CosThMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst5->Draw("cthl>>h5_CosThMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
     
    // Push back the scaling
    double int1_cth_mu_cc1pi0, int2_cth_mu_cc1pi0, int3_cth_mu_cc1pi0, int4_cth_mu_cc1pi0, int5_cth_mu_cc1pi0;

    int1_cth_mu_cc1pi0 = h1_CosThMu_area_cc1pi0->Integral();
    int2_cth_mu_cc1pi0 = h2_CosThMu_area_cc1pi0->Integral();
    int3_cth_mu_cc1pi0 = h3_CosThMu_area_cc1pi0->Integral();
    int4_cth_mu_cc1pi0 = h4_CosThMu_area_cc1pi0->Integral();
    int5_cth_mu_cc1pi0 = h5_CosThMu_area_cc1pi0->Integral();
    
    int_cth_mu_cc1pi0.push_back( 1 / int1_cth_mu_cc1pi0 );
    int_cth_mu_cc1pi0.push_back( 1 / int2_cth_mu_cc1pi0 );
    int_cth_mu_cc1pi0.push_back( 1 / int3_cth_mu_cc1pi0 );
    int_cth_mu_cc1pi0.push_back( 1 / int4_cth_mu_cc1pi0 );
    int_cth_mu_cc1pi0.push_back( 1 / int5_cth_mu_cc1pi0 );
    
    // Push them onto the vector
    CosThMu_area_cc1pi0_stack.push_back(h1_CosThMu_area_cc1pi0);
    CosThMu_area_cc1pi0_stack.push_back(h2_CosThMu_area_cc1pi0);
    CosThMu_area_cc1pi0_stack.push_back(h3_CosThMu_area_cc1pi0);
    CosThMu_area_cc1pi0_stack.push_back(h4_CosThMu_area_cc1pi0);
    CosThMu_area_cc1pi0_stack.push_back(h5_CosThMu_area_cc1pi0);

    HistStacker(CosThMu_area_cc1pi0_stack, leg_entries, int_cth_mu_cc1pi0, "CosThetaMu Area Normalised, CC1#pi^{0}", "CosThMu_CC1Pi0_Area_Plot.png", "cos(#theta_{#mu})", "SBND events, area normalised");

    delete h1_CosThMu_area_cc1pi0;
    delete h2_CosThMu_area_cc1pi0;
    delete h3_CosThMu_area_cc1pi0;
    delete h4_CosThMu_area_cc1pi0;
    delete h5_CosThMu_area_cc1pi0;
    // -------------------------------------------------------------------------
    //                          PMu Absolute norm cc0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_abs_cc0pi_stack;

    // The Ev histograms
    TH1D *h1_PMu_abs_cc0pi = new TH1D("h1_PMu_abs_cc0pi","P_{#mu} Absolute Normalisation, CC0#pi",40,0,1.5);
    TH1D *h2_PMu_abs_cc0pi = new TH1D("h2_PMu_abs_cc0pi","P_{#mu} Absolute Normalisation, CC0#pi",40,0,1.5);
    TH1D *h3_PMu_abs_cc0pi = new TH1D("h3_PMu_abs_cc0pi","P_{#mu} Absolute Normalisation, CC0#pi",40,0,1.5);
    TH1D *h4_PMu_abs_cc0pi = new TH1D("h4_PMu_abs_cc0pi","P_{#mu} Absolute Normalisation, CC0#pi",40,0,1.5);
    TH1D *h5_PMu_abs_cc0pi = new TH1D("h5_PMu_abs_cc0pi","P_{#mu} Absolute Normalisation, CC0#pi",40,0,1.5);

    gst1->Draw("pl>>h1_PMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst2->Draw("pl>>h2_PMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst3->Draw("pl>>h3_PMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst4->Draw("pl>>h4_PMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    gst5->Draw("pl>>h5_PMu_abs_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 )  ");
    
    // Push them onto the vector
    PMu_abs_cc0pi_stack.push_back(h1_PMu_abs_cc0pi);
    PMu_abs_cc0pi_stack.push_back(h2_PMu_abs_cc0pi);
    PMu_abs_cc0pi_stack.push_back(h3_PMu_abs_cc0pi);
    PMu_abs_cc0pi_stack.push_back(h4_PMu_abs_cc0pi);
    PMu_abs_cc0pi_stack.push_back(h5_PMu_abs_cc0pi);

    HistStacker(PMu_abs_cc0pi_stack, leg_entries, norms, "PMu Absolute Normalisation, CC0#pi", "PMu_cc0pi_Plot.png", "P_{#mu}", "Number of SBND events");

    delete h1_PMu_abs_cc0pi;
    delete h2_PMu_abs_cc0pi;
    delete h3_PMu_abs_cc0pi;
    delete h4_PMu_abs_cc0pi;
    delete h5_PMu_abs_cc0pi;
    // -------------------------------------------------------------------------
    //                          PMu Absolute norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_abs_cc1pip_stack;

    // The Ev histograms
    TH1D *h1_PMu_abs_cc1pip = new TH1D("h1_PMu_abs_cc1pip","P_{#mu} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h2_PMu_abs_cc1pip = new TH1D("h2_PMu_abs_cc1pip","P_{#mu} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h3_PMu_abs_cc1pip = new TH1D("h3_PMu_abs_cc1pip","P_{#mu} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h4_PMu_abs_cc1pip = new TH1D("h4_PMu_abs_cc1pip","P_{#mu} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);
    TH1D *h5_PMu_abs_cc1pip = new TH1D("h5_PMu_abs_cc1pip","P_{#mu} Absolute Normalisation, CC1#pi^{+}",40,0,1.5);

    gst1->Draw("pl>>h1_PMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst2->Draw("pl>>h2_PMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst3->Draw("pl>>h3_PMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst4->Draw("pl>>h4_PMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst5->Draw("pl>>h5_PMu_abs_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    
    // Push them onto the vector
    PMu_abs_cc1pip_stack.push_back(h1_PMu_abs_cc1pip);
    PMu_abs_cc1pip_stack.push_back(h2_PMu_abs_cc1pip);
    PMu_abs_cc1pip_stack.push_back(h3_PMu_abs_cc1pip);
    PMu_abs_cc1pip_stack.push_back(h4_PMu_abs_cc1pip);
    PMu_abs_cc1pip_stack.push_back(h5_PMu_abs_cc1pip);

    HistStacker(PMu_abs_cc1pip_stack, leg_entries, norms, "PMu Absolute Normalisation, CC1#pi^{+}", "PMu_CC1PiP_Plot.png", "P_{#mu}", "Number of SBND events");

    delete h1_PMu_abs_cc1pip;
    delete h2_PMu_abs_cc1pip;
    delete h3_PMu_abs_cc1pip;
    delete h4_PMu_abs_cc1pip;
    delete h5_PMu_abs_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          PMu Absolute norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_abs_cc1pi0_stack;

    // The Ev histograms
    TH1D *h1_PMu_abs_cc1pi0 = new TH1D("h1_PMu_abs_cc1pi0","P_{#mu} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h2_PMu_abs_cc1pi0 = new TH1D("h2_PMu_abs_cc1pi0","P_{#mu} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h3_PMu_abs_cc1pi0 = new TH1D("h3_PMu_abs_cc1pi0","P_{#mu} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h4_PMu_abs_cc1pi0 = new TH1D("h4_PMu_abs_cc1pi0","P_{#mu} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);
    TH1D *h5_PMu_abs_cc1pi0 = new TH1D("h5_PMu_abs_cc1pi0","P_{#mu} Absolute Normalisation, CC1#pi^{0}",40,0,1.5);

    gst1->Draw("pl>>h1_PMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst2->Draw("pl>>h2_PMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst3->Draw("pl>>h3_PMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst4->Draw("pl>>h4_PMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst5->Draw("pl>>h5_PMu_abs_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    
    // Push them onto the vector
    PMu_abs_cc1pi0_stack.push_back(h1_PMu_abs_cc1pi0);
    PMu_abs_cc1pi0_stack.push_back(h2_PMu_abs_cc1pi0);
    PMu_abs_cc1pi0_stack.push_back(h3_PMu_abs_cc1pi0);
    PMu_abs_cc1pi0_stack.push_back(h4_PMu_abs_cc1pi0);
    PMu_abs_cc1pi0_stack.push_back(h5_PMu_abs_cc1pi0);

    HistStacker(PMu_abs_cc1pi0_stack, leg_entries, norms, "PMu Absolute Normalisation, CC1#pi^{0}", "PMu_CC1Pi0_Plot.png", "P_{#mu}", "Number of SBND events");

    delete h1_PMu_abs_cc1pi0;
    delete h2_PMu_abs_cc1pi0;
    delete h3_PMu_abs_cc1pi0;
    delete h4_PMu_abs_cc1pi0;
    delete h5_PMu_abs_cc1pi0;
    
    // -------------------------------------------------------------------------
    //                          PMu Area norm cc0pi
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_area_cc0pi_stack;
    vector< double > int_p_mu_cc0pi;

    // The Ev histograms
    TH1D *h1_PMu_area_cc0pi = new TH1D("h1_PMu_area_cc0pi","P_{#mu} Area Normalised, CC0#pi",40,0,1.5);
    TH1D *h2_PMu_area_cc0pi = new TH1D("h2_PMu_area_cc0pi","P_{#mu} Area Normalised, CC0#pi",40,0,1.5);
    TH1D *h3_PMu_area_cc0pi = new TH1D("h3_PMu_area_cc0pi","P_{#mu} Area Normalised, CC0#pi",40,0,1.5);
    TH1D *h4_PMu_area_cc0pi = new TH1D("h4_PMu_area_cc0pi","P_{#mu} Area Normalised, CC0#pi",40,0,1.5);
    TH1D *h5_PMu_area_cc0pi = new TH1D("h5_PMu_area_cc0pi","P_{#mu} Area Normalised, CC0#pi",40,0,1.5);

    gst1->Draw("pl>>h1_PMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst2->Draw("pl>>h2_PMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst3->Draw("pl>>h3_PMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst4->Draw("pl>>h4_PMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    gst4->Draw("pl>>h5_PMu_area_cc0pi", "cc && ( nfpi0 + nfpim + nfpip == 0 ) ");
    
    // Push back the scaling
    double int1_p_mu_cc0pi, int2_p_mu_cc0pi, int3_p_mu_cc0pi, int4_p_mu_cc0pi, int5_p_mu_cc0pi;

    int1_p_mu_cc0pi = h1_PMu_area_cc0pi->Integral();
    int2_p_mu_cc0pi = h2_PMu_area_cc0pi->Integral();
    int3_p_mu_cc0pi = h3_PMu_area_cc0pi->Integral();
    int4_p_mu_cc0pi = h4_PMu_area_cc0pi->Integral();
    int5_p_mu_cc0pi = h5_PMu_area_cc0pi->Integral();
    
    int_p_mu_cc0pi.push_back( 1 / int1_p_mu_cc0pi );
    int_p_mu_cc0pi.push_back( 1 / int2_p_mu_cc0pi );
    int_p_mu_cc0pi.push_back( 1 / int3_p_mu_cc0pi );
    int_p_mu_cc0pi.push_back( 1 / int4_p_mu_cc0pi );
    int_p_mu_cc0pi.push_back( 1 / int5_p_mu_cc0pi );
    
    // Push them onto the vector
    PMu_area_cc0pi_stack.push_back(h1_PMu_area_cc0pi);
    PMu_area_cc0pi_stack.push_back(h2_PMu_area_cc0pi);
    PMu_area_cc0pi_stack.push_back(h3_PMu_area_cc0pi);
    PMu_area_cc0pi_stack.push_back(h4_PMu_area_cc0pi);
    PMu_area_cc0pi_stack.push_back(h5_PMu_area_cc0pi);

    HistStacker(PMu_area_cc0pi_stack, leg_entries, int_p_mu_cc0pi, "PMu Area Normalised, CC0#pi", "PMu_cc0pi_Area_Plot.png", "P_{#mu}", "SBND events, area normalised");

    delete h1_PMu_area_cc0pi;
    delete h2_PMu_area_cc0pi;
    delete h3_PMu_area_cc0pi;
    delete h4_PMu_area_cc0pi;
    delete h5_PMu_area_cc0pi;
    
    // -------------------------------------------------------------------------
    //                      PMu Area norm CC1pip
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_area_cc1pip_stack;
    vector< double > int_p_mu_cc1pip;

    // The Ev histograms
    TH1D *h1_PMu_area_cc1pip = new TH1D("h1_PMu_area_cc1pip","P_{#mu} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h2_PMu_area_cc1pip = new TH1D("h2_PMu_area_cc1pip","P_{#mu} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h3_PMu_area_cc1pip = new TH1D("h3_PMu_area_cc1pip","P_{#mu} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h4_PMu_area_cc1pip = new TH1D("h4_PMu_area_cc1pip","P_{#mu} Area Normalised, CC1#pi^{+}",40,0,1.5);
    TH1D *h5_PMu_area_cc1pip = new TH1D("h5_PMu_area_cc1pip","P_{#mu} Area Normalised, CC1#pi^{+}",40,0,1.5);

    gst1->Draw("pl>>h1_PMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )  ");
    gst2->Draw("pl>>h2_PMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst3->Draw("pl>>h3_PMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst4->Draw("pl>>h4_PMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
    gst5->Draw("pl>>h5_PMu_area_cc1pip", "cc && nfpip == 1 && ( nfpim + nfpi0 == 0 )   ");
     
    // Push back the scaling
    double int1_p_mu_cc1pip, int2_p_mu_cc1pip, int3_p_mu_cc1pip, int4_p_mu_cc1pip, int5_p_mu_cc1pip;

    int1_p_mu_cc1pip = h1_PMu_area_cc1pip->Integral();
    int2_p_mu_cc1pip = h2_PMu_area_cc1pip->Integral();
    int3_p_mu_cc1pip = h3_PMu_area_cc1pip->Integral();
    int4_p_mu_cc1pip = h4_PMu_area_cc1pip->Integral();
    int5_p_mu_cc1pip = h5_PMu_area_cc1pip->Integral();
    
    int_p_mu_cc1pip.push_back( 1 / int1_p_mu_cc1pip );
    int_p_mu_cc1pip.push_back( 1 / int2_p_mu_cc1pip );
    int_p_mu_cc1pip.push_back( 1 / int3_p_mu_cc1pip );
    int_p_mu_cc1pip.push_back( 1 / int4_p_mu_cc1pip );
    int_p_mu_cc1pip.push_back( 1 / int5_p_mu_cc1pip );
    
    // Push them onto the vector
    PMu_area_cc1pip_stack.push_back(h1_PMu_area_cc1pip);
    PMu_area_cc1pip_stack.push_back(h2_PMu_area_cc1pip);
    PMu_area_cc1pip_stack.push_back(h3_PMu_area_cc1pip);
    PMu_area_cc1pip_stack.push_back(h4_PMu_area_cc1pip);
    PMu_area_cc1pip_stack.push_back(h5_PMu_area_cc1pip);

    HistStacker(PMu_area_cc1pip_stack, leg_entries, int_p_mu_cc1pip, "PMu Area Normalised, CC1#pi^{+}", "PMu_CC1PiP_Area_Plot.png", "P_{#mu}", "SBND events, area normalised");

    delete h1_PMu_area_cc1pip;
    delete h2_PMu_area_cc1pip;
    delete h3_PMu_area_cc1pip;
    delete h4_PMu_area_cc1pip;
    delete h5_PMu_area_cc1pip;
    
    // -------------------------------------------------------------------------
    //                          PMu Area norm CC1pi0
    // -------------------------------------------------------------------------
    
    // Vector of histograms and strings to pass to the stacking function
    vector< TH1D* > PMu_area_cc1pi0_stack;
    vector< double > int_p_mu_cc1pi0;

    // The Ev histograms
    TH1D *h1_PMu_area_cc1pi0 = new TH1D("h1_PMu_area_cc1pi0","P_{#mu} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h2_PMu_area_cc1pi0 = new TH1D("h2_PMu_area_cc1pi0","P_{#mu} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h3_PMu_area_cc1pi0 = new TH1D("h3_PMu_area_cc1pi0","P_{#mu} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h4_PMu_area_cc1pi0 = new TH1D("h4_PMu_area_cc1pi0","P_{#mu} Area Normalised, CC1#pi^{0}",40,0,1.5);
    TH1D *h5_PMu_area_cc1pi0 = new TH1D("h5_PMu_area_cc1pi0","P_{#mu} Area Normalised, CC1#pi^{0}",40,0,1.5);
    
    gst1->Draw("pl>>h1_PMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst2->Draw("pl>>h2_PMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst3->Draw("pl>>h3_PMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst4->Draw("pl>>h4_PMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
    gst5->Draw("pl>>h5_PMu_area_cc1pi0", "cc && nfpi0 == 1 && ( nfpim + nfpip == 0 )  ");
     
    // Push back the scaling
    double int1_p_mu_cc1pi0, int2_p_mu_cc1pi0, int3_p_mu_cc1pi0, int4_p_mu_cc1pi0, int5_p_mu_cc1pi0;

    int1_p_mu_cc1pi0 = h1_PMu_area_cc1pi0->Integral();
    int2_p_mu_cc1pi0 = h2_PMu_area_cc1pi0->Integral();
    int3_p_mu_cc1pi0 = h3_PMu_area_cc1pi0->Integral();
    int4_p_mu_cc1pi0 = h4_PMu_area_cc1pi0->Integral();
    int5_p_mu_cc1pi0 = h5_PMu_area_cc1pi0->Integral();
    
    int_p_mu_cc1pi0.push_back( 1 / int1_p_mu_cc1pi0 );
    int_p_mu_cc1pi0.push_back( 1 / int2_p_mu_cc1pi0 );
    int_p_mu_cc1pi0.push_back( 1 / int3_p_mu_cc1pi0 );
    int_p_mu_cc1pi0.push_back( 1 / int4_p_mu_cc1pi0 );
    int_p_mu_cc1pi0.push_back( 1 / int5_p_mu_cc1pi0 );
    
    // Push them onto the vector
    PMu_area_cc1pi0_stack.push_back(h1_PMu_area_cc1pi0);
    PMu_area_cc1pi0_stack.push_back(h2_PMu_area_cc1pi0);
    PMu_area_cc1pi0_stack.push_back(h3_PMu_area_cc1pi0);
    PMu_area_cc1pi0_stack.push_back(h4_PMu_area_cc1pi0);
    PMu_area_cc1pi0_stack.push_back(h5_PMu_area_cc1pi0);

    HistStacker(PMu_area_cc1pi0_stack, leg_entries, int_p_mu_cc1pi0, "PMu Area Normalised, CC1#pi^{0}", "PMu_CC1Pi0_Area_Plot.png", "P_{#mu}", "SBND events, area normalised");

    delete h1_PMu_area_cc1pi0;
    delete h2_PMu_area_cc1pi0;
    delete h3_PMu_area_cc1pi0;
    delete h4_PMu_area_cc1pi0;
    delete h5_PMu_area_cc1pi0;
    
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

        if( i != 4){
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
    
    // ----------------------------------------------------------------------
    //                                Draw
    // ----------------------------------------------------------------------


    gStyle->SetEndErrorSize(3);
    gStyle->SetErrorX(0);
*/
    // Draw the histograms
    //hists[0]->Draw( "e1" );
    hists[0]->Draw();
    for ( int i = 1; i < n_hists; ++i ) {        
     
        // For now, don't draw G17_01b
        if( i != 4 ){
            // Draw the histograms 
            //hists[i]->Draw( "samee1" );
            hists[i]->Draw( "same" );
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
             + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0 ){
         
                // Get the values needed
                e   = b_mu_e->GetLeaf("El")->GetValue();
                p   = b_mu_p->GetLeaf("pl")->GetValue();
                cth = b_theta->GetLeaf("cthl")->GetValue(); 
            
                reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( e - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 
        
                reco_mc = TMath::Abs( reco - double(b_nu_e->GetLeaf("Ev")->GetValue()) );

                // Make the vectors of reconstructed and reconstructed-MC energy for CC0pi
                MC_reco_E_CC.push_back(reco_mc);
                reco_E_CC.push_back(reco);

        }
        // For NC0pi
        else if( b_nc->GetLeaf("nc")->GetValue() != 0 
                 && b_nfpip->GetLeaf("nfpip")->GetValue()
                  + b_nfpim->GetLeaf("nfpim")->GetValue()
                  + b_nfpi0->GetLeaf("nfpi0")->GetValue() == 0 ){
              
                // Get the values
                e   = b_mu_e->GetLeaf("El")->GetValue();
                p   = b_mu_p->GetLeaf("pl")->GetValue();
                cth = b_theta->GetLeaf("cthl")->GetValue(); 
            
                reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( e - ( 1 / ( 2 * m_n) ) * m_mu * m_mu  ); 
        
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
                 vector< double > &n_cc_fsi,
                 vector< double > &n_nc_fsi,
                 vector< double > &n_cc_mc_fsi,
                 vector< double > &n_nc_mc_fsi ){

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
    n_cc_mc_fsi.push_back (ncc0pi);
    n_cc_mc_fsi.push_back (ncc1pip);
    n_cc_mc_fsi.push_back (ncc1pim);
    n_cc_mc_fsi.push_back (ncc1pi0);
    n_cc_mc_fsi.push_back (ncc2pip);
    n_cc_mc_fsi.push_back (ncc2pim);
    n_cc_mc_fsi.push_back (ncc2pi0);
    n_cc_mc_fsi.push_back (nccpippim);
    n_cc_mc_fsi.push_back (nccpippi0);
    n_cc_mc_fsi.push_back (nccpimpi0);
    n_cc_mc_fsi.push_back (ncc3pi);
    n_cc_mc_fsi.push_back (ncccoh);

    n_nc_mc_fsi.push_back (nnc0pi);
    n_nc_mc_fsi.push_back (nnc1pip);
    n_nc_mc_fsi.push_back (nnc1pim);
    n_nc_mc_fsi.push_back (nnc1pi0);
    n_nc_mc_fsi.push_back (nnc2pip);
    n_nc_mc_fsi.push_back (nnc2pim);
    n_nc_mc_fsi.push_back (nnc2pi0);
    n_nc_mc_fsi.push_back (nncpippim);
    n_nc_mc_fsi.push_back (nncpippi0);
    n_nc_mc_fsi.push_back (nncpimpi0);
    n_nc_mc_fsi.push_back (nnc3pi);
    n_nc_mc_fsi.push_back (nnccoh);

    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pip) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pim) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc1pi0) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc2pip) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc2pim) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc2pi0) );
    n_cc_fsi.push_back ( TMath::Floor( norm * nccpippim) );
    n_cc_fsi.push_back ( TMath::Floor( norm * nccpippi0) );
    n_cc_fsi.push_back ( TMath::Floor( norm * nccpimpi0) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc3pi) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncccoh) );

    n_nc_fsi.push_back ( TMath::Floor( norm * nnc0pi) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc1pip) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc1pim) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc1pi0) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc2pip) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc2pim) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc2pi0) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nncpippim) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nncpippi0) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nncpimpi0) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnc3pi) );
    n_nc_fsi.push_back ( TMath::Floor( norm * nnccoh) );

    file << " -----------SBND----------" << endl; 
    file << " ------------------------- " << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi * norm << endl;
    file << setprecision(5) << " CC1Pi+   : " << ncc1pip * norm  << endl;
    file << setprecision(5) << " CC1Pi-   : " << ncc1pim * norm  << endl;
    file << setprecision(5) << " CC1Pi0   : " << ncc1pi0 * norm  << endl;
    file << setprecision(5) << " CC2Pi+   : " << ncc2pip * norm  << endl;
    file << setprecision(5) << " CC2Pi-   : " << ncc2pim * norm  << endl;
    file << setprecision(5) << " CC2Pi0   : " << ncc2pi0 * norm  << endl;
    file << setprecision(5) << " CCPi+Pi- : " << nccpippim * norm  << endl;
    file << setprecision(5) << " CCPi+Pi0 : " << nccpippi0 * norm  << endl;
    file << setprecision(5) << " CCPi-Pi0 : " << nccpimpi0 * norm  << endl;
    file << setprecision(5) << " CC>3Pi   : " << ncc3pi * norm  << endl;
    file << setprecision(5) << " CCCOH    : " << ncccoh * norm  << endl;
    file << " ------------------------- " << endl;
    file << setprecision(5) << " NC0Pi    : " << nnc0pi * norm  << endl;
    file << setprecision(5) << " NC1Pi+   : " << nnc1pip * norm  << endl;
    file << setprecision(5) << " NC1Pi-   : " << nnc1pim * norm  << endl;
    file << setprecision(5) << " NC1Pi0   : " << nnc1pi0 * norm  << endl;
    file << setprecision(5) << " NC2Pi+   : " << nnc2pip * norm  << endl;
    file << setprecision(5) << " NC2Pi-   : " << nnc2pim * norm  << endl;
    file << setprecision(5) << " NC2Pi0   : " << nnc2pi0 * norm  << endl;
    file << setprecision(5) << " NCPi+Pi- : " << nncpippim * norm  << endl;
    file << setprecision(5) << " NCPi+Pi0 : " << nncpippi0 * norm  << endl;
    file << setprecision(5) << " NCPi-Pi0 : " << nncpimpi0 * norm  << endl;
    file << setprecision(5) << " NC>3Pi   : " << nnc3pi * norm  << endl;
    file << setprecision(5) << " NCCOH    : " << nnccoh * norm  << endl;
    file << " ------------------------- " << endl;
    file << " ------------MC-----------" << endl; 
    file << " ------------------------- " << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi << endl;
    file << setprecision(5) << " CC1Pi+   : " << ncc1pip << endl;
    file << setprecision(5) << " CC1Pi-   : " << ncc1pim << endl;
    file << setprecision(5) << " CC1Pi0   : " << ncc1pi0 << endl;
    file << setprecision(5) << " CC2Pi+   : " << ncc2pip << endl;
    file << setprecision(5) << " CC2Pi-   : " << ncc2pim << endl;
    file << setprecision(5) << " CC2Pi0   : " << ncc2pi0 << endl;
    file << setprecision(5) << " CCPi+Pi- : " << nccpippim << endl;
    file << setprecision(5) << " CCPi+Pi0 : " << nccpippi0 << endl;
    file << setprecision(5) << " CCPi-Pi0 : " << nccpimpi0 << endl;
    file << setprecision(5) << " CC>3Pi   : " << ncc3pi << endl;
    file << setprecision(5) << " CCCOH    : " << ncccoh << endl;
    file << " ------------------------- " << endl;
    file << setprecision(5) << " NC0Pi    : " << nnc0pi << endl;
    file << setprecision(5) << " NC1Pi+   : " << nnc1pip << endl;
    file << setprecision(5) << " NC1Pi-   : " << nnc1pim << endl;
    file << setprecision(5) << " NC1Pi0   : " << nnc1pi0 << endl;
    file << setprecision(5) << " NC2Pi+   : " << nnc2pip << endl;
    file << setprecision(5) << " NC2Pi-   : " << nnc2pim << endl;
    file << setprecision(5) << " NC2Pi0   : " << nnc2pi0 << endl;
    file << setprecision(5) << " NCPi+Pi- : " << nncpippim << endl;
    file << setprecision(5) << " NCPi+Pi0 : " << nncpippi0 << endl;
    file << setprecision(5) << " NCPi-Pi0 : " << nncpimpi0 << endl;
    file << setprecision(5) << " NC>3Pi   : " << nnc3pi << endl;
    file << setprecision(5) << " NCCOH    : " << nnccoh << endl;
    file << " ------------------------- " << endl;
    file << " ------------------------- " << endl;
}

// -------------------------------------------------------------------------
//                      Make final state tables
// -------------------------------------------------------------------------
void MakeTable( const m_outer &n_cc_vect,
                const m_outer &n_nc_vect,
                const vector< string > interactions,
                ostream &file ){
    
    // Iterators 
    typedef m_outer::const_iterator map_it;

    // Make a table from an input map - of - maps to print nice things
    // Number of columns and rows to be made
    int n_models, n_interactions;
    n_models = n_cc_vect.size();
    n_interactions = interactions.size();
    
    // Get the model names
    vector< string > m_names;

    for (map_it it = n_cc_vect.begin(); it != n_cc_vect.end(); it++ ){
        m_names.push_back(it->first);
    }

    // Begin the tabular environment in the LaTeX output file for n_predictions columns
    file << "\\begin{longtable}{| l || * {" << n_models << "}{c | } }" << endl;
    file << "\\hline" << endl;
 
    // Fill the first line with "Configurations" for all the right-hand columns
    file << "  & \\multicolumn{ " << n_models << " }{c|}{ \\textbf{ Model Configurations } } \\\\" << endl;
    file << "\\hline" << endl;
 
 
    // Fill a line of the table with the prediction names
    file << " \\textbf{ Hadronic Final State } & " ;
    for( int i = 0; i < n_models - 1; ++i ){
        file << "\\rotatebox{90}{ \\textbf{ " << m_names[i] << " } } & ";
    }
    file << "\\rotatebox{90}{ \\textbf{ " << m_names[n_models - 1] << " } } \\\\" << endl;
    file << " \\hline " << endl;
 
    file << " \\multicolumn{ " << n_models + 1 << " }{ | c | }{ \\textit{ Charged Current } } \\\\ " << endl;
    file << " \\hline " << endl;

    // -------------------------------------------------------------------
    //      Loop over the outer maps and fill the rest of the table                                   
    // -------------------------------------------------------------------
   
    for( int i = 0; i < n_interactions; ++i ){

        file << interactions[i] << " & ";

        for( map_it it = n_cc_vect.begin(), it1 = --n_cc_vect.end(); it != it1; ++it ){

            file << setprecision(5) <<  it->second[i] << " & ";
        }
        map_it it_cc = --n_cc_vect.end();
        file << setprecision(5) <<  it_cc->second[i] << " \\\\ " << endl;
    } 
   
    file << " \\hline " << endl;
    
    file << " \\multicolumn{ " << n_models + 1 << " }{ | c | }{ \\textit{ Neutral Current } } \\\\ " << endl;
    file << " \\hline " << endl;

    for( int i = 0; i < n_interactions; ++i ){
    
        file << interactions[i] << " & ";

        for( map_it it2 = n_nc_vect.begin(), it3 = --n_nc_vect.end(); it2 != it3; ++it2 ){
        
            file << setprecision(5) <<  it2->second[i] << " & ";
        }
         
        map_it it_nc = --n_nc_vect.end();
        file << setprecision(5) <<  it_nc->second[i] << " \\\\ " << endl;
    } 
    file << " \\hline " << endl;
    
    // -------------------------------------------------------------------
    //                    End the table                                   
    // -------------------------------------------------------------------
 
    // End the tabular environment
    file << "\\end{longtable}" << endl;
    

}

