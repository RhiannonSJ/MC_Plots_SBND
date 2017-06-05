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

#include "function_defs.h"

using namespace std; 

int table() {
    
    TFile f_flux("/hepstore/rjones/Exercises/Fluxes/sbn_FHC_flux_hist.root");
    if (f_flux.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " SBND flux file is open " << endl;
    }
    
    // -------------------------------------------------------------------------
    //                         Open the cross section files
    // -------------------------------------------------------------------------
    
    
    // G16_01b
    TFile f1_xsec("/hepstore/rjones/Exercises/Flavours/G16_01b/xsec_files/total_xsec.root");
    if (f1_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01b xsec file is open " << endl;
    }
    
    // G16_02a
    TFile f2_xsec("/hepstore/rjones/Exercises/Flavours/G16_02a/xsec_files/total_xsec.root");
    if (f2_xsec.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_02a xsec file is open " << endl;
        cout << "===========================" << endl;
    }
    
    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open G16_01b
    TFile f1("/hepstore/rjones/Exercises/Flavours/G16_01b/sbnd/1M/gntp.10000.gst.root");
    if (f1.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_01b event file 1 is open " << endl;
    }

    // Open G16_02a
    TFile f2("/hepstore/rjones/Exercises/Flavours/G16_02a/sbnd/1M/gntp.10000.gst.root");
    if (f2.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " G16_02a event file 1 is open " << endl;
    }
    // -------------------------------------------------------------------------
    //                     Calculate the reco energy difference
    // -------------------------------------------------------------------------
    
    // Get the trees we want from the root files
    TTree *gst1 = (TTree*) f1.Get("gst");
    TTree *gst2 = (TTree*) f2.Get("gst");
    // -------------------------------------------------------------------------
    //                          Get the normalisations
    // -------------------------------------------------------------------------
    
    vector< double > norms;
    norms.push_back( Norm(1000000, f1_xsec, f_flux) );
    norms.push_back( Norm(1000000, f2_xsec, f_flux) );
    
    // -------------------------------------------------------------------------
    //           Charged and neutral current maps to put into table
    // -------------------------------------------------------------------------
    m_outer cc_proc_model_ints;
    m_outer nc_proc_model_ints;
    m_outer cc_reco_model_ints;
    m_outer nc_reco_model_ints;
    
    vector< double > n_cc_fsi_1;
    vector< double > n_cc_fsi_2;
    
    vector< double > n_nc_fsi_1;
    vector< double > n_nc_fsi_2;
    
    vector< double > n_cc_proc_1;
    vector< double > n_cc_proc_2;
    
    vector< double > n_nc_proc_1;
    vector< double > n_nc_proc_2;
    
    // txt file to compare with TeX output
    ofstream file_n;
    file_n.open("n_interactions.txt");

    file_n << " G16_01b " << endl;

    FSINumbers( gst1, file_n, norms[0], n_cc_proc_1, n_nc_proc_1, n_cc_fsi_1, n_nc_fsi_1 );
    
    cc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_cc_fsi_1 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_nc_fsi_1 ) );
    cc_proc_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_cc_proc_1 ) );
    nc_proc_model_ints.insert( pair< string, vector< double > >( "G17\\_01b", n_nc_proc_1 ) );
    
    file_n << " G16_02a " << endl;
    
    FSINumbers( gst2, file_n, norms[1],n_cc_proc_2, n_nc_proc_2, n_cc_fsi_2, n_nc_fsi_2 );
    
    cc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_02a", n_cc_fsi_2 ) );
    nc_reco_model_ints.insert( pair< string, vector< double > >( "G17\\_02a", n_nc_fsi_2 ) );
    cc_proc_model_ints.insert( pair< string, vector< double > >( "G17\\_02a", n_cc_proc_2 ) );
    nc_proc_model_ints.insert( pair< string, vector< double > >( "G17\\_02a", n_nc_proc_2 ) );
    
    // Make vector of final hadronic states
    // Using LaTeX syntax
    vector< string > FHS;
    vector< string > PP;

    FHS.push_back( "Inclusive" );
    FHS.push_back( "0 \\( \\pi \\)" );
    FHS.push_back( "\\hspace{.3cm} 0 \\( \\pi \\) 0p" );
    FHS.push_back( "\\hspace{.3cm} 0 \\( \\pi \\) 1p" );
    FHS.push_back( "\\hspace{.3cm} 0 \\( \\pi \\) 2p" );
    FHS.push_back( "\\hspace{.3cm} 0 \\( \\pi \\) 3p" );
    FHS.push_back( "\\hspace{.3cm} 0 \\( \\pi > \\) 3p" );
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
    
    PP.push_back( "Quasi-elastic" );
    PP.push_back( "Meson Exchange Current" );
    PP.push_back( "Deep Inelastic Scattering" );
    PP.push_back( "Coherent" );

    ofstream file_reco;
    file_reco.open( "FSI_Reco_Table.tex" );
    MakeTable( cc_reco_model_ints, nc_reco_model_ints, cc_proc_model_ints, nc_proc_model_ints, FHS, PP, file_reco );
    
    return 0;
}

