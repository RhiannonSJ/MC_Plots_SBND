/*
 * The functions used throughout the plotting of MC info
 *
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : June 2017
 *
*/

#include "function_defs.h"

using namespace std; 

// -------------------------------------------------------------------------
//                      normalisation function
// -------------------------------------------------------------------------
double Norm(int n_events, TFile &xsec_file, TFile &flux_file){
    
    // The constant values used throughout
    double Na         = 6.022e23;                 // Avogadro
    double A_Ar       = 0.039948;                 // Argon molar mass, kg/mol
    double rho_Ar     = 1390;                     // Liquid argon density, kg/m^3
    double M_fid_mu   = rho_Ar * 55;              // 55 m^3 fiducial volume for nu_mu * denisty ~ 77000
    double POT_sbnd   = 6.6e20;                   // POT
    double n_sbnd_mu  = ( M_fid_mu * Na ) / A_Ar; // Number of target particles in muon fiducial volume in argon
    double scale      = 1e-48;                    // Sort units from flux * xsec: (10^-34m^2/50MeV/10e6POT/m^2) -> (20*10^-6*10^-34)/GeV/POT -> 2e-39  
   

    // Get the flux histograms
    TH1D *h_flux        = (TH1D*) flux_file.Get("h_numu_110m");

    // Get the numu x_sec graphs
    TGraph *cc          = (TGraph*) xsec_file.Get("nu_mu_Ar40/tot_cc");
    TGraph *nc          = (TGraph*) xsec_file.Get("nu_mu_Ar40/tot_nc");
   
    TH1D * h_mu_prod    = new TH1D(*h_flux);

    int nx      = h_flux->GetNbinsX();

    // Get the product of flux and CC xsec + product of flux and NC xsec
    for ( int i = 0; i <= nx; ++i ){
        double mu_centre     = h_flux->GetBinCenter(i);
        double mu_content    = h_flux->GetBinContent(i);
        double mu_xsec_cc    = cc->Eval(mu_centre);
        double mu_xsec_nc    = nc->Eval(mu_centre);
        
        double mu_product    = (mu_xsec_cc * mu_content) + (mu_xsec_nc * mu_content);
        h_mu_prod->SetBinContent(i,mu_product);
    }
    
    // Integrate the products of flux and xsec
    double mu_prod    =  h_mu_prod->Integral();
   
    // Total number of numu, numubar, nue, nuebar events in SBND
    double n_mu    = mu_prod    * n_sbnd_mu * POT_sbnd * scale;
    
    // Normalisation of total number of events with the 1,000,000 sample
    return ( n_mu    / double(n_events) );
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

    // Draw the histograms
    hists[0]->Draw();
    for ( int i = 1; i < n_hists; ++i ) {        
     
        // For now, don't draw G17_01b
        if( i != 4 ){
            // Draw the histograms 
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
//             hist stacking function with statistical errors
// -------------------------------------------------------------------------
void ErrHistStacker ( vector< TH1D* >   &hists,
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

    // Draw the histograms
    hists[0]->Draw( "e1" );
    hists[0]->Draw();
    for ( int i = 1; i < n_hists; ++i ) {        
     
        // For now, don't draw G17_01b
        if( i != 4 ){
            // Draw the histograms 
            hists[i]->Draw( "samee1" );
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
                 vector< double > &n_cc_proc,
                 vector< double > &n_nc_proc,
                 vector< double > &n_cc_fsi,
                 vector< double > &n_nc_fsi ){

    // Firstly, get out the trees we want to look at
    // Need both cc and nc with varying number of outgoing pions
    //      - cc : charged current FSI
    //      - nc : neutral current FSI
    
    // Set the branch addresses for these leaves
    TBranch *b_pdg = event_tree->GetBranch("pdgf");
    TBranch *b_cc  = event_tree->GetBranch("cc");
    TBranch *b_nc  = event_tree->GetBranch("nc");
    TBranch *b_coh = event_tree->GetBranch("coh");
    TBranch *b_qel = event_tree->GetBranch("qel");
    TBranch *b_dis = event_tree->GetBranch("dis");
    TBranch *b_mec = event_tree->GetBranch("mec");
    TBranch *b_pi0 = event_tree->GetBranch("nfpi0");
    TBranch *b_pip = event_tree->GetBranch("nfpip");
    TBranch *b_pim = event_tree->GetBranch("nfpim");
    TBranch *b_p   = event_tree->GetBranch("nfp");

    // Charged current counters
    int ncc       = 0;
    int ncc0pi    = 0;
    int ncc0pi0p  = 0;
    int ncc0pi1p  = 0;
    int ncc0pi2p  = 0;
    int ncc0pi3p  = 0;
    int ncc0pimp  = 0;
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
    int nccqel    = 0;
    int nccdis    = 0;
    int nccmec    = 0;

    // Neutral current counters
    int nnc       = 0;
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
    int nncqel    = 0;
    int nncdis    = 0;
    int nncmec    = 0;

    // Get the number of events which contain final state muons
    int n_values = event_tree->GetEntries(); // Number of entries to loop over
    
    // Loop and count for various conditions
    for( int i = 0; i < n_values; ++i){
    
        // Get the current entry
        event_tree->GetEntry(i);
    
        // Charged current
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 ){ 
            ++ncc;
        }

        // Charged current
        // CC0Pi
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 ){
            
            ++ncc0pi;
        }

        // Charged current
        // CC0Pi 0p
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 
             && b_p->GetLeaf("nfp")->GetValue() == 0 ){
            
            ++ncc0pi0p;
        }

        // Charged current
        // CC0Pi 1p
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 
             && b_p->GetLeaf("nfp")->GetValue() == 1 ){
            
            ++ncc0pi1p;
        }

        // Charged current
        // CC0Pi 2p
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 
             && b_p->GetLeaf("nfp")->GetValue() == 2 ){
            
            ++ncc0pi2p;
        }

        // Charged current
        // CC0Pi 3p
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 
             && b_p->GetLeaf("nfp")->GetValue() == 3 ){
            
            ++ncc0pi3p;
        }

        // Charged current
        // CC0Pi > 3p
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_pi0->GetLeaf("nfpi0")->GetValue() == 0 
             && b_pip->GetLeaf("nfpip")->GetValue() == 0 
             && b_pim->GetLeaf("nfpim")->GetValue() == 0 
             && b_p->GetLeaf("nfp")->GetValue() > 3 ){
            
            ++ncc0pimp;
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
       
        // CCQE
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_qel->GetLeaf("qel")->GetValue() !=0 ){
            
            ++nccqel;
        }
        // CCDIS
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_dis->GetLeaf("dis")->GetValue() !=0 ){
            
            ++nccdis;
        }
        // CCMEC
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 
             && b_mec->GetLeaf("mec")->GetValue() !=0 ){
            
            ++nccmec;
        }
        //----------------------------------------------------------
        // Neutral current
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 ){ 
            ++nnc;
        }

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
        
        // NCQE
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_qel->GetLeaf("qel")->GetValue() !=0 ){
            
            ++nncqel;
        }
        
        // NCDIS
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_dis->GetLeaf("dis")->GetValue() !=0 ){
            
            ++nncdis;
        }
        
        // NCMEC
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 
             && b_mec->GetLeaf("mec")->GetValue() !=0 ){
            
            ++nncmec;
        }
    }


    // Now fill the map with the values
    // Do this in such a way that the key can be printed straight into a table
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi0p) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi1p) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi2p) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pi3p) );
    n_cc_fsi.push_back ( TMath::Floor( norm * ncc0pimp) );
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

    n_cc_proc.push_back ( TMath::Floor( norm * nccqel) );
    n_cc_proc.push_back ( TMath::Floor( norm * nccmec) );
    n_cc_proc.push_back ( TMath::Floor( norm * nccdis) );
    n_cc_proc.push_back ( TMath::Floor( norm * ncccoh) );

    n_nc_fsi.push_back ( TMath::Floor( norm * nnc) );
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
    
    n_nc_proc.push_back ( TMath::Floor( norm * nncqel) );
    n_nc_proc.push_back ( TMath::Floor( norm * nccmec) );
    n_nc_proc.push_back ( TMath::Floor( norm * nccdis) );
    n_nc_proc.push_back ( TMath::Floor( norm * ncccoh) );

    file << " -----------SBND----------" << endl; 
    file << " ------------------------- " << endl;
    file << setprecision(5) << " CC Inc   : " << ncc * norm << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi * norm << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi * norm << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi * norm << endl;
    file << setprecision(5) << " CC0Pi    : " << ncc0pi * norm << endl;
    file << setprecision(5) << " CC0Pi 0p : " << ncc0pi0p * norm << endl;
    file << setprecision(5) << " CC0Pi 1p : " << ncc0pi1p * norm << endl;
    file << setprecision(5) << " CC0Pi 2p : " << ncc0pi2p * norm << endl;
    file << setprecision(5) << " CC0Pi 3p : " << ncc0pi3p * norm << endl;
    file << setprecision(5) << " CC0Pi mp : " << ncc0pimp * norm << endl;
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
    file << setprecision(5) << " NC Inc   : " << nnc * norm  << endl;
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
}

// -------------------------------------------------------------------------
//                      Make final state tables
// -------------------------------------------------------------------------
void MakeTable( const m_outer &n_cc_vect,
                const m_outer &n_nc_vect,
                const m_outer &n_cc_proc_vect,
                const m_outer &n_nc_proc_vect,
                const vector< string > interactions,
                const vector< string > processes,
                ostream &file ){
    
    // Iterators 
    typedef m_outer::const_iterator map_it;

    // Make a table from an input map - of - maps to print nice things
    // Number of columns and rows to be made
    int n_models, n_interactions, n_processes;
    n_models       = n_cc_vect.size();
    n_interactions = interactions.size();
    n_processes    = processes.size();
    
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
    file << " & " ;
    for( int i = 0; i < n_models - 1; ++i ){
        file << "\\rotatebox{90}{ \\textbf{ " << m_names[i] << " } } & ";
    }
    file << "\\rotatebox{90}{ \\textbf{ " << m_names[n_models - 1] << " } } \\\\" << endl;
    file << " \\hline " << endl;
 
    // -------------------------------------------------------------------
    //      Loop over the outer maps and fill the rest of the table                                   
    // -------------------------------------------------------------------
   
    file << " \\multicolumn{ " << n_models + 1 << " }{ | c | }{ \\textit{ Charged Current } } \\\\ " << endl;
    
    file << " \\hline " << endl;
    
    file << " \\textbf{Hadronic Final State } & " ;
    
    for ( int i = 0; i < n_models - 1; ++i ){
        file << " & " << endl;
    }
    
    file << " \\\\ " << endl;

    for( int i = 0; i < n_interactions; ++i ){

        file << interactions[i] << " & ";

        for( map_it it = n_cc_vect.begin(), it1 = --n_cc_vect.end(); it != it1; ++it ){

            file << setprecision(5) <<  it->second[i] << " & ";
        }
        map_it it_cc = --n_cc_vect.end();
        file << setprecision(5) <<  it_cc->second[i] << " \\\\ " << endl;
    } 
   
    file << " \\textbf{Physical Process } &" << endl;
    for ( int i = 0; i < n_models - 1; ++i ){
        file << " & " << endl;
    }
    file << " \\\\ " << endl;


    for( int i = 0; i < n_processes; ++i ){

        file << processes[i] << " & ";

        for( map_it it = n_cc_proc_vect.begin(), it1 = --n_cc_proc_vect.end(); it != it1; ++it ){

            file << setprecision(5) <<  it->second[i] << " & ";
        }
        map_it it_cc = --n_cc_proc_vect.end();
        file << setprecision(5) <<  it_cc->second[i] << " \\\\ " << endl;
    } 
   
    file << " \\hline " << endl;
    
    file << " \\multicolumn{ " << n_models + 1 << " }{ | c | }{ \\textit{ Neutral Current } } \\\\ " << endl;
    
    file << " \\hline " << endl;
    
    file << " \\textbf{Hadronic Final State } & " ;
    
    for ( int i = 0; i < n_models - 1; ++i ){
        file << " & " << endl;
    }
    
    file << " \\\\ " << endl;


    for( int i = 0; i < n_interactions; ++i ){

        if ( i < 2 ){
            file << interactions[i] << " & ";

            for( map_it it2 = n_nc_vect.begin(), it3 = --n_nc_vect.end(); it2 != it3; ++it2 ){
        
                file << setprecision(5) <<  it2->second[i] << " & ";
            }
         
            map_it it_nc = --n_nc_vect.end();
            file << setprecision(5) <<  it_nc->second[i] << " \\\\ " << endl;
    
        }
        if ( i > 6 ){
            file << interactions[i] << " & ";

            for( map_it it2 = n_nc_vect.begin(), it3 = --n_nc_vect.end(); it2 != it3; ++it2 ){
        
                file << setprecision(5) <<  it2->second[i - 5] << " & ";
            }
         
            map_it it_nc = --n_nc_vect.end();
            file << setprecision(5) <<  it_nc->second[i - 5] << " \\\\ " << endl;
    
        }
    }
    
    
    file << " \\textbf{Physical Process } &" << endl;
    for ( int i = 0; i < n_models - 1; ++i ){
        file << " & " << endl;
    }
    file << " \\\\ " << endl;


    for( int i = 0; i < n_processes; ++i ){

        file << processes[i] << " & ";

        for( map_it it = n_nc_proc_vect.begin(), it1 = --n_nc_proc_vect.end(); it != it1; ++it ){

            file << setprecision(5) <<  it->second[i] << " & ";
        }
        map_it it_nc = --n_nc_proc_vect.end();
        file << setprecision(5) <<  it_nc->second[i] << " \\\\ " << endl;
    } 
   
    file << " \\hline " << endl;
    
    // -------------------------------------------------------------------
    //                    End the table                                   
    // -------------------------------------------------------------------
 
    // End the tabular environment
    file << "\\end{longtable}" << endl;
    

}

