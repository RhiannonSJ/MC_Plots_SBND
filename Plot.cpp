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
double Plot::Norm(int n_events, TFile &xsec_file, TFile &flux_file){
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
void Plot::HistStacker ( vector< TH1D* >   &hists,
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
void Plot::RecoNuE( TTree *event_tree,
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
void Plot::FSPNumbers( TTree *event_tree,
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
void Plot::FSINumbers( TTree *event_tree,
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
void Plot::MakeTable( const m_outer &n_cc_vect,
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

