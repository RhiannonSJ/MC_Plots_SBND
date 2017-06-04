/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : February 2017
 *
 *
*/

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TObjArray.h"

// -------------------------------------------------------------------------
//                          Typedefs 
// -------------------------------------------------------------------------

typedef map< string, vector < double > > m_outer; 
typedef map< string, double > m_map;

// -------------------------------------------------------------------------
// Normalisation:
// For each model configuration, the distributions need to be normalised 
//      The flux_i will be the same for all
//      The cross sections will vary by model configuration
//      The ratio 
//          (targetdist_MB / targetdist_SBND)^2 = (450 / 110 )^2 = 16.7 
//          this accounts for the different distances of MiniBooNE and SBND 
//          from the target when using the MiniBooNE flux
//
// The arguments passed to the normalisation function are:
//      n_events   : The number of events we have generated
//                   This needs to be taken from the files
//      x_sec file : The cross section file for the relevant model config
//      flux  file : The MB (eventually SBND) flux file
// -------------------------------------------------------------------------

double Norm ( int   n_events, 
              TFile &xsec_file, 
              TFile &flux_file );

// -------------------------------------------------------------------------
// Stacking the histograms:
// The arguments passed to the stacking function are:
//      hists            : vector of TH1Ds
//      leg_entries      : vector of strings for legend entries
//      norm             : vector of normalisation values
//      x_axis, y_axis   : axis labels
//      title            : histogram title
//      file_name        : name of the file to save the plot to
//      min, mix, n_bins : the min and max histogram x bins and the number of bins
//
// It will assign a different colour to each histogram it has to stack and 
// will set the style and min and max bin positions according to the 
// smallest min and largest max bin in the vector
// It will also take the number of bins from one of the histograms - since 
// they should all be the same, or similar
// -------------------------------------------------------------------------
void HistStacker ( vector< TH1D* >   &hists,
                   vector< string >  &leg_entries,
                   vector< double >  &norm,
                   const char* title,
                   const char* file_name,
                   const char* x_axis,
                   const char* y_axis );

// -------------------------------------------------------------------------
// Calculating the reconstructed energy of the neutrinos to compare with the 
// MC value
// -------------------------------------------------------------------------
void RecoNuE( TTree *event_tree,
              vector< double > &reco_E_CC, 
              vector< double > &reco_E_NC,
              vector< double > &MC_reco_E_CC, 
              vector< double > &MC_reco_E_NC );

// -------------------------------------------------------------------------
// Make a map of the number of final state particles in each model
// configuration
// -------------------------------------------------------------------------
void FSPNumbers( TTree *event_tree,
                 m_map &n_fsp );

// -------------------------------------------------------------------------
// Make a map of the number of different final state interactions 
// in each model configuration
// -------------------------------------------------------------------------
void FSINumbers( TTree *event_tree,
                 ostream &file,
                 double norm,
                 vector< double > &n_cc_fsi,
                 vector< double > &n_nc_fsi,
                 vector< double > &n_cc_mc_fsi,
                 vector< double > &n_nc_mc_fsi );

// -------------------------------------------------------------------------
// Make a table to compare the number of different final state particles
// and interactions in each model configuration
// -------------------------------------------------------------------------
void MakeTable( const m_outer &n_cc_vect,
                const m_outer &n_nc_vect,
                const vector< string > interactions,
                ostream &file );

// -------------------------------------------------------------------------
