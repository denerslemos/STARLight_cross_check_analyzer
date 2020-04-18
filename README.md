# STARLight_cross_check_analyzer

Analyzer created to make cross-check histograms for STARLight MC since the Forest code is not yet available for CMSSW11X. I've included the .cc and .py used to make the histograms.

The root output file includes:

- **hiHF_hist**: 1D histogram with transverse energy sum of HF tower (HF+ + HF-) using centrality collection;
- **etHFtowerSum_hist**: 1D histogram with transverse energy sum of HF tower (HF+ + HF-) using tower collection;
- **hiHF_vs_etHFtowerSum**: 2D histogram comparing hiHF_hist and etHFtowerSum_hist: Just to check;
- **hiNpix_hist**: 1D histogram with the number of pixel hits (using centrality collection);
- **hiNtracks_hist**: 1D histogram with the number of tracks; (using centrality collection);
- **hiHF_vs_hiNpix**: 2D histogram comparing hiHF_hist and hiNpix_hist;
- **hiHF_vs_hiNtracks**: 2D histogram comparing hiHF_hist and hiNtracks_hist;
- **hiHF_plus_hist**: 1D histogram with transverse energy sum of HF tower only in HF+ side (using centrality collection);
- **hiHF_minus_hist**: 1D histogram with transverse energy sum of HF tower only in HF- side (using centrality collection);
- **hiHF_plus_vs_minus**: 2D histogram comparing hiHF_plus_hist and hiHF_minus_hist;
- **trk_pt**: 1D histogram with transverse momentum ($p_{T}$) of the tracks in generalTracks;
- **trk_eta**: 1D histogram with pseudorapidity (&eta;) of the tracks in generalTracks;
- **trk_phi**: 1D histogram with azhimutal distribution (&phi;) of the tracks in generalTracks;
- **tower_pt**: 1D histogram with transverse momentum ($p_{T}$) of the HF in tower collection;
- **tower_eta**: 1D histogram with pseudorapidity (&eta;) of the HF in tower collection;
- **tower_phi**: 1D histogram with azhimutal distribution (&phi;) of the HF in tower collection;


 h<sub>&theta;</sub>(x)



