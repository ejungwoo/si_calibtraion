root -q -b -l make_histograms.C            # output = data/start*.hist.root
root -q -b -l run_energy_correction.C      # output = data/start*.c0.root        containing c0 correction parameters (tree)
root -q -b -l run_slope_correction.C       # output = data/start*.c1.root        containing c1 correction parameters (tree)
root -q -b -l make_c1_histograms.C         # output = data/start*.hist_c1.root   containing c0, c1 corrected histograms
root -q -b -l run_ballistic_correction.C   # output = data/start*.c2.root        containing c2 correction parameters (tree)
root -q -b -l make_c2_histograms.C         # output = data/start*.hist_c2.root   containing c2 corrected histograms
root summary.C               # to check result
