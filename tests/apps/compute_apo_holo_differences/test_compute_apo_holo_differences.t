Standard use case for app on TCRs
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --select-entities tcr \
  > -o test_tcr_apo_holo_fw_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-5 test_tcr_apo_holo_fw_align.csv > test_entries
  $ cut -d, -f1-5 $TESTDIR/reference/tcr_apo_holo_fw_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f6 test_tcr_apo_holo_fw_align.csv | sed 1d > test_values
  $ cut -d, -f6 $TESTDIR/reference/tcr_apo_holo_fw_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Now aligning TCR loops before computing
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --select-entities tcr \
  > --align-entities \
  > -o test_tcr_apo_holo_loop_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-5 test_tcr_apo_holo_loop_align.csv > test_entries
  $ cut -d, -f1-5 $TESTDIR/reference/tcr_apo_holo_loop_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f6 test_tcr_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f6 $TESTDIR/reference/tcr_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Looking at Per-Residue Changes in TCR loops
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --log-level error \
  > --select-entities tcr \
  > --align-entities \
  > --per-residue \
  > --per-residue-measurements rmsd ca_distance chi_angle_change com_distance \
  > -o test_tcr_per_res_apo_holo_loop_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-8 test_tcr_per_res_apo_holo_loop_align.csv > test_entries
  $ cut -d, -f1-8 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f9 test_tcr_per_res_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f9 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f10 test_tcr_per_res_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f10 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f11 test_tcr_per_res_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f11 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f12 test_tcr_per_res_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f12 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Running only on Residue COM
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --log-level error \
  > --select-entities tcr \
  > --align-entities \
  > --per-residue \
  > --per-residue-measurements com_distance \
  > -o test_tcr_per_res_apo_holo_loop_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-8 test_tcr_per_res_apo_holo_loop_align.csv > test_entries
  $ cut -d, -f1-8 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f9 test_tcr_per_res_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f12 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Including anchors in the calculations
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --log-level error \
  > --select-entities tcr \
  > --align-entities \
  > --per-residue \
  > --per-residue-measurements rmsd ca_distance chi_angle_change com_distance \
  > --num-anchors 5 \
  > -o test_tcr_per_res_apo_holo_loop_align_anchors.csv \
  > $TESTDIR/data

  $ cut -d, -f1-8 test_tcr_per_res_apo_holo_loop_align_anchors.csv > test_entries
  $ cut -d, -f1-8 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align_anchors.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f9 test_tcr_per_res_apo_holo_loop_align_anchors.csv | sed 1d > test_values
  $ cut -d, -f9 $TESTDIR/reference/tcr_per_res_apo_holo_loop_align_anchors.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Computing D-scores
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --log-level error \
  > --select-entities tcr \
  > --align-entities \
  > --per-residue \
  > --per-residue-measurements d_score \
  > --num-anchors 6 \
  > -o test_tcr_per_res_apo_holo_d_score.csv \
  > $TESTDIR/data

  $ cut -d, -f1-8 test_tcr_per_res_apo_holo_d_score.csv > test_entries
  $ cut -d, -f1-8 $TESTDIR/reference/tcr_per_res_apo_holo_d_score.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f9 test_tcr_per_res_apo_holo_d_score.csv | sed 1d > test_values
  $ cut -d, -f9 $TESTDIR/reference/tcr_per_res_apo_holo_d_score.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Now on the MHC side
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --select-entities pmhc \
  > -o test_pmhc_apo_holo.csv \
  > $TESTDIR/data

  $ cut -d, -f1-4 test_pmhc_apo_holo.csv > test_entries
  $ cut -d, -f1-4 $TESTDIR/reference/pmhc_apo_holo.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f5 test_pmhc_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f5 $TESTDIR/reference/pmhc_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

... and adding TCR contact information
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --select-entities pmhc \
  > --pmhc-tcr-contact-residues 1058 1059 1061 1061A 1062 1063 1065 1066 1068 1069 1070 1072 1072A 1073 1076 1077 1080 18 58 59 61 62 65 66 68 69 70 72 73 75 76 79 80 \
  > -o test_pmhc_tcr_contact_apo_holo.csv \
  > $TESTDIR/data

  $ cut -d, -f1-5 test_pmhc_tcr_contact_apo_holo.csv > test_entries
  $ cut -d, -f1-5 $TESTDIR/reference/pmhc_tcr_contact_apo_holo.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f6 test_pmhc_tcr_contact_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f6 $TESTDIR/reference/pmhc_tcr_contact_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Looking at per-residue changes of pMHCs
  $ python -m tcr_pmhc_interface_analysis.apps.compute_apo_holo_differences \
  > --log-level error \
  > --select-entities pmhc \
  > --per-residue \
  > --per-residue-measurements rmsd ca_distance chi_angle_change com_distance \
  > -o test_pmhc_per_res_apo_holo.csv \
  > $TESTDIR/data

  $ cut -d, -f1-7 test_pmhc_per_res_apo_holo.csv > test_entries
  $ cut -d, -f1-7 $TESTDIR/reference/pmhc_per_res_apo_holo.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f8 test_pmhc_per_res_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f8 $TESTDIR/reference/pmhc_per_res_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f9 test_pmhc_per_res_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f9 $TESTDIR/reference/pmhc_per_res_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f10 test_pmhc_per_res_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f10 $TESTDIR/reference/pmhc_per_res_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f11 test_pmhc_per_res_apo_holo.csv | sed 1d > test_values
  $ cut -d, -f11 $TESTDIR/reference/pmhc_per_res_apo_holo.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"