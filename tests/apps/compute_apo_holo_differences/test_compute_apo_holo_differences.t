Standard use case for app
  $ python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences \
  > -o test_tcr_apo_holo_fw_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-5 test_tcr_apo_holo_fw_align.csv > test_entries
  $ cut -d, -f1-5 $TESTDIR/reference/tcr_apo_holo_fw_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f6 test_tcr_apo_holo_fw_align.csv | sed 1d > test_values
  $ cut -d, -f6 $TESTDIR/reference/tcr_apo_holo_fw_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Now aligning loops before computing
  $ python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences \
  > --align-loops \
  > -o test_tcr_apo_holo_loop_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-5 test_tcr_apo_holo_loop_align.csv > test_entries
  $ cut -d, -f1-5 $TESTDIR/reference/tcr_apo_holo_loop_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f6 test_tcr_apo_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f6 $TESTDIR/reference/tcr_apo_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

Looking at Per-Residue Changes
  $ python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences \
  > --log-level error \
  > --align-loops \
  > --per-residue \
  > -o test_tcr_apo_per_res_holo_loop_align.csv \
  > $TESTDIR/data

  $ cut -d, -f1-8 test_tcr_apo_per_res_holo_loop_align.csv > test_entries
  $ cut -d, -f1-8 $TESTDIR/reference/tcr_apo_per_res_holo_loop_align.csv > reference_entries
  $ diff test_entries reference_entries

  $ cut -d, -f9 test_tcr_apo_per_res_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f9 $TESTDIR/reference/tcr_apo_per_res_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f10 test_tcr_apo_per_res_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f10 $TESTDIR/reference/tcr_apo_per_res_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f11 test_tcr_apo_per_res_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f11 $TESTDIR/reference/tcr_apo_per_res_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"

  $ cut -d, -f12 test_tcr_apo_per_res_holo_loop_align.csv | sed 1d > test_values
  $ cut -d, -f12 $TESTDIR/reference/tcr_apo_per_res_holo_loop_align.csv | sed 1d > reference_values
  $ python -c "import numpy as np; test_vals = np.loadtxt('test_values'); ref_vals = np.loadtxt('reference_values'); np.testing.assert_array_almost_equal(test_vals, ref_vals)"