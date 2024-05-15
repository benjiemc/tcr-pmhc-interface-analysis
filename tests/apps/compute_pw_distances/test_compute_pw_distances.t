Test defaults
  $ python -m tcr_pmhc_interface_analysis.apps.compute_pw_distances --log-level error -o test $TESTDIR/data

  $ diff $TESTDIR/reference/structure_names.txt test/structure_names.txt

  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr1_alpha_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr1_alpha_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"
  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr2_alpha_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr2_alpha_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"
  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr3_alpha_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr3_alpha_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"
  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr1_beta_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr1_beta_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"
  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr2_beta_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr2_beta_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"
  $ python -c "import numpy as np; import os; \
  > test_dir = os.environ['TESTDIR']; \
  > test_vals = np.loadtxt('test/cdr3_beta_distance_matrix.txt'); \
  > ref_vals = np.loadtxt(f'{test_dir}/reference/cdr3_beta_distance_matrix.txt'); \
  > np.testing.assert_array_almost_equal(test_vals, ref_vals)"