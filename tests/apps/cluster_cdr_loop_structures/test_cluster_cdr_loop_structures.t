Test app.
  $ python -m tcr_pmhc_interface_analysis.apps.cluster_cdr_loop_structures \
  > -o test.csv \
  > $TESTDIR/data/structure_names.txt \
  > $TESTDIR/data/*_distance_matrix.txt

  $ diff test.csv $TESTDIR/reference/clusters.csv
