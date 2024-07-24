Test app.
  $ python -m tcr_pmhc_interface_analysis.apps.sample_ots \
  > --seed 123 \
  > -n 5 \
  > -o test.csv \
  > $TESTDIR/data

  $ diff test.csv $TESTDIR/reference/sample.csv
