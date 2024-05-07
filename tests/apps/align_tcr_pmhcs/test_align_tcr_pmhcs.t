Align TCR-pMHC based on holo structures.
  $ python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs -o test $TESTDIR/data

  $ diff test/apo_holo_summary.csv $TESTDIR/reference/apo_holo_summary.csv

  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1kgc_D-E_tcr.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/1kgc_D-E_tcr.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1m05_A-B-E_pmhc.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/1m05_A-B-E_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1m05_C-D-F_pmhc.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/1m05_C-D-F_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1mi5_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/1mi5_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/3sko_A-B-C_pmhc.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/3sko_A-B-C_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/3x13_A-B-C_pmhc.pdb $TESTDIR/reference/1mi5_D-E-C-A-B_tcr_pmhc/3x13_A-B-C_pmhc.pdb