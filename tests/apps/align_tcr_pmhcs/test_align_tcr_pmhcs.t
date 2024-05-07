Align TCR-pMHC based on holo structures.
  $ python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs -o test $TESTDIR/data/apo-holo

  $ diff test/apo_holo_summary.csv $TESTDIR/reference/apo-holo/apo_holo_summary.csv

  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1kgc_D-E_tcr.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/1kgc_D-E_tcr.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1m05_A-B-E_pmhc.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/1m05_A-B-E_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1m05_C-D-F_pmhc.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/1m05_C-D-F_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/1mi5_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/1mi5_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/3sko_A-B-C_pmhc.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/3sko_A-B-C_pmhc.pdb
  $ diff test/1mi5_D-E-C-A-B_tcr_pmhc/3x13_A-B-C_pmhc.pdb $TESTDIR/reference/apo-holo/1mi5_D-E-C-A-B_tcr_pmhc/3x13_A-B-C_pmhc.pdb

Align TCR-pMHC holo structures on CDR or pMHC
  $ python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs --log-level error --only-holo -o test-holo $TESTDIR/data/only-holo

  $ diff test-holo/apo_holo_summary.csv $TESTDIR/reference/only-holo/apo_holo_summary.csv

  $ diff test-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3vxu_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3vxu_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3w0w_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3w0w_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3vxu_I-J-H-F-G_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-GTYNQGGKLI-MNHEY-SMNVEV-ASSGASHEQY/3vxu_I-J-H-F-G_tcr_pmhc.pdb
  $ diff test-holo/hla_a_02_01_LLFGFPVYV/3qfj_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/hla_a_02_01_LLFGFPVYV/3qfj_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/hla_a_02_01_LLFGFPVYV/3d3v_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/hla_a_02_01_LLFGFPVYV/3d3v_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/hla_a_02_01_LLFGFPVYV/3d39_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/hla_a_02_01_LLFGFPVYV/3d39_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3qfj_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3qfj_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3d3v_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3d3v_D-E-C-A-B_tcr_pmhc.pdb
  $ diff test-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3d39_D-E-C-A-B_tcr_pmhc.pdb $TESTDIR/reference/only-holo/DRGSQS-IYSNGD-AVTTDSWGKLQ-MNHEY-SVGAGI-ASRPGLAGGRPEQY/3d39_D-E-C-A-B_tcr_pmhc.pdb