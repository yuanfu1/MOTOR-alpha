set(
  fppFiles
  "test_loadtxt_qp.fypp"
  "test_savetxt_qp.fypp"
)
fypp_f90("${fyppFlags}" "${fppFiles}" outFiles)

ADDTEST(loadtxt)
ADDTEST(savetxt)

ADDTEST(loadtxt_qp)
ADDTEST(savetxt_qp)
set_tests_properties(loadtxt_qp PROPERTIES LABELS quadruple_precision)
set_tests_properties(savetxt_qp PROPERTIES LABELS quadruple_precision)

ADDTEST(getline)
ADDTEST(npy)
ADDTEST(open)
ADDTEST(parse_mode)
