#!/bin/bash
# test the GUI
# unitary test 
check_exit_code()
{
exit_code=$1
test_script=$2
if [ ${exit_code} -ne 0 ]; then
   echo ">>>>>> ERROR ${test_script}"
else
   echo ">>>>>> OK ${test_script}"
fi  
}
want_coverage="yes"
if [ "$1" == "nocoverage" ] ; then
echo "do not use coverage"
want_coverage="no"
fi

datetest=`date +%Y%m%d`

for var in RTTOV_GUI_PREFIX RTTOV_GUI_PROFILE_DIR RTTOV_GUI_COEFF_DIR RTTOV_GUI_EMISS_DIR
do
  echo $var
  envvar=`echo $var`
  if [ -z "${envvar}" ] ; then
    echo "environment variable $var not set"
    exit 1
  fi
done
# define the WRK_DIR for the test
echo "RTTOV_GUI_WRK_DIR=$RTTOV_GUI_WRK_DIR"
if [ -z "${RTTOV_GUI_WRK_DIR}" ] ; then
  echo "environment variable RTTOV_GUI_WRK_DIR not set"
  exit 1
else 
   # make a specific working directory for the test suite (in order not to disturb interactive usage of the GUI)
   mkdir -p $RTTOV_GUI_WRK_DIR/test_wrk_dir
fi
export RTTOV_GUI_WRK_DIR=$RTTOV_GUI_WRK_DIR/test_wrk_dir


log_file=${RTTOV_GUI_WRK_DIR}/test_${datetest}.log
echo "test_gui $datetest" > ${log_file}
echo "RTTOV_GUI_WRK_DIR=$RTTOV_GUI_WRK_DIR" >> ${log_file}

cd $RTTOV_GUI_PREFIX
coverage_report_file=coverage.report.log
 
has_coverage="no"
command -v coverage  > /dev/null 2>&1 &&  has_coverage="yes"
echo "has_coverage=${has_coverage}"

for mytest in test_options test_project test_karchive test_kmatrix test_run test_aerosols test_atlas test_clouds test_full test_nlte test_ozone test_runpc test_updates12_2
do
  echo "testing ${mytest}" | tee -a ${log_file}
if [ "${has_coverage}" == "yes" -a "${want_coverage}" == "yes" ]
then
  rm -f .coverage
  coverage run -m test.${mytest}   >> ${log_file} 2>&1 
  check_exit_code $? ${mytest} | tee -a ${log_file}
  echo "testing ${mytest}" >> ${coverage_report_file}
  coverage report -m >> ${coverage_report_file} 2>&1
else
  python -m test.${mytest} >> ${log_file} 2>&1
  check_exit_code $? ${mytest} | tee -a ${log_file}
fi
  echo | tee -a ${log_file}
done


