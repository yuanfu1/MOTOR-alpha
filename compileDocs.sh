source pathEnv.sh

# Build the doxygen api docs.
cd ${ETC_DIR}/doxygen
doxygen MOTOR.DA.Doxyfile

echo $?
if [ $? -ne 0 ]; then
    echo "Failed to build the doxygen api docs."
    exit 1
fi

# Build the sphnix docs.
echo ${ETC_DIR}/sphinx
cd ${ETC_DIR}/sphinx
make html

if [ $? -ne 0 ]; then
    echo "Failed to build the sphinx api docs."
    exit 2
fi

# Copy to build file
rm -rf ${BIN_DIR}/doxygen_html
rm -rf ${BIN_DIR}/sphinx_html
cp -rf ${ETC_DIR}/sphinx/build/html ${BIN_DIR}/sphinx_html
cp -rf ${ETC_DIR}/doxygen/html ${BIN_DIR}/sphinx_html/doxygen_html

# Output
echo =====================================================
echo Build docs of MOTOR-DA successful!
echo API docs: ${BIN_DIR}/sphinx_html/doxygen_html/index.html
echo sphnix docs: ${BIN_DIR}/sphinx_html/index.html
