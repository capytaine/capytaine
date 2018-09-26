"%PYTHON%" setup.py build_ext
"%PYTHON%" setup.py build_ext --inplace
"%PYTHON%" setup.py install
if errorlevel 1 exit 1

