#
# Do not try to run needs to be sourced
#   . 00.SETUP.sh
#
# MACS2 2.2.9.1 and IDR 2.0.3 both build under the python 3.10 on the
# IRIS PATH, so the venv uses whatever python3 is current.
#
UNAME=$(uname)

python3 -m venv venv

. venv/bin/activate
pip install --upgrade pip
pip install numpy==1.23.0
pip install matplotlib==3.9.0

cd code/idr
pip install scipy==1.13.1

# IRIS
if [ $UNAME == "Linux" ]; then
  python3 setup.py install
fi

# Mac
if [ $UNAME == "Darwin" ]; then
  echo pip install wheel
  echo pip install --no-build-isolation .
fi

cd ../..
pip install MACS2==2.2.9.1

deactivate
