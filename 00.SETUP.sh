#
# Do not try to run needs to be sourced
#   . 00.SETUP.sh

#
UNAME=$(uname)
HOST=$(hostname)

if [ $HOST == "terra" ]; then
  /juno/work/bic/socci/opt/common/CentOS_7/python/python-3.9.7/bin/python3 -m venv venv
else
  python3 -m venv venv
fi

. venv/bin/activate
pip install --upgrade pip
pip install numpy==1.23.0
pip install matplotlib==3.9.0

cd code/idr
pip install scipy==1.13.1

# TERRA
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

