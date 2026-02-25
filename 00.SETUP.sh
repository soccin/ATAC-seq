/juno/work/bic/socci/opt/common/CentOS_7/python/python-3.9.7/bin/python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install numpy==1.23.0
pip install matplotlib==3.9.0
pip install MACS2==2.2.9.1
cd code/idr
pip install scipy==1.13.1
python3 setup.py install
cd ../..
deactivate

