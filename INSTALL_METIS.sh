mkdir -p metis
wget -c http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
tar -zxvf metis-5.1.0.tar.gz 
cd metis-5.1.0
make config prefix=$(pwd)/../metis
make
make install
rm -rf metis-5.1.0.tar.gz