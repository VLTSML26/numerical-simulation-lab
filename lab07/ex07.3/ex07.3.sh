# Preparations
./clean.sh
make
# SOLID 
echo "==== SOLID ===="
cp ../input_NVE/input.solid ../input_NVE/input.dat
cp ../input_NVE/config.fcc ../input_NVE/config.0
echo "=== start"
./NVE.exe start solid
echo "=== measure"
cp ../input_NVE/config.final ../input_NVE/config.0
cp ../input_NVE/old.final ../input_NVE/old.0
./NVE.exe measure solid
# LIQUID
echo "==== LIQUID ===="
cp ../input_NVE/input.liquid ../input_NVE/input.dat
cp ../input_NVE/config.fcc ../input_NVE/config.0
echo "=== start"
./NVE.exe start liquid
echo "=== measure"
cp ../input_NVE/config.final ../input_NVE/config.0
cp ../input_NVE/old.final ../input_NVE/old.0
./NVE.exe measure liquid
# GAS
echo "==== GAS ===="
cp ../input_NVE/input.gas ../input_NVE/input.dat
cp ../input_NVE/config.fcc ../input_NVE/config.0
echo "=== start"
./NVE.exe start gas
echo "=== measure"
cp ../input_NVE/config.final ../input_NVE/config.0
cp ../input_NVE/old.final ../input_NVE/old.0
./NVE.exe measure gas
make clean
