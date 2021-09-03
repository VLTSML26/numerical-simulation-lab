./clean.sh
make
echo "==== SOLID"
cp ../input_NVT/input.solid ../input_NVT/input.dat
cp ../input_NVT/config.fcc ../input_NVT/config.0
./NVT.exe
echo "==== LIQUID"
cp ../input_NVT/input.liquid ../input_NVT/input.dat
cp ../input_NVT/config.fcc ../input_NVT/config.0
./NVT.exe
echo "==== GAS"
cp ../input_NVT/input.gas ../input_NVT/input.dat
cp ../input_NVT/config.fcc ../input_NVT/config.0
./NVT.exe
rm -rf *.o
rm -rf *.exe
