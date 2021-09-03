# Preparations
./clean.sh
make
# Solid
echo "===== SOLID"
cp ../input/input.solid ../input/input.dat
cp ../input/config.fcc ../input/config.0
echo "=== start"
./NVE.exe start
for i in {1..5}
do
	echo "=== equilibrate $i"
	cp ../input/config.final ../input/config.0
	cp ../input/old.final ../input/old.0
	./NVE.exe equilibrate
done
echo "=== measure"
cp ../input/config.final ../input/config.0
cp ../input/old.final ../input/old.0
./NVE.exe measure
# Liquid
echo "===== LIQUID"
cp ../input/input.liquid ../input/input.dat
cp ../input/config.fcc ../input/config.0
echo "=== start"
./NVE.exe start
for i in {1..5}
do
	echo "=== equilibrate $i"
	cp ../input/config.final ../input/config.0
	cp ../input/old.final ../input/old.0
	./NVE.exe equilibrate
done
echo "=== measure"
cp ../input/config.final ../input/config.0
cp ../input/old.final ../input/old.0
./NVE.exe measure
# Gas
echo "===== GAS"
cp ../input/input.gas ../input/input.dat
cp ../input/config.fcc ../input/config.0
echo "=== start"
./NVE.exe start
for i in {1..5}
do
	echo "=== equilibrate $i"
	cp ../input/config.final ../input/config.0
	cp ../input/old.final ../input/old.0
	./NVE.exe equilibrate
done
echo "=== measure"
cp ../input/config.final ../input/config.0
cp ../input/old.final ../input/old.0
./NVE.exe measure
echo "measure.out columns: 0-step 1-epot 2-ekin 3-temp 4-etot"
echo "measure.out rows: phases are appended in order solid-liquid-gas"
make clean
