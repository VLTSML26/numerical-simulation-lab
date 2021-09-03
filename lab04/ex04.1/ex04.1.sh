# it is necessary to clean up before starting
# the simulation since we append data (ios::app)
./clean.sh
make
# copy initial configuration and config-data
cp ../input/input.solid ../input/input.dat
cp ../input/config.fcc ../input/config.0
# 1 start, 8 equilibrations and 1 measure
echo "=== start"
./ex04.1.exe start
for i in {1..8}
do
	echo "=== equilibrate $i"
	cp ../input/config.final ../input/config.0
	cp ../input/old.final ../input/old.0
	./ex04.1.exe equilibrate
done
echo "=== measure"
cp ../input/config.final ../input/config.0
cp ../input/old.final ../input/old.0
./ex04.1.exe measure
make clean
