./clean.sh
cp config.fcc config.0
make
echo "=== start"
./MolDyn.exe start
for i in {1..5}
do
	echo "=== equilibrate $i"
	cp config.final config.0
	cp old.final old.0
	./MolDyn.exe equilibrate
done
for j in {1..1}
do
	echo "=== measure $j"
	cp config.final config.0
	cp old.final old.0
	./MolDyn.exe measure
done
echo "measure.out: 0-step 1-temp 2-etot 3-ekin 4-epot"
