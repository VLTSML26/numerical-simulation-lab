./clean.sh
g++ -o MolDyn.exe main.cpp
cp config.fcc config.0
echo "start"
./MolDyn.exe start
for i in {1..1}
do
	echo "equilibrate $i"
	cp config.final config.0
	cp old.final old.0
	./MolDyn.exe equilibrate
done
for j in {1..1}
do
	echo "measure $j"
	cp config.final config.0
	cp old.final old.0
	./MolDyn.exe measure
done

