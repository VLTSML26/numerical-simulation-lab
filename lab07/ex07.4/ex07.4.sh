./clean.sh
make
echo "======== NVT SIMULATIONS:"
for phase in solid liquid gas
do
	echo $phase
	cp ../input_NVT/input.$phase ../input_NVT/input.dat
	cp ../input_NVT/config.fcc ../input_NVT/config.0
	./NVT.exe $phase
done
echo "======== NVE SIMULATIONS:"
for phase in solid liquid gas
do
	echo $phase
	cp ../input_NVE/input.$phase ../input_NVE/input.dat
	cp ../input_NVE/config.fcc ../input_NVE/config.0
	echo "== start"
	./NVE.exe start $phase
	for i in {1..5}
	do
		echo "== equilibrate $i"
		cp ../input_NVE/config.final ../input_NVE/config.0
		cp ../input_NVE/old.final ../input_NVE/old.0
		./NVE.exe equilibrate $phase
	done
	echo "== measure"
	cp ../input_NVE/config.final ../input_NVE/config.0
	cp ../input_NVE/old.final ../input_NVE/old.0
	./NVE.exe measure $phase
done
echo "measure.out columns: 0-step 1-epot 2-ekin 3-temp 4-etot"
echo "measure.out rows: phases are appended in order solid-liquid-gas"
make clean
