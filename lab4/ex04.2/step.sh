for i in {1..10}
do
	echo "step $i"
	cp config.final config.0
	./MolDyn_NVE.exe
done
