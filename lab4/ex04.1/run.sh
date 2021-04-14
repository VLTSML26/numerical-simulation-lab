cp config.fcc config.0
./MolDyn true
for i in {2..4}
do
	echo "step $i"
	cp config.final config.0
	cp old.final old.0
	./MolDyn false
done
