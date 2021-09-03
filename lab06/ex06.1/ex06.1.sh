./clean.sh
make
echo "======== NO EXTERNAL FIELD"
echo "===== Metropolis"
./Ising.exe metropolis 20 0
rm input/config.final
echo "===== Gibbs"
./Ising.exe gibbs 20 0
rm input/config.final
echo "======== EXTERNAL FIELD h = 0.02"
echo "===== Metropolis"
./Ising.exe metropolis 20 0.02
rm input/config.final
echo "===== Gibbs"
./Ising.exe gibbs 20 0.02
make clean
