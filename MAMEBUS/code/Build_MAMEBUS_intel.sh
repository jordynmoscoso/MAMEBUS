source /opt/intel/bin/compilervars.sh intel64
icc -parallel -par-report3 -O3 -c rktvd.c defs.c
icc -parallel -par-report3 -O3 mamebus.c -o mamebus.exe -Xlinker rktvd.o defs.o
rm *.o
