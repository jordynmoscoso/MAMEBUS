source /opt/intel/bin/compilervars.sh intel64
icc -parallel -par-report3 -O3 -c ab.c defs.c
icc -parallel -par-report3 -O3 mamebus.c -o mamebus.exe -Xlinker ab.o defs.o
rm *.o
