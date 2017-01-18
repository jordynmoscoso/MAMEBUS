gcc -O3 -ffast-math -c rktvd.c defs.c
gcc -O3 -ffast-math -lm mamebus.c -o mamebus.exe -Xlinker rktvd.o defs.o
rm *.o
