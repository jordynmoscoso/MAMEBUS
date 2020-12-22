gcc -O3 -ffast-math -c ab.c defs.c
gcc -O3 -ffast-math -lm mamebus.c -o mamebus.exe -Xlinker ab.o defs.o
rm *.o
