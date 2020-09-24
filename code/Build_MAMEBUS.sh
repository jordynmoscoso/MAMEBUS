gcc -O3 -ffast-math -c ab3.c defs.c
gcc -O3 -ffast-math -lm mamebus.c -o mamebus.exe -Xlinker ab3.o defs.o
rm *.o
