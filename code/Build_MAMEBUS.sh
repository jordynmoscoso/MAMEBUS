gcc -O3 -ffast-math -c defs.c
gcc -O3 -ffast-math -lm mamebus.c -o mamebus.exe -Xlinker defs.o
rm *.o
