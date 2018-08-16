// Header file for defs.c
#ifndef _DEFS_H_
#define _DEFS_H_

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "nsode.h"

// Max name length for files containing input parameter arrays
#define MAX_PARAMETER_FILENAME_LENGTH 256

// Time-integration method identifiers
#define TIMESTEPPING_RKTVD1 0   // First-order total-variation-diminishing Runge-Kutta
#define TIMESTEPPING_RKTVD2 1   // Second-order total-variation-diminishing Runge-Kutta
#define TIMESTEPPING_RKTVD3 2   // Third-order total-variation-diminishing Runge-Kutta

// Biogeochemical method identifiers
#define BGC_NONE 0              // No biogeochemistry
#define BGC_NITRATEONLY 1       // Single nitrate model
#define BGC_NPZD 2              // Size structured NPZD model

// Tracer advection scheme identifiers
#define ADVECTION_CENTERED 0    // Fluxes based on mean phi at cell boundaries
#define ADVECTION_KT00 1        // Kurganov-Tadmor scheme for conservation laws (see K&T 2000)

// Momentum scheme identifiers
#define MOMENTUM_NONE 0         // No explicit momentum time stepping, mean velocities have prescribed boundary layer structure only
#define MOMENTUM_TTTW 1         // Momentum evolves under time-dependend turbulent thermal wind approximation

// Handy constants
#define _PI  3.14159265358979323846
#define _2PI 6.28318530717958647692

// Efficient way to square and cube numbers
#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// Alternatives for simple powers
#define POW2(x)  ((x)*(x))
#define POW3(x)  ((x)*(x)*(x))
#define POW4(x)  ((x)*(x)*(x)*(x))
#define POW5(x)  ((x)*(x)*(x)*(x)*(x))
#define POW6(x)  ((x)*(x)*(x)*(x)*(x)*(x))
#define POW7(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW8(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW9(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW10(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))

#define MATALLOC3(mat3,Nx,Ny,Nz)                                \
    mat3 = matalloc3(Nx,Ny,Nz);                                 \
    if (mat3 == NULL)                                           \
    {                                                           \
      fprintf(stderr,"ERROR: Unable to allocate memory\r\n");   \
      return 0;                                                 \
    }

#define MATALLOC(mat,Nx,Ny)                                     \
    mat = matalloc(Nx,Ny);                                      \
    if (mat == NULL)                                            \
    {                                                           \
      fprintf(stderr,"ERROR: Unable to allocate memory\r\n");   \
      return 0;                                                 \
    }

#define VECALLOC(vec,N)                                         \
    vec = vecalloc(N);                                          \
    if (vec == NULL)                                            \
    {                                                           \
      fprintf(stderr,"ERROR: Unable to allocate memory\r\n");   \
      return 0;                                                 \
    }

#define CHECKALLOC(ptr,size)                                    \
    ptr = malloc(size);                                         \
    if (ptr == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\r\n"); \
        return 0;                                               \
    }

// Input parameter data
typedef struct paramdata
{
    char * name;
    char * type;
    void * pvar;
    bool read;
}
paramdata;

// Vertical stretching function
real stretch_ROMS (real sigma, real h_c, real theta_s, real theta_b, real h_b);

// To define input parameters
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read);
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm);

// Data I/O
void printMatrix (FILE * outfile, real ** mat, uint m, uint n);
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm);
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm);
bool readVector (char * fname, real * vec, uint len, FILE * errstrm);

// Wrappers for allocating/freeing arrays
real *** matalloc3 (uint m1, uint m2, uint m3);
void matfree3 (real *** mat);
real ** matalloc (uint m, uint n);
void matfree (real ** mat);
real * vecalloc (uint m);
void vecfree (real * vec);
void vec2mat (real * vec, real *** pmat, uint m, uint n);
void vec2mat3 (real * vec, real **** pmat3, uint Nr, uint Nc, uint Ns);

// Functions for limiting
real max3 (real v1, real v2, real v3);
real min3 (real v1, real v2, real v3);
real minmod (real v1, real v2, real v3);
real limMin (real val, real minVal);

// Kahan summation
real kahanSum (real * vec, int N);

// Thomas algorithm
void thomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n);

#endif // #defined _DEFS_H_
