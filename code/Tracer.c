 /**
  * Tracer.c
  *
  * Andrew Stewart
  *
  * Advects and mixes Nd isotopes using the streamfunction
  * output from Overturning.c.
  *
  */
#include <time.h>
#include <math.h>

#include "rktvd.h"
#include "ot_defs.h"


#define NPARAMS 19


// Current iteration data
real ** phi_w;

// Time derivative of current iteration
real ** dt_phi_w = NULL;

// Work arrays for KT scheme
real ** HHy = NULL; // Advective (hyperbolic) fluxes
real ** HHz = NULL;
real ** PPy = NULL; // Diffusive (parabolic) fluxes
real ** PPz = NULL; 
real ** dphi_dy = NULL;
real ** dphi_dz = NULL;
real ** phi_yp = NULL;
real ** phi_ym = NULL;
real ** phi_zp = NULL;
real ** phi_zm = NULL;

// Implicit solver arrays
real ** impl_A = NULL;
real ** impl_B = NULL;
real ** impl_C = NULL;
real ** impl_D = NULL;

// Parameters
real Ly = 0;
real Lz = 0;
real kaph0 = 1000;
real kapv0 = 0;

// Kurganov-Tadmor minmod-limiting parameter
real sigma = 1.4; 

// Parameter arrays
real ** psi_r = NULL; // Residual streamfunction
real ** v_r = NULL; // Residual streamfunction
real ** w_r = NULL; // Residual streamfunction
real ** kapv_psi = NULL; // Eddy vertical diffusivity at psi-gridpoints
real ** kapv_w = NULL; // Eddy vertical diffusivity at w-gridpoints
real ** kaph_psi = NULL; // Eddy lateral diffusivity at psi-gridpoints
real ** kaph_v = NULL; // Eddy lateral diffusivity at v-gridpoints
real ** kaph_w = NULL; // Eddy lateral diffusivity at w-gridpoints
real ** ss_psi = NULL; // Neutral slope at psi-gridpoints
real ** ss_v = NULL; // Neutral slope at v-gridpoints
real ** ss_w = NULL; // Neutral slope at w-gridpoints
real ** ssq_w = NULL; // Squared neutral slope at w-gridpoints
real ** phi_relax = NULL; // Tracer relaxation values
real ** T_relax = NULL; // Relaxation times
real ** wSink = NULL; // Particle sinking velocity
real ** pdRatio = NULL; // Particle to dissolved tracer ratio

// Grid size - number of grid boxes in the domain
int Ny = 0;
int Nz = 0;
int Nphi = 0;

// Grid spacings
real dy = 0;    
real dz = 0;
real _dy = 0;
real _dz = 0;
real _2dy = 0;
real _2dz = 0;
real dysq = 0;
real dzsq = 0;
real dt_dzsq = 0;

// CFL-limited time step
real cfl_min_step = 0;

// Name of the program (for error messages)
char * progname = NULL;

// Time-stepping method
const uint method_t = METHOD_RKTVD2;

// Spatial discretisation
const uint method_s = METHOD_KT;


/**
 * tderiv
 *
 * Calculates the time derivatives at all phi gridpoints. Returns the
 * largest time step that keeps both the advective and diffusive parts
 * of the buoyancy advection equation stable.
 *
 */
real tderiv (const real t, const real * data, real * dt_data, const uint numvars)
{
    // Looping variables
    int j,k;
    real ypos,zpos;

    // Construct output matrix from the dt_data vector
    memset(dt_data,0,numvars*sizeof(real));
    vec2mat(dt_data,&dt_phi_w,Ny,Nz);

    // Copy phi from the 'data' array to the work array with space for ghost points
    for (j = 0; j < Ny; j ++)
    {
        memcpy(phi_w[j+1]+1,data+j*Nz,Nz*sizeof(real));
    }

    // Some useful reminders for interpreting the code:
    // phi_w[j][k] corresponds to phi_{j,k}
    // psi_r[j][k] corresponds to psi_r{j+1/2,k+1/2}
    // FFy[j][k] corresponds to F^(y)_{j+1/2,k+1}
    // FFz[j][k] corresponds to F^(z)_{j+1,k+1/2}
    // kap_w[j][k] corresponds to kappa_{j+1,k+1/2}
    // etc.



    //////////////////////////////////////
    ///// BEGIN SETTING GHOST POINTS /////
    //////////////////////////////////////

    // Top/bottom boundaries
    for (j = 1; j <= Ny; j ++)
    {     
        phi_w[j][0] = 2*phi_w[j][1]-phi_w[j][2];
        phi_w[j][Nz+1] = 2*phi_w[j][Nz]-phi_w[j][Nz-1];
    }

    // South/north boundaries
    for (k = 1; k <= Nz; k ++)
    {     
        phi_w[0][k] = 2*phi_w[1][k]-phi_w[2][k];
        phi_w[Ny+1][k] = 2*phi_w[Ny][k]-phi_w[Ny-1][k];
    }

    // Corners (note linear extrapolation to corners is consistent in either direction)
    phi_w[0][0] = 2*phi_w[1][0]-phi_w[2][0];
    phi_w[0][Nz+1] = 2*phi_w[1][Nz+1]-phi_w[2][Nz+1];
    phi_w[Ny+1][0] = 2*phi_w[Ny][0]-phi_w[Ny-1][0];
    phi_w[Ny+1][Nz+1] = 2*phi_w[Ny][Nz+1]-phi_w[Ny-1][Nz+1];  

    ////////////////////////////////////
    ///// END SETTING GHOST POINTS /////
    ////////////////////////////////////
 


    //////////////////////////////////
    ///// BEGIN INTEGRATION CODE /////
    //////////////////////////////////
    
    // Arrays used by the Kurganov-Tadmor scheme
    if (method_s == METHOD_KT)
    {
        
#pragma parallel
        
        // Determine limited slopes at cell centres
        for (j = 1; j < Ny+1; j ++)
        {
            for (k = 1; k < Nz+1; k ++)
            {
                dphi_dy[j][k] = minmod( sigma * (phi_w[j+1][k]-phi_w[j][k]) * _dy,
                                        (phi_w[j+1][k]-phi_w[j-1][k]) * _2dy,
                                        sigma * (phi_w[j][k]-phi_w[j-1][k]) * _dy  );

                dphi_dz[j][k] = minmod( sigma * (phi_w[j][k+1]-phi_w[j][k]) * _dz,
                                        (phi_w[j][k+1]-phi_w[j][k-1]) * _2dz,
                                        sigma * (phi_w[j][k]-phi_w[j][k-1]) * _dz  );
            }
        }
        
#pragma parallel

        // Interpolate phi to cell y-faces
        for (j = 0; j < Ny+1; j ++)
        {
            for (k = 0; k < Nz; k ++)
            {
                phi_ym[j][k] = phi_w[j][k+1] + 0.5*dy*dphi_dy[j][k+1];
                phi_yp[j][k] = phi_w[j+1][k+1] - 0.5*dy*dphi_dy[j+1][k+1];
            }
        }
        
#pragma parallel
        
        // Interpolate phi to cell z-faces
        for (j = 0; j < Ny; j ++)
        {
            for (k = 0; k < Nz+1; k ++)
            {
                phi_zm[j][k] = phi_w[j+1][k] + 0.5*dz*dphi_dz[j+1][k];
                phi_zp[j][k] = phi_w[j+1][k+1] - 0.5*dz*dphi_dz[j+1][k+1];
            }
        }
    }

#pragma parallel
    
    // Calculate y-fluxes
    for (k = 0; k < Nz; k ++)
    {
        // No flux across boundaries
        HHy[0][k] = 0;
        HHy[Ny][k] = 0;
        PPy[0][k] = 0;
        PPy[Ny][k] = 0;

        for (j = 1; j < Ny; j ++)
        {      
            // Flux calculation depends on spatial discretisation method
            switch (method_s)
            {
                case METHOD_BASIC:
                {
                    HHy[j][k] = v_r[j][k] * 0.5*(phi_w[j][k+1]+phi_w[j+1][k+1]);

                    PPy[j][k] = kaph_v[j][k] * (phi_w[j+1][k+1]-phi_w[j][k+1]) * _dy
                                + 0.5*kaph_v[j][k]*ss_v[j][k] 
                                  * (phi_w[j+1][k+2]-phi_w[j+1][k]+phi_w[j][k+2]-phi_w[j][k]) * _2dz;

                    break;
                }
                case METHOD_KT:
                {
                    HHy[j][k] = 0.5*( v_r[j][k]*(phi_ym[j][k]+phi_yp[j][k]) 
                                - fabs(v_r[j][k])*(phi_yp[j][k]-phi_ym[j][k]) );

                    PPy[j][k] = kaph_v[j][k] * (phi_w[j+1][k+1]-phi_w[j][k+1]) * _dy
                                + 0.5*kaph_v[j][k]*ss_v[j][k] * (dphi_dz[j+1][k+1]+dphi_dz[j][k+1]);

                    break;
                }
                default:
                {
                    fprintf(stderr,"ERROR: Unknown spatial discretisation specified\n");
                    break;
                }              
            }    
        }        
    }

#pragma parallel
    
    // Calculate z-fluxes
    for (j = 0; j < Ny; j ++)
    {
        // No flux across top boundary
        HHz[j][Nz] = 0;
        PPz[j][Nz] = 0;
        
        // Allow advective flux across bottom boundary - here w_r is only nonzero due to particulate sinking flux
        PPz[j][0] = 0;
        // Flux calculation depends on spatial discretisation method
        switch (method_s)
        {
            case METHOD_BASIC:
            {
                HHz[j][0] = w_r[j][k] * 0.5*(phi_w[j+1][0]+phi_w[j+1][1]);
                break;
            }
            case METHOD_KT:
            {
                HHz[j][0] = 0.5*( w_r[j][0]*(phi_zm[j][0]+phi_zp[j][0])
                                 - fabs(w_r[j][0])*(phi_zp[j][0]-phi_zm[j][0]) );
                break;
            }
            default:
            {
                fprintf(stderr,"ERROR: Unknown spatial discretisation specified\n");
                break;
            }
        }
        
        // Interior gridpoints
        for (k = 1; k < Nz; k ++)
        {
            // Flux calculation depends on spatial discretisation method
            switch (method_s)
            {
                case METHOD_BASIC:
                {
                    HHz[j][k] = w_r[j][k] * 0.5*(phi_w[j+1][k]+phi_w[j+1][k+1]);

                    PPz[j][k] = (kaph_w[j][k]*ssq_w[j][k]) 
                                    * (phi_w[j+1][k+1]-phi_w[j+1][k]) * _dz
                                + 0.5*kaph_w[j][k]*ss_w[j][k] 
                                    * (phi_w[j+2][k+1]-phi_w[j][k+1]+phi_w[j+2][k]-phi_w[j][k]) * _2dy;

                    break;
                }
                case METHOD_KT:
                {
                    HHz[j][k] = 0.5*( w_r[j][k]*(phi_zm[j][k]+phi_zp[j][k]) 
                                - fabs(w_r[j][k])*(phi_zp[j][k]-phi_zm[j][k]) );

                    PPz[j][k] = (kaph_w[j][k]*ssq_w[j][k]) 
                                    * (phi_w[j+1][k+1]-phi_w[j+1][k]) * _dz
                                + 0.5*kaph_w[j][k]*ss_w[j][k] 
                                    * (dphi_dy[j+1][k+1]+dphi_dy[j+1][k]);

                    break;
                }
                default:
                {
                    fprintf(stderr,"ERROR: Unknown spatial discretisation specified\n");
                    break;
                }              
            }                        
        }
    }
    
#pragma parallel

    // Calculate time derivative of phi
    for (j = 0; j < Ny; j ++)
    {
        ypos = (j+0.5)*dy;

        for (k = 0; k < Nz; k ++)
        {
            zpos = -Lz + (k+0.5)*dz;

            // Add hyperbolic and parabolic fluxes, plus linear decay
            dt_phi_w[j][k] = ( (HHy[j][k] - HHy[j+1][k]) * _dy
                             + (HHz[j][k] - HHz[j][k+1]) * _dz 
                             + (PPy[j+1][k] - PPy[j][k]) * _dy
                             + (PPz[j][k+1] - PPz[j][k]) * _dz ); 

            // Negative relaxation time means no relaxation
            if (T_relax[j][k] >= 0)
            {           
                // Instantaneous relaxation: phi cannot change
                if (T_relax[j][k] == 0.0)
                {
                    dt_phi_w[j][k] = 0.0;
                }
                // Otherwise relax phi towards phi_relax with timescale T_relax
                else 
                {
                    dt_phi_w[j][k] -= (phi_w[j+1][k+1]-phi_relax[j][k]) / T_relax[j][k];
                }            
            }
        }
    }

    ////////////////////////////////
    ///// END INTEGRATION CODE /////
    ////////////////////////////////


    return cfl_min_step; 
}

/**
 *  printUsage
 *
 *  Prints an information message describing correct usage of the
 *  program in terms of required parameters.
 *
 */
void printUsage (void)
{
    printf
    (
        "USAGE: %s <infile> <outfile> \n"
        "  \n"
        "  where <infile> is the name of the input parameter file, and\n"
        "  <outfile> is the name of the output data file.\n"
        "  \n"
        "  The input file must specify a series of input parameters\n"
        "  in name-value form, separated by whitespace, i.e.\n"
        "  <name> <value> [<name> <value> [...]]\n"
        "  \n"
        "  name               value\n"
        "  \n"
        "  Ny                 Number of grid cells in the y-direction. Must be >0.\n"
        "  Nz                 Number of grid cells in the z-direction. Must be >0.\n"
        "  Ly                 Length of the domain. Must be > 0.\n"
        "  Lz                 Domain height. Must be >0.\n"
        "  cflFrac            CFL number. The time step dt will be chosen at each\n"
        "                     iteration such that dt=cfl*dt_max, where dt_max is the\n"
        "                     estimated maximum stable time step. Must be >0.\n"
        "  targetRes          Target L2 error between successive iterations. The\n"
        "                     calculation will stop when this condition is met.\n"
        "  maxTime            Maximum time for the integration, in\n"
        "                     case the target residual is not achieved.\n"
        "  monitorFrequency   Frequency with which to write out iteration data. If\n"
        "                     set to zero, only the final data will be written.\n"
        "                     Otherwise must be >0.\n"
        "  kaph0              Lateral diffusivity. If kaphFile is not specified,\n"
        "                     then the coefficient is set uniformly\n"
        "                     to kaph0 across the entire domain. \n"
        "                     Optional, default is 1000, must be > 0.\n"
        "  kapv0              Diapycnal diffusivity.\n"
        "                     Optional, default is 0, must be >= 0.\n"
        "  kaphFile           File containing an Ny+1 x Nz+1 matrix (column-major order)\n"
        "                     of lateral diffusion coefficients at cell\n"
        "                     corners. Optional, all elements must be > 0, default is\n"
        "                     kaph0 everywhere.\n"
        "  kapvFile           File containing an Ny+1 x Nz+1 matrix (column-major order)\n"
        "                     of diapycnal diffusion coefficients at cell\n"
        "                     corners. Optional, all elements must be >= 0, default is\n"
        "                     kapv0 everywhere.\n"
        "  slopeFile          File containing an Ny+1 x Nz+1 matrix (column-major order)\n"
        "                     of isopycnal slopes at cell corners.\n"
        "  psirFile           File containing an Ny+1 x Nz+1 matrix (column-major order)\n"
        "                     of residual streamfunction values at cell corners.\n"
        "  relaxTracerFile    File containing an Ny x Nz matrix (column-major order) of\n"
        "                     tracer relaxation concentrations. Optional, default is 1\n"
        "                     everywhere.\n"
        "  relaxTimeFile      File containing an Ny x Nz matrix (column-major order) of\n"
        "                     tracer relaxation timescales. Negative values imply no\n"
        "                     relaxation. Zero values imply instantaneous relaxation.\n"
        "                     Optional, default is -1 (no relaxation) everywhere.\n"
        "  initTracerFile     File containing an Ny x Nz matrix (column-major order)\n"
        "                     of initial tracer concentrations. Optional, default is\n"
        "                     0 everywhere.\n"
        "  wSinkFile          File containing an Ny x Nz+1 matrix (column-major order)\n"
        "                     of particle sinking velocities. Optional, default is\n"
        "                     0 everywhere.\n"
        "  pdRatioFile        File containing an Ny x Nz+1 matrix (column-major order)\n"
        "                     of particle to dissolved tracer concentration ratios.\n"
        "                     Optional, default is 0 everywhere.\n"
        " \n"
        ,
        progname
    ); 
} 



/**
 * main
 *
 * Program entry point. Initialises all of the required memory,
 * reads in input parameters, performs iterations of the required
 * numerical method, and finally cleans up.
 *
 */
int main (int argc, char ** argv)
{
    // Time parameters
    real tmin = 0.0;        // Start time
    real tmax = 0;          // End time
    real t = 0;             // Time counter
    real dt = 0;            // Time step
    real dt_s = 0;          // Data save time step
    real t_next = 0;        // Next save time
    real cflFrac = 1;       // CFL number
    real wc = 0;            // Time-interpolation weights
    real wn = 0;
 
    // Convergence parameters
    real targetRes = 0;
    real res = 0;
    uint nIters = 0;

    // File descriptor for output file
    FILE * outfile = NULL;

    // Storage array for the tracer
    real ** phi = NULL;
  
    // For time-interpolation
    real ** phi_interp = NULL;

    // Pointer to output vector - used for implicit vertical diffusion
    real ** phi_outmat = NULL;

    // Work arrays for time derivatives - expressed as vectors rather 
    // than matrices for compatibility with time integration methods
    real * phi_in = NULL; // Input to time-integration method
    real * phi_out = NULL; // Output from time-integration method
    real * phi_buf = NULL; // Buffer for time-integration method and time-interpolation

    // Stores data required for parsing input parameters
    paramdata params[NPARAMS];

    // Filename holders for input parameter arrays
    char kaphFile[MAX_PARAMETER_FILENAME_LENGTH];
    char kapvFile[MAX_PARAMETER_FILENAME_LENGTH];
    char slopeFile[MAX_PARAMETER_FILENAME_LENGTH];
    char psirFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
    char initTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char wSinkFile[MAX_PARAMETER_FILENAME_LENGTH];
    char pdRatioFile[MAX_PARAMETER_FILENAME_LENGTH];

    // Looping variables
    int i = 0;
    int j = 0;
    int k = 0;
 
    // For setting initial conditions
    real zpos = 0;

    // To store current time
    time_t now;
    FILE * tfile = NULL;

    // For calculating CFLs
    real v_max = 0;
    real w_max = 0;   
    real ydiff = 0;
    real ydiff_max = 0;
    real zdiff = 0;
    real zdiff_max = 0;
    real cfl_v = 0; // y-advection CFL
    real cfl_w = 0; // z-advection CFL
    real cfl_y = 0; // y-diffusive CFL
    real cfl_z = 0; // z-diffusive CFL

    // Define input parameter data
    setParam(params,0,"Ny","%u",&Ny,false);
    setParam(params,1,"Nz","%u",&Nz,false);
    setParam(params,2,"Ly","%lf",&Ly,false);
    setParam(params,3,"Lz","%lf",&Lz,false);
    setParam(params,4,"cflFrac","%lf",&cflFrac,false);
    setParam(params,5,"targetRes","%le",&targetRes,false);
    setParam(params,6,"maxTime","%lf",&tmax,false);
    setParam(params,7,"monitorFrequency","%lf",&dt_s,true);
    setParam(params,8,"kaph0","%lf",&kaph0,true);
    setParam(params,9,"kapv0","%lf",&kapv0,true);
    setParam(params,10,"kaphFile","%s",kaphFile,true);
    setParam(params,11,"kapvFile","%s",kapvFile,true);
    setParam(params,12,"slopeFile","%s",slopeFile,false);
    setParam(params,13,"psirFile","%s",psirFile,false);
    setParam(params,14,"relaxTracerFile","%s",relaxTracerFile,true);
    setParam(params,15,"relaxTimeFile","%s",relaxTimeFile,true);
    setParam(params,16,"initTracerFile","%s",initTracerFile,true);
    setParam(params,17,"wSinkFile","%s",wSinkFile,true);
    setParam(params,18,"pdRatioFile","%s",pdRatioFile,true);
	
    // Default file name parameters - zero-length strings
    kaphFile[0] = '\0';
    kapvFile[0] = '\0';
    slopeFile[0] = '\0';
    psirFile[0] = '\0';
    relaxTracerFile[0] = '\0';
    relaxTimeFile[0] = '\0';
    initTracerFile[0] = '\0';

    // First program argument always carries the program name
    progname = argv[0];

    // Check that file names have been specified
    if (argc < 3)
    {
        fprintf(stderr,"ERROR: Not enough input file names supplied\n");
        printUsage();
        return 0;
    }

    // Calculate elapsed time and write to a file
    tfile = fopen("time.txt","w");
    if (tfile != NULL)
    {
        time(&now);
        fprintf(tfile,"Program started at %s\n", ctime(&now));
        fflush(tfile);
    }    

    //////////////////////////////////////////
    ///// BEGIN READING INPUT PARAMETERS /////
    //////////////////////////////////////////

    // Attempt to read input parameter data. Errors will be printed
    // to stderr, and will result in 'false' being returned.
    if (!readParams(argv[1],params,NPARAMS,stderr))
    {
        printUsage();
        return 0;
    }

    // Check that required parameters take legal values
    if ( (Ny <= 0) || 
         (Nz <= 0) || 
         (Ly <= 0.0) || 
         (Lz <= 0.0) || 
         (targetRes <= 0.0) ||
         (dt_s < 0.0) ||
         (tmax <= 0.0) || 
         (cflFrac <= 0.0) ||
         (kaph0 < 0.0) ||
         (kapv0 < 0.0) ||
         (strlen(kaphFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(kapvFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(slopeFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(psirFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(relaxTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(relaxTimeFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
         (strlen(initTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) )
    {
        fprintf(stderr,"ERROR: Invalid input parameter values\n");
        printUsage();
        return 0;
    }

    // Calculate total work array size
    Nphi = Ny*Nz;
    
    // Calculate grid spacings
    dy = Ly/Ny;
    dz = Lz/Nz;
    _dy = 1/dy;
    _dz = 1/dz;
    _2dy = 1/(2*dy);
    _2dz = 1/(2*dz);
    dysq = dy*dy;
    dzsq = dz*dz;

    // If no save interval is specified, only write final data
    if (dt_s == 0.0)
    {
        dt_s = 2*tmax;
    }
 
    // Current time and next save point
    t = tmin;
    t_next = tmin + dt_s;

    // Convergence residual
    res = 10*targetRes;

    ////////////////////////////////////////
    ///// END READING INPUT PARAMETERS /////
    ////////////////////////////////////////

    

    ///////////////////////////////////
    ///// BEGIN MEMORY ALLOCATION /////
    ///////////////////////////////////

    // Allocate the 'phi_in' vector and then reference 
    // it as a matrix to create 'phi'
    VECALLOC(phi_in,Nphi);
    VECALLOC(phi_out,Nphi);
    VECALLOC(phi_buf,Nphi);
    vec2mat(phi_in,&phi,Ny,Nz);
    vec2mat(phi_buf,&phi_interp,Ny,Nz);
    vec2mat(phi_out,&phi_outmat,Ny,Nz);
  
    // Allocate parameter arrays
    MATALLOC(kapv_psi,Ny+1,Nz+1);
    MATALLOC(kapv_w,Ny,Nz+1);
    MATALLOC(kaph_psi,Ny+1,Nz+1);
    MATALLOC(kaph_v,Ny+1,Nz);
    MATALLOC(kaph_w,Ny,Nz+1);
    MATALLOC(ss_psi,Ny+1,Nz+1);
    MATALLOC(ss_v,Ny+1,Nz);
    MATALLOC(ss_w,Ny,Nz+1);
    MATALLOC(ssq_w,Ny,Nz+1);
    MATALLOC(psi_r,Ny+1,Nz+1);
    MATALLOC(v_r,Ny+1,Nz);
    MATALLOC(w_r,Ny,Nz+1);
    MATALLOC(phi_relax,Ny,Nz);
    MATALLOC(T_relax,Ny,Nz);
    MATALLOC(wSink,Ny,Nz+1);
    MATALLOC(pdRatio,Ny,Nz+1);

    // These pointers will be used to create a matrix from the input vector in 'tderiv'
    CHECKALLOC(dt_phi_w,Ny*sizeof(real *));

    // Work array for 'tderiv' function, including space for ghost points
    MATALLOC(phi_w,Ny+2,Nz+2);

    // Fluxes at cell faces
    MATALLOC(HHy,Ny+1,Nz);
    MATALLOC(HHz,Ny,Nz+1);
    MATALLOC(PPy,Ny+1,Nz);
    MATALLOC(PPz,Ny,Nz+1);  
    
    // Work arrays for estimating phi at cell faces
    MATALLOC(dphi_dy,Ny+2,Nz+2);
    MATALLOC(dphi_dz,Ny+2,Nz+2);
    MATALLOC(phi_yp,Ny+1,Nz);
    MATALLOC(phi_ym,Ny+1,Nz);
    MATALLOC(phi_zp,Ny,Nz+1);
    MATALLOC(phi_zm,Ny,Nz+1);
    
    // Implicit solver arrays
    MATALLOC(impl_A,Ny,Nz);
    MATALLOC(impl_B,Ny,Nz);
    MATALLOC(impl_C,Ny,Nz);
    MATALLOC(impl_D,Ny,Nz);

    /////////////////////////////////
    ///// END MEMORY ALLOCATION /////
    /////////////////////////////////


    

    /////////////////////////////////////
    ///// BEGIN PARAMETER DEFAULTS  /////
    /////////////////////////////////////

#pragma parallel
    
    // Ny+1 x Nz+1 arrays 
    for (j = 0; j < Ny+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            // Uniform diffusivities
            kaph_psi[j][k] = kaph0; 
            kapv_psi[j][k] = kapv0; 
         }
    }
    
#pragma parallel

    // Ny x Nz arrays 
    for (j = 0; j < Ny; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Set phi to zero initially by default
            phi[j][k] = 0;
  
            // Tracer relaxation values are all 1 by default
            phi_relax[j][k] = 1;

            // No relaxation by default
            T_relax[j][k] = -1;
         }
    }

    //////////////////////////////////
    ///// END PARAMETER DEFAULTS /////
    //////////////////////////////////



    ////////////////////////////////////////
    ///// BEGIN READING PARAMETER DATA /////
    ////////////////////////////////////////
  
    // Read input matrices and vectors
    if ( ( (strlen(kaphFile) > 0)       &&  !readMatrix(kaphFile,kaph_psi,Ny+1,Nz+1,stderr) )  ||
         ( (strlen(kapvFile) > 0)       &&  !readMatrix(kapvFile,kapv_psi,Ny+1,Nz+1,stderr) )  ||
         ( (strlen(slopeFile) > 0)      &&  !readMatrix(slopeFile,ss_psi,Ny+1,Nz+1,stderr) )  ||
         ( (strlen(psirFile) > 0)       &&  !readMatrix(psirFile,psi_r,Ny+1,Nz+1,stderr) )  ||
         ( (strlen(relaxTracerFile) > 0) && !readMatrix(relaxTracerFile,phi_relax,Ny,Nz,stderr) )  ||
         ( (strlen(relaxTimeFile) > 0) &&   !readMatrix(relaxTimeFile,T_relax,Ny,Nz,stderr) )  ||
         ( (strlen(initTracerFile) > 0) &&  !readMatrix(initTracerFile,phi,Ny,Nz,stderr) ) ||
         ( (strlen(wSinkFile) > 0) &&       !readMatrix(wSinkFile,wSink,Ny,Nz+1,stderr) ) ||
         ( (strlen(pdRatioFile) > 0) &&     !readMatrix(pdRatioFile,pdRatio,Ny,Nz,stderr) )  )
    {
        printUsage();
        return 0;
    }
  
#pragma parallel
    
    // Diffusivity must be positive
    for (j = 0; j < Ny+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            if (kaph_psi[j][k] < 0)
            {
                fprintf(stderr,"kaphFile may contain only values that are >=0.");
                printUsage();
                return 0;
            }
            if (kapv_psi[j][k] < 0)
            {
                fprintf(stderr,"kapvFile may contain only values that are >=0.");
                printUsage();
                return 0;
            }
        }
    }

    // Streamfunction must vanish on boundaries
    for (j = 0; j < Ny+1; j ++)
    {
        if ((psi_r[j][0] != 0) || (psi_r[j][Nz] != 0))
        {
            fprintf(stderr,"Residual streamfunction must vanish at domain boundaries.");
            printUsage();
            return 0;
        }
    }
    for (k = 0; j < Nz+1; k ++)
    {
        if ((psi_r[0][k] != 0) || (psi_r[Ny][k] != 0))
        {
            fprintf(stderr,"Residual streamfunction must vanish at domain boundaries.");
            printUsage();
            return 0;
        }
    }    

#pragma parallel
    
    // If T_relax==0 then phi must be fixed at its relaxation value 
    for (j = 0; j < Ny; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            if (T_relax[j][k] == 0.0) 
            {                 
                phi[j][k] = phi_relax[j][k];      
            }
        }
    }
     
    //////////////////////////////////////
    ///// END READING PARAMETER DATA /////
    //////////////////////////////////////

    

    //////////////////////////////////////////////////////////////////
    ///// BEGIN CALCULATING PARAMETERS ON ALTERNATIVE GRIDPOINTS /////
    //////////////////////////////////////////////////////////////////
    
#pragma parallel
    
    // Parameters on v-gridpoints
    for (j = 0; j < Ny+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            kaph_v[j][k] = 0.5*(kaph_psi[j][k+1] + kaph_psi[j][k]);
            ss_v[j][k] = 0.5*(ss_psi[j][k+1] + ss_psi[j][k]);
            v_r[j][k] = (psi_r[j][k] - psi_r[j][k+1])*_dz;
        }
    }

#pragma parallel
    
    // Parameters on w-gridpoints
    for (j = 0; j < Ny; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            kaph_w[j][k] = 0.5*(kaph_psi[j+1][k] + kaph_psi[j][k]);
            kapv_w[j][k] = 0.5*(kapv_psi[j+1][k] + kapv_psi[j][k]);
            ss_w[j][k] = 0.5*(ss_psi[j+1][k] + ss_psi[j][k]);
            ssq_w[j][k] = ss_w[j][k]*ss_w[j][k];
            w_r[j][k] = (psi_r[j+1][k] - psi_r[j][k])*_dy + wSink[j][k]*pdRatio[j][k]/(1+pdRatio[j][k]);
        }
    }
    
    ////////////////////////////////////////////////////////////////
    ///// END CALCULATING PARAMETERS ON ALTERNATIVE GRIDPOINTS /////
    ////////////////////////////////////////////////////////////////

      

    //////////////////////////////////
    ///// BEGIN CALCULATING CFLS /////
    //////////////////////////////////

#pragma parallel
    
    // For y-fluxes  
    for (j = 0; j <= Ny; j ++)
    {
        for (k = 0; k < Nz; k ++) 
        { 
            // Max advecting velocity
            if (fabs(v_r[j][k]) > v_max)
            {
                v_max = fabs(v_r[j][k]);
            }

            /// Max effective diffusivity
            ydiff = kaph_v[j][k];
            if (ydiff > ydiff_max)
            {
                ydiff_max = ydiff;
            }
        }
    }
 
#pragma parallel
    
    // For z-fluxes
    for (j = 0; j < Ny; j ++)
    {
        for (k = 0; k <= Nz; k ++) 
        { 
            // Max advecting velocity
            if (fabs(w_r[j][k]) > w_max)
            {
                w_max = fabs(w_r[j][k]);
            }

            // Max effective diffusivity
            zdiff = kaph_w[j][k]*ssq_w[j][k];
            if (zdiff > zdiff_max)
            {
                zdiff_max = zdiff;
            }
        }
    }

    // Calculate CFL criteria
    cfl_v = 0.5*dy/v_max;
    cfl_w = 0.5*dz/w_max;
    cfl_y = 0.5*dysq/ydiff_max;
    cfl_z = 0.5*dzsq/zdiff_max;

    // Actual CFL-limted time step
    cfl_min_step = fmin(fmin(cfl_v,cfl_w),fmin(cfl_y,cfl_z));
 
    ////////////////////////////////
    ///// END CALCULATING CFLS /////
    ////////////////////////////////



    // Open the output file
    outfile = fopen(argv[2],"w");
    if (outfile == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file name\n");
        printUsage();
        return 0;
    }

    // Write out the initial data if a save interval is specified
    if (dt_s > 0.0)
    {
        fprintf(outfile,"%e\n",t);
        printMatrix(outfile,phi,Ny,Nz);
    }

    // Write initial time to time file
    if (tfile != NULL)
    {
        fprintf(tfile,"%e ",t);
        fflush(tfile);
    }

    // Numerical time-integration loop - keep going while the residual exceeds 
    // the target and the current time does not exceed the max time
    while ((res > targetRes) && (t < tmax))
    {
        // Perform a single numerical time-step
        switch (method_t)
        {
            case METHOD_RKTVD1:
            {
                dt = rktvd1(&t,phi_in,phi_out,cflFrac,Nphi,&tderiv);
                break;
            }
            case METHOD_RKTVD2:
            {
                dt = rktvd2(&t,phi_in,phi_out,phi_buf,cflFrac,Nphi,&tderiv);
                break;
            }
            case METHOD_RKTVD3:
            {
                dt = rktvd3(&t,phi_in,phi_out,phi_buf,cflFrac,Nphi,&tderiv);
                break;
            }
            default:
            {
                fprintf(stderr,"ERROR: Unknown time-integration method\n");
                break;
            }
        }

#pragma parallel
        
        // Implicit vertical diffusion
        dt_dzsq = dt*_dz*_dz;
        for (j = 0; j < Ny; j ++)
        {            
            for (k = 1; k < Nz-1; k ++)
            {
                impl_A[j][k] = - kapv_w[j][k] * dt_dzsq;
                impl_B[j][k] = 1 + (kapv_w[j][k]+kapv_w[j][k+1]) * dt_dzsq;
                impl_C[j][k] = - kapv_w[j][k+1] * dt_dzsq;
                impl_D[j][k] = phi_outmat[j][k];
            }
           
            impl_A[j][0] = 0;            
            impl_C[j][0] = - kapv_w[j][1] * dt_dzsq;
            impl_B[j][0] = 1 - impl_C[j][0];
            impl_D[j][0] = phi_outmat[j][0];

            impl_A[j][Nz-1] = - kapv_w[j][Nz-1] * dt_dzsq;
            impl_B[j][Nz-1] = 1 - impl_A[j][Nz-1];
            impl_C[j][Nz-1] = 0;                
            impl_D[j][Nz-1] = phi_outmat[j][Nz-1];
                        
            thomas(impl_A[j],impl_B[j],impl_C[j],impl_D[j],phi_outmat[j],Nz);
        }

        // If the time step has taken us past a save point (or multiple 
        // save points), interpolate and write out the data
        while ((dt_s > 0) && (t >= t_next))
        {
            wc = (t-t_next) / dt;
            wn = 1 - wc;

            // Interpolate values at t_next into the phi_buf array
            for (i = 0; i < Nphi; i ++)
            {
                phi_buf[i] = wc*phi_in[i] + wn*phi_out[i];
            }

            // Write out the interpolated data from the phi_in array
            fprintf(outfile,"%le\n",t_next);
            printMatrix(outfile,phi_interp,Ny,Nz); // Recall phi_interp is a pointer to phi_buf

            if (tfile != NULL)
            {
                fprintf(tfile,"%e ",t_next);
                fflush(tfile);
            }

            // Update the next save time
            t_next += dt_s;
            
            // Printing out the residual can be quite a useful way of 
            // keeping track of the computation's progress
            printf("%e\n",res);            
            fflush(stdout);
        }

        // Calculate the residual as an L2-norm between adjacent iterations
        res = 0;
        for (i = 0; i < Nphi; i ++)
        {
            res += SQUARE((phi_in[i]-phi_out[i])/dt);
        }
        res = sqrt(res/Nphi);

        // Copy the next iteration from vars_out back to vars,
        // ready for the next time step
        memcpy(phi_in,phi_out,Nphi*sizeof(real));

        // Increment iteration count
        nIters += 1;
    }

    // NaN residual clearly indicates a problem
    if (isnan(res))
    {        
        fprintf(stderr,"ERROR: Computation blew up: residual==NaN\n");
        printUsage();
        return 0;
    }

    // Write out the final data
    fprintf(outfile,"%e\n",t);
    printMatrix(outfile,phi,Ny,Nz);

    // Print completion time
    if (tfile != NULL)
    {
        time(&now);
        fprintf(tfile,"\nProgram completed at %s\n", ctime(&now));
        fprintf(tfile,"\nres = %e, target res = %e\n", res, targetRes);
        fclose(tfile);
    }



    ///// BEGIN cleanup /////  

    fclose(outfile);

    vecfree(phi_in);
    vecfree(phi_out);
    vecfree(phi_buf);
    free(phi);
    free(phi_interp);

    free(dt_phi_w);
    matfree(phi_w);

    matfree(kapv_psi);
    matfree(kapv_w);
    matfree(kaph_psi);
    matfree(kaph_v);
    matfree(kaph_w);
    matfree(ss_psi);
    matfree(ss_v);
    matfree(ss_w);
    matfree(ssq_w);
    matfree(psi_r);
    matfree(v_r);
    matfree(w_r);
    matfree(phi_relax);
    matfree(T_relax);
    matfree(wSink);
    matfree(pdRatio);

    matfree(HHy);
    matfree(HHz);
    matfree(PPy);
    matfree(PPz);

    matfree(dphi_dy);
    matfree(dphi_dz);
    matfree(phi_ym);
    matfree(phi_yp);    
    matfree(phi_zm);
    matfree(phi_zp);

    matfree(impl_A);
    matfree(impl_B);
    matfree(impl_C);
    matfree(impl_D);
    
    ///// END cleanup /////


    return 0;
}
