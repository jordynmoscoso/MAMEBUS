/**
 * mamebus.c
 *
 * Andrew Stewart
 * Jordyn Moscoso
 *
 * Core code file for the Meridionally-Averaged Model of Eastern Boundary Upwelling Systems (MAMEBUS).
 * Integrates residual-mean buoyancy and tracer equations in an Eastern Boundary Current-like domain.
 *
 */
#include <time.h>

#include "defs.h"
#include "ab.h"



// Total number of input parameters - must match the number of parameters defined in main()
#define NPARAMS 54


//////////////////////////////////
///// BEGIN GLOBAL VARIABLES /////
//////////////////////////////////

// Grid size
int Ntracs = 0;
int Nx = 0;
int Nz = 0;
int Ntot = 0;

// Physical parameters
real Lx = 0;
real Ly = 0;
real Lz = 0;
real Kconv0 = 10;
real rho0 = 1025;
real f0 = 1e-4;
bool use_sml = true;
bool use_bbl = true;
real Hsml = 50.0;
real Hbbl = 50.0;
int tlength = 0;
real r_bbl = 0; // drag coefficient
real Kdia0 = 1e-5;
real Ksml = 1e-1;
real Kbbl = 1e-1;
real nu_h = 0;
real nu_v = 0;

// Biogeochemical parameters
uint bgcModel = 0;
int MP = 0;
int MZ = 0;
int MD = 1;                 // Always have one small detrital group.
int MN = 1;                 // Always have one nutrient
int npbgc = 0;              // Counts number of biogeochemical parameters
int nallo = 0;              // Value for allometric rates
int idxAllo = 0;            // Number of allometric rates
real Nint = 0;              // placeholder for the total domain nitrate

// Scaling Constants
real day = 86400;                 // Seconds in a day
real year = 31449600;             // Seconds in a year (52 7-day weeks)
real grav = 9.8;                  // m/s^2 (gravity)
real tref = 20;                   // Reference surface temperature
real area = 0;
real alpha = 1e-4;                // thermal expansion coefficient

// Parameter arrays
real *** phi_init = NULL;     // Initial condition
real *** phi_north = NULL;    // Northern profiles of tracers
real *** phi_south = NULL;    // Southern profiles of tracers
real ** phi_flux = NULL;
real * hb_in = NULL;          // Ocean depth
real * tau = NULL;            // Surface wind stress
real ** Kgm_psi_ref = NULL;   // Reference GM diffusivity
real ** Kiso_psi_ref = NULL;  // Reference isopycnal diffusivity
real ** Kdia_psi_ref = NULL;  // Reference diapycnal diffusivity
real *** phi_relax = NULL;    // Tracer relaxation values
real *** T_relax = NULL;      // Tracer relaxation time scale

// BGC parameter arrays
real ** bgcRates = NULL;
real ** theta_p = NULL;       // predator-prey interaction for Z on P grazing
real ** theta_z = NULL;       // predator-prey interaction for Z on Z grazing
real * bgc_params = NULL;     // Vector containing biogeochemical parameters from bgc_setup
real * lpvec = NULL;
real * lzvec = NULL;

// Numerical parameters
real KT00_sigma = 1;         // Kurganov-Tadmor minmod-limiting parameter
bool limSlopes = true;    // Use Cox slope limiting
real Smax = 0.1;        // Max isopycnal slope
const int idx_uvel = 0;   // Index of zonal momentum variable in list of tracers
const int idx_vvel = 1;   // Index of meridional momentum variable in list of tracers
const int idx_buoy = 2;   // Index of buoyancy variable in list of tracers
const int idx_nitrate = 3;
const int max_det = 3;

// Sigma-coordinate grid parameters
real h_c = 1e16;
real theta_s = 0;
real theta_b = 0;

// Grids and grid spacings
real ds = 0;
real dx = 0;
real _dx = 0;
real _2dx = 0;
real dxsq = 0;
real minVal = 1e-12;      // Minimum values for harmonic averages
real cff = 0;             // dummy variable to hold numerators and coefficients
real * sigma_phi = NULL;  // Stretched coordinate grids
real * sigma_psi = NULL;
real ** _dz_phi = NULL;   // Grid spacings
real ** _2dz_phi = NULL;
real ** _dzsq_psi = NULL;
real ** _dz_w = NULL;
real ** _dzsq_w = NULL;
real ** dz_u = NULL;
real ** _dz_u = NULL;
real ** dA_psi = NULL;
real * hb_phi = NULL;     // Ocean depths and bottom slope
real * hb_psi = NULL;
real * sb_phi = NULL;
real * sb_psi = NULL;
real ** ZZ_phi = NULL;    // Depth grids
real ** ZZ_psi = NULL;
real ** ZZ_u = NULL;
real ** ZZ_w = NULL;


// Debugging parameter
bool debug = false;
bool nanfound = false;

// Name of the program (for error messages)
char * progname = NULL;

// Time-stepping scheme
uint timeSteppingScheme = TIMESTEPPING_AB3;

// Tracer advection scheme
uint advectionScheme = ADVECTION_KT00;

// Momentum scheme
uint momentumScheme = MOMENTUM_TTTW;

// Pressure scheme
uint pressureScheme = PRESSURE_CUBIC;

// Output filenames
static const char OUTN_ZZ_PHI[] = "ZZ_PHI";
static const char OUTN_ZZ_PSI[] = "ZZ_PSI";
static const char OUTN_ZZ_U[] = "ZZ_U";
static const char OUTN_ZZ_W[] = "ZZ_W";
static const char OUTN_PSIM[] = "PSIM";
static const char OUTN_PSIE[] = "PSIE";
static const char OUTN_PSIR[] = "PSIR";
static const char OUTN_TRAC[] = "TRAC";

// Debugging
static const char OUTN_DBDX_CUBIC[] = "DB_DX_CUBIC";
static const char OUTN_DBDX_LINEAR[] = "DB_DX_LINEAR";

// Work arrays for calculating time derivatives
real * phi_wrk_V = NULL;      // Current iteration data - work array
real *** phi_wrk = NULL;
real *** dphi_dt_wrk = NULL;  // Time derivative of current iteration - work array
real ** HHx = NULL;           // Hyperbolic tracer fluxes
real ** HHz = NULL;
real ** PPx = NULL;           // Parabolic tracer fluxes
real ** PPz = NULL;
real ** psi_r = NULL;         // Residual streamfunction
real ** u_r = NULL;
real ** w_r = NULL;
real ** psi_m = NULL;         // Eulerian-mean streamfunction
real ** psi_e = NULL;         // Eddy streamfunction
real ** dphi_dx = NULL;       // Tracer gradients at cell centers
real ** dphi_dz = NULL;
real ** phi_xp = NULL;        // Extrapolated tracer values on cell faces
real ** phi_xm = NULL;
real ** phi_zp = NULL;
real ** phi_zm = NULL;
real ** Sgm_psi = NULL;       // Isopycnal slope (for residual streamfunction)
real ** Siso_u = NULL;
real ** Siso_w = NULL;
real * impl_A = NULL;         // Implicit diffusion solver arrays
real * impl_B = NULL;
real * impl_C = NULL;
real * impl_D = NULL;
real ** wsink_wrk = NULL;
real ** Kgm_psi = NULL;       // GM diffusivity
real ** Kgm_u = NULL;
real ** Kgm_w = NULL;
real ** Kiso_u = NULL;        // Isopycnal diffusivity
real ** Kiso_w = NULL;
real ** Kdia_w = NULL;        // Diapycnal diffusivity
real ** BPa = NULL;           // Baroclinic Pressure
real ** BPx = NULL;           // Zonal barolinic Pressure gradient
real ** BPy = NULL;           // Meridional baroclinic pressure gradient
real ** BBy = NULL;           // Buoyancy
real ** Nbuoy = NULL;
real ** rho_north = NULL;     // density at the northern edge of the domain
real ** rho_south = NULL;     // density at the southen edge of the domain
real _Ly = 0;

// Pointers for pressure calculation scheme.
real ** drx = NULL;     // Stores elementary differences
real ** drz = NULL;
real ** dzx = NULL;
real ** dzz = NULL;
real ** hrx = NULL;     // Stores hyperbolic differences
real ** hrz = NULL;
real ** hzx = NULL;
real ** hzz = NULL;
real ** rhos = NULL;    // Density holder
real ** P = NULL;       // Pressure
real ** FC = NULL;      // integrated density
real ** FX = NULL;


// Boundary layer work arrays
uint * k_sml = NULL;
real * wn_sml = NULL;
real * wp_sml = NULL;
uint * k_bbl = NULL;
real * wn_bbl = NULL;
real * wp_bbl = NULL;
real ** db_dx = NULL;
real ** db_dz = NULL;
real ** db_dx_wrk = NULL;
real ** db_dz_wrk = NULL;
real ** db_dx_lin = NULL;
real ** bot_nflux = NULL;


// work arrays for phytoplankton and zooplankton time tendencies
real ** hetMat = NULL;                      // heterotrophic grazing matrix between prey zooplankton and predator zooplankton
real ** grazMat = NULL;                      // grazing matrix between phytoplankton and zooplankton
real * GPvec = NULL;                     // loss of biomass from phytoplankton due to grazing
real * GZvec = NULL;                     // growth of biomass to zooplankton from grazing
real * PTvec = NULL;                     // total phytoplankton avaialble to zooplankton vector
real * ZTvec = NULL;                     // total zooplankton avaiable to predator zooplankton
real * UPvec = NULL;                     // uptake vector
real * HMvec = NULL;                     // loss due to heterotrophic grazing
real * HPvec = NULL;                     // growth due to heterotrophic grazing
real * MPvec = NULL;                     // mortality of phytoplankton
real * MZvec = NULL;                     // mortality of zooplankton

////////////////////////////////
///// END GLOBAL VARIABLES /////
////////////////////////////////










/**
 * calcPsim
 *
 * Calculates the mean overturning streamfunction from the mean momentum equation.
 * The result is stored in the psi_m matrix.
 *
 */
void calcPsim (const real t, real ** uvel, real ** psi_m)
{
    int j,k;
    real z;
    
    // Zero streamfunction at top/bottom boundaries
    for (j = 0; j < Nx+1; j ++)
    {
        psi_m[j][0] = 0;
        psi_m[j][Nz] = 0;
    }
    // Zero streamfunction at east/west boundaries
    for (k = 0; k < Nz+1; k ++)
    {
        psi_m[0][k] = 0;
        psi_m[Nx][k] = 0;
    }
    
    // If we are not evolving momentum explicitly then we prescribe simple
    // uniform Ekman velocities in the SML and BBL
    if (momentumScheme == MOMENTUM_NONE)
    {
        
        
#pragma parallel
        
        // Current implementation: uniform horizontal velocity in SML and BBL
        for (j = 1; j < Nx; j ++)
        {
            for (k = 1; k < Nz; k ++)
            {
                z = ZZ_psi[j][k];
                
                psi_m[j][k] = tau[j]/(rho0*f0);
                
                if (use_sml && (z > -Hsml))
                {
                    psi_m[j][k] *= (-z/Hsml);
                }
                if (use_bbl && (z < -hb_psi[j] + Hbbl))
                {
                    psi_m[j][k] *= ((z+hb_psi[j])/Hbbl);
                }
            }
        }
    }
    // Otherwise just integrate u-velocity to get the mean streamfunction
    else
    {
        for (j = 1; j < Nx; j ++)
        {
            for (k = 1; k < Nz; k ++)
            {
                psi_m[j][k] = psi_m[j][k-1] - dz_u[j][k-1]*uvel[j][k-1];
            }
        }
    }
    
}










/**
 *
 * calcKiso
 *
 * Calculates the Redi (1981) isopycnal diffusivity (Kiso).
 *
 */
void calcKiso (const real t, real ** buoy, real ** Kiso_u, real ** Kiso_w)
{
    int j,k;
    
    // Current implementation just keeps Kiso constant
    
#pragma parallel
    
    // Diffusivity on u-gridpoints
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            Kiso_u[j][k] = 0.5*(Kiso_psi_ref[j][k+1] + Kiso_psi_ref[j][k]);
        }
    }
    
#pragma parallel
    
    // Diffusivity on w-gridpoints
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            Kiso_w[j][k] = 0.5*(Kiso_psi_ref[j+1][k] + Kiso_psi_ref[j][k]);
        }
    }
    
}












/**
 *
 * calcKdia
 *
 * Calculates the diapycnal diffusivity (Kdia).
 *
 */
void calcKdia (const real t, real ** buoy, real ** Kdia_w)
{
    int j,k;
    
    // Current implementation takes constant Kdia and adds
    // a large diffusivity where the water column is statically
    // unstable
    
#pragma parallel
    
    // Parameters on w-gridpoints
    for (j = 0; j < Nx; j ++)
    {
        for (k = 1; k < Nz; k ++)
        {
            Kdia_w[j][k] = 0.5*(Kdia_psi_ref[j+1][k] + Kdia_psi_ref[j][k]);
            if (buoy[j][k] < buoy[j][k-1])
            {
                Kdia_w[j][k] += Kconv0;
            }
        }
        
        Kdia_w[j][0] = 0.5*(Kdia_psi_ref[j+1][0] + Kdia_psi_ref[j][0]);
        Kdia_w[j][Nz] = 0.5*(Kdia_psi_ref[j+1][Nz] + Kdia_psi_ref[j][Nz]);
    }
    
}


















/**
 * surfStructFun
 *
 * Surface structure function for the modified Ferrari et al. (2008) boundary layer parameterization. Aligns the slope used
 * to calculate the eddy streamfunction and symmetric diffusion tensor with the ocean surface as the surface is approached.
 *
 * z is the current vertical position
 * h_sml is the surface mixed layer thickness
 * _lambda_sml is the reciprocal of the vertical eddy lengthscale at the SML base, and defines the vertical derivative of G at the SML base
 *
 *
 */
real surfStructFun (real z, real h_sml, real _lambda)
{
    real G = 1.0;
    
    if (z > -h_sml)
    {
//        G = -z/h_sml;
        G = -(1+h_sml*_lambda)*SQUARE(z/h_sml) - (2+h_sml*_lambda)*(z/h_sml);
        
    }
    
    return G;
    
}


/**
 * botStructFun
 *
 * Bottom structure function for the modified Ferrari et al. (2008) boundary layer parameterization. Reduces the effective
 * slope used to calculate the eddy streamfunction to zero at the ocean bed.
 *
 */
real botStructFun (real z, real h_bbl, real h_b, real lambda)
{
    real G = 1.0;

    if (z < h_b + h_bbl)
    {
        G = ((z+h_b)/h_bbl)*(2 - ((z+h_b)/h_bbl));
    }
    
    return G;
}


/**
 * slopeStructFun
 *
 * Slope structure function for the modified Ferrari et al. (2008) boundary layer parameterization. Aligns the effective
 * slope for the symmetric diffusion tensor with the ocean bed as the ocean bed is approached.
 *
 */
real slopeStructFun (real z, real h_bbl, real h_b)
{
    real G = 1.0;
    
    if (z < h_b + h_bbl)
    {
        G = (z+h_b)/h_bbl;
    }
    
    return G;
}

















/*
 * Calculate the pressure gradient using a linear pressure interpolation
 * this calculation is the linear version of the  "Sigma coordinate
 * primitive form" in the Shchepetkin and McWilliams (2003) paper.
 *
 */

void calcPressure(const real t, real ** buoy)
{
    int j,k;
    real z;
    
    // Working variables
    real zeta = 0;
    
    real OneFifth = 1.0/5.0;
    real OneTwelfth = 1.0/12.0;

    real cff1 = 0;
    real cff2 = 0;
    
    
    // Calculate the density field
    for (j = 0; j < Nx; j++)
    {
        for (k = 0; k < Nz; k++)
        {
            rhos[j][k] = (1-alpha*(buoy[j][k] - tref));
        }
    }
    

    
    // Vertical Calculation
    // Compute elementary differences
    for (j = 0; j < Nx; j++)
    {
        for (k = 0; k < Nz-1; k++)
        {
            // starts with +1/2 (0) and ends with Nz-3/2 (Nz-2)
            drz[j][k] = rhos[j][k+1] - rhos[j][k];
        }
    }

    
    // Calculate the harmonic averages:
    for (j = 0; j < Nx; j++)
    {
        for (k = 1; k < Nz-1; k++)
        {
            cff = drz[j][k-1]*drz[j][k];
            if (cff > minVal)
            {
                hrz[j][k] = 2*cff/(drz[j][k-1] + drz[j][k]);
            }
            else
            {
                hrz[j][k] = 0;
            }
        }

        // Set the boundary conditions
        hrz[j][0] = 1.5*(rhos[j][1] - rhos[j][0]) - 0.5*hrz[j][1];
        hrz[j][Nz-1] = 1.5*(rhos[j][Nz-1] - rhos[j][Nz-2]) - 0.5*hrz[j][Nz-2];
    }


    // Calculate the pressure at the surface
    for (j = 0; j < Nx; j++)
    {
        cff = ZZ_w[j][Nz] - ZZ_phi[j][Nz-1];
        P[j][Nz-1] = grav*cff*( rhos[j][Nz-1]
                         + 0.5 * cff * ( rhos[j][Nz-1] - rhos[j][Nz-2] )/( ZZ_phi[j][Nz-1] - ZZ_phi[j][Nz-2] ) );
    }


    // Calculate the pressure by vertically integrating
    for (j = 0; j < Nx; j++)
    {
        for (k = Nz-2; k >= 0; k --)
        {
            P[j][k] = P[j][k+1] + grav * 0.5 * ( (rhos[j][k+1] + rhos[j][k])*( ZZ_phi[j][k+1] - ZZ_phi[j][k] )
                        - OneFifth * ( ( hrz[j][k+1] - hrz[j][k] )*( ZZ_phi[j][k+1] - ZZ_phi[j][k]
                            - OneTwelfth*( hzz[j][k+1] + hzz[j][k] ) )
                        - ( hzz[j][k+1] - hzz[j][k] )*( rhos[j][k+1] - rhos[j][k]
                            - OneTwelfth* ( hrz[j][k+1] + hrz[j][k] ) ) ) );
        }
    }


    // Horizontal Calculations
    // Elementary Differences
    for (j = 0; j < Nx-1; j++)
    {
        for (k = 0; k < Nz; k++)
        {
            drx[j][k] = rhos[j+1][k] - rhos[j][k];
        }
    }

    // Calculate the harmonic average:
    for (k = 0; k < Nz; k++)
    {
        for (j = 1; j < Nx-1; j++)
        {
            cff = drx[j-1][k]*drx[j][k];
            if (cff > minVal)
            {
                hrx[j][k] = 2*cff/(drx[j-1][k] + drx[j][k]);
            }
            else
            {
                hrx[j][k] = 0;
            }
        }

        hrx[0][k] = 1.5*(rhos[1][k] - rhos[0][k]) - 0.5*hrx[1][k];
        hrx[Nx-1][k] = 1.5*(rhos[Nx-1][k] - rhos[Nx-2][k]) - 0.5*hrx[Nx-2][k];
    }


    // Calculate the gradient along the horizontal component of the line integral
    for (k = 0; k < Nz; k++)
    {
        for (j = 1; j < Nx; j++)
        {
            FC[j][k] = 0.5 * ( (rhos[j][k] + rhos[j-1][k])*(ZZ_phi[j][k] - ZZ_phi[j-1][k])
                              - OneFifth * ( ( hrx[j][k] - hrx[j-1][k] )*( ZZ_phi[j][k] - ZZ_phi[j-1][k]
                                        - OneTwelfth * ( hzx[j][k] + hzx[j-1][k] ) )
                                    - (hzx[j][k] - hzx[j-1][k])*( rhos[j][k] - rhos[j-1][k]
                                        - OneTwelfth * ( hrx[j][k] - hrx[j-1][k] ) ) ) );
        }
    }


    // Calculate the zonal pressure gradient
    // BPx is stored on the u grid points, meaning that FC starts on the 1/2 (between phi_0 and phi_1) points
    // And the left boundary BPx[0][k] occurs on the left of phi_0 so on the -1/2 grid points.
    for (j = 1; j < Nx; j++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Set the boundary to zero
            BPx[0][k] = 0;
            BPx[Nx][k] = 0;

            BPx[j][k] = (grav * FC[j][k] - (P[j-1][k] - P[j][k]))*_dx;
        }
    }
    

}
















/**
 *
 * calcSlopes
 *
 * Calculates isopycnal slopes everywhere. The result is stored in
 * Sgm_psi (psi-gridpoints), Siso_u (u-gridpoints) and Siso_w (w-gridpoints).
 *
 * In the top and bottom boundary layers, the slopes are replaced
 * by effective slopes
 *   Sgm_psi -> S_e
 *   Siso_u,Siso_w -> S'_e
 * where S_e approaches zero at the boundary (so that the GM streamfunction
 * vanishes there) and S'_e approaches the top or bottom slope at the boundary,
 * ensuring that the eddy residual tracer flux is oriented parallel to the
 * boundary.
 *
 */
void calcSlopes (     const real        t,
                 real **           buoy,
                 real **           Sgm_psi,
                 real **           Siso_u,
                 real **           Siso_w)
{
    int j,k;
    real db_dz_sml, db_dz_bbl;
    real G_sml, G_bbl, G_slope;
    real z;
    real d2b_dz2;
    real _lambda_sml, lambda_bbl;
    real dzm = 0;
    real dzp = 0;
    real dpm = 0;
    real dpp = 0;
    
    
    // Working variables
    real zeta = 0;
    
    real OneFifth = 1.0/5.0;
    real OneTwelfth = 1.0/12.0;
    real cff1 = 0;
    real cff2 = 0;
    
#pragma parallel

    switch(pressureScheme)
    {
        // If the pressure scheme is the linear pressure gradient, then calculate the buoyancy gradient here
        case PRESSURE_LINEAR:
        {
            // Calculate horizontal and vertical buoyancy gradients in (x,z) space
            for (j = 1; j < Nx; j ++)
            {
                for (k = 1; k < Nz; k ++)
                {
                    db_dz[j][k] = 0.5 * ( (buoy[j][k]-buoy[j][k-1])*_dz_w[j][k] + (buoy[j-1][k]-buoy[j-1][k-1])*_dz_w[j-1][k] );
                    db_dx[j][k] = 0.5 * ( (buoy[j][k]-buoy[j-1][k])*_dx + (buoy[j][k-1]-buoy[j-1][k-1])*_dx );
                    db_dx[j][k] -= (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx * db_dz[j][k];
//                    fprintf(stderr,"j = %d, k = %d, db_dx = %le, db_dz = %le \n",j,k,db_dx[j][k],db_dz[j][k]);
                }
            }
            
            break;
        }
        case PRESSURE_CUBIC: // Calculate the buoyancy gradient using the Shchepetkin & McWilliams (2003) scheme.
        {
//Linear
            for (j = 1; j < Nx; j ++)
            {
                for (k = 1; k < Nz; k ++)
                {
                    db_dz[j][k] = 0.5 * ( (buoy[j][k]-buoy[j][k-1])*_dz_w[j][k] + (buoy[j-1][k]-buoy[j-1][k-1])*_dz_w[j-1][k] );
                    db_dx_lin[j][k] = 0.5 * ( (buoy[j][k]-buoy[j-1][k])*_dx + (buoy[j][k-1]-buoy[j-1][k-1])*_dx );
                    db_dx_lin[j][k] -= (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx * db_dz[j][k];
//                    fprintf(stderr,"j = %d, k = %d, db_dx = %le, db_dz = %le \n",j,k,db_dx[j][k],db_dz[j][k]);
                }
            }
            
                        
            // Vertical Calculation
            // Compute elementary differences
            for (j = 0; j < Nx; j++)
            {
                for (k = 0; k < Nz-1; k++)
                {
                    // starts with +1/2 (0) and ends with Nz-3/2 (Nz-2)
                    drz[j][k] = buoy[j][k+1] - buoy[j][k];
                }
            }


            // Calculate the harmonic averages:
            for (j = 0; j < Nx; j++)
            {
                for (k = 1; k < Nz-1; k++)
                {
                    cff = drz[j][k-1]*drz[j][k];
                    if (cff > minVal)
                    {
                        hrz[j][k] = 2*cff/(drz[j][k-1] + drz[j][k]);
                    }
                    else
                    {
                        hrz[j][k] = 0;
                    }
                }

                // Set the boundary conditions
                hrz[j][0] = 1.5*(buoy[j][1] - buoy[j][0]) - 0.5*hrz[j][1];
                hrz[j][Nz-1] = 1.5*(buoy[j][Nz-1] - buoy[j][Nz-2]) - 0.5*hrz[j][Nz-2];
            }


            // Calculate the component of the integral at the surface
            for (j = 0; j < Nx; j++)
            {
                cff = ZZ_w[j][Nz] - ZZ_phi[j][Nz-1];
                FX[j][Nz-1] = cff*(buoy[j][Nz-1]
                                 + 0.5 * cff * ( buoy[j][Nz-1] - buoy[j][Nz-2] )/( ZZ_phi[j][Nz-1] - ZZ_phi[j][Nz-2] ) );
            }


            // Calculate the vertical component of the gradient
            // At the moment the correction term is causing problems, so we set it to zero.
            cff2 = OneFifth;
            for (j = 0; j < Nx; j++)
            {
                for (k = Nz-2; k >= 1; k --)
                {
                    FX[j][k] = 0.5 * ( (buoy[j][k] + buoy[j][k-1])*( ZZ_phi[j][k] - ZZ_phi[j][k-1] )
                                - cff2 * ( ( hrz[j][k] - hrz[j][k-1] )*( ZZ_phi[j][k] - ZZ_phi[j][k-1]
                                    - OneTwelfth*( hzz[j][k] + hzz[j][k-1] ) )
                                - ( hzz[j][k] - hzz[j][k-1] )*( buoy[j][k] - buoy[j][k-1]
                                    - OneTwelfth* ( hrz[j][k] + hrz[j][k-1] ) ) ) );
                }
            }


            // Horizontal Calculations
            // Elementary Differences
            for (j = 0; j < Nx-1; j++)
            {
                for (k = 0; k < Nz; k++)
                {
                    drx[j][k] = buoy[j+1][k] - buoy[j][k];
                }
            }

            // Calculate the harmonic average:
            for (k = 0; k < Nz; k++)
            {
                for (j = 1; j < Nx-1; j++)
                {
                    cff = drx[j-1][k]*drx[j][k];
                    if (cff > minVal)
                    {
                        hrx[j][k] = 2*cff/(drx[j-1][k] + drx[j][k]);
                    }
                    else
                    {
                        hrx[j][k] = 0;
                    }
                }

                hrx[0][k] = 1.5*(rhos[1][k] - rhos[0][k]) - 0.5*hrx[1][k];
                hrx[Nx-1][k] = 1.5*(rhos[Nx-1][k] - rhos[Nx-2][k]) - 0.5*hrx[Nx-2][k];

                // Test extending the boundary only
//                hrx[0][k] = hrx[1][k];
//                hrx[Nx-1][k] = hrx[Nx-2][k];
            }


            // Calculate the gradient along the horizontal component of the line integral
            for (k = 0; k < Nz; k++)
            {
                for (j = 1; j < Nx; j++)
                {
                    FC[j][k] = 0.5 * ( (buoy[j][k] + buoy[j-1][k])*(ZZ_phi[j][k] - ZZ_phi[j-1][k])
                                      - OneFifth * ( ( hrx[j][k] - hrx[j-1][k] )*( ZZ_phi[j][k] - ZZ_phi[j-1][k]
                                                - OneTwelfth * ( hzx[j][k] + hzx[j-1][k] ) )
                                            - (hzx[j][k] - hzx[j-1][k])*( buoy[j][k] - buoy[j-1][k]
                                                - OneTwelfth * ( hrx[j][k] - hrx[j-1][k] ) ) ) );
                }
            }



            // Calculate the bouyancy gradient
            for (j = 1; j < Nx; j++)
            {
                for (k = 1; k < Nz; k ++)
                {
                    cff = 1.0/dA_psi[j][k];
                    if (k == Nz-1)
                    {
                        // Use a linear calculation for the surface.
                        db_dx[j][k] = 0.5 * ( (buoy[j][k]-buoy[j-1][k])*_dx + (buoy[j][k-1]-buoy[j-1][k-1])*_dx );
                        db_dx[j][k] -= (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx * db_dz[j][k];
                    }
                    else
                    {
                        // Calculate the buoyancy gradient using the cubic spline interpolation
                        db_dx[j][k] = -(FC[j][k] - FC[j][k-1] + FX[j-1][k] - FX[j][k])*cff;
                    }
                }
            }



            // Set the boundary to zero
            // these should never be used, but defined just in case
            for (k = 0; k < Nz; k++)
            {
                db_dx[0][k] = 0;
                db_dx[Nx][k] = 0;
            }
            

            for (j = 1; j < Nx; j++)
            {
                for (k = 1; k < Nz; k++)
                {
                    cff2 = 0.5 * ( (buoy[j][k]-buoy[j-1][k])*_dx + (buoy[j][k-1]-buoy[j-1][k-1])*_dx );
                    cff2 -= (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx * db_dz[j][k];
                    db_dx_wrk[j][k] = cff2;
                }
            }
            
            
            break;
        }
        default:
        {
            printf("Error: PressureScheme is undefined \n");
            break;
        }
            
    }

    
    // Calculate the effective isopycnal slope S_e everywhere
    for (j = 1; j < Nx; j ++)
    {
        db_dz[j][0] = 0; // N.B. THESE SHOULD NEVER BE USED
        db_dx[j][0] = 0; // Here we just set then so that they are defined
        db_dz[j][Nz] = 0;
        db_dx[j][Nz] = 0;
        
        // Calculate vertical gradients at SML base and BBL top
        if (use_sml)
        {
            // Buoyancy gradient at SML base
            db_dz_sml = db_dz[j][k_sml[j]]*wp_sml[j] + db_dz[j][k_sml[j]+1]*wn_sml[j];
            
            // Calculate reciprocal of vertical eddy length scale at SML base
            d2b_dz2 = (db_dz[j][k_sml[j]+1]-db_dz[j][k_sml[j]]) / (ZZ_psi[j][k_sml[j]+1]-ZZ_psi[j][k_sml[j]]);
            _lambda_sml = - d2b_dz2 / db_dz_sml;
        }
        if (use_bbl)
        {
            // Buoyancy gradient at BBL top
            db_dz_bbl = db_dz[j][k_bbl[j]]*wn_bbl[j] + db_dz[j][k_bbl[j]-1]*wp_bbl[j];
            lambda_bbl = 0;
        }
        
        // Construct effective isopycnal slope
        for (k = 1; k < Nz; k ++)
        {
            z = ZZ_psi[j][k];
            
            if (use_sml && (z > -Hsml))
            {
                G_sml = surfStructFun(z,Hsml,_lambda_sml);
                Sgm_psi[j][k] = - G_sml * db_dx[j][k] / db_dz_sml;
            }
            else if (use_bbl && (z < -hb_psi[j] + Hbbl))
            {
                G_bbl = botStructFun(z,Hbbl,hb_psi[j],lambda_bbl);
                Sgm_psi[j][k] = - G_bbl * db_dx[j][k] / db_dz_bbl;
            }
            else
            {
                // Calculate slope on psi-gridpoints
                Sgm_psi[j][k] = - db_dx[j][k] / db_dz[j][k];
            }
                // Cox slope-limiting.
                if (limSlopes)
                {
                    if (Sgm_psi[j][k] > Smax)
                    {
                        Sgm_psi[j][k] = Smax;
                    }
                    if (Sgm_psi[j][k] < -Smax)
                    {
                        Sgm_psi[j][k] = -Smax;
                    }
                }

        }
    }
    
    // Isopycnal slope on the boundaries
    // N.B. THESE SHOULD NEVER BE USED
    // except to set the slope on u and w points below
    for (k = 1; k < Nz; k ++)
    {
        Sgm_psi[0][k] = Sgm_psi[1][k];
        Sgm_psi[Nx][k] = Sgm_psi[Nx-1][k];
    }
    for (j = 1; j < Nx; j ++)
    {
        Sgm_psi[j][0] = Sgm_psi[j][1];
        Sgm_psi[j][Nz] = Sgm_psi[j][Nz-1];
    }
    Sgm_psi[0][0] = Sgm_psi[1][1];
    Sgm_psi[0][Nz] = Sgm_psi[1][Nz-1];
    Sgm_psi[Nx][0] = Sgm_psi[Nx-1][1];
    Sgm_psi[Nx][Nz] = Sgm_psi[Nx-1][Nz-1];
    
#pragma parallel
    
    // Calculate modified effective slope S'_e on u-gridpoints
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Interpolate S_e
            Siso_u[j][k] = 0.5*(Sgm_psi[j][k+1] + Sgm_psi[j][k]);
            
            // Augment to calculate modified effective slope S'_e
            if (use_bbl)
            {
                z = ZZ_u[j][k];
                if (z < -hb_psi[j] + Hbbl)
                {
                    G_slope = slopeStructFun(z,Hbbl,hb_psi[j]);
                    Siso_u[j][k] += (1-G_slope)*sb_psi[j];
                }
            }
        }
    }
    
#pragma parallel
    
    // Calculate modified effective slope S'_e on w-gridpoints
    for (k = 0; k < Nz+1; k ++)
    {
        for (j = 0; j < Nx; j ++)
        {
            // Interpolate S_e
            Siso_w[j][k] = 0.5*(Sgm_psi[j+1][k] + Sgm_psi[j][k]);
            
            // Augment to calculate modified effective slope S'_e
            if (use_bbl)
            {
                z = ZZ_w[j][k];
                if (z < -hb_phi[j] + Hbbl)
                {
                    G_slope = slopeStructFun(z,Hbbl,hb_phi[j]);
                    Siso_w[j][k] += (1-G_slope)*sb_phi[j];
                }
            }
        }
    }
    
    
    
}















/**
 *
 * calcKgm
 *
 * Calculates the Gent-McWilliams (1990) buoyancy diffusivity (Kgm).
 *
 */
void calcKgm (const real t, real ** buoy, real ** Kgm_psi, real ** Kgm_u, real ** Kgm_w)
{
    int j,k;
    real ucon = 0.015;
    real ld = 0;
    real n_col = 0;
    real cff = 0;
    real db_dz_zint = 0, db_dx_zint = 0, srho_zavg = 0, delta = 0, Kfac = 0; // Eddy parameterization parameters
    real alpha_K = 0.05; // Kappa tuning parameter
    real Kmin = 100; // Minimum GM diffusivity
#pragma parallel
    
    
    
    // Loop over model grid and compute Kgm
    for (j = 0; j < Nx+1; j ++)
    {
        /*
        // Integrate buoyancy gradients and compute slope parameter
        db_dz_zint = 0;
        db_dx_zint = 0;
        for (k = 1; k < Nz; k ++)
        {
            db_dx_zint = db_dx_zint + db_dx[j][k]*(ZZ_u[j][k]-ZZ_u[j][k-1]);
            db_dz_zint = db_dz_zint + db_dz[j][k]*(ZZ_u[j][k]-ZZ_u[j][k-1]);
        }
        srho_zavg = - db_dx_zint / db_dz_zint;
        delta = sb_psi[j]/srho_zavg;
        
        //alpha_K = 0.05;
        //Kfac = (1 + 0.5*sqrt( SQUARE(1-fabs(delta)) + 4*SQUARE(alpha)*SQUARE(fabs(delta)) )
        //          - 0.5*sqrt( SQUARE(1+fabs(delta)) + 4*SQUARE(alpha)*SQUARE(fabs(delta)) ) );
        Kfac = (2.5/1000) * ( fabs(delta) + 1/(0.05*(fabs(delta)+0.05)) );
        */
        
        for (k = 0; k < Nz+1; k ++)
        {
            //Kgm_psi[j][k] = Kgm_psi_ref[j][k]*Kfac;
            Kgm_psi[j][k] = Kgm_psi_ref[j][k]; // Constant Kgm
            //Kgm_psi[j][k] = fmax(Kgm_psi[j][k],Kmin);
        }
    }
    
#pragma parallel
    
    // Diffusivity on u-gridpoints
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            Kgm_u[j][k] = 0.5*(Kgm_psi[j][k+1] + Kgm_psi[j][k]);
        }
    }
    
#pragma parallel
    
    // Diffusivity on w-gridpoints
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            Kgm_w[j][k] = 0.5*(Kgm_psi[j+1][k] + Kgm_psi[j][k]);
        }
    }
    
}



























/**
 * calcPsie
 *
 * Calculates the eddy streamfunction from the mean buoyancy field.
 *
 */
void calcPsie (const real t, real ** buoy, real ** Kgm_psi, real ** Sgm_psi, real ** psi_e)
{
    int j,k;
    
    // Zero streamfunction at top/bottom boundaries
    for (j = 0; j < Nx+1; j ++)
    {
        psi_e[j][0] = 0;
        psi_e[j][Nz] = 0;
    }
    // Zero streamfunction at east/west boundaries
    for (k = 0; k < Nz+1; k ++)
    {
        psi_e[0][k] = 0;
        psi_e[Nx][k] = 0;
    }
    
#pragma parallel
    
    // Calculate psi_e everywhere
    for (j = 1; j < Nx; j ++)
    {
        for (k = 1; k < Nz; k ++)
        {
            psi_e[j][k] = Kgm_psi[j][k]*Sgm_psi[j][k];
        }
    }
    
}




















/**
 * calcPsir
 *
 * Calculates the residual streamfunction from the mean momentum equation
 * and mean buoyancy field.
 *
 */
void calcPsir (real ** psi_m, real ** psi_e, real ** psi_r, real ** u_r, real ** w_r)
{
    int j,k;
    
#pragma parallel
    
    // Calculate residual streamfunction - this step is really just for code clarity
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            psi_r[j][k] = psi_e[j][k] + psi_m[j][k];
        }
    }
    
#pragma parallel
    
    for (k = 0; k < Nz; k ++)
    {
        for (j = 0; j < Nx+1; j ++)
        {
            // x-velocity
            u_r[j][k] = - (psi_r[j][k+1]-psi_r[j][k]) * _dz_u[j][k];
        }
    }
    
#pragma parallel
    
    // Calculate z-fluxes
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            // z-velocity
            w_r[j][k] = (psi_r[j+1][k]-psi_r[j][k]) * _dx;
        }
    }
}
















/**
 * tderiv_relax
 *
 * Calculates the time tendency of each tracer due to tracer relaxation.
 *
 */
void tderiv_relax (const real t, real *** phi, real *** dphi_dt)
{
    // Looping variables
    int i,j,k;
    
    for (i = 0; i < Ntracs; i ++)
    {
#pragma parallel
        
        // Calculate time derivative of phi
        for (j = 0; j < Nx; j ++)
        {
            for (k = 0; k < Nz; k ++)
            {
                // Negative relaxation time means no relaxation
                if (T_relax[i][j][k] >= 0)
                {
                    // Instantaneous relaxation: phi cannot change
                    // N.B. This is enforced again after implicit diffusion is applied
                    if (T_relax[i][j][k] == 0.0)
                    {
                        dphi_dt[i][j][k] = 0.0;
                    }
                    // Otherwise relax phi towards phi_relax with timescale T_relax
                    else
                    {
                        dphi_dt[i][j][k] -= (phi[i][j][k]-phi_relax[i][j][k]) / T_relax[i][j][k];
                    }
                }
            }
        }
        
        
        // Calculate the surface fluxes
        for (j = 0; j < Nx; j++)
        {
            dphi_dt[i][j][Nz-1] += _dz_phi[j][Nz-1]*phi_flux[i][j];
        }
        
    }
    
}








void npzd(const real t, const int j, const int k, real *** phi, real *** dphi_dt)

{
    int jj = 0;
    int kk = 0;
    int ip = 0;
    int iz = 0;
    int jz = 0;
    
    // TO DO
    // SIZE DIFFUSION
    
    
    // load in the biogeochemical variables
    real qsw = bgc_params[0];          // maximum surface irradiance
    real kw = bgc_params[1];           // light attunation of water
    real kc = bgc_params[2];           // light attunation of phytoplankton
    real tref = bgc_params[3];         // reference temperature for temp-dependent uptake
    real tr = bgc_params[4];           // coefficient of temperature limitaiton function
    real mpfrac = bgc_params[5];       // mortality of phytoplankton
    real pdia = bgc_params[6];         // size diffusion parameter
    real eff = bgc_params[7];          // efficiency of Z on P grazing
    real seff = bgc_params[8];         // efficiency of Z on Z grazing
    real kp = bgc_params[9];           // nutrient saturation parameter for grazing
    real dxg = bgc_params[10];         // grazing profile witdh
    
    // The following are in units of 1/day
    real mztwo = bgc_params[11];       // density dependent zooplankton grazing
    real rmax = bgc_params[12];        // remineralization rate
    real wsink = bgc_params[13];       // sinking speed of detritus
    int idxUptake = bgc_params[14];    // index for uptake in bgcRates vector
    int idxSat = bgc_params[15];       // index for phytoplankton nutrient saturation in bgcRates vector
    int idxGrazing = bgc_params[16];   // index for grazing in bgcRates vector
    int idxPredPrey = bgc_params[17];  // index for predator prey length ratio in bgcRates vector
    
    // calculate the indexes for phytoplankton
    int idx_phyto = idx_nitrate+1;
    int idx_zoo = idx_nitrate+MP+1;
    int idx_detritus = idx_nitrate+MP+MZ+1;
    
    // placeholders for the physical variables
    real I0 = qsw*0.45;
    real IR = 0;
    real Isurf = qsw*0.45;
    real ifrac = 0;                     // fraction of light dependent uptake
    real tfrac = 0;                     // fraction of temperature dependent uptake in the cell
    real kpar = 0;
    real T = phi[idx_buoy][j][k];       // temperatre in the cell
    real dz = 0;
    
    // placeholders for the biogeochemcial variables
    real N = phi[idx_nitrate][j][k];
    real P = 0;
    real ZO = 0;
    real D = phi[idx_detritus][j][k];
    real ptot = 0;
    real ztot = 0;
    real umax = 0;
    real gmax = 0;
    real kn = 0;
    real mmentis = 0;
    real uptake = 0;
    real uptot = 0;
    real remin = 0;
    real sf_flux = 0;
    
    real grazing = 0;
    real pf = 0;
    real gp = 0;
    real gz = 0;
    real gptot = 0;
    real gztot = 0;
    real gmess = 0;
    
    real pmort = 0;
    real zmort = 0;
    real dmort = 0;
    real Dsink = 0;
    
    real fi = 0;
    real fo = 0;

    // calculate the light in the kth cell
    for (kk = Nz-1; kk >= k; kk --) // integrate from the top down
    {
        ptot = 0;
        
        // sum the total phytoplankton in the cell
        for (ip = 0; ip < MP; ip ++)
        {
            ptot += phi[idx_phyto+ip][j][kk];
//            ptot = 0;
        }

        // calculate the coefficient of light attenuation
        dz = -(ZZ_w[j][kk+1] - ZZ_w[j][kk]);
        kpar = kw + ptot*kc;
        IR = I0/(1 - (dz*kpar) ); // integrate vertically
        I0 = IR;
    }

    // temperature and light dependence
    tfrac = exp(tr*(T - tref));
    ifrac = I0/sqrt( pow(I0,2) + pow(Isurf,2) );

//    fprintf(stderr,"ZZ = %f tfrac = %.3e, ifrac %.3e \n",ZZ_phi[j][k],tfrac,ifrac);


    /*
     /
     / START OF THE BIOGEOCHEMCIAL MODEL EQUATIONS
     /
     */

    //
    // Build the nutrient uptake profile
    //
    dmort = 0;
    uptot = 0;
    for (ip = 0; ip < MP; ip++)
    {
        kn = bgcRates[ip][idxSat];
        mmentis = N/(N + kn);
        umax = bgcRates[ip][idxUptake];
        P = phi[idx_nitrate+1+ip][j][k];
        uptake = (umax/day)*tfrac*ifrac*mmentis*P;
        uptake = uptake*(1-exp(-N)); // limits excess uptake
        uptot += uptake;
        UPvec[ip] = uptake;
        
        
        pmort = mpfrac*(umax/day)*P;
        MPvec[ip] = pmort;
        dmort += pmort;
    }
    
    //
    // Build grazing of phytoplankton by zooplankton
    //
    for (iz = 0; iz < MZ; iz ++)
    {
        pf = 0;
        for (ip = 0; ip < MP; ip++)
        {
            P = phi[idx_phyto+ip][j][k];
            pf += theta_p[ip][iz]*P;        // calculate the total bioavailable phytoplankton
        }
        PTvec[iz] = pf;
    }
    
    // calculate the grazing profile
    for (iz = 0; iz < MZ; iz ++)
    {
        gmax = bgcRates[iz][idxGrazing];
        ZO = phi[idx_zoo+iz][j][k];
        for (ip = 0; ip < MP; ip++)
        {
            P = phi[idx_phyto+ip][j][k];
            grazing = (gmax/day)*PTvec[ip]/( PTvec[ip]*P + kp  );
            grazing = grazing*(1-exp(-PTvec[ip]));
            grazing = grazing*theta_p[ip][iz];
            grazMat[ip][iz] = grazing*P*ZO;
//            fprintf(stderr,"ip = %d iz = %d, Gmat = %.3e \n",ip,iz,theta_p[ip][iz]);
        }
    }
    
    // Sum all the grazing calculated previously
    gptot = 0;
    for (ip = 0; ip < MP; ip++)
    {
    gp = 0;
        for (iz = 0; iz < MZ; iz ++)
        {
            gp += grazMat[ip][iz];
    
        }
        GPvec[ip] = gp;
        gptot += gp;
    }
    
    
    // sum growth for Z
    gztot = 0;
    for (iz = 0; iz < MZ; iz ++)
    {
    gz = 0;
        for (ip = 0; ip < MP; ip++)
        {
            gz += eff*grazMat[ip][iz];
            
        }
        GZvec[iz] = gz;
        gztot += gz;
    }
    
    // detritus
    gmess = gptot - gztot; // enforced to conserve nutrients because of numerical error.
    remin = (rmax/day)*D;
    
    // zooplankton mortality
    for (iz = 0; iz < MZ; iz ++)
    {
        ztot += phi[idx_zoo+iz][j][k];
    }
    
    // calculate the density dependent mortality
    for (iz = 0; iz < MZ; iz ++)
    {
        ZO = phi[idx_zoo+iz][j][k];
        zmort = (mztwo/day)*ztot*ZO;
        MZvec[iz] = zmort;
        dmort += zmort;
    }
    
    // Calculate sinking fluxes
       if (k == Nz-1) // This is the surface, so there is no flux through the surface.
       {
           fi = sf_flux;
           fo = 0.5*fabs(wsink/day)*(phi[idx_detritus][j][k] + phi[idx_detritus][j][k-1]);
       }
       else if (k == 0) // there is no flux out of the domain
       {
           
       }
       else
       {
           fi = 0.5*(wsink/day)*(phi[idx_detritus][j][k+1] + phi[idx_detritus][j][k]);
           fo = 0.5*fabs(wsink/day)*(phi[idx_detritus][j][k] + phi[idx_detritus][j][k-1]);
       }
       
       Dsink = (fi-fo)*_dz_phi[j][k];
    
    
    // update all time tendencies
    dphi_dt[idx_nitrate][j][k] += -uptot + remin;
    for (ip = 0; ip < MP; ip++)
    {
        if (phi[idx_phyto+ip][j][k] < 1e-8) // set a minimum value for phytoplankton
        {
            phi[idx_phyto+ip][j][k] = 1e-8;
        }
        else // if this is above the value, then update the time tendencies
        {
            dphi_dt[idx_phyto+ip][j][k] += UPvec[ip] - GPvec[ip] - MPvec[ip];
        }
        
    }
    for (iz = 0; iz < MZ; iz ++)
    {
        if (phi[idx_zoo+iz][j][k] < 1e-9)  // set a minimum value for zooplankton (an order of magnitude smaller than P)
        {
            phi[idx_zoo+iz][j][k] = 1e-9;
        }
        else // update time tendencies if it is
        {
            dphi_dt[idx_zoo+iz][j][k] += GZvec[iz] - MZvec[iz];
        }
    }
    dphi_dt[idx_detritus][j][k] += dmort + gmess - remin + Dsink;

}















/**
 * doSizeDiffusion
 *
 * Calculates the size diffusion for phytoplankton and zooplankton
 *
 */
void doSizeDiffusion(const real t, const int j, const int k, real *** phi, real *** dphi_dt)
{
//    int ip, iz;
//    real pdia = bgc_params[6];         // size diffusion parameter
//    real dxp = 0;                      // Difference in size for plus one
//    real dxm = 0;                      // Difference in size for minus one
//    real dp = 0;                       // Coefficient for plus one
//    real dc = 0;                       // Coefficient for center
//    real dm = 0;                       // Coefficient for minus one
//    real PP = 0;                       // P plus one
//    real PC = 0;                       // P center
//    real PM = 0;                       // P minuns one
//
//    // calculate the indexes for phytoplankton
//    int idx_phyto = idx_nitrate+1;
//    int idx_zoo = idx_nitrate+MP+1;
//    int idx_detritus = idx_nitrate+MP+MZ+1;
//
//
//    // Calculate the effect of size diffusion
//    for (ip = 0; ip < MP; ip++)
//    {
//        if (ip == 0) // boundary conditions
//        {
//            PP = phi[idx_phyto+ip+1][j][k];
//            dphi_dt[idx_phyto+ip][j][k] += PP;
//        }
//        else if (ip == MP-1)
//        {
//            PM = phi[idx_phyto+ip-1][j][k];
//            dphi_dt[idx_phyto+ip][j][k] += PM;
//        }
//        else
//        {
//            PP = phi[idx_phyto+ip+1][j][k];
//            PC = phi[idx_phyto+ip][j][k];
//            PM = phi[idx_phyto+ip-1][j][k];
//
//            // calculate the distances in size space for size diffusion
//            dxp = log10(lpvec[ip+1]) - log10(lpvec[ip]);
//            dxm = log10(lpvec[ip]) - log10(lpvec[ip-1]);
//
//            dm = pdia/(pow(dxm,2));
//            dc = -pdia*( 1/(dxp+dxm) + 1/(pow(dxm,2)) );
//            dp = pdia/(dxp*dxm);
//
//            dphi_dt[idx_phyto+ip][j][k] += (dp*PP + dc*PC + dm*PM);
//        }
//    }
//
//
//
    
    
    
}




















/**
 * tderiv_bgc
 *
 * Calculates the time tendency of any tracer phi due to biogeochemical processes.
 *
 */
void tderiv_bgc (const real t, real *** phi, real *** dphi_dt)
{
    int i,j,k; // Physical coutners
    

    if (bgcModel == BGC_NONE) // only need to check if there is no bgc
    {
        
    }
    else if (bgcModel == BGC_NPZD)  // both bgc models use the same base model equations
    {
        for (j = 0; j < Nx; j++)
        {
            for (k = Nz-1; k >= 0; k--)
            {
                npzd(t,j,k,phi,dphi_dt);        // calculate the npzd model
            }
        }
    }
    else
    {
        for (j = 0; j < Nx; j++)
        {
            for (k = Nz-1; k >= 0; k--)
            {
                npzd(t,j,k,phi,dphi_dt);        // calculate ROEM
            }
        }
    }
    
}










/**
 * tderiv_mom
 *
 * Calculates the time tendency of momentum tracers.
 *
 */

void tderiv_mom (const real t, real *** phi, real *** dphi_dt)
{
    // Looping variables
    int i, j, k;
    
    // Pointers to velocity and buoyancy matrices
    real ** uvel = NULL;
    real ** vvel = NULL;
    real ** buoy = NULL;
    real ** du_dt = NULL;
    real ** dv_dt = NULL;
    
    // Pointers to velocity and buoyancy matrices
    uvel = phi[idx_uvel];
    vvel = phi[idx_vvel];
    buoy = phi[idx_buoy];
    du_dt = dphi_dt[idx_uvel];
    dv_dt = dphi_dt[idx_vvel];
    
    calcPressure(t,buoy);
    
    // Calculate the tendency due to the coriolis force and add to the pressure term
    for (j = 1; j < Nx; j++)
    {
        for (k = 0; k < Nz; k++)
        {
            du_dt[j][k] = f0*vvel[j][k] - BPx[j][k];
            dv_dt[j][k] = -f0*uvel[j][k] - BPy[j][k];
            
            if (j == 1)
            {
                du_dt[j][k] += nu_h*(uvel[j+1][k]-uvel[j][k])/dxsq;
                dv_dt[j][k] += nu_h*(vvel[j+1][k]-vvel[j][k])/dxsq;
            }
            else if (j == Nx-1)
            {
                
                du_dt[j][k] += nu_h*(uvel[j-1][k]-uvel[j][k])/dxsq;
                dv_dt[j][k] += nu_h*(vvel[j-1][k]-vvel[j][k])/dxsq;
            }
            else
            {
                du_dt[j][k] += nu_h*(uvel[j+1][k]-2*uvel[j][k]+uvel[j-1][k])/dxsq;
                dv_dt[j][k] += nu_h*(vvel[j+1][k]-2*vvel[j][k]+vvel[j-1][k])/dxsq;
            }
            
            if (k == 0)
            {
                du_dt[j][k] += nu_v * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz)
                                    * (uvel[j][k+1]-uvel[j][k])/((ZZ_u[j][k+1]-ZZ_u[j][k])*(ZZ_psi[j][k+1]-ZZ_psi[j][k]));
                dv_dt[j][k] += nu_v * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz)
                                    * (vvel[j][k+1]-vvel[j][k])/((ZZ_u[j][k+1]-ZZ_u[j][k])*(ZZ_psi[j][k+1]-ZZ_psi[j][k]));
            }
            else if (k == Nz-1)
            {
                du_dt[j][k] += nu_v * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz)
                                    * (uvel[j][k-1]-uvel[j][k])/((ZZ_u[j][k]-ZZ_u[j][k-1])*(ZZ_psi[j][k+1]-ZZ_psi[j][k]));
                dv_dt[j][k] += nu_v * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz)
                                    * (vvel[j][k-1]-vvel[j][k])/((ZZ_u[j][k]-ZZ_u[j][k-1])*(ZZ_psi[j][k+1]-ZZ_psi[j][k]));
            }
            else
            {
                du_dt[j][k] += nu_v
                            * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz) // Scale with grid
                            * (
                                    (uvel[j][k+1]-uvel[j][k])/(ZZ_u[j][k+1]-ZZ_u[j][k]) + (uvel[j][k-1]-uvel[j][k])/(ZZ_u[j][k]-ZZ_u[j][k-1])
                               )
                              / (ZZ_psi[j][k+1]-ZZ_psi[j][k]);
                dv_dt[j][k] += nu_v
                            * SQUARE(ZZ_psi[j][k+1]-ZZ_psi[j][k])/SQUARE(Lz/Nz) // Scale with grid
                            * (
                                    (vvel[j][k+1]-vvel[j][k])/(ZZ_u[j][k+1]-ZZ_u[j][k]) + (vvel[j][k-1]-vvel[j][k])/(ZZ_u[j][k]-ZZ_u[j][k-1])
                               )
                                / (ZZ_psi[j][k+1]-ZZ_psi[j][k]);
            }
            
        }
    
    }
    
    // Add surface/bottom momentum fluxes
    for (j = 1; j < Nx; j++)
    {
        du_dt[j][0] -= r_bbl*uvel[j][0]/(ZZ_psi[j][1] - ZZ_psi[j][0]); // bottom momentum flux
        dv_dt[j][0] -= r_bbl*vvel[j][0]/(ZZ_psi[j][1] - ZZ_psi[j][0]);
        
        
        dv_dt[j][Nz-1] += tau[j]/(rho0*(ZZ_psi[j][Nz] - ZZ_psi[j][Nz-1])); // surface momentum flux due to wind forcing
        
    }
    
    for (k = 0; k < Nz; k++)
    {
        // u,v tendencies are zero at the western edge of the domain
        du_dt[0][k] = 0;
        dv_dt[0][k] = 0;
    }
    
    
}











/**
 * do_adv_diff
 *
 * Calculates the time tendency of any tracer phi due to advection/diffusion.
 *
 */
void do_adv_diff (  const real    t,
                  real **       phi,
                  real **       phi_north,
                  real **       phi_south,
                  real **       dphi_dt,
                  bool          is_buoy,
                  real **       u_r,
                  real **       vvel,
                  real **       w_r,
                  real **       Kiso_u,
                  real **       Kiso_w,
                  real **       Siso_u,
                  real **       Siso_w)
{
    /// Looping variables
    int j,k;
    
    // True isopycnal slope, relative to sigma coordinates
    real Siso_true = 0;
    real z;
    
    // meridional velocity
    real v_r = 0; // meridional velocity holder
    
    /////////////////////////////////////////////////////
    ///// Arrays used by the Kurganov-Tadmor scheme /////
    /////////////////////////////////////////////////////
    
    if (advectionScheme == ADVECTION_KT00)
    {
        // Near-boundary derivatives
        for (j = 0; j < Nx; j ++)
        {
            dphi_dz[j][0] = (phi[j][1]-phi[j][0]) * _dz_w[j][1];
            dphi_dz[j][Nz-1] = (phi[j][Nz-1]-phi[j][Nz-2]) * _dz_w[j][Nz-1];
        }
        for (k = 0; k < Nz; k ++)
        {
            dphi_dx[0][k] = (phi[1][k]-phi[0][k]) * _dx;
            dphi_dx[Nx-1][k] = (phi[Nx-1][k]-phi[Nx-2][k]) * _dx;
        }
        
#pragma parallel
        
        // Determine limited y-slopes at cell centres
        for (j = 1; j < Nx-1; j ++)
        {
            for (k = 0; k < Nz; k ++)
            {
                dphi_dx[j][k] = minmod( KT00_sigma * (phi[j+1][k]-phi[j][k]) * _dx,
                                       (phi[j+1][k]-phi[j-1][k]) * _2dx,
                                       KT00_sigma * (phi[j][k]-phi[j-1][k]) * _dx  );
            }
        }
        
#pragma parallel
        
        // Determine limited z-slopes at cell centres
        for (j = 0; j < Nx; j ++)
        {
            for (k = 1; k < Nz-1; k ++)
            {
                dphi_dz[j][k] = minmod( KT00_sigma * (phi[j][k+1]-phi[j][k]) * _dz_w[j][k+1],
                                       (phi[j][k+1]-phi[j][k-1]) * _2dz_phi[j][k],
                                       KT00_sigma * (phi[j][k]-phi[j][k-1]) * _dz_w[j][k]  );
            }
        }
        
#pragma parallel
        
        // Interpolate phi to cell y-faces
        for (j = 1; j < Nx; j ++)
        {
            for (k = 0; k < Nz; k ++)
            {
                phi_xm[j][k] = phi[j-1][k] + 0.5*dx*dphi_dx[j-1][k];
                phi_xp[j][k] = phi[j][k] - 0.5*dx*dphi_dx[j][k];
            }
        }
        
#pragma parallel
        
        // Interpolate phi to cell z-faces
        for (j = 0; j < Nx; j ++)
        {
            for (k = 1; k < Nz; k ++)
            {
                phi_zm[j][k] = phi[j][k-1] + (ZZ_w[j][k]-ZZ_phi[j][k-1])*dphi_dz[j][k-1];
                phi_zp[j][k] = phi[j][k] + (ZZ_w[j][k]-ZZ_phi[j][k])*dphi_dz[j][k];
            }
        }
    }
    
    
    
    /////////////////////////////////
    ///// Calculation of fluxes /////
    /////////////////////////////////
    
#pragma parallel
    
    // Calculate y-fluxes
    for (k = 0; k < Nz; k ++)
    {
        // No flux across boundaries
        HHx[0][k] = 0;
        HHx[Nx][k] = 0;
        PPx[0][k] = 0;
        PPx[Nx][k] = 0;
        
        for (j = 1; j < Nx; j ++)
        {
            // Current depth
            z = ZZ_u[j][k];
            
            // True isopycnal slope defines direction of isopycnal mixing in sigma coordinates
            Siso_true = Siso_u[j][k] - (ZZ_phi[j][k]-ZZ_phi[j-1][k])*_dx;
            
            // Flux calculation depends on spatial discretisation method
            switch (advectionScheme)
            {
                case ADVECTION_CENTERED:
                {
                    // Hyperbolic fluxes
                    HHx[j][k] = u_r[j][k] * 0.5*(phi[j][k]+phi[j+1][k]);
                    
                    // Parabolic fluxes
                    if ( is_buoy && (!use_sml || (z < -Hsml)) && (!use_bbl || (z > -hb_psi[j]+Hbbl)) )
                    {
                        PPx[j][k] = 0;
                    }
                    else
                    {
                        PPx[j][k] = Kiso_u[j][k] * (phi[j][k]-phi[j-1][k]) * _dx;
                        
                        if (k == 0)
                        {
                            PPx[j][k] += Kiso_u[j][k] * Siso_true
                            * 0.5 * ( (phi[j][k+1]-phi[j][k])*_dz_w[j][k+1] + (phi[j-1][k+1]-phi[j-1][k])*_dz_w[j-1][k+1] );
                        }
                        else if (k == Nz-1)
                        {
                            PPx[j][k] += Kiso_u[j][k] * Siso_true
                            * 0.5 * ( (phi[j][k]-phi[j][k-1])*_dz_w[j][k] + (phi[j-1][k]-phi[j-1][k-1])*_dz_w[j-1][k] );
                        }
                        else
                        {
                            PPx[j][k] += Kiso_u[j][k] * Siso_true
                            * 0.5 * ( (phi[j][k+1]-phi[j][k-1])*_2dz_phi[j][k] + (phi[j-1][k+1]-phi[j-1][k-1])*_2dz_phi[j-1][k] );
                        }
                    }
                    
                    break;
                }
                case ADVECTION_KT00:
                {
                    // Hyperbolic fluxes
                    HHx[j][k] = 0.5*( u_r[j][k]*(phi_xm[j][k]+phi_xp[j][k])
                                     - fabs(u_r[j][k])*(phi_xp[j][k]-phi_xm[j][k]) );
                    
                    // Parabolic fluxes
                    if ( is_buoy && (!use_sml || (z < -Hsml)) && (!use_bbl || (z > -hb_psi[j]+Hbbl)) )
                    {
                        PPx[j][k] = 0;
                    }
                    else
                    {
                        PPx[j][k] = Kiso_u[j][k] * (
                                                    (phi[j][k]-phi[j-1][k]) * _dx
                                                    + Siso_true * 0.5*(dphi_dz[j][k]+dphi_dz[j-1][k]) );
                    }
                    
                    // Boundary layer fluxes
                    
                    break;
                }
                default:
                {
                    fprintf(stderr,"ERROR: Unknown spatial discretisation specified\n");
                    break;
                }
            }
            
            // Using sigma coordinates introduces a factor of hb in the fluxes. This is divided
            // out after the flux divergence is taken
            HHx[j][k] *= dz_u[j][k];
            PPx[j][k] *= dz_u[j][k];
        }
    }
    
#pragma parallel
    
    // Calculate z-fluxes
    for (j = 0; j < Nx; j ++)
    {
        // No flux across boundaries
        HHz[j][0] = 0;
        HHz[j][Nz] = 0;
        PPz[j][0] = 0;
        PPz[j][Nz] = 0;
        
        // Interior gridpoints
        for (k = 1; k < Nz; k ++)
        {
            // Current depth
            z = ZZ_w[j][k];
            
            // True isopycnal slope defines direction of isopycnal mixing in sigma coordinates
            Siso_true = Siso_w[j][k] - (ZZ_psi[j+1][k]-ZZ_psi[j][k])*_dx;
            
            // Flux calculation depends on spatial discretisation method
            switch (advectionScheme)
            {
                case ADVECTION_CENTERED:
                {
                    HHz[j][k] = w_r[j][k] * 0.5*(phi[j][k-1]+phi[j][k]);
                    
                    if ( is_buoy && (!use_sml || (z < -Hsml)) && (!use_bbl || (z > -hb_phi[j]+Hbbl)) )
                    {
                        PPz[j][k] = 0;
                    }
                    else
                    {
                        PPz[j][k] = (Kiso_w[j][k]*Siso_true*Siso_true) * (phi[j][k]-phi[j][k-1]) * _dz_w[j][k];
                        
                        if (j == 0)
                        {
                            PPz[j][k] += Kiso_w[j][k]*Siso_true
                            * 0.5*(phi[j+1][k]-phi[j][k]+phi[j+1][k-1]-phi[j][k-1]) * _dx;
                        }
                        else if (j == Nx-1)
                        {
                            PPz[j][k] += Kiso_w[j][k]*Siso_true
                            * 0.5*(phi[j][k]-phi[j-1][k]+phi[j][k-1]-phi[j-1][k-1]) * _dx;
                        }
                        else
                        {
                            PPz[j][k] += Kiso_w[j][k]*Siso_true
                            * 0.5*(phi[j+1][k]-phi[j-1][k]+phi[j+1][k-1]-phi[j-1][k-1]) * _2dx;
                        }
                    }
                    
                    break;
                }
                case ADVECTION_KT00:
                {
                    HHz[j][k] = 0.5*( w_r[j][k]*(phi_zm[j][k]+phi_zp[j][k])
                                     - fabs(w_r[j][k])*(phi_zp[j][k]-phi_zm[j][k]) );
                    
                    // Parabolic fluxes
                    if ( is_buoy && (!use_sml || (z < -Hsml)) && (!use_bbl || (z > -hb_phi[j]+Hbbl)) )
                    {
                        PPz[j][k] = 0;
                    }
                    else
                    {
                        PPz[j][k] = Kiso_w[j][k] * Siso_true *
                        (
                         Siso_true * (phi[j][k]-phi[j][k-1]) * _dz_w[j][k]
                         + 0.5*(dphi_dx[j][k]+dphi_dx[j][k-1])
                         );
                    }
                    
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
    
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ///// Calculate the meridional component of advection /////
    //////////////////////////////////////////////////////////////////////////////////////////
    
//    for (j = 0; j < Nx; j++)
//    {
//        for (k = 0; k < Nz; k++)
//        {
//            v_r = 0.5*(vvel[j][k] + vvel[j+1][k]);
//            if (v_r < 0)
//            {
//                dphi_dt[j][k] += _Ly*v_r*(phi[j][k] - phi_north[j][k]);
//            }
//            else if (v_r > 0)
//            {
//                dphi_dt[j][k] += _Ly*v_r*(phi[j][k] - phi_south[j][k]);
//            }
//        }
//    }
    
    
    /////////////////////////////////////////////
    ////// Calculate time derivative of phi /////
    /////////////////////////////////////////////
    
#pragma parallel
    
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Add fluxes
            dphi_dt[j][k] += ( - (HHx[j+1][k] - HHx[j][k]) * _dx * _dz_phi[j][k] // Hyperbolic fluxes
                              - (HHz[j][k+1] - HHz[j][k]) * _dz_phi[j][k]
                              + (PPx[j+1][k] - PPx[j][k]) * _dx * _dz_phi[j][k] // Parabolic fluxes
                              + (PPz[j][k+1] - PPz[j][k]) * _dz_phi[j][k] );
        }
    }
    
}











/**
 * tderiv_adv_diff
 *
 * Calculates the time tendency of all tracers due to advection/diffusion.
 * Returns the largest time step that keeps both the advective and diffusive
 * parts of the tracer advection equation stable.
 *
 */
real tderiv_adv_diff (const real t, real *** phi, real *** dphi_dt)
{
    // Looping variables
    int i, j, k;
    
    // Pointers to velocity and buoyancy matrices
    real ** uvel = NULL;
    real ** vvel = NULL;
    real ** buoy = NULL;
    bool is_buoy = false;
    
    // True isopycnal slope, relative to sigma coordinates
    real Sgm_true = 0;
    real Siso_true = 0;
    
    // For CFL calculations
    real cfl_dt = 0;
    real cfl_u = 0;
    real cfl_w = 0;
    real cfl_y = 0;
    real cfl_z = 0;;
    real cfl_igw = 0;
    real cfl_sink = 0;
    real cfl_bgc = 0;
    real u_max = 0;
    real w_dz = 0;
    real w_dz_max = 0;
    real xdiff = 0;
    real xdiff_max = 0;
    real zdiff_dzsq = 0;
    real zdiff_dzsq_max = 0;
    real cfl_phys = 0;
    real n_col = 0;
    real n_max = 0;
    real wsink = 0;
    
    // Slope calculations
    real siso = 0;
    real sgm = 0;
    real z = 0;
    real sink_dz = 0;
    real sink_dz_max = 0;
    real umax = 0;
    real wmax = 0;
    real uval = 0;
    real siso_max = 0;
    real sgm_max = 0;
    
    ////////////////////////////////////////
    ///// BEGIN CALCULATING TENDENCIES /////
    ////////////////////////////////////////
    
    // Pointers to velocity and buoyancy matrices
    uvel = phi[idx_uvel];
    vvel = phi[idx_vvel];
    buoy = phi[idx_buoy];
    
    // Calculate isopycnal slopes
    calcSlopes(t,buoy,Sgm_psi,Siso_u,Siso_w);
    
    // Calculate Gent-McWilliams diffusivity Kgm
    calcKgm(t,buoy,Kgm_psi,Kgm_u,Kgm_w);
    
    // Calculate isopycnal diffusivity Kiso
    calcKiso(t,buoy,Kiso_u,Kiso_w);
        
    // Calculate mean and eddy components of overturning streamfunction
    calcPsim(t,uvel,psi_m);
    calcPsie(t,buoy,Kgm_psi,Sgm_psi,psi_e);
    
    // Compute residual streamfunction and velocities
    calcPsir(psi_m,psi_e,psi_r,u_r,w_r);
    
    // Loop over tracers and calculate advective/diffusive tendency
    for (i = 0; i < Ntracs; i ++)
    {
        // No advection/diffusion for momentum - their tendency is handled separately in tderiv_mom
        if ((i == idx_uvel) || (i == idx_vvel))
        {
            continue;
        }
        
        // Advection/diffusion depends on whether the tracer is the buoyancy variable
        if (i == idx_buoy)
        {
//            is_buoy = true;
            is_buoy = false; // testing the buoyancy mixing along isopycnals
        }
        else
        {
            is_buoy = false;
        }
        
        do_adv_diff(t,phi[i],phi_north[i],phi_south[i],dphi_dt[i],is_buoy,u_r,vvel,w_r,Kiso_u,Kiso_w,Siso_u,Siso_w);
    }
    
    //////////////////////////////////////
    ///// END CALCULATING TENDENCIES /////
    //////////////////////////////////////
    
    
    
    //////////////////////////////////
    ///// BEGIN CALCULATING CFLS /////
    //////////////////////////////////
    
#pragma parallel
    
    // For x-fluxes
    for (j = 0; j <= Nx; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Max advecting velocity
            uval = fabs(u_r[j][k]);
            if (uval > u_max)
            {
                u_max = uval;
            }
            
            /// Max effective diffusivity
            xdiff = fmax(Kiso_u[j][k],Kgm_u[j][k]);
            if (xdiff > xdiff_max)
            {
                xdiff_max = xdiff;
            }
        }
    }
    
    
#pragma parallel
    
    // For z-fluxes
    // N.B. Carefully chosen indices so that we can calculate diffusive CFLs
    // associated with Kgm and Kiso (almost) everywhere.
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k <= Nz; k ++)
        {
            z = ZZ_u[j][k];
            // Max advecting velocity
            w_dz = fabs(w_r[j][k]) * _dz_w[j][k];
            if (w_dz > w_dz_max)
            {
                w_dz_max = w_dz;
            }
            
//            // Calculate CFL conditions for sinking speeds
            sink_dz = wsink * _dz_w[j][k];
            if (sink_dz > sink_dz_max)
            {
                sink_dz_max = sink_dz;

            }
            
            // Max effective diffusivity
            // A.L.S.: For reasons I don't fully understand, the GM eddy advection
            // stability follows a diffusive scaling rather than an advective one.
            // We therefore take the maximum of Kgm and Kiso here. This CFL constraint
            // (i.e. pseudo-vertical diffusivity) tends to be the limiting one.
            
            // "true" slopes in x/sigma coordinates, i.e. slopes relative to grid coordinates
            Sgm_true = j>0 ? Sgm_psi[j][k] - (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx : 0; // Basically just skips j==0 case
            Siso_true = Siso_w[j][k] - (ZZ_psi[j+1][k]-ZZ_psi[j][k])*_dx;

            
            // N.B. Indexing here is slightly off for the GM component, but good enough for a CFL estimate
            zdiff_dzsq = fmax( Kiso_w[j][k]*SQUARE(Siso_true)*_dzsq_w[j][k], Kgm_psi[j][k]*SQUARE(Sgm_true)*_dzsq_psi[j][k] );
            if (zdiff_dzsq > zdiff_dzsq_max)
                {
                    zdiff_dzsq_max = zdiff_dzsq;
            }
            
            siso = fabs(Sgm_true);
            sgm = fabs(Siso_true);
            if (sgm > sgm_max)
            {
                sgm_max = sgm;
            }
            if (siso > sgm_max)
            {
                siso_max = siso;
            }
        }
    }
    
    for (j = 1; j < Nx; j++)
    {
        for (k = 1; k < Nz; k++)
        {
            Nbuoy[j][k] = sqrt(0.5*alpha*grav*( (buoy[j][k] - buoy[j][k-1])*_dz_w[j][k] + (buoy[j-1][k] - buoy[j-1][k-1])*_dz_w[j-1][k] )  ) ;
        }
    }
    
    
    // Integrate the N term
    for (j = 0; j < Nx; j++)
    {
        n_col = 0;
        for (k = Nz-1; k > 0; k--)
        {
            n_col += 0.5*(Nbuoy[j][k] + Nbuoy[j][k-1])/_dz_phi[j][k];
        }
        
        if (n_max < n_col)
        {
            n_max = n_col/_PI;
        }
    }
    
    
    
    // Calculate CFL criteria
    cfl_u = 0.5*dx/u_max;
    cfl_w = 0.5/w_dz_max;
    cfl_y = 0.5*dxsq/xdiff_max;
    cfl_z = 0.5/zdiff_dzsq_max;
    cfl_igw = 0.5*dx/n_max;
    cfl_sink = 0.5*dx/sink_dz_max;

    // Actual CFL-limted time step
    cfl_phys = fmin(fmin(cfl_u,cfl_w),fmin(cfl_y,cfl_z));
    cfl_phys = fmin(fmin(cfl_sink,cfl_igw),cfl_phys);
    cfl_dt = cfl_phys;
    
    if ( cfl_phys > 60*60*24 )
    {
        fprintf(stderr,"!!: Large timestep detected (> 1 day), printing cfl_dt \n");
        fprintf(stderr,"cfl_u = %f, cfl_w = %f, cfl_y = %f \n", cfl_u, cfl_w, cfl_y);
        fprintf(stderr,"cfl_z = %f, cfl_sink = %f, cfl_igw = %f \n", cfl_z, cfl_sink, cfl_igw);
    }
    
    ////////////////////////////////
    ///// END CALCULATING CFLS /////
    ////////////////////////////////
    
    
    return cfl_dt;
}




















/**
 * tderiv
 *
 * Calculates the time derivatives at all buoyancy/tracer gridpoints. Returns the
 * largest time step that keeps both the advective and diffusive parts
 * of the buoyancy/tracer advection-diffusion equations stable.
 *
 */
real tderiv (const real t, const real * data, real * dt_data, const uint numvars)
{
    // CFL timesteps
    real cfl_dt = 0;
    
    // Construct output matrix from the dt_data vector
    memset(dt_data,0,numvars*sizeof(real));
    vec2mat3(dt_data,&dphi_dt_wrk,Ntracs,Nx,Nz);
    
    // Copy tracer data from the 'data' array to the work array
    memcpy(phi_wrk_V,data,Ntot*sizeof(real));
    
    // Calculate momentum tendencies
    if (momentumScheme != MOMENTUM_NONE)
    {
        tderiv_mom (t, phi_wrk, dphi_dt_wrk);
    }
    
    // Calculate tracer tendencies due to advection/diffusion
    cfl_dt = tderiv_adv_diff (t, phi_wrk, dphi_dt_wrk);

    // Calculate tracer tendencies due to biogeochemistry
    tderiv_bgc (t, phi_wrk, dphi_dt_wrk);
    
    // Calculate tracer tendencies due to relaxation
    tderiv_relax (t, phi_wrk, dphi_dt_wrk);
    
    return cfl_dt;
}













/**
 *
 * do_impl_diff
 *
 * Computes diapycnal diffusivity and performs implicit vertical diffusion of tracers.
 *
 */
void do_impl_diff (real t, real dt, real *** phi)
{
    int i,j,k;
    real ** buoy;
    
    // Calculate diapycnal mixing coefficient
    buoy = phi[idx_buoy];
    calcKdia(t,buoy,Kdia_w);
    
#pragma parallel
    
    // Solve for each column
    for (j = 0; j < Nx; j ++)
    {
        // Solve for each tracer
        for (i = 0; i < Ntracs; i ++)
        {
            // Implicit vertical diffusion coefficients
            for (k = 0; k < Nz; k ++)
            {
                impl_A[k] = k == 0 ? 0 : - Kdia_w[j][k] * dt * _dz_w[j][k] * _dz_phi[j][k];
                impl_C[k] = k == Nz-1 ? 0 : - Kdia_w[j][k+1] * dt * _dz_w[j][k+1] * _dz_phi[j][k];
                impl_B[k] = 1 - impl_A[k] - impl_C[k];
                impl_D[k] = phi[i][j][k];
            }
            
            // Solve algebraically
            thomas(impl_A,impl_B,impl_C,impl_D,phi[i][j],Nz);
        }
    }
    
}






/**
 *
 * do_pressure_correct
 *
 * Applies zonal pressure gradient correction to enforce non-divergence of the zonal/vertical mean flow.
 * This amounts to simply subtracting the depth-averaged velocity from each u-point.
 *
 */
void do_pressure_correct (real *** phi)
{
    int j,k;
    real ** uvel = NULL;
    real * udz = NULL;
    real u_avg = 0;
    
    // Extract zonal velocity matrix
    uvel = phi[idx_uvel];
    
    // We'll abuse the zonal velocity at the western edge of the domain,
    // which is just zero everywhere, to serve as a work vector for
    // calculating depth-integrated velocities
    udz = uvel[0];
    
#pragma parallel
    
    // Loop through columns
    for (j = 1; j < Nx; j ++)
    {
        
        // Implicit vertical diffusion coefficients
        for (k = 0; k < Nz; k ++)
        {
            udz[k] = uvel[j][k] * dz_u[j][k];
        }
        
        // Calculate depth-averaged velocity via Kahan summation to ensure numerical accuracy
        u_avg = kahanSum(udz,Nz) / hb_psi[j];
        
        // Subtract average velocity from each u-point
        for (k = 0; k < Nz; k ++)
        {
            uvel[j][k] -= u_avg;
        }
        
    }
    
    // Zero zonal velocity at the western wall
    for (k = 0; k < Nz; k ++)
    {
        uvel[0][k] = 0;
    }
    
}













/**
 *
 * constructOutputName
 *
 * Constructs an output file path from the destination directory (outdir),
 * variable name (varname, defined at the top of this file), tracer number (i),
 * and output index (n). The resulting file path is stored in outfile. If i<0
 * then no tracer number is specified in the output file. If n<0 then no
 * output index is specified.
 *
 */
void constructOutputName (char * outdir, const char * varname, int i, int n, char * outfile)
{
    char nstr[MAX_PARAMETER_FILENAME_LENGTH];
    char istr[MAX_PARAMETER_FILENAME_LENGTH];
    
    // Create a string for the current iteration number
    sprintf(nstr,"%d",n);
    
    // Create a string for the current isopycnal layer
    if (i >= 0)
    {
        sprintf(istr,"%d",i);
    }
    else
    {
        sprintf(istr,"");
    }
    
    // Construct output file path
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,varname);
    strcat(outfile,istr);
    if (n >= 0)
    {
        strcat(outfile,"_n=");
        strcat(outfile,nstr);
    }
    strcat(outfile,".dat");
}













/**
 *
 * writeOutputFile
 *
 * Convenience function to write a 2D matrix to a specified output file.
 *
 */
bool writeOutputFile (char * outfile, real ** mat, uint m, uint n)
{
    FILE * outfd = NULL;
    
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,mat,m,n);
    fclose(outfd);
    return true;
}
















/**
 *
 * writeModelGrid
 *
 * Writes model vertical grid to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
bool writeModelGrid (real ** ZZ_phi, real ** ZZ_psi, real ** ZZ_u, real ** ZZ_w, char * outdir)
{
    char outfile[MAX_PARAMETER_FILENAME_LENGTH];
    
    // Write depths on phi-points
    constructOutputName(outdir,OUTN_ZZ_PHI,-1,-1,outfile);
    if (!writeOutputFile(outfile,ZZ_phi,Nx,Nz)) return false;
    
    // Write depths on psi-points
    constructOutputName(outdir,OUTN_ZZ_PSI,-1,-1,outfile);
    if (!writeOutputFile(outfile,ZZ_psi,Nx,Nz)) return false;
    
    // Write depths on u-points
    constructOutputName(outdir,OUTN_ZZ_U,-1,-1,outfile);
    if (!writeOutputFile(outfile,ZZ_u,Nx,Nz)) return false;
    
    // Write depths on w-points
    constructOutputName(outdir,OUTN_ZZ_W,-1,-1,outfile);
    if (!writeOutputFile(outfile,ZZ_w,Nx,Nz)) return false;
    
    return true;
}











/**
 *
 * writeModelState
 *
 * Writes model variables to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
bool writeModelState (const int t, const int n, real *** phi, char * outdir)
{
    int i = 0;
    real ** buoy = NULL;
    char outfile[MAX_PARAMETER_FILENAME_LENGTH];
    char nstr[MAX_PARAMETER_FILENAME_LENGTH];
    char istr[MAX_PARAMETER_FILENAME_LENGTH];
    FILE * outfd = NULL;
    real ** uvel = NULL;
    real ** vvel = NULL;
    
    // Calculate residual streamfunction
    uvel = phi[idx_uvel];
    vvel = phi[idx_vvel];
    buoy = phi[idx_buoy];
    calcKgm(t,buoy,Kgm_psi,Kgm_u,Kgm_w);
    calcSlopes(t,buoy,Sgm_psi,Siso_u,Siso_w);
    calcPsim(t,uvel,psi_m);
    calcPsie(t,buoy,Kgm_psi,Sgm_psi,psi_e);
    calcPsir(psi_m,psi_e,psi_r,u_r,w_r);
    
    // Loop over tracers
    for (i = 0; i < Ntracs; i ++)
    {
        // Write iteration data for this tracer
        constructOutputName(outdir,OUTN_TRAC,i,n,outfile);
        if (!writeOutputFile(outfile,phi[i],Nx,Nz)) return false;
    }
    
    // Write iteration data for the residual streamfunction
    constructOutputName(outdir,OUTN_PSIR,-1,n,outfile);
    if (!writeOutputFile(outfile,psi_r,Nx+1,Nz+1)) return false;
    
    // Write iteration data for the mean streamfunction
    constructOutputName(outdir,OUTN_PSIM,-1,n,outfile);
    if (!writeOutputFile(outfile,psi_m,Nx+1,Nz+1)) return false;
    
    // Write iteration data for the eddy streamfunction
    constructOutputName(outdir,OUTN_PSIE,-1,n,outfile);
    if (!writeOutputFile(outfile,psi_e,Nx+1,Nz+1)) return false;
    
    // Debugging
    if (debug && t > 0.1)
    {
        // Write iteration data for the file of interest
        constructOutputName(outdir,OUTN_DBDX_CUBIC,-1,n,outfile);
        if (!writeOutputFile(outfile,db_dx,Nx+1,Nz+1)) return false; // cubic buoyancy gradient
        
        
        // Write iteration data for the file of interest
        constructOutputName(outdir,OUTN_DBDX_LINEAR,-1,n,outfile);
        if (!writeOutputFile(outfile,db_dx_lin,Nx+1,Nz+1)) return false; // cubic buoyancy gradient
        
        // Loop over tracers
        for (i = 0; i < Ntracs; i ++)
        {
            // Write iteration data for this tracer
            constructOutputName(outdir,OUTN_TRAC,i,n,outfile);
            if (!writeOutputFile(outfile,phi[i],Nx,Nz)) return false;
        }
        
        // Write iteration data for the residual streamfunction
        constructOutputName(outdir,OUTN_PSIR,-1,n,outfile);
        if (!writeOutputFile(outfile,psi_r,Nx+1,Nz+1)) return false;
        
        // Write iteration data for the mean streamfunction
        constructOutputName(outdir,OUTN_PSIM,-1,n,outfile);
        if (!writeOutputFile(outfile,psi_m,Nx+1,Nz+1)) return false;
        
        // Write iteration data for the eddy streamfunction
        constructOutputName(outdir,OUTN_PSIE,-1,n,outfile);
        if (!writeOutputFile(outfile,psi_e,Nx+1,Nz+1)) return false;
    }
    
    return true;
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
     "USAGE: %s <infile> <outdir> \n"
     "  \n"
     "  where <infile> is the name of the input parameter file, and\n"
     "  <outdir> is the directory into which output files should be\n"
     "  written.\n"
     "  \n"
     "  The input file must specify a series of input parameters\n"
     "  in name-value form, separated by whitespace, i.e.\n"
     "  <name> <value> [<name> <value> [...]]\n"
     "  \n"
     "  All matrix parameters should be formatted in column-major\n"
     "  order.\n"
     "  \n"
     "  name                  value\n"
     "  \n"
     "  Ntracs                Number of tracer variables in the simulation. Must be >3,\n"
     "                        as the first 3 tracers are reserved for zonal velocity,\n"
     "                        meridional velocity, and buoyancy.\n"
     "  Nx                    Number of grid cells in the y-direction. Must be >0.\n"
     "  Nz                    Number of grid cells in the z-direction. Must be >0.\n"
     "  Lx                    Domain width. Must be >0.\n"
     "  Lz                    Domain height. Must be >0.\n"
     "  \n"
     "  cflFrac               CFL number. The time step dt will be chosen at each\n"
     "                        iteration such that dt=cfl*dt_max, where dt_max is the\n"
     "                        estimated maximum stable time step. Must be >0.\n"
     "  startTime             Time at which integration should start.\n"
     "                        Optional. If unset or negative then the run is initialized\n"
     "                        from the time of the pickup file, if restart is specified. If\n"
     "                        startTime is set and >=0  then it overrides the pickup time.\n"
     "  endTime               End time for the integration, in\n"
     "                        case the target residual is not achieved.\n"
     "  monitorFrequency      Frequency with which to write out iteration data. If\n"
     "                        set to zero, only the final data will be written.\n"
     "                        Otherwise must be >0.\n"
     "                        Optional, default is 0, must be >= 0.\n"
     "  restart               Set to 0 to initialize the simulation from t=0\n"
     "                        using input from the initFile parameter. Set to\n"
     "                        1 to pick up from a previous run using output files\n"
     "                        specified by startIdx as input files for this run.\n"
     "                        Optional - default is 0.\n"
     "  startIdx              Specifies the index of the output files to use as\n"
     "                        input files for this run. The output index counts\n"
     "                        the number of output files written in any given\n"
     "                        simulation. If the startTime parameter is not specified\n"
     "                        then the simulation start time is calculated as\n"
     "                        startIdx*monitorFrequency. Optional - default is 0.\n"
     "  checkConvergence      Set true to check whether the model has converged to a steady\n"
     "                        state after each time step. Convergence criteria for each\n"
     "                        are set via the targetResFile parameter. Default is 0\n"
     "                        (no convergence checking).\n"
     
     "  \n"
     "  f0                    Coriolis parameter.\n"
     "                        Optional, default is 10^{-4} rad/s, must be =/= 0.\n"
     "  rho0                  Reference density.\n"
     "                        Optional, default is 1000 kg/m^3, must be > 0.\n"
     "  Kconv0                Convective adjustment diffusivity.\n"
     "                        Optional, default is 10 m^2/s, must be > 0.\n"
     "  Hsml                  Depth of the surface mixed layer. Must be >=0.\n"
     "                        If equal to 0 then no SML is imposed\n"
     "  Hbbl                  Depth of the bottom boundary layer. Must be >=0.\n"
     "                        If equal to 0 then no BBL is imposed\n"
     "  r_bbl                 Drag coefficient in the bottom boundary layer \n"
     "  Ksml                  Max vertical diffusivity in surface mixed layer.\n"
     "                        Optional, default is 0.1 m^2/s, must be > 0.\n"
     "  Kbbl                  Max vertical diffusivity in bottom boundary layer.\n"
     "                        Optional, default is 0.1 m^2/s, must be > 0.\n"
     "  nu_h                  Horizontal (iso-sigma) viscosity. Default is 0.\n"
     "  nu_v                  Vertial viscosity. Actual viscosity is scaled with\n"
     "                        local vertical grid spacing. Specified value is applied\n"
     "                        to grid cells of thickness H/Lz. Default is 0.\n"
     "  \n"
     "  h_c                   Depth parameter controlling the range of depths over\n"
     "                        which the vertical coordinate the coordinate is\n"
     "                        approximately aligned with geopotentials.\n"
     "                        Must be >0. Default is 1e16 (pure sigma coordinate).\n"
     "  theta_s               Surface sigma coordinate stretching parameter,\n"
     "                        defined `meaningfully' between 0 and 10.\n"
     "                        Default is 0 (no surface stretching).\n"
     "  theta_b               Bottom sigma coordinate stretching parameter,\n"
     "                        defined `meaningfully' between 0 and 4.\n"
     "  \n"
     "  timeStepppingScheme   Selects time stepping scheme. See defs.h for\n"
     "                        identifiers. Default is TIMESTEPPING_RKTVD2.\n"
     "  advectionScheme       Selects advection scheme. See defs.h for\n"
     "                        identifiers. Default is ADVECTION_KT00.\n"
     "  momentumScheme        Selects momentum scheme. See defs.h for\n"
     "                        identifiers. Default is MOMENTUM_NONE.\n"
     "  pressureScheme        Selects the pressure scheme. See defs.h for \n"
     "                        identifiers. Default is PRESSURE_CUBIC. \n"
     "  bgcModel              Selects biogeochemical model type. See SWdefs.h for\n"
     "                        numerical identifiers. Default is BGC_NONE.\n"
     "  KT00_sigma            Kurganov-Tadmor (2000) minmod parameter.\n"
     "                        Should be >=1 and <=2. Default is 1.\n"
     "  \n"
     "  targetResFile         File containing an Ntracs x 1 vector of target L2\n"
     "                        errors between successive iterations. The\n"
     "                        calculation will stop when this condition is met\n"
     "                        for all tracers.\n"
     "  initFile              File containing an Ntracs x Nx x Nz array of initial\n"
     "                        tracer values over the whole domain.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  northTracerFile       File containing an Ntracs x Nx x Nz array of initial\n"
     "                        tracer values over the northern edge of the domain.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  southTracerFile       File containing an Ntracs x Nx x Nz array of initial\n"
     "                        tracer values over the southern edge of the domain.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  fluxFile              File containing an Ntracs x Nx array of tracer values\n"
     "                        across the surface of the domain.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  topogFile             File containing an Nx+1 x 1 array of ocean depths\n"
     "                        at grid cell corners. All elements must be > 0.\n"
     "                        Optional, default is Lz everywhere.\n"
     "  tauFile               File containing an Nx+1 x tlength array of wind stresses\n"
     "                        at grid cell corners.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  KgmFile               File containing an Nx+1 x Nz+1 array of GM eddy\n"
     "                        diffusivities at grid cell corners. All\n"
     "                        elements must be >=0.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  KisoFile              File containing an Nx+1 x Nz+1 array of isopycnal\n"
     "                        diffusivities at grid cell corners. All\n"
     "                        elements must be >=0.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  relaxTracerFile       File containing an Ntracers x Nx x Nz array of tracer\n"
     "                        relaxation target values over the whole domain.\n"
     "                        Optional, default is 0 everywhere.\n"
     "  relaxTimeFile         File containing an Ntracers x Nx x Nz array of tracer\n"
     "                        relaxation time scales over the whole domain. Negative \n"
     "                        values imply no relaxation. Zero values imply\n"
     "                        instantaneous relaxation. Optional, default is -1\n"
     "                        (no relaxation) everywhere.\n"
     " \n"
     " \n"
     ,
     progname
     );
}












/**
 * main
 *
 * Program entry point. Initialises all of the required memory,
 * reads in input parameters, performs time integration loop,
 * and finally cleans up.
 *
 */
int main (int argc, char ** argv)
{
    
    
    
    
    ///////////////////////////////////////////
    ///// BEGIN MAIN VARIABLE DEFINITIONS /////
    ///////////////////////////////////////////
    
    // Time parameters
    real tmin = -1;        // Start time
    real tmax = 0;          // End time
    real t = 0;             // Time counter
    real dt = 0;            // Time step
    real dt_s = 0;          // Data save time step
    real t_next = 0;        // Next save time
    real cflFrac = 1;       // CFL number
    real wc = 0;            // Time-interpolation weights
    real wn = 0;
    
    // Iteration parameters
    bool restart = false;   // Restart flag
    real n_s = 0;           // Save counter
    uint n0 = 0;            // Initial save index
    
    // Convergence parameters
    real * targetRes = 0;
    bool checkConvergence = false;
    bool targetReached = false;
    real * res = 0;
    uint nIters = 0;
    
    // Output directory
    char * outdir = NULL;
    
    // Looping variables
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;
    int jp = 0;
    int jz = 0;
    
    // Pointer to buoyancy matrix
    real ** buoy = NULL;
    
    // Work arrays for time derivatives - expressed as vectors rather
    // than matrices for compatibility with time integration methods
    real * phi_in_V = NULL; // Input to time-integration method
    real * phi_out_V = NULL; // Output from time-integration method
    real * phi_buf_V = NULL; // Buffer for time-integration method (needed for RK2 or RK3)
    real * phi_int_V = NULL; // Buffer for time-interpolation
    real *** phi_in = NULL;
    real *** phi_out = NULL;
    real *** phi_int = NULL;
    
    // Work arrays, variables for AB methods, and variables for
    // the shoelace method to calculate the area of a grid cell.
    real * dt_vars = NULL;
    real * dt_vars_1 = NULL;
    real * dt_vars_2 = NULL;
    real h = 0;
    real h1 = 0;
    real h2 = 0;
    real xj = 0;
    real xjp1 = 0;
    
    // BGC parameters
    real wsink = 0;
    
    // Stores data required for parsing input parameters
    paramdata params[NPARAMS];
    int paramcntr = 0;
    
    // Filename holders for input parameter arrays
    char targetResFile[MAX_PARAMETER_FILENAME_LENGTH];
    char initFile[MAX_PARAMETER_FILENAME_LENGTH];
    char northTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char southTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char fluxFile[MAX_PARAMETER_FILENAME_LENGTH];
    char topogFile[MAX_PARAMETER_FILENAME_LENGTH];
    char tauFile[MAX_PARAMETER_FILENAME_LENGTH];
    char bgcFile[MAX_PARAMETER_FILENAME_LENGTH];
    char bgcRatesFile[MAX_PARAMETER_FILENAME_LENGTH];
    char lpFile[MAX_PARAMETER_FILENAME_LENGTH];
    char lzFile[MAX_PARAMETER_FILENAME_LENGTH];
    char KgmFile[MAX_PARAMETER_FILENAME_LENGTH];
    char KisoFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
    
    // To store current time
    time_t now;
    FILE * tfile = NULL;
    
    /////////////////////////////////////////
    ///// END MAIN VARIABLE DEFINITIONS /////
    /////////////////////////////////////////
    
    
    
    ///////////////////////////////////////////
    ///// BEGIN DEFINING INPUT PARAMETERS /////
    ///////////////////////////////////////////
    
    // Domain dimensions and grid sizes (5)
    setParam(params,paramcntr++,"Ntracs","%u",&Ntracs,false);
    setParam(params,paramcntr++,"Nx","%u",&Nx,false);
    setParam(params,paramcntr++,"Nz","%u",&Nz,false);
    setParam(params,paramcntr++,"Lx","%lf",&Lx,false);
    setParam(params,paramcntr++,"Ly","%lf",&Ly,false);
    setParam(params,paramcntr++,"Lz","%lf",&Lz,false);
    
    // Time stepping/output parameters (7)
    setParam(params,paramcntr++,"cflFrac","%lf",&cflFrac,false);
    setParam(params,paramcntr++,"startTime","%lf",&tmin,true);
    setParam(params,paramcntr++,"endTime","%lf",&tmax,false);
    setParam(params,paramcntr++,"monitorFrequency","%lf",&dt_s,true);
    setParam(params,paramcntr++,"restart","%d",&restart,true);
    setParam(params,paramcntr++,"startIdx","%u",&n0,true);
    setParam(params,paramcntr++,"checkConvergence","%d",&checkConvergence,true);
    
    // Physical constants (11)
    setParam(params,paramcntr++,"rho0","%lf",&rho0,true);
    setParam(params,paramcntr++,"f0","%lf",&f0,true);
    setParam(params,paramcntr++,"Kconv","%lf",&Kconv0,true);
    setParam(params,paramcntr++,"Hsml","%lf",&Hsml,true);
    setParam(params,paramcntr++,"Hbbl","%lf",&Hbbl,true);
    setParam(params,paramcntr++,"r_bbl","%lf",&r_bbl,true);
    setParam(params,paramcntr++,"Kdia0","%lf",&Kdia0,true);
    setParam(params,paramcntr++,"Ksml","%lf",&Ksml,true);
    setParam(params,paramcntr++,"Kbbl","%lf",&Kbbl,true);
    setParam(params,paramcntr++,"nu_h","%lf",&nu_h,true);
    setParam(params,paramcntr++,"nu_v","%lf",&nu_v,true);
    
    // Sigma-coordinate parameters (3)
    setParam(params,paramcntr++,"h_c","%le",&h_c,true);
    setParam(params,paramcntr++,"theta_s","%lf",&theta_s,true);
    setParam(params,paramcntr++,"theta_b","%lf",&theta_b,true);
    
    // Scheme selectors (6)
    setParam(params,paramcntr++,"timeSteppingScheme","%u",&timeSteppingScheme,true);
    setParam(params,paramcntr++,"advectionScheme","%u",&advectionScheme,true);
    setParam(params,paramcntr++,"momentumScheme","%u",&momentumScheme,true);
    setParam(params,paramcntr++,"pressureScheme","%u",&pressureScheme,true);
    setParam(params,paramcntr++,"bgcModel","%u",&bgcModel,true);
    setParam(params,paramcntr++,"KT00_sigma","%lf",&KT00_sigma,true);
    
    // Biogeochemical parameter inputs (5)
    setParam(params,paramcntr++,"MP","%u",&MP,true); // number of phytoplankton size classes
    setParam(params,paramcntr++,"MZ","%u",&MZ,true); // number of zooplankton size classes
    setParam(params,paramcntr++,"npbgc","%u",&npbgc,true); // number of bgc parameters
    setParam(params,paramcntr++,"nallo","%u",&nallo,true); // number of bgc parameters
    setParam(params,paramcntr++,"idxAllo","%u",&idxAllo,true); // number of bgc parameters
    
    // Input file names (15)
    setParam(params,paramcntr++,"targetResFile","%s",&targetResFile,true);
    setParam(params,paramcntr++,"initFile","%s",initFile,false);
    setParam(params,paramcntr++,"northTracerFile","%s",northTracerFile,false);
    setParam(params,paramcntr++,"southTracerFile","%s",southTracerFile,false);
    setParam(params,paramcntr++,"fluxFile","%s",fluxFile,false);
    setParam(params,paramcntr++,"topogFile","%s",topogFile,true);
    setParam(params,paramcntr++,"tlength","%u",&tlength,false);
    setParam(params,paramcntr++,"tauFile","%s",tauFile,true);
    setParam(params,paramcntr++,"bgcFile","%s",bgcFile,true);
    setParam(params,paramcntr++,"bgcRatesFile","%s",bgcRatesFile,true);
    setParam(params,paramcntr++,"lpFile","%s",lpFile,true);
    setParam(params,paramcntr++,"lzFile","%s",lzFile,true);
    setParam(params,paramcntr++,"KgmFile","%s",KgmFile,true);
    setParam(params,paramcntr++,"KisoFile","%s",KisoFile,true);
    setParam(params,paramcntr++,"relaxTracerFile","%s",relaxTracerFile,true);
    setParam(params,paramcntr++,"relaxTimeFile","%s",relaxTimeFile,true);
    
    // Default file name parameters - zero-length strings
    targetResFile[0] = '\0';
    initFile[0] = '\0';
    northTracerFile[0] = '\0';
    southTracerFile[0] = '\0';
    fluxFile[0] = '\0';
    topogFile[0] = '\0';
    tauFile[0] = '\0';
    bgcFile[0] = '\0';
    bgcRatesFile[0] = '\0';
    lpFile[0] = '\0';
    lzFile[0] = '\0';
    KgmFile[0] = '\0';
    KisoFile[0] = '\0';
    relaxTracerFile[0] = '\0';
    relaxTimeFile[0] = '\0';
    
    /////////////////////////////////////////
    ///// END DEFINING INPUT PARAMETERS /////
    /////////////////////////////////////////
    
    
    
    //////////////////////////////////////////
    ///// BEGIN READING INPUT PARAMETERS /////
    //////////////////////////////////////////
    
    // First program argument always carries the program name
    progname = argv[0];
    
    // Check that file names have been specified
    if (argc < 3)
    {
        fprintf(stderr,"ERROR: Not enough input file names supplied\n");
        printUsage();
        return 0;
    }
    
    // Output directory
    outdir = argv[2];
    
    // Calculate elapsed time and write to a file
    tfile = fopen("time.txt","w");
    if (tfile != NULL)
    {
        time(&now);
        fprintf(tfile,"Program started at %s\n", ctime(&now));
        fflush(tfile);
    }
    
    // Attempt to read input parameter data. Errors will be printed
    // to stderr, and will result in 'false' being returned.
    if (!readParams(argv[1],params,NPARAMS,stderr))
    {
        printUsage();
        return 0;
    }
    
    // Check that required parameters take legal values
    if (  (Ntracs < 3) ||
        (Nx <= 0) ||
        (Nz <= 0) ||
        (Lx <= 0) ||
        (Lz <= 0) ||
        (dt_s < 0.0) ||
        (tmax <= 0.0) ||
        (cflFrac <= 0.0) ||
        (rho0 <= 0.0) ||
        (f0 == 0.0) ||
        (h_c <= 0.0) ||
        (strlen(targetResFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(initFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(northTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(southTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(fluxFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(topogFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(tauFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(KgmFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(KisoFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(relaxTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(relaxTimeFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(lpFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(lzFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(bgcRatesFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(bgcFile) > MAX_PARAMETER_FILENAME_LENGTH)           )
    {
        fprintf(stderr,"ERROR: Invalid input parameter values\n");
        printUsage();
        return 0;
    }
    
    // Calculate total work array size
    Ntot = Ntracs*Nx*Nz;
    
    // Calculate grid spacings
    ds = 1.0/Nz;
    dx = Lx/Nx;
    _dx = 1/dx;
    _2dx = 1/(2*dx);
    dxsq = dx*dx;
    
    // SML/BBL flags
    use_sml = Hsml > 0;
    use_bbl = Hbbl > 0;
    
    // Set the integration start time
    // Note: If dt_s is not specified then this will just be zero
    if (tmin < 0)
    {
        if (restart)
        {
            tmin = n0*dt_s;
            if (tmin >= tmax)
            {
                fprintf(stderr,"Integration start time (=startIdx*monitorFrequency) exceeds integration end time (endTime)");
                printUsage();
                return 0;
            }
        }
        else
        {
            tmin = 0;
        }
    }
    else
    {
        if (tmin >= tmax)
        {
            fprintf(stderr,"Integration start time (startTime) exceeds integration end time (endTime)");
            printUsage();
            return 0;
        }
    }
    t = tmin;
    
    // Current time and next save point
    if (dt_s == 0.0)
    {
        t_next = tmin + 2*tmax; // If no save interval is specified, only write final data
    }
    else
    {
        t_next = tmin + dt_s;
    }
    
    // If restarting then initialize output count accordingly
    if (restart)
    {
        n_s = n0;
    }
    else
    {
        n_s = 0;
    }
    
    ////////////////////////////////////////
    ///// END READING INPUT PARAMETERS /////
    ////////////////////////////////////////
    
    
    
    ///////////////////////////////////
    ///// BEGIN MEMORY ALLOCATION /////
    ///////////////////////////////////
    
    // Allocate tracer time stepping data  as 1D vectors, then reference
    // them using 3D matrices
    VECALLOC(phi_in_V,Ntot);
    VECALLOC(phi_out_V,Ntot);
    VECALLOC(phi_buf_V,Ntot);
    VECALLOC(phi_int_V,Ntot);
    vec2mat3(phi_in_V,&phi_in,Ntracs,Nx,Nz);
    vec2mat3(phi_out_V,&phi_out,Ntracs,Nx,Nz);
    vec2mat3(phi_int_V,&phi_int,Ntracs,Nx,Nz);
    
    // Allocate arrays for AB methods
    VECALLOC(dt_vars,Ntot);
    VECALLOC(dt_vars_1,Ntot);
    VECALLOC(dt_vars_2,Ntot);
    
    // Allocate parameter arrays
    VECALLOC(res,Ntracs);
    VECALLOC(targetRes,Ntracs);
    MATALLOC3(phi_init,Ntracs,Nx,Nz);
    MATALLOC3(phi_north,Ntracs,Nx,Nz);
    MATALLOC3(phi_south,Ntracs,Nx,Nz);
    MATALLOC(phi_flux,Ntracs,Nx);
    VECALLOC(hb_in,Nx+2);
    VECALLOC(tau,Nx+1);
    MATALLOC(Kgm_psi_ref,Nx+1,Nz+1);
    MATALLOC(Kiso_psi_ref,Nx+1,Nz+1);
    MATALLOC(Kdia_psi_ref,Nx+1,Nz+1);
    MATALLOC3(phi_relax,Ntracs,Nx,Nz);
    MATALLOC3(T_relax,Ntracs,Nx,Nz);
    
    // BGC parameter arrays
    MATALLOC(bgcRates,nallo,idxAllo);
    MATALLOC(theta_p,MP,MZ);
    MATALLOC(theta_z,MZ,MZ);
    MATALLOC(grazMat,MP,MZ);
    MATALLOC(hetMat,MP,MZ);
    VECALLOC(bgc_params,npbgc);
    VECALLOC(lpvec,MP);
    VECALLOC(lzvec,MZ);
    VECALLOC(UPvec,MP);
    VECALLOC(GPvec,MP);
    VECALLOC(GZvec,MZ);
    VECALLOC(PTvec,MZ);
    VECALLOC(ZTvec,MZ);
    VECALLOC(HMvec,MZ);
    VECALLOC(HPvec,MZ);
    VECALLOC(MPvec,MP);
    VECALLOC(MZvec,MZ);
    
    
    // Grid vectors
    VECALLOC(sigma_phi,Nz);
    VECALLOC(sigma_psi,Nz+1);
    MATALLOC(_dz_phi,Nx,Nz);
    MATALLOC(_2dz_phi,Nx,Nz);
    MATALLOC(_dzsq_psi,Nx+1,Nz+1);
    MATALLOC(_dz_w,Nx,Nz+1);
    MATALLOC(_dzsq_w,Nx,Nz+1);
    MATALLOC(dz_u,Nx+1,Nz);
    MATALLOC(_dz_u,Nx+1,Nz);
    MATALLOC(dA_psi,Nx,Nz);
    VECALLOC(hb_phi,Nx);
    VECALLOC(hb_psi,Nx+1);
    VECALLOC(sb_phi,Nx);
    VECALLOC(sb_psi,Nx+1);
    MATALLOC(ZZ_psi,Nx+1,Nz+1);
    MATALLOC(ZZ_phi,Nx,Nz);
    MATALLOC(ZZ_u,Nx+1,Nz);
    MATALLOC(ZZ_w,Nx,Nz+1);
    
    // Work arrays for 'tderiv' function
    VECALLOC(phi_wrk_V,Ntot);
    vec2mat3(phi_wrk_V,&phi_wrk,Ntracs,Nx,Nz);
    vec2mat3(phi_wrk_V,&dphi_dt_wrk,Ntracs,Nx,Nz); // N.B. dphi_dt_wrk will be reassigned when tderiv is first called
    
    MATALLOC(HHx,Nx+1,Nz);
    MATALLOC(HHz,Nx,Nz+1);
    MATALLOC(PPx,Nx+1,Nz);
    MATALLOC(PPz,Nx,Nz+1);
    MATALLOC(psi_r,Nx+1,Nz+1);
    MATALLOC(u_r,Nx+1,Nz);
    MATALLOC(w_r,Nx,Nz+1);
    MATALLOC(psi_m,Nx+1,Nz+1);
    MATALLOC(psi_e,Nx+1,Nz+1);
    MATALLOC(dphi_dx,Nx,Nz);
    MATALLOC(dphi_dz,Nx,Nz);
    MATALLOC(phi_xp,Nx+1,Nz);
    MATALLOC(phi_xm,Nx+1,Nz);
    MATALLOC(phi_zp,Nx,Nz+1);
    MATALLOC(phi_zm,Nx,Nz+1);
    MATALLOC(Sgm_psi,Nx+1,Nz+1);
    MATALLOC(Siso_u,Nx+1,Nz);
    MATALLOC(Siso_w,Nx,Nz+1);
    VECALLOC(impl_A,Nz);
    VECALLOC(impl_B,Nz);
    VECALLOC(impl_C,Nz);
    VECALLOC(impl_D,Nz);
    MATALLOC(wsink_wrk,Nx,Nz+1);
    MATALLOC(Kgm_psi,Nx+1,Nz+1);
    MATALLOC(Kgm_u,Nx+1,Nz);
    MATALLOC(Kgm_w,Nx,Nz+1);
    MATALLOC(Kiso_u,Nx+1,Nz);
    MATALLOC(Kiso_w,Nx,Nz+1);
    MATALLOC(Kdia_w,Nx,Nz+1);
    MATALLOC(BPa,Nx,Nz+1);
    MATALLOC(BPx,Nx+1,Nz);
    MATALLOC(BPy,Nx,Nz);
    MATALLOC(BBy,Nx,Nz);
    MATALLOC(Nbuoy,Nx,Nz);
    
    // Pressure calculation scheme
    MATALLOC(rhos,Nx,Nz);
    MATALLOC(drz,Nx,Nz+1);
    MATALLOC(drx,Nx+1,Nz);
    MATALLOC(dzz,Nx,Nz+1);
    MATALLOC(dzx,Nx+1,Nz);
    MATALLOC(hrx,Nx+1,Nz);
    MATALLOC(hrz,Nx,Nz+1);
    MATALLOC(hzx,Nx+1,Nz);
    MATALLOC(hzz,Nx,Nz+1);
    MATALLOC(P,Nx,Nz);
    MATALLOC(FX,Nx,Nz+1);
    MATALLOC(FC,Nx+1,Nz);
    MATALLOC(rho_north,Nx,Nz);
    MATALLOC(rho_south,Nx,Nz);
    
    
    // Boundary layer work arrays
    k_sml = malloc((Nx+1)*sizeof(uint));
    VECALLOC(wn_sml,Nx+1);
    VECALLOC(wp_sml,Nx+1);
    k_bbl = malloc((Nx+1)*sizeof(uint));
    VECALLOC(wn_bbl,Nx+1);
    VECALLOC(wp_bbl,Nx+1);
    MATALLOC(db_dx,Nx+1,Nz+1);
    MATALLOC(db_dz,Nx+1,Nz+1);
    MATALLOC(db_dx_lin,Nx+1,Nz+1);
    MATALLOC(db_dx_wrk,Nx+1,Nz+1);
    MATALLOC(db_dz_wrk,Nx+1,Nz+1);
    MATALLOC(bot_nflux,max_det,Nx);
    
    /////////////////////////////////
    ///// END MEMORY ALLOCATION /////
    /////////////////////////////////
    
    
    
    ////////////////////////////////////////
    ///// BEGIN READING PARAMETER DATA /////
    ////////////////////////////////////////
    
    // Read input matrices and vectors
    if (  ( (strlen(targetResFile) > 0)   &&   !readVector(targetResFile,targetRes,Ntracs,stderr) ) ||
        ( (strlen(topogFile) > 0)         &&   !readVector(topogFile,hb_in,Nx+2,stderr) ) ||
        ( (strlen(tauFile) > 0)           &&   !readVector(tauFile,tau,Nx+1,stderr) ) ||
        ( (strlen(bgcFile) > 0)           &&   !readVector(bgcFile,bgc_params,npbgc,stderr) ) ||
        ( (strlen(lpFile) > 0)            &&   !readVector(lpFile,lpvec,MP,stderr) ) ||
        ( (strlen(lzFile) > 0)            &&   !readVector(lzFile,lzvec,MZ,stderr) ) ||
        ( (strlen(bgcRatesFile) > 0)      &&   !readMatrix(bgcRatesFile,bgcRates,nallo,idxAllo,stderr) ) ||
        ( (strlen(KgmFile) > 0)           &&   !readMatrix(KgmFile,Kgm_psi_ref,Nx+1,Nz+1,stderr) ) ||
        ( (strlen(KisoFile) > 0)          &&   !readMatrix(KisoFile,Kiso_psi_ref,Nx+1,Nz+1,stderr) )  ||
        ( (strlen(relaxTracerFile) > 0)   &&   !readMatrix3(relaxTracerFile,phi_relax,Ntracs,Nx,Nz,stderr) )  ||
        ( (strlen(northTracerFile) > 0)   &&   !readMatrix3(northTracerFile,phi_north,Ntracs,Nx,Nz,stderr) )  ||
        ( (strlen(southTracerFile) > 0)   &&   !readMatrix3(southTracerFile,phi_south,Ntracs,Nx,Nz,stderr) )  ||
        ( (strlen(fluxFile) > 0)          &&   !readMatrix(fluxFile,phi_flux,Ntracs,Nx,stderr) )  ||
        ( (strlen(relaxTimeFile) > 0)     &&   !readMatrix3(relaxTimeFile,T_relax,Ntracs,Nx,Nz,stderr) )  )
    {
        printUsage();
        return 0;
    }
    
    // If a restart is specified then read from the output of a previous run
    if (restart)
    {
        // Check initial state
        for (i = 0; i < Ntracs; i ++)
        {
            // Overwrite initFile parameter as we won't be using it
            constructOutputName(outdir,OUTN_TRAC,i,n0,initFile);
            if (!readMatrix(initFile,phi_init[i],Nx,Nz,stderr))
            {
                fprintf(stderr,"Unable to read pickup files.\n");
                printUsage();
                return 0;
            }
        }
    }
    // Read from initialization files
    else
    {
        if ( (strlen(initFile) > 0)  &&  !readMatrix3(initFile,phi_init,Ntracs,Nx,Nz,stderr) )
        {
            fprintf(stderr,"Unable to read initialization files.\n");
            printUsage();
            return 0;
        }
    }
    
    // Default initial condition is zero tracer everywhere
    if (!restart && strlen(initFile)==0)
    {
        memset(*(*phi_init),0,(Ntot)*sizeof(real));
    }
    
    // Default initial condition is zero tracer everywhere
    if (!restart && strlen(northTracerFile)==0)
    {
        memset(*(*phi_north),0,(Ntot)*sizeof(real));
    }
    
    // Default initial condition is zero tracer everywhere
    if (!restart && strlen(southTracerFile)==0)
    {
        memset(*(*phi_south),0,(Ntot)*sizeof(real));
    }
    
    // Default initial condition is zero flux at the surface
    if (!restart && strlen(fluxFile)==0)
    {
        memset(*phi_flux,0,(Ntracs*Nx)*sizeof(real));
    }
    
    // Default ocean depth is constant everywhere
    if (strlen(topogFile)==0)
    {
        for (j = 0; j < Nx+2; j ++)
        {
            hb_in[j] = Lz;
        }
    }
    
    // Default wind stress is zero everywhere
    if (strlen(tauFile)==0)
    {
        memset(tau,0,(tlength*Nx+1)*sizeof(real));
    }
    
    // Default Kgm is zero everywhere
    if (strlen(KgmFile)==0)
    {
        memset(*Kgm_psi_ref,0,(Nx*Nz)*sizeof(real));
    }
    
    // Default Kiso is zero everywhere
    if (strlen(KisoFile)==0)
    {
        memset(*Kiso_psi_ref,0,(Nx*Nz)*sizeof(real));
    }
    
    // Default tracer relaxation target is zero everywhere
    if (strlen(relaxTracerFile)==0)
    {
        memset(*(*phi_relax),0,(Ntot)*sizeof(real));
    }
    
    // Default tracer relaxation is no relaxation everywhere
    if (strlen(relaxTimeFile)==0)
    {
        for (i = 0; i < Ntracs; i ++)
        {
            for (j = 0; j < Nx; j ++)
            {
                for (k = 0; k < Nz; k ++)
                {
                    T_relax[i][j][k] = -1;
                }
            }
        }
    }
    
    // Target residuals must be positive
    if (checkConvergence)
    {
        for (i = 0; i < Ntracs; i ++)
        {
            if (targetRes[i]<=0)
            {
                fprintf(stderr,"targetResFile may contain only values that are >0.");
                printUsage();
                return 0;
            }
        }
    }
    
#pragma parallel
    
    // Diffusivities must be positive
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            if (Kgm_psi_ref[j][k] < 0)
            {
                fprintf(stderr,"KgmFile may contain only values that are >0.");
                printUsage();
                return 0;
            }
            if (Kiso_psi_ref[j][k] < 0)
            {
                fprintf(stderr,"KisoFile may contain only values that are >0.");
                printUsage();
                return 0;
            }
            if (Kdia_psi_ref[j][k] < 0)
            {
                fprintf(stderr,"Kdia may contain only values that are >0.");
                printUsage();
                return 0;
            }
        }
    }
    
    // Topographic depths must be positive
    for (j = 0; j < Nx+2; j ++)
    {
        if (hb_in[j]<=0)
        {
            fprintf(stderr,"topogFile may contain only values that are >0.");
            printUsage();
            return 0;
        }
    }
    
    //////////////////////////////////////
    ///// END READING PARAMETER DATA /////
    //////////////////////////////////////
    
    
    
    
    ////////////////////////////
    ///// BEGIN GRID SETUP /////
    ////////////////////////////
    
    // Sigma coordinates
    for (k = 0; k < Nz; k ++)
    {
        sigma_phi[k] = -1 + (k+0.5)*ds;
    }
    for (k = 0; k < Nz+1; k ++)
    {
        sigma_psi[k] = -1 + k*ds;
    }
    
    // Water column depths on phi-points
    for (j = 0; j < Nx; j ++)
    {
        // Topographic depth on phi-gridpoints
        hb_phi[j] = hb_in[j+1]; // Note +1 offset because hb_in has length Nx+2
        
        // Topographic slope on phi-gridpoints
        sb_phi[j] =  - (hb_in[j+2]-hb_in[j]) * _2dx; // Note +1 offset because hb_in has length Nx+2
    }
    
    // Water column depths on psi-points
    for (j = 0; j < Nx+1; j ++)
    {
        // Topographic depth on psi-gridpoints
        hb_psi[j] = 0.5*(hb_in[j]+hb_in[j+1]);
        
        // Topographic slope on psi-gridpoints
        sb_psi[j] =  - (hb_in[j+1]-hb_in[j]) * _dx;
    }
    
    // True depths at phi-points
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            ZZ_phi[j][k] = stretch_ROMS(sigma_phi[k],h_c,theta_s,theta_b,hb_phi[j]);
        }
    }
    
    // True depths at psi-points
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            ZZ_psi[j][k] = stretch_ROMS(sigma_psi[k],h_c,theta_s,theta_b,hb_psi[j]);
        }
    }
    
    // True depths at u-points
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            ZZ_u[j][k] = stretch_ROMS(sigma_phi[k],h_c,theta_s,theta_b,hb_psi[j]);
        }
    }
    
    // True depths at w-points
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            ZZ_w[j][k] = stretch_ROMS(sigma_psi[k],h_c,theta_s,theta_b,hb_phi[j]);
        }
    }
    
    // Grid spacings on phi-points
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            _dz_phi[j][k] = 1/(ZZ_w[j][k+1]-ZZ_w[j][k]);
            
            if ((k == 0) || (k == Nz-1))
            {
                _2dz_phi[j][k] = 0.5*_dz_phi[j][k];       // Should never be used
            }
            else
            {
                _2dz_phi[j][k] = 1/(ZZ_phi[j][k+1] - ZZ_phi[j][k-1]);
            }
        }
    }
    
    // Grid spacings on u-points
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            dz_u[j][k] = (ZZ_psi[j][k+1]-ZZ_psi[j][k]);
            _dz_u[j][k] = 1/dz_u[j][k];
        }
    }
    
    // Grid spacings on w-points
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            if (k == 0)
            {
                _dz_w[j][k] = 1/(2*(ZZ_phi[j][k]-ZZ_w[j][k]));       // Should never be used
            }
            else if (k == Nz)
            {
                _dz_w[j][k] = 1/(2*(ZZ_w[j][k]-ZZ_phi[j][k-1]));       // Should never be used
            }
            else
            {
                _dz_w[j][k] = 1/(ZZ_phi[j][k]-ZZ_phi[j][k-1]);
            }
            
            _dzsq_w[j][k] = SQUARE(_dz_w[j][k]);
        }
    }
    
    // Grid spacings on psi-points
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            if (k == 0)
            {
                _dzsq_psi[j][k] = SQUARE( 1/(2*(ZZ_u[j][k]-ZZ_psi[j][k])) );       // Should never be used
            }
            else if (k == Nz)
            {
                _dzsq_psi[j][k] = SQUARE( 1/(2*(ZZ_psi[j][k]-ZZ_u[j][k-1])) );       // Should never be used
            }
            else
            {
                _dzsq_psi[j][k] = SQUARE( 1/(ZZ_u[j][k]-ZZ_u[j][k-1]) );
            }
        }
    }
    
    // Grid areas on the psi-points (with verticies on the phi points)
    // Calculate the area of the square using the shoelace method
    for (j = 1; j < Nx; j++)
    {
        for (k = 1; k < Nz; k++)
        {
            dA_psi[j][k] = 0.5*dx*fabs( ZZ_phi[j-1][k-1] + ZZ_phi[j][k-1] - ZZ_phi[j][k] - ZZ_phi[j-1][k]);
        }
    }
    
    // Vertical Calculation
    // Compute elementary differences
    for (j = 0; j < Nx; j++)
    {
        for (k = 0; k < Nz-1; k++)
        {
            // starts with +1/2 (0) and ends with Nz-3/2 (Nz-2)
            dzz[j][k] = ZZ_phi[j][k+1] - ZZ_phi[j][k];
        }
    }

    
    // Calculate the harmonic averages:
    for (j = 0; j < Nx; j++)
    {
        for (k = 1; k < Nz-1; k++)
        {
            cff = dzz[j][k-1]*dzz[j][k];
            if (cff > minVal)
            {
                hzz[j][k] = 2*cff/(dzz[j][k-1] + dzz[j][k]);
            }
            else
            {
                hzz[j][k] = 0;
            }
        }

        // Set the boundary conditions
        hzz[j][0] = 1.5*(ZZ_phi[j][1] - ZZ_phi[j][0]) - 0.5*hzz[j][1];
        hzz[j][Nz-1] = 1.5*(ZZ_phi[j][Nz-1] - ZZ_phi[j][Nz-2]) - 0.5*hzz[j][Nz-2];
    }
    
    // Horizontal Calculations
    // Elementary Differences
    for (j = 0; j < Nx-1; j++)
    {
        for (k = 0; k < Nz; k++)
        {
            dzx[j][k] = ZZ_phi[j+1][k] - ZZ_phi[j][k];
        }
    }

    // Calculate the harmonic average:
    for (k = 0; k < Nz; k++)
    {
        for (j = 1; j < Nx-1; j++)
        {
            cff = dzx[j-1][k]*dzx[j][k];
            if (cff > minVal)
            {
                hzx[j][k] = 2*cff/(dzx[j-1][k] + dzx[j][k]);
            }
            else
            {
                hzx[j][k] = 0;
            }
        }

        hzx[0][k] = 1.5*(ZZ_phi[1][k] - ZZ_phi[0][k]) - 0.5*hzx[1][k];
        hzx[Nx-1][k] = 1.5*(ZZ_phi[Nx-1][k] - ZZ_phi[Nx-2][k]) - 0.5*hzx[Nx-2][k];
    }
    
    
    
    
    // Calculate the meridional pressure gradients
    if (Ly < 0 || Ly == 0)
    {
        _Ly = 0; // check to see if the domain is infinite or undefined
        for (j = 0; j < Nx; j++)
        {
            for (k = 0; k < Nz-1; k++)
            {
                BPy[j][k] = 0;
            }
        }
    }
    else
    {
        _Ly = 1/Ly; // set the gradients to zero since they haven't been extensively tested
        
        // Calculate the density gradient of the northern and southern
        // edges of the domin
        for (j = 0; j < Nx; j++)
        {
            for (k = 0; k < Nz; k++)
            {
                rho_north[j][k] = (1-alpha*(phi_north[idx_buoy][j][k] - tref));
                rho_south[j][k] = (1-alpha*(phi_south[idx_buoy][j][k] - tref));
            }
        }
        
        // Calculate the pressure at the (Nz-1)-th grid point
        for (j = 0; j < Nx; j++)
        {
            cff = ZZ_psi[j][Nz] - ZZ_u[j][Nz-1];
            BPy[j][Nz-1] = _Ly*grav*cff*( rho_north[j][Nz-1]
                                     + 0.5 * cff * ( rho_north[j][Nz-1] - rho_north[j][Nz-2] )/( ZZ_u[j][Nz-1] - ZZ_u[j][Nz-2] )
                                     -  (rho_south[j][Nz-1]
                                     + 0.5 * cff * ( rho_south[j][Nz-1] - rho_south[j][Nz-2] )/( ZZ_u[j][Nz-1] - ZZ_u[j][Nz-2] ) )  );
            BPy[j][Nz-1] = _Ly*(1/rho0)*(rho_north[j][Nz-1] - rho_south[j][Nz-1]);
        }

        // Calculate the pressure by vertically integrating
        for (j = 0; j < Nx; j++)
        {
            for (k = Nz-1; k >= 0; k --)
            {
               BPy[j][k] += BPy[j][k+1] + _Ly*(grav/rho0)*( 0.5*(rho_north[j][k+1] + rho_north[j][k])  - 0.5*(rho_south[j][k+1] + rho_south[j][k]) )*(ZZ_u[j][k+1] - ZZ_u[j][k]);
            }
        }
        
        for (k = 0; k < Nz; k++)
        {
            BPy[0][k] = 0; // set the pressure gradient on the western wall to zero
        }
        
        
    }
        
    
    // Calculate diapyncal diffusivity file from ksml and kbbl
    for (j = 0; j < Nx; j++)
    {
        for (k = 0; k < Nz; k++)
        {
            Kdia_psi_ref[j][k] = Kdia0;
            
            // Calculate the diapyncal profile in the surface boundary layer
            if (ZZ_psi[j][k] > -Hsml)
            {
                Kdia_psi_ref[j][k] = Kdia0 + Ksml*(27/4)*( (-ZZ_psi[j][k])/Hsml * pow(1 + ZZ_psi[j][k]/Hsml,2) );
            }
            
            // Calculate the diapyncal profile in the bottom boundary layer
            if (ZZ_psi[j][k] < -hb_psi[j] + Hbbl)
            {
                Kdia_psi_ref[j][k] = Kdia0 + Kbbl*(27/4)*(  (ZZ_psi[j][k] + hb_psi[j])/Hbbl * pow(1 - (ZZ_psi[j][k] + hb_psi[j])/Hbbl ,2)  );
            }
        }
    }
    

    // Calculate the predator-prey interaction matrices
    // Z on P grazing
    int idxPredPrey = bgc_params[17];  // read in index for predator prey length ratio in bgcRates vector
    real dxg = bgc_params[10];         // read in grazing profile witdh
    for (jp = 0; jp < MP; jp++)
    {
        for (jz = 0; jz < MZ; jz++)
        {
            theta_p[jp][jz] = exp( -1*pow( (log10(lpvec[jp]) - log10(bgcRates[jz][idxPredPrey]))/dxg,2) );
        }
    }
    // Z on Z grazing
    for (jp = 0; jp < MZ; jp++) // use the same indexing variable for zooplankton so there's no variable loading
    {
        for (jz = 0; jz < MZ; jz++)
        {
            theta_z[jp][jz] = exp( -1*pow( (log10(lzvec[jp]) - log10(bgcRates[jz][idxPredPrey]))/dxg,2) );
        }
    }
    
    
    
    // Write model grid to output files for post-simulation analysis
    if (!writeModelGrid(ZZ_phi,ZZ_psi,ZZ_u,ZZ_w,outdir))
    {
        fprintf(stderr,"Unable to write model grid");
        printUsage();
        return 0;
    }
    
    // For calculating buoyancy gradients at SML base
    if (use_sml)
    {
        for (j = 0; j < Nx+1; j ++)
        {
            k_sml[j] = 0; // Default - bottom of the ocean
            for (k = Nz-1; k >= 0; k --)
            {
                if (ZZ_psi[j][k] < -Hsml)
                {
                    k_sml[j] = k;
                    break;
                }
            }
            // If the SML is shallower than the first vertical grid point,
            // set the index and weights to calculate db/dz at the first interior
            // psi-gridpoint
            if (k_sml[j] == Nz-1)
            {
                wp_sml[j] = 1;
                wn_sml[j] = 0;
            }
            // If the SML extends to the very bottom of the ocean, set the linear
            // interpolation weights to calculate db/dz at the first interior
            // psi-gridpoint
            else if (k_sml[j] == 0)
            {
                wp_sml[j] = 0;
                wn_sml[j] = 1;
            }
            // Default linear interpolation weights for db/dz at the SML base
            else
            {
                wn_sml[j] = (-Hsml-ZZ_psi[j][k_sml[j]]) / (ZZ_psi[j][k_sml[j]+1]-ZZ_psi[j][k_sml[j]]);
                wp_sml[j] = 1.0 - wn_sml[j];
            }
        }
    }
    
    // For calculating buoyancy gradients at BBL top
    if (use_bbl)
    {
        for (j = 0; j < Nx+1; j ++)
        {
            k_bbl[j] = Nz; // Default - top of the ocean
            for (k = 1; k <= Nz; k ++)
            {
                if (ZZ_psi[j][k] > -hb_psi[j] + Hbbl)
                {
                    k_bbl[j] = k;
                    break;
                }
            }
            // If the BBL occupies the full water column depth
            // set the index and weights to calculate db/dz at the first interior
            // psi-gridpoint
            if (k_bbl[j] == Nz)
            {
                wp_bbl[j] = 1;
                wn_bbl[j] = 0;
            }
            // If the BBL is shallower than the bottom grid cell, set the linear
            // interpolation weights to calculate db/dz at the first interior
            // psi-gridpoint
            else if (k_bbl[j] == 1)
            {
                wp_bbl[j] = 0;
                wn_bbl[j] = 1;
            }
            // Default linear interpolation weights for db/dz at the BBL top
            else
            {
                wp_bbl[j] = (ZZ_psi[j][k_bbl[j]] - (-hb_psi[j]+Hbbl)) / (ZZ_psi[j][k_bbl[j]] - ZZ_psi[j][k_bbl[j]-1]);
                wn_bbl[j] = 1.0 - wp_bbl[j];
            }
        }
    }

    
    
    
    //////////////////////////
    ///// END GRID SETUP /////
    //////////////////////////
    
    
    
    ////////////////////////////////////
    ///// BEGIN INITIAL CONDITIONS /////
    ////////////////////////////////////
    
#pragma parallel
    
    // Set the initial tracer fields
    for (i = 0; i < Ntracs; i ++)
    {
        for (j = 0; j < Nx; j ++)
        {
            for (k = 0; k < Nz; k ++)
            {
                phi_in[i][j][k] = phi_init[i][j][k];
                if (T_relax[i][j][k] == 0.0)
                {
                    phi_in[i][j][k] = phi_relax[i][j][k];
                }
            }
        }
    }
    
    
    // Initialize convergence residuals
    if (checkConvergence)
    {
        for (i = 0; i < Ntracs; i ++)
        {
            res[i] = 10*targetRes[i];
        }
    }
    
    // Write out the initial data if a save interval is specified
    if (dt_s > 0.0)
    {
        if (!writeModelState(t,0,phi_in,outdir))
        {
            fprintf(stderr,"Unable to write model initial state");
            printUsage();
            return 0;
        }
    }
    
    // Write initial time to time file
    if (tfile != NULL)
    {
        fprintf(tfile,"%e ",t);
        fflush(tfile);
    }
    
    
    //////////////////////////////////
    ///// END INITIAL CONDITIONS /////
    //////////////////////////////////
    
    
    
    
    ///////////////////////////////////////
    ///// BEGIN TIME INTEGRATION LOOP /////
    ///////////////////////////////////////
    
    clock_t tic = clock();
    
    fprintf(stderr,"There are %d tracers in this model \n",Ntracs);
    
    
    // Numerical time-integration loop - keep going while the residual exceeds
    // the target and the current time does not exceed the max time
    while (!targetReached && (t < tmax))
    {

        // Step 1: Perform a single numerical time-step for all physically explicit terms in the equations
        switch (timeSteppingScheme)
        {
            case TIMESTEPPING_AB1:
            {
                dt = ab1(&t,phi_in_V,phi_out_V,dt_vars,cflFrac,Ntot,&tderiv);
                break;
            }
            case TIMESTEPPING_AB2:
            {
                // calculate the first few timesteps with lower order schemes
                if (nIters == 0)
                {
                    
                    dt = ab1(&t,phi_in_V,phi_out_V,dt_vars,cflFrac,Ntot,&tderiv);
                    // Save h1 for the next time step.
                    h1 = dt;
                    
                    
                    // copy over data for the next timestep
                    memcpy(dt_vars_1,dt_vars,Ntot*sizeof(real));
                }
                else
                {
                    dt = ab2(&t,phi_in_V,phi_out_V,dt_vars,dt_vars_1,cflFrac,h1,Ntot,&tderiv);
                    // Save h1 for the next time step.
                    h1 = dt;
                }
                break;
            }
            case TIMESTEPPING_AB3:
            {
                
                // calculate the first few timesteps with lower order schemes
                if (nIters == 0)
                {
                    dt = ab1(&t,phi_in_V,phi_out_V,dt_vars,cflFrac,Ntot,&tderiv);
                    // Save h1 for the next time step.
                    h1 = dt;
                    // copy over data for the next timestep
                    memcpy(dt_vars_1,dt_vars,Ntot*sizeof(real));
                    
                    // save this data for the n == 3 time step
                    memcpy(dt_vars_2,dt_vars,Ntot*sizeof(real));
                    h2 = dt; // time step
                }
                else if (nIters == 1)
                {
                    dt = ab2(&t,phi_in_V,phi_out_V,dt_vars,dt_vars_1,cflFrac,h1,Ntot,&tderiv);
                    // Save h1 for the next time step.
                    h1 = dt;
                }
                else  // Once we have the other two iterations, we can use the third order step
                {
                    dt = ab3(&t,phi_in_V,phi_out_V,dt_vars,dt_vars_1,dt_vars_2,cflFrac,h1,h2,Ntot,&tderiv);
                    // save for next time step
                    h2 = h1;
                    h1 = dt;
                }
                
                break;
            }
            default:
            {
                fprintf(stderr,"ERROR: Unknown time-integration method\n");
                break;
            }
        }
        
        
        

        // Step 2: Add implicit vertical diffusion and remineralization
        do_impl_diff(t,dt,phi_out);
        
        
        // Step 3 (optional): apply zonal barotropic pressure gradient correction
        do_pressure_correct(phi_out);
        
        
#pragma parallel
        
        // Step 4: Enforce zero tendency where relaxation time is zero
        for (i = 0; i < Ntracs; i ++)
        {
            for (j = 0; j < Nx; j ++)
            {
                for (k = 0; k < Nz; k ++)
                {
                    if (T_relax[i][j][k] == 0.0)
                    {
                        phi_out[i][j][k] = phi_relax[i][j][k];
                    }
                }
            }
        }
        
        
        
        
#pragma parallel
        
        // Step 5 (optional): Calculate the residuals as an L2-norm between adjacent iterations
        if (checkConvergence)
        {
            targetReached = true;
            for (i = 0; i < Ntracs; i ++)
            {
                res[i] = 0;
                for (j = 0; j < Nx; j ++)
                {
                    for (k = 0; k < Nz; k ++)
                    {
                        res[i] += SQUARE((phi_in[i][j][k]-phi_out[i][j][k])/dt);
                    }
                }
                res[i] = sqrt(res[i]/(Nx*Nz));
                
                // Check whether we've reached ALL target residuals for the various tracers
                if (res[i] > targetRes[i])
                {
                    targetReached = false;
                }
                
                // NaN residual clearly indicates a problem
                if (isnan(res[i]))
                {
                    fprintf(stderr,"ERROR: Computation blew up: residual==NaN\n");
                    printUsage();
                    return 0;
                }
            }
        }
        
        
        
        // Step 6: If the time step has taken us past a save point (or multiple
        // save points), interpolate and write out the data
        while ((dt_s > 0) && (t >= t_next))
        {
            wc = (t-t_next) / dt;
            wn = 1 - wc;
            
#pragma parallel
            
            // Interpolate values at t_next into the phi_int matrix
            for (i = 0; i < Ntracs; i ++)
            {
                for (j = 0; j < Nx; j ++)
                {
                    for (k = 0; k < Nz; k ++)
                    {
                        phi_int[i][j][k] = wc*phi_in[i][j][k] + wn*phi_out[i][j][k];
                    }
                }
            }
            
            // Write out the interpolated data from the phi_int matrix
            n_s += 1;
            if (!writeModelState(t_next,n_s,phi_int,outdir)) // Recall phi_int is a pointer to phi_int_V
            {
                fprintf(stderr,"Unable to write model state");
                printUsage();
                return 0;
            }
            
            // Write current time to time file
            if (tfile != NULL)
            {
                fprintf(tfile,"%e ",t_next);
                fflush(tfile);
            }
            
            // Update the next save time
            t_next += dt_s;
            
            // Printing out the buoyancy residual can be quite a useful way of
            // keeping track of the computation's progress
            printf("%e, time = %f \n",res[idx_buoy],t/day);
            fflush(stdout);
        }
        
        
        // this write the model state after a single timestep, and stops the calculation
        if (debug == true)
        {
            if (!writeModelState(t,n_s,phi_in,outdir))
            {
                fprintf(stderr,"Unable to write model state");
                printUsage();
                return 0;
            }
            
            targetReached = true;
        }
        
        
        if (t > tmax || dt > 60*60*24)
        {
            fprintf(stderr,"ERROR: maximum time exceeded OR cfl_dt > 1 day \n");
            printUsage();
            return 0;
        }
        
        // Copy the next iteration from phi_out back to phi_in,
        // ready for the next time step
        memcpy(phi_in_V,phi_out_V,Ntot*sizeof(real));
        
        // Increment iteration count
        nIters += 1;
    }
    
    ////////////////////////////////
    ///// END TIME INTEGRATION /////
    ////////////////////////////////
    
    
    clock_t toc = clock();
    fprintf(stderr,"Elapsed: %f seconds \n", (double)(toc-tic)/CLOCKS_PER_SEC);
    
    
    /////////////////////////
    ///// BEGIN CLEANUP /////
    /////////////////////////
    
    // Write out the final model state
    if (dt_s == 0)
    {
        n_s += 1;
        if (!writeModelState(t,n_s,phi_in,outdir))
        {
            fprintf(stderr,"Unable to write model final state");
            printUsage();
            return 0;
        }
    }
    
    // Print completion time in time file
    if (tfile != NULL)
    {
        time(&now);
        fprintf(tfile,"\nProgram completed at %s\n", ctime(&now));
        for (i = 0; i < Ntracs; i ++)
        {
            fprintf(tfile,"\nres %d = %e, target res %d = %e\n", i, res[i], i, targetRes[i]);
        }
        fclose(tfile);
    }
    
    ///////////////////////
    ///// END CLEANUP /////
    ///////////////////////
    
    
    
    
    return 0;
}

