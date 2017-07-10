/**
 * mamebus.c
 *
 * Andrew Stewart
 *
 * Core code file for the Meridionally-Averaged Model of Eastern Boundary Upwelling Systems (MAMEBUS).
 * Integrates residual-mean buoyancy and tracer equations in an Eastern Boundary Current-like domain.
 *
 */
#include <time.h>

#include "rktvd.h"
#include "defs.h"



// Total number of input parameters - must match the number of parameters defined in main()
#define NPARAMS 25




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
real Lz = 0;
real Kconv0 = 10;
real rho0 = 1e3;
real f0 = 1e-4;
bool use_sml = true;
bool use_bbl = true;
real Hsml = 50.0;
real Hbbl = 50.0;
real tot_depth = 3000.0;  // Depth of the ocean

// Parameter arrays
real *** phi_init = NULL;     // Initial condition
real * hb_in = NULL;          // Ocean depth
real * tau = NULL;            // Surface wind stress
real ** Kgm_psi_ref = NULL;   // Reference GM diffusivity
real ** Kiso_psi_ref = NULL;  // Reference isopycnal diffusivity
real ** Kdia_psi_ref = NULL;  // Reference diapycnal diffusivity
real *** phi_relax = NULL;    // Tracer relaxation values
real *** T_relax = NULL;      // Tracer relaxation time scale

// Numerical parameters
real sigma = 1.4;         // Kurganov-Tadmor minmod-limiting parameter
bool limSlopes = false;    // Use Cox slope limiting
real Smax = 0.1;        // Max isopycnal slope
const int idx_buoy = 0;   // Index of buoyancy variable in list of tracers

// Grid parameters
real h_c = 1e16;
real theta_s = 0;
real theta_b = 0;

// Grid spacings
real ds = 0;
real dx = 0;
real _dx = 0;
real _2dx = 0;
real dxsq = 0;
real * sigma_phi = NULL;  // Stretched coordinate grids
real * sigma_psi = NULL;
real ** _dz_phi = NULL;   // Grid spacings
real ** _2dz_phi = NULL;
real ** _dzsq_psi = NULL;
real ** _dz_w = NULL;
real ** _dzsq_w = NULL;
real ** dz_u = NULL;
real ** _dz_u = NULL;
real * hb_phi = NULL;     // Ocean depths and bottom slope
real * hb_psi = NULL;
real * sb_phi = NULL;
real * sb_psi = NULL;
real ** ZZ_phi = NULL;    // Depth grids
real ** ZZ_psi = NULL;
real ** ZZ_u = NULL;
real ** ZZ_w = NULL;

// Name of the program (for error messages)
char * progname = NULL;

// Time-stepping method
const uint method_t = METHOD_RKTVD2;

// Spatial discretisation
const uint method_s = METHOD_KT;

// Output filenames
static const char OUTN_ZZ_PHI[] = "ZZ_PHI";
static const char OUTN_ZZ_PSI[] = "ZZ_PSI";
static const char OUTN_ZZ_U[] = "ZZ_U";
static const char OUTN_ZZ_W[] = "ZZ_W";
static const char OUTN_PSIM[] = "PSIM";
static const char OUTN_PSIE[] = "PSIE";
static const char OUTN_PSIR[] = "PSIR";
static const char OUTN_TRAC[] = "TRAC";

// Work arrays for calculating time derivatives
real * phi_wrk_V = NULL;      // Current iteration data - work array
real *** phi_wrk = NULL;
real *** dphi_dt_wrk = NULL;  // Time derivative of current iteration - work array
real ** HHx = NULL;           // Hyperbolic tracer fluxes
real ** HHz = NULL;
real ** PPx = NULL;           // Parabolic tracer fluxes
real ** PPz = NULL;
real ** FFx = NULL;           // Diabatic surface/bottom boundary layer fluxes
real ** FFz = NULL;
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
real ** Kgm_psi = NULL;       // GM diffusivity
real ** Kgm_u = NULL;
real ** Kgm_w = NULL;
real ** Kiso_u = NULL;        // Isopycnal diffusivity
real ** Kiso_w = NULL;
real ** Kdia_w = NULL;        // Diapycnal diffusivity

// Boundary layer work arrays
uint * k_sml = NULL;
real * wn_sml = NULL;
real * wp_sml = NULL;
uint * k_bbl = NULL;
real * wn_bbl = NULL;
real * wp_bbl = NULL;
real * db_dx = NULL;
real * db_dz = NULL;

////////////////////////////////
///// END GLOBAL VARIABLES /////
////////////////////////////////





// TODO ened to add surface bottom BLs
// TODO need to add meridional pressure gradient term here
// TODO solve diffusion equation to calculate u and v through water column


/**
 * calcPsim
 *
 * Calculates the mean overturning streamfunction from the mean momentum equation.
 * The result is stored in the psi_m matrix.
 *
 */
void calcPsim (const real t, real ** buoy, real ** psi_m)
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
    
#pragma parallel
    
    // Current implementation just keeps Kgm constant
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            Kgm_psi[j][k] = Kgm_psi_ref[j][k];
        }
    }
    
#pragma parallel
    
    // Diffusivity on u-gridpoints
    for (j = 0; j < Nx+1; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            Kgm_u[j][k] = 0.5*(Kgm_psi_ref[j][k+1] + Kgm_psi_ref[j][k]);
        }
    }
    
#pragma parallel
    
    // Diffusivity on w-gridpoints
    for (j = 0; j < Nx; j ++)
    {
        for (k = 0; k < Nz+1; k ++)
        {
            Kgm_w[j][k] = 0.5*(Kgm_psi_ref[j+1][k] + Kgm_psi_ref[j][k]);
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
















// TODO upgrade these to Ferrari et al form

/**
 * surfStructFun
 *
 * Surface structure function for the modified Ferrari et al. (2008) boundary layer parameterization. Aligns the slope used
 * to calculate the eddy streamfunction and symmetric diffusion tensor with the ocean surface as the surface is approached.
 *
 * z is the current vertical position
 * h_sml is the surface mixed layer thickness
 * _lambda is the reciprocal of the vertical eddy lengthscale at the SML base, and defines the vertical derivative of G at the SML base
 *
 */
real surfStructFun (real z, real h_sml, real _lambda)
{
    real G = 1.0;
    
    if (z > -h_sml)
    {
        //G = -z/h_sml;
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
        G = (z+h_b)/h_bbl;
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




















// TODO need to handle the case in which -Hsml < -h_b+Hbbl, i.e. in which BLs overlap

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
    real _lambda_sml, _lambda_bbl;
    
#pragma parallel
    
    // Calculate the effective isopycnal slope S_e everywhere
    for (j = 1; j < Nx; j ++)
    {
        // Calculate horizontal and vertical buoyancy gradients in (x,z) space
        for (k = 1; k < Nz; k ++)
        {
            db_dz[k] = 0.5 * ( (buoy[j][k]-buoy[j][k-1])*_dz_w[j][k] + (buoy[j-1][k]-buoy[j-1][k-1])*_dz_w[j-1][k] );
            db_dx[k] = 0.5 * ( (buoy[j][k]-buoy[j-1][k])*_dx + (buoy[j][k-1]-buoy[j-1][k-1])*_dx );
            db_dx[k] -= (ZZ_w[j][k]-ZZ_w[j-1][k])*_dx * db_dz[k];
        }
        db_dz[0] = 0; // N.B. THESE SHOULD NEVER BE USED
        db_dx[0] = 0; // Here we just set then so that they are defined
        db_dz[Nz] = 0;
        db_dx[Nz] = 0;
        
        // Calculate vertical gradients at SML base and BBL top
        if (use_sml)
        {
            // Buoyancy gradient at SML base
            db_dz_sml = db_dz[k_sml[j]]*wp_sml[j] + db_dz[k_sml[j]+1]*wn_sml[j];
            
            // Calculate reciprocal of vertical eddy length scale at SML base
            d2b_dz2 = (db_dz[k_sml[j]+1]-db_dz[k_sml[j]]) / (ZZ_psi[j][k_sml[j]+1]-ZZ_psi[j][k_sml[j]]);
            _lambda_sml = - d2b_dz2 / db_dz_sml;
        }
        if (use_bbl)
        {
            // Buoyancy gradient at BBL top
            db_dz_bbl = db_dz[k_bbl[j]]*wn_bbl[j] + db_dz[k_bbl[j]-1]*wp_bbl[j];
            
            // TODO
            _lambda_bbl = 0;
        }
        
        // Construct effective isopycnal slope
        for (k = 1; k < Nz; k ++)
        {
            z = ZZ_psi[j][k];
            
            if (use_sml && (z > -Hsml))
            {
                // TODO calculate lambda
                G_sml = surfStructFun(z,Hsml,_lambda_sml);
                Sgm_psi[j][k] = - G_sml * db_dx[k] / db_dz_sml;
            }
            else if (use_bbl && (z < -hb_psi[j] + Hbbl))
            {
                G_bbl = botStructFun(z,Hbbl,hb_psi[j],_lambda_bbl);
                Sgm_psi[j][k] = - G_bbl * db_dx[k] / db_dz_sml;
            }
            else
            {
                // Calculate slope on psi-gridpoints
                Sgm_psi[j][k] = - db_dx[k] / db_dz[k];
                
                // TODO replace with exponential taper - see Griffies (2004)
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
        
    }
    
}








/*
 
 WORKING HERE
 
 */


/**
 * tderiv_bgc
 *
 * Calculates the time tendency of any tracer phi due to biogeochemical processes.
 *
 */
void tderiv_bgc (const real t, real *** phi, real *** dphi_dt)
{
    int i,j,k;

    // Parameters
    real irr_0 = 340;                   // W/m^2 (Can eventually include a seasonal amplitude)
    real irr_scaleheight = 30;          // m
    real a = 0.6;                       // 1/d
    real b = 1.066;                     // From Sarmiento & Gruber (2006)
    real c = 1;                         // deg C
    real f = 0.5;                       // fraction of exported material
    real monod = 0;                     // nutrient limitation term
    real alpha = 0.025;                 // d^-1 (W/m^2)^-1
    
    // Variables
    real temp_flux                      // holder for remin value
    real flux                           // flux of remineralized nutrient
    real t_uptake                       // temperature dependent uptake rate
    real remin                          // remineralization
    real remin_in_box                   // value of remineralization in box
    real irradiance                     // amount of light that penetrates from the surface
    real lk                             // light saturation constant
    real l_uptake                       // light dependent uptake rate
    real uptake                         // uptake in each grid box
    real T                              // Placeholder for Temperature
    real N                              // Placeholder for Nitrate
    real scale_height                   // scale height from 50 to 200 for remineralization (surface to floor)
    real dz =                           // vertical grid spacing placeholder
    
    // Build temperature dependent uptake
    for (j = 0; j < Nx; j++)
    {
        for (k = Nz; k > 0; k--)
        {
            // Temperature and Irradiance
            T = phi[0][j][k];    // temperature in grid box
            t_uptake = a * pow(b, c * T);    // temperature dependent uptake rate
            
            irradiance = irr_0 * exp(-ZZ_u[j][k] / irr_scaleheight);  // value of irradiance
            lk = t_uptake / alpha;
            l_uptake = irradiance / sqrt(POW2(irradiance) + POW2(lk));
            
            // Uptake and Remineralizaiton
            N = phi[2][j][k];
            uptake = (l_uptake * t_uptake) * POW2(N) / (N + monod);
            remin_in_box = (1-f) * uptake;
            
            // Scale Height for Remin
            scale_height = (150/tot_depth) * z + 50;
            dz = 1/_dz_w[j][k];
            
            phi = (temp_phi - (dz * f) * uptake) / ( 1 + dz/scale_height );
            remin = - flux / scale_height;
            
            temp_phi = flux;             // Move to the next vertical grid cell
            
            // Return any extra flux of nutrients to bottom grid cell
            if (k == 0)
            {
                remin -= phi / dz;
                temp_phi = 0;           // Reset flux of nutrients at the surface to zero
            }
            
            dphi_dt[2][j][k] += remin + remin_in_box - uptake;
        }
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
                  real **       dphi_dt,
                  bool          is_buoy,
                  real **       u_r,
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
    
    /////////////////////////////////////////////////////
    ///// Arrays used by the Kurganov-Tadmor scheme /////
    /////////////////////////////////////////////////////
    
    if (method_s == METHOD_KT)
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
                dphi_dx[j][k] = minmod( sigma * (phi[j+1][k]-phi[j][k]) * _dx,
                                       (phi[j+1][k]-phi[j-1][k]) * _2dx,
                                       sigma * (phi[j][k]-phi[j-1][k]) * _dx  );
            }
        }
        
#pragma parallel
        
        // Determine limited z-slopes at cell centres
        for (j = 0; j < Nx; j ++)
        {
            for (k = 1; k < Nz-1; k ++)
            {
                dphi_dz[j][k] = minmod( sigma * (phi[j][k+1]-phi[j][k]) * _dz_w[j][k+1],
                                       (phi[j][k+1]-phi[j][k-1]) * _2dz_phi[j][k],
                                       sigma * (phi[j][k]-phi[j][k-1]) * _dz_w[j][k]  );
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
        FFx[0][k] = 0;
        FFx[Nx][k] = 0;
        
        for (j = 1; j < Nx; j ++)
        {
            // Current depth
            z = ZZ_u[j][k];
            
            // True isopycnal slope defines direction of isopycnal mixing in sigma coordinates
            Siso_true = Siso_u[j][k] - (ZZ_phi[j][k]-ZZ_phi[j-1][k])*_dx;
            
            // Flux calculation depends on spatial discretisation method
            switch (method_s)
            {
                case METHOD_CENTERED:
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
                case METHOD_KT:
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
                                                    + (Siso_u[j][k] + ZZ_u[j][k]*sb_psi[j]/hb_psi[j]) * 0.5*(dphi_dz[j][k]+dphi_dz[j-1][k]) );
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
        FFz[j][0] = 0;
        FFz[j][Nz] = 0;
        
        // Interior gridpoints
        for (k = 1; k < Nz; k ++)
        {
            // Current depth
            z = ZZ_w[j][k];
            
            // True isopycnal slope defines direction of isopycnal mixing in sigma coordinates
            Siso_true = Siso_w[j][k] - (ZZ_psi[j+1][k]-ZZ_psi[j][k])*_dx;
            
            // Flux calculation depends on spatial discretisation method
            switch (method_s)
            {
                case METHOD_CENTERED:
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
                case METHOD_KT:
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
                              + (PPz[j][k+1] - PPz[j][k]) * _dz_phi[j][k]
                              - (FFx[j+1][k] - FFx[j][k]) * _dx * _dz_phi[j][k] // Diabatic BL fluxes
                              - (FFz[j][k+1] - FFz[j][k]) * _dz_phi[j][k] );
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
    
    // Pointer to buoyancy matrix
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
    real u_max = 0;
    real w_dz = 0;
    real w_dz_max = 0;
    real xdiff = 0;
    real xdiff_max = 0;
    real zdiff_dzsq = 0;
    real zdiff_dzsq_max = 0;
    
    
    
    ////////////////////////////////////////
    ///// BEGIN CALCULATING TENDENCIES /////
    ////////////////////////////////////////
    
    // Pointer to buoyancy matrix
    buoy = phi[idx_buoy];
    
    // Calculate Gent-McWilliams diffusivity Kgm
    calcKgm(t,buoy,Kgm_psi,Kgm_u,Kgm_w);
    
    // Calculate isopycnal diffusivity Kiso
    calcKiso(t,buoy,Kiso_u,Kiso_w);
    
    // Calculate isopycnal slopes
    calcSlopes(t,buoy,Sgm_psi,Siso_u,Siso_w);
    
    // Calculate mean and eddy components of overturning streamfunction
    calcPsim(t,buoy,psi_m);
    calcPsie(t,buoy,Kgm_psi,Sgm_psi,psi_e);
    
    // Compute residual streamfunction and velocities
    calcPsir(psi_m,psi_e,psi_r,u_r,w_r);
    
    // Loop over tracers and calculate advective/diffusive tendency
    for (i = 0; i < Ntracs; i ++)
    {
        // Advection/diffusion depends on whether the tracer is the buoyancy variable
        if (i == idx_buoy)
        {
            is_buoy = true;
        }
        else
        {
            is_buoy = false;
        }
        
        do_adv_diff(t,phi[i],dphi_dt[i],is_buoy,u_r,w_r,Kiso_u,Kiso_w,Siso_u,Siso_w);
    }
    
    // Add tendency due to biogeochemistry
    tderiv_bgc(t,phi,dphi_dt);
    
    // Add tendency due to relaxation
    tderiv_relax(t,phi,dphi_dt);
    
    //////////////////////////////////////
    ///// END CALCULATING TENDENCIES /////
    //////////////////////////////////////
    
    
    
    //////////////////////////////////
    ///// BEGIN CALCULATING CFLS /////
    //////////////////////////////////
    
#pragma parallel
    
    // For y-fluxes
    for (j = 0; j <= Nx; j ++)
    {
        for (k = 0; k < Nz; k ++)
        {
            // Max advecting velocity
            if (fabs(u_r[j][k]) > u_max)
            {
                u_max = fabs(u_r[j][k]);
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
            // Max advecting velocity
            w_dz = fabs(w_r[j][k]) * _dz_w[j][k];
            if (w_dz > w_dz_max)
            {
                w_dz_max = w_dz;
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
        }
    }
    
    
    
    // Calculate CFL criteria
    cfl_u = 0.5*dx/u_max;
    cfl_w = 0.5/w_dz_max;
    cfl_y = 0.5*dxsq/xdiff_max;
    cfl_z = 0.5/zdiff_dzsq_max;
    
    // Actual CFL-limted time step
    cfl_dt = fmin(fmin(cfl_u,cfl_w),fmin(cfl_y,cfl_z));
    
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
 * writeModelGrid
 *
 * Writes model vertical grid to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
bool writeModelGrid (real ** ZZ_phi, real ** ZZ_psi, real ** ZZ_u, real ** ZZ_w, char * outdir)
{
    char outfile[MAX_PARAMETER_FILENAME_LENGTH];
    FILE * outfd = NULL;
    
    // Write depths on phi-points
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_ZZ_PHI);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,ZZ_phi,Nx,Nz);
    fclose(outfd);
    
    // Write depths on psi-points
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_ZZ_PSI);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,ZZ_psi,Nx+1,Nz+1);
    fclose(outfd);
    
    // Write depths on u-points
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_ZZ_U);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,ZZ_u,Nx+1,Nz);
    fclose(outfd);
    
    // Write depths on w-points
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_ZZ_W);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,ZZ_w,Nx,Nz+1);
    fclose(outfd);
    
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
    
    // Create a string for the current iteration number
    sprintf(nstr,"%d",n);
    
    // Calculate residual streamfunction
    buoy = phi[idx_buoy];
    calcKgm(t,buoy,Kgm_psi,Kgm_u,Kgm_w);
    calcSlopes(t,buoy,Sgm_psi,Siso_u,Siso_w);
    calcPsim(t,buoy,psi_m);
    calcPsie(t,buoy,Kgm_psi,Sgm_psi,psi_e);
    calcPsir(psi_m,psi_e,psi_r,u_r,w_r);
    
    // Loop over tracers
    for (i = 0; i < Ntracs; i ++)
    {
        // Create a string for the current tracer number
        sprintf(istr,"%d",i);
        
        // Write iteration data for this tracer
        strcpy(outfile,outdir);
        strcat(outfile,"/");
        strcat(outfile,OUTN_TRAC);
        strcat(outfile,istr);
        strcat(outfile,"_n=");
        strcat(outfile,nstr);
        strcat(outfile,".dat");
        outfd = fopen(outfile,"w");
        if (outfd == NULL)
        {
            fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
            return false;
        }
        printMatrix(outfd,phi[i],Nx,Nz);
        fclose(outfd);
    }
    
    // Write iteration data for the residual streamfunction
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_PSIR);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,psi_r,Nx+1,Nz+1);
    fclose(outfd);
    
    // Write iteration data for the mean streamfunction
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_PSIM);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,psi_m,Nx+1,Nz+1);
    fclose(outfd);
    
    // Write iteration data for the eddy streamfunction
    strcpy(outfile,outdir);
    strcat(outfile,"/");
    strcat(outfile,OUTN_PSIE);
    strcat(outfile,"_n=");
    strcat(outfile,nstr);
    strcat(outfile,".dat");
    outfd = fopen(outfile,"w");
    if (outfd == NULL)
    {
        fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
        return false;
    }
    printMatrix(outfd,psi_e,Nx+1,Nz+1);
    fclose(outfd);
    
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
     "  name               value\n"
     "  \n"
     "  Ntracs             Number of tracer variables in the simulation. Must be >0.\n"
     "  Nx                 Number of grid cells in the y-direction. Must be >0.\n"
     "  Nz                 Number of grid cells in the z-direction. Must be >0.\n"
     "  Lx                 Domain length. Must be > 0.\n"
     "  Lz                 Domain height. Must be >0.\n"
     "  cflFrac            CFL number. The time step dt will be chosen at each\n"
     "                     iteration such that dt=cfl*dt_max, where dt_max is the\n"
     "                     estimated maximum stable time step. Must be >0.\n"
     "  maxTime            Maximum time for the integration, in\n"
     "                     case the target residual is not achieved.\n"
     "  monitorFrequency   Frequency with which to write out iteration data. If\n"
     "                     set to zero, only the final data will be written.\n"
     "                     Otherwise must be >0.\n"
     "                     Optional, default is 0, must be >= 0.\n"
     "  f0                 Reference Coriolis parameter.\n"
     "                     Optional, default is 10^{-4} rad/s, must be =/= 0.\n"
     "  h_c                Depth parameter controlling the range of depths over\n"
     "                     which the vertical coordinate the coordinate is\n"
     "                     approximately aligned with geopotentials.\n"
     "                     Must be >0. Default is 1e16 (pure sigma coordinate).\n"
     "  theta_s            Surface sigma coordinate stretching parameter,\n"
     "                     defined `meaningfully' between 0 and 10.\n"
     "                     Default is 0 (no surface stretching).\n"
     "  theta_b            Bottom sigma coordinate stretching parameter,\n"
     "                     defined `meaningfully' between 0 and 4.\n"
     "  Hsml               Depth of the surface mixed layer. Must be >=0.\n"
     "                     If equal to 0 then no SML is imposed\n"
     "  Hbbl               Depth of the bottom boundary layer. Must be >=0.\n"
     "                     If equal to 0 then no BBL is imposed\n"
     "  targetResFile      File containing an Ntracs x 1 vector of target L2\n"
     "                     errors between successive iterations. The\n"
     "                     calculation will stop when this condition is met\n"
     "                     for all tracers.\n"
     "  initFile           File containing an Ntracs x Nx x Nz array of initial\n"
     "                     tracer values over the whole domain.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  topogFile          File containing an Nx+1 x 1 array of ocean depths\n"
     "                     at grid cell corners. All elements must be > 0.\n"
     "                     Optional, default is Lz everywhere.\n"
     "  tauFile            File containing an Nx+1 x 1 array of wind stresses\n"
     "                     at grid cell corners.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  KgmFile            File containing an Nx+1 x Nz+1 array of GM eddy\n"
     "                     diffusivities at grid cell corners. All\n"
     "                     elements must be >=0.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  KisoFile           File containing an Nx+1 x Nz+1 array of isopycnal\n"
     "                     diffusivities at grid cell corners. All\n"
     "                     elements must be >=0.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  KdiaFile           File containing an Nx+1 x Nz+1 matrix of diapycnal\n"
     "                     diffusivities at grid cell corners. All\n"
     "                     elements must be >= 0.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  relaxTracerFile    File containing an Ntracers x Nx x Nz array of tracer\n"
     "                     relaxation target values over the whole domain.\n"
     "                     Optional, default is 0 everywhere.\n"
     "  relaxTimeFile      File containing an Ntracers x Nx x Nz array of tracer\n"
     "                     relaxation time scales over the whole domain. Negative \n"
     "                     values imply no relaxation. Zero values imply\n"
     "                     instantaneous relaxation. Optional, default is -1\n"
     "                     (no relaxation) everywhere.\n"
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
    real n_s = 0;           // Number of saves
    real t_next = 0;        // Next save time
    real cflFrac = 1;       // CFL number
    real wc = 0;            // Time-interpolation weights
    real wn = 0;
    
    // Convergence parameters
    real * targetRes = 0;
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
    
    // Pointer to buoyancy matrix
    real ** buoy = NULL;
    
    // Work arrays for time derivatives - expressed as vectors rather
    // than matrices for compatibility with time integration methods
    real * phi_in_V = NULL; // Input to time-integration method
    real * phi_out_V = NULL; // Output from time-integration method
    real * phi_buf_V = NULL; // Buffer for time-integration method
    real * phi_int_V = NULL; // Buffer for time-interpolation
    real *** phi_in = NULL;
    real *** phi_out = NULL;
    real *** phi_int = NULL;
    
    // Stores data required for parsing input parameters
    paramdata params[NPARAMS];
    
    // Filename holders for input parameter arrays
    char targetResFile[MAX_PARAMETER_FILENAME_LENGTH];
    char initFile[MAX_PARAMETER_FILENAME_LENGTH];
    char topogFile[MAX_PARAMETER_FILENAME_LENGTH];
    char tauFile[MAX_PARAMETER_FILENAME_LENGTH];
    char KgmFile[MAX_PARAMETER_FILENAME_LENGTH];
    char KisoFile[MAX_PARAMETER_FILENAME_LENGTH];
    char KdiaFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTracerFile[MAX_PARAMETER_FILENAME_LENGTH];
    char relaxTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
    
    // To store current time
    time_t now;
    FILE * tfile = NULL;
    
    // Define input parameter data
    int paramcntr = 0;
    setParam(params,paramcntr++,"Ntracs","%u",&Ntracs,false);
    setParam(params,paramcntr++,"Nx","%u",&Nx,false);
    setParam(params,paramcntr++,"Nz","%u",&Nz,false);
    setParam(params,paramcntr++,"Lx","%lf",&Lx,false);
    setParam(params,paramcntr++,"Lz","%lf",&Lz,false);
    setParam(params,paramcntr++,"cflFrac","%lf",&cflFrac,false);
    setParam(params,paramcntr++,"maxTime","%lf",&tmax,false);
    setParam(params,paramcntr++,"monitorFrequency","%lf",&dt_s,true);
    setParam(params,paramcntr++,"rho0","%lf",&rho0,true);
    setParam(params,paramcntr++,"f0","%lf",&f0,true);
    setParam(params,paramcntr++,"Kconv","%lf",&Kconv0,true);
    setParam(params,paramcntr++,"h_c","%le",&h_c,true);
    setParam(params,paramcntr++,"theta_s","%lf",&theta_s,true);
    setParam(params,paramcntr++,"theta_b","%lf",&theta_b,true);
    setParam(params,paramcntr++,"Hsml","%lf",&Hsml,true);
    setParam(params,paramcntr++,"Hbbl","%lf",&Hbbl,true);
    setParam(params,paramcntr++,"targetResFile","%s",&targetResFile,false);
    setParam(params,paramcntr++,"initFile","%s",initFile,false);
    setParam(params,paramcntr++,"topogFile","%s",topogFile,true);
    setParam(params,paramcntr++,"tauFile","%s",tauFile,true);
    setParam(params,paramcntr++,"KgmFile","%s",KgmFile,true);
    setParam(params,paramcntr++,"KisoFile","%s",KisoFile,true);
    setParam(params,paramcntr++,"KdiaFile","%s",KdiaFile,true);
    setParam(params,paramcntr++,"relaxTracerFile","%s",relaxTracerFile,true);
    setParam(params,paramcntr++,"relaxTimeFile","%s",relaxTimeFile,true);
    
    // Default file name parameters - zero-length strings
    targetResFile[0] = '\0';
    initFile[0] = '\0';
    topogFile[0] = '\0';
    tauFile[0] = '\0';
    KgmFile[0] = '\0';
    KisoFile[0] = '\0';
    KdiaFile[0] = '\0';
    relaxTracerFile[0] = '\0';
    relaxTimeFile[0] = '\0';
    
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
    if (  (Ntracs <= 0) ||
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
        (strlen(topogFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(tauFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(KgmFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(KisoFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(KdiaFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(relaxTracerFile) > MAX_PARAMETER_FILENAME_LENGTH) ||
        (strlen(relaxTimeFile) > MAX_PARAMETER_FILENAME_LENGTH) )
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
    
    // Current time and next save point
    t = tmin;
    if (dt_s == 0.0)
    {
        t_next = 2*tmax; // If no save interval is specified, only write final data
    }
    else
    {
        t_next = tmin + dt_s;
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
    
    // Allocate parameter arrays
    VECALLOC(res,Ntracs);
    VECALLOC(targetRes,Ntracs);
    MATALLOC3(phi_init,Ntracs,Nx,Nz);
    VECALLOC(hb_in,Nx+2);
    VECALLOC(tau,Nx+1);
    MATALLOC(Kgm_psi_ref,Nx+1,Nz+1);
    MATALLOC(Kiso_psi_ref,Nx+1,Nz+1);
    MATALLOC(Kdia_psi_ref,Nx+1,Nz+1);
    MATALLOC3(phi_relax,Ntracs,Nx,Nz);
    MATALLOC3(T_relax,Ntracs,Nx,Nz);
    
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
    MATALLOC(FFx,Nx+1,Nz);
    MATALLOC(FFz,Nx,Nz+1);
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
    MATALLOC(Kgm_psi,Nx+1,Nz+1);
    MATALLOC(Kgm_u,Nx+1,Nz);
    MATALLOC(Kgm_w,Nx,Nz+1);
    MATALLOC(Kiso_u,Nx+1,Nz);
    MATALLOC(Kiso_w,Nx,Nz+1);
    MATALLOC(Kdia_w,Nx,Nz+1);
    
    // Boundary layer work arrays
    k_sml = malloc((Nx+1)*sizeof(uint));
    VECALLOC(wn_sml,Nx+1);
    VECALLOC(wp_sml,Nx+1);
    k_bbl = malloc((Nx+1)*sizeof(uint));
    VECALLOC(wn_bbl,Nx+1);
    VECALLOC(wp_bbl,Nx+1);
    VECALLOC(db_dx,Nz+1);
    VECALLOC(db_dz,Nz+1);
    
    /////////////////////////////////
    ///// END MEMORY ALLOCATION /////
    /////////////////////////////////
    
    
    
    
    ////////////////////////////////////////
    ///// BEGIN READING PARAMETER DATA /////
    ////////////////////////////////////////
    
    // Read input matrices and vectors
    if (  ( (strlen(targetResFile) > 0)   &&  !readVector(targetResFile,targetRes,Ntracs,stderr) ) ||
        ( (strlen(initFile) > 0)        &&  !readMatrix3(initFile,phi_init,Ntracs,Nx,Nz,stderr) ) ||
        ( (strlen(topogFile) > 0)       &&  !readVector(topogFile,hb_in,Nx+2,stderr) ) ||
        ( (strlen(tauFile) > 0)         &&  !readVector(tauFile,tau,Nx+1,stderr) ) ||
        ( (strlen(KgmFile) > 0)         &&  !readMatrix(KgmFile,Kgm_psi_ref,Nx+1,Nz+1,stderr) ) ||
        ( (strlen(KisoFile) > 0)        &&  !readMatrix(KisoFile,Kiso_psi_ref,Nx+1,Nz+1,stderr) )  ||
        ( (strlen(KdiaFile) > 0)        &&  !readMatrix(KdiaFile,Kdia_psi_ref,Nx+1,Nz+1,stderr) )  ||
        ( (strlen(relaxTracerFile) > 0) &&  !readMatrix3(relaxTracerFile,phi_relax,Ntracs,Nx,Nz,stderr) )  ||
        ( (strlen(relaxTimeFile) > 0)   &&  !readMatrix3(relaxTimeFile,T_relax,Ntracs,Nx,Nz,stderr) )  )
    {
        printUsage();
        return 0;
    }
    
    
    // Default initial condition is zero tracer everywhere
    if (strlen(initFile)==0)
    {
        memset(*(*phi_init),0,(Ntot)*sizeof(real));
    }
    
    // Default ocean depth is zero everywhere
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
        memset(tau,0,(Nx+1)*sizeof(real));
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
    
    // Default Kdia is zero everywhere
    if (strlen(KdiaFile)==0)
    {
        memset(*Kdia_psi_ref,0,(Nx*Nz)*sizeof(real));
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
    for (i = 0; i < Ntracs; i ++)
    {
        if (targetRes[i]<=0)
        {
            fprintf(stderr,"targetResFile may contain only values that are >0.");
            printUsage();
            return 0;
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
                fprintf(stderr,"KdiaFile may contain only values that are >0.");
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
    
    // Convergence residuals
    for (i = 0; i < Ntracs; i ++)
    {
        res[i] = 10*targetRes[i];
    }
    
    //////////////////////////////////
    ///// END INITIAL CONDITIONS /////
    //////////////////////////////////
    
    
    
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
    
    // Numerical time-integration loop - keep going while the residual exceeds
    // the target and the current time does not exceed the max time
    while (!targetReached && (t < tmax))
    {
        // Step 1: Perform a single numerical time-step for all explicit terms in the equations
        switch (method_t)
        {
            case METHOD_RKTVD1:
            {
                dt = rktvd1(&t,phi_in_V,phi_out_V,cflFrac,Ntot,&tderiv);
                break;
            }
            case METHOD_RKTVD2:
            {
                dt = rktvd2(&t,phi_in_V,phi_out_V,phi_buf_V,cflFrac,Ntot,&tderiv);
                break;
            }
            case METHOD_RKTVD3:
            {
                dt = rktvd3(&t,phi_in_V,phi_out_V,phi_buf_V,cflFrac,Ntot,&tderiv);
                break;
            }
            default:
            {
                fprintf(stderr,"ERROR: Unknown time-integration method\n");
                break;
            }
        }
        
        // Step 2: Add implicit vertical diffusion
        do_impl_diff(t,dt,phi_out);
        
#pragma parallel
        
        // Step 3: Enforce zero tendency where relaxation time is zero
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
        
        // Step 4: Calculate the residuals as an L2-norm between adjacent iterations
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
        
        // Step 5: If the time step has taken us past a save point (or multiple
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
            printf("%e\n",res[idx_buoy]);
            fflush(stdout);
        }
        
        // Copy the next iteration from phi_out back to phi_in,
        // ready for the next time step
        memcpy(phi_in_V,phi_out_V,Ntot*sizeof(real));
        
        // Increment iteration count
        nIters += 1;
    }
    
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
    
    
    
    return 0;
}
