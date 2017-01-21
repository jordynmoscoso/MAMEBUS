/**
 * defs.c
 *
 * Andrew Stewart
 *
 * Contains a series of subroutines used in mamebus.c.
 *
 */

#include "defs.h"


/**
 * stretch_ROMS
 *
 * Implements the topographic stretching function of A. Shchepetkin (2010). 
 *
 * sigma is the fractional vertical stretching coordinate, ranging from -1 to 0.
 * h_c is a (positive) depth parameter controlling the range of depths over which 
 *    the coordinates are approximately aligned with geopotentials.
 * theta_s is the surface stretching parameter, defined "meaningfully" between 0 and 10
 * theta_b is the bottom stretching parameter, defined "meaningfully" between 0 and 4
 * h_b is the thickness of the fluid column
 *
 */
real stretch_ROMS (real sigma, real h_c, real theta_s, real theta_b, real h_b)
{
  real z, C;
  
  // Surface refinement function
  if (theta_s > 0)
  {
    C = (1 - cosh(theta_s*sigma)) / (cosh(theta_s) - 1);
  }
  else
  {
    C = - sigma*sigma;
  }
  
  // Augment surface refinement function to include bottom refinement
  if (theta_b > 0)
  {
    C = (exp(theta_b*C) - 1) / (1 - exp(-theta_b));
  }
  
  // Calculate actual depth
  z = h_b * (h_c*sigma + h_b*C) / (h_c + h_b);
  
  return z;
}

/**
 * printMatrix
 *
 * Convenience method to write a matrix of data (mat, dimensions m by n)
 * to a data file (outfile).
 *
 */
void printMatrix (FILE * outfile, real ** mat, uint m, uint n)
{
    uint i,j;

    // Write array of data
    for (j = 0; j < n; j ++)
    {
        for (i = 0; i < m; i ++)
        {
           fprintf(outfile,"%e ",mat[i][j]);
        }
        fprintf(outfile,"\r\n");
    }
    fflush(outfile);
}


/**
 *  matalloc3
 *
 *  Convenience method to allocate a real matrix with dimensions
 *  (Nr,Nc,Ns). Returns a pointer to an array of pointers that point 
 *  to an array of pointers that point to the beginning of each row 
 *  in the matrix.
 *
 */
real *** matalloc3 (uint Nr, uint Nc, uint Ns)
{
  real * slices = malloc(Nr*Nc*Ns*sizeof(real));
  real ** cols = malloc(Nr*Nc*sizeof(real *));
  real *** rows = malloc(Nr*sizeof(real **));
  
  uint i, j;
  
  if ((slices == NULL) || (cols == NULL) || (rows == NULL))
  {
    return NULL;
  }
  
  memset(slices,0,Nr*Nc*Ns*sizeof(real));
  
  for (i = 0; i < Nr; i ++)
  {
    rows[i] = cols + Nc*i;
    for (j = 0; j < Nc; j ++)
    {
      cols[i*Nc+j] = slices + Ns*j + Nc*Ns*i;
    }
  }
  
  return rows;
}


/**
 *  matfree3
 *
 *  Frees memory used by a matrix allocated with matalloc3.
 *
 */
void matfree3 (real *** mat3)
{
  free(*(*mat3));
  free(*mat3);
  free(mat3);
}


/**
 *  matalloc
 *
 *  Convenience method to allocate a real matrix with dimensions
 *  (m,n). Returns a pointer to an array of pointers that point to the
 *  beginning of each row in the matrix.
 *
 */
real ** matalloc (uint m, uint n)
{
    real * data = malloc(m*n*sizeof(real));
    real ** mat = malloc(m*sizeof(real *));
    uint i;

    if ((data == NULL) || (mat == NULL))
    {
        return NULL;
    }

    memset(data,0,m*n*sizeof(real));
    mat[0] = data;

    for (i = 1; i < m; i ++)
    {
        mat[i] = mat[i-1] + n;
    }

    return mat;
}


/**
 *  matfree
 *
 *  Frees memory used by a matrix allocated with matalloc.
 *
 */
void matfree (real ** mat)
{
    free(*mat);
    free(mat);
}


/**
 *  vecalloc
 *
 *  Convenience method to allocate a real vector of length m.
 *  Returns a pointer to the first element in the allocated vector.
 *
 */
real * vecalloc (uint m)
{
    real * vec = malloc(m*sizeof(real));
    memset(vec,0,m*sizeof(real));
    return vec;
}


/**
 *  vecfree
 *
 *  Frees memory used by a vector allocated with vecalloc.
 *
 */
void vecfree (real * vec)
{
    free(vec);
}


/**
 *  vec2mat
 *
 * Creates a matrix from a vector of length m*n. If the matrix pointed
 * to by pmat (*pmat) is NULL, a vector of pointers is created that must
 * later be freed.
 */
void vec2mat (real * vec, real *** pmat, uint m, uint n)
{
  int i;
  
  if (*pmat == NULL)
  {
    *pmat = malloc(m*sizeof(real *));
  }
  
  (*pmat)[0] = vec;
  
  for (i = 1; i < m; i ++)
  {
    (*pmat)[i] = (*pmat)[i-1] + n;
  }
}


/**
 *  vec2mat3
 *
 * Creates a 3-dimensional matrix from a vector of length Nr*Nc*Ns. If the matrix pointed
 * to by pmat3 (*pmat3) is NULL, a vector of pointers is created that must
 * later be freed.
 */
void vec2mat3 (real * vec, real **** pmat3, uint Nr, uint Nc, uint Ns)
{
  int i,j;
  real *** rows = NULL;
  real ** cols = NULL;
  real * slices = NULL;
  
  // Change nomenclature for clarity
  slices = vec;
  
  // If the matrix doesn't exist then create one
  if (*pmat3 == NULL)
  {
    cols = malloc(Nr*Nc*sizeof(real *));
    rows = malloc(Nr*sizeof(real **));
  }
  else
  {
    rows = *pmat3;
    cols = *rows;
  }
  
  // Can't create a matrix from nothing
  if ((slices == NULL) || (cols == NULL) || (rows == NULL))
  {
    if (pmat3 != NULL)
    {
      *pmat3 = NULL;
    }
    return;
  }
  
  // Assign row and column pointers
  for (i = 0; i < Nr; i ++)
  {
    rows[i] = cols + Nc*i;
    for (j = 0; j < Nc; j ++)
    {
      cols[i*Nc+j] = slices + Ns*j + Nc*Ns*i;
    }
  }
  
  *pmat3 = rows;
}


/**
 * setParam
 *
 * Convenience method to set the values in the paramdata struct at the specified
 * index (idx) in an array (params). The method sets the parameter name (name)
 * parameter type (type, corresponding to the format string required for fscanf),
 * the pointer to the variable that will hold the parameter value (pvar) and
 * whether the parameter should be marked as having been read (read).
 */
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read)
{
    params[idx].name = name;
    params[idx].type = type;
    params[idx].pvar = pvar;
    params[idx].read = read;
}


/**
 * readParams
 * 
 * Reads data into the specified array of paramdata structures
 * 'params' of length 'nparams' from the file specified by
 * 'infname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 */
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm)
{
    FILE * infile = NULL;
    bool readok = true;
    char inbuf[256];
    int i = 0;

    // Basic error checking
    if (infname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readParams'\r\n");
        }
        return false;
    }
    if (params == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL parameter information supplied to 'readParams'\r\n");
        }
        return false;
    }

    // Open the input file and check that it can be read 
    infile = fopen(infname,"r");
    if (infile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not find input parameter file\r\n");
        }
        return false;
    }

    // Scan input file for parameter values
    while (fscanf(infile,"%s",inbuf) != EOF)
    {
        // Check whether the string read matches any parameter names
        for (i = 0; i < nparams; i ++)
        {            
            if (strcmp(inbuf,params[i].name) == 0)
            {
                if (fscanf(infile,params[i].type,params[i].pvar) == EOF)
                {        
                    if (errstrm != NULL)
                    {
                        fprintf(errstrm,"ERROR: Could not read parameter %s\r\n",inbuf);
                    }
                    readok = false;
                }
                else
                {                    
                    params[i].read = true;
                }
                break;
            }
        }


        // Parameter name doesn't match any of those in the list
        if (i == nparams)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Parameter not recognised: %s\r\n",inbuf);
            }
            readok = false;
        }
    }


    // Check that required parameter values have been specified
    for (i = 0; i < nparams; i ++)
    {
        if (!params[i].read)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Required parameter not specified: %s\r\n",params[i].name);
            }
            readok = false;
        }
    }

    // We don't need the input file any more
    fclose(infile);
    return readok;
}


/**
 * readVector
 *
 * Reads 'len' real values into the vector 'vec' from the file
 * specified by 'fname'. Any error will result in 'false' being
 * returned, and an error message being printed to the specified
 * FILE/stream 'errstrm'.
 *
 */
bool readVector (char * fname, real * vec, uint len, FILE * errstrm)
{
    // Attempt to open the topography file
    FILE * pfile = NULL;
    int i = 0;    
    bool readok = true;

    // Basic error checking
    if (fname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readVector'\r\n");
        }
        return false;
    }
    if (vec == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL vector supplied to 'readVector'\r\n");
        }
        return false;
    }
 
    // Attempt to open the file
    pfile = fopen(fname,"r");

    // Check that the file exists
    if (pfile == NULL)
    {
        if (errstrm != NULL)
        {
	    fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
        }
        return false;
    }

    // Read data from the file
    for (i = 0; i < len; i ++)
    {
        if (fscanf(pfile,"%le",&vec[i]) == EOF)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",i,len,fname);
            }
            readok = false;
            break;
       }
    }

    // Close the file and exit
    fclose(pfile);
    return readok;
}


/**
 * readMatrix
 *
 * Reads 'm'x'n' real values into the matrix 'mat' in column-major
 * order (i.e. reads mat[1][1], mat[2][1], mat[3][1], ... , mat[1][2],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm)
{
  // Attempt to open the topography file
  FILE * pfile = NULL;
  int i = 0;
  int j = 0;
  bool readok = true;
  
  // Basic error checking
  if (fname == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix'\r\n");
    }
    return false;
  }
  if (mat == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix'\r\n");
    }
    return false;
  }
  for (i = 0; i < m; i ++)
  {
    if (mat[i] == NULL)
    {
      if (errstrm != NULL)
      {
        fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix'\r\n");
      }
      return false;
    }
  }
  
  // Attempt to open the file
  pfile = fopen(fname,"r");
  
  // Check that the file exists
  if (pfile == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
    }
    return false;
  }
  
  // Read in 'len' real values from the file
  for (j = 0; j < n; j ++)
  {
    for (i = 0; i < m; i ++)
    {
      if (fscanf(pfile,"%le",&mat[i][j]) == EOF)
      {
        if (errstrm != NULL)
        {
          fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",j*m+i,m*n,fname);
        }
        readok = false;
        break;
      }
    }
  }
  
  // Close the file and exit
  fclose(pfile);
  return readok;
}


/**
 * readMatrix3
 *
 * Reads 'Nr'x'Nc'x'Ns' real values into the matrix 'mat3' in column-major
 * order (i.e. reads mat[1][1][1], mat[2][1][1], mat[3][1][1], ... , mat[1][2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm)
{
  FILE * pfile = NULL;
  int i = 0;
  int j = 0;
  int k = 0;
  bool readok = true;
  
  // Basic error checking
  if (fname == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  if (mat3 == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  for (i = 0; i < Nr; i ++)
  {
    if (mat3[i] == NULL)
    {
      if (errstrm != NULL)
      {
        fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
      }
      return false;
    }
    for (j = 0; j < Nc; j ++)
    {
      if (mat3[i][j] == NULL)
      {
        if (errstrm != NULL)
        {
          fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix3'\r\n");
        }
        return false;
      }
    }
  }
  
  // Attempt to open the file
  pfile = fopen(fname,"r");
  
  // Check that the file exists
  if (pfile == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
    }
    return false;
  }
  
  // Read data from the file
  for (k = 0; k < Ns; k ++)
  {
    for (j = 0; j < Nc; j ++)
    {
      for (i = 0; i < Nr; i ++)
      {
        if (fscanf(pfile,"%le",&mat3[i][j][k]) == EOF)
        {
          if (errstrm != NULL)
          {
            fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",k*Nr*Nc+j*Nr+i,Nr*Nc*Ns,fname);
          }
          readok = false;
          break;
        }
      }
    }
  }
  
  // Close the file and exit
  fclose(pfile);
  return readok;
}


/**
 * max3
 *
 * Calculates the maximum of three real values. Not a very pretty implementation,
 * but always employs the minimal number of comparisons (2) to calculate the max.
 *
 */
real max3 (real v1, real v2, real v3)
{
    if (v1 > v2)
    {
        if (v1 > v3)
        {
            return v1;
        }
        else
        {
            return v3;
        }
    }
    else
    {
        if (v2 > v3)
        {
            return v2;
        }
        else
        {
            return v3;
        }
    }
}


/**
 * min3
 *
 * Calculates the minimum of three real values. Not a very pretty implementation,
 * but always employs the minimal number of comparisons (2) to calculate the min.
 *
 */
real min3 (real v1, real v2, real v3)
{
    if (v1 < v2)
    {
        if (v1 < v3)
        {
            return v1;
        }
        else
        {
            return v3;
        }
    }
    else
    {
        if (v2 < v3)
        {
            return v2;
        }
        else
        {
            return v3;
        }
    }
}


/**
 * minmod
 *
 * Implements the minmod function of Kurganov and Tadmor (2000) for three
 * quantities.
 *
 */
real minmod (real v1, real v2, real v3)
{
    real test = v1*v2*v3;
    
    if (v1 < 0 && v2 < 0 && v3 < 0)
    {
        return max3(v1,v2,v3);
    }
    else if (v1 > 0 && v2 > 0 && v3 > 0)
    {
        return min3(v1,v2,v3);
    }
    else
    {
        return 0;
    }
}


/**
 * limMin
 *
 * Returns val, or, if the size of val is smaller than minVal,
 * returns sign(val)*minVal.
 *
 */
real limMin (real val, real minVal)
{
    real val_lim = val;
    
    if ((val < minVal) && (val >= 0))  
    {
        val_lim = minVal;
    }
    if ((val > -minVal) && (val < 0))  
    {
        val_lim = -minVal;
    }

    return val_lim;
}


/**
 * thomas
 *
 * Performs the thomas algorithm. Here a, b and c are the sub-diagonal, diagonal
 * and super-diagonal in the tridiagonal coefficient matrix, and d is the vector
 * on the right hand side. The solution is filled into x.
 * Warning: will modify c and d!
 */
void thomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n)
{
  real id;
  int i;
    
	// Modify the coefficients
	c[0] /= b[0];	// Division by zero risk
	d[0] /= b[0];	// Division by zero would imply a singular matrix
	for (i = 1; i < n; i ++)
  {
		id = 1 / (b[i] - c[i-1]*a[i]);     // Division by zero risk
		c[i] *= id;	                       // Last value calculated is redundant
		d[i] = (d[i] - d[i-1]*a[i]) * id;
	}
    
	// Now back substitute
	x[n - 1] = d[n - 1];
	for (i = n-2; i >= 0; i --)
  {
		x[i] = d[i] - c[i] * x[i + 1];
  }
}
