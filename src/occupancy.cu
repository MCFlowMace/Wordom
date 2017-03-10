#include <math.h>
#include <stdio.h> 
#include <stdlib.h>

__global__ void normalize(float a[3])
{
  float  b;

  b = sqrtf((float)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] /= b;
  a[1] /= b;
  a[2] /= b;
}



__global__ void cross(float a[3], float b[3], float c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

__global__ void dot(float a[3], float b[3],float* c)
{
  //return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  *c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}



/*
 * setup_rotation() 
 *
 *      given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
__global__ void setup_rotation(const float *ref_xlist,
                    const float *mov_xlist, 
                    int& n_list,
                    float R[3][3],
                    float* E0)
{
  int i, j, n;

  /* initialize */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++) 
      R[i][j] = 0.0;
  *E0 = 0.0;

  for (n=0; n<n_list; n++) 
  {
    /* 
     * E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) 
     */
    for (i=0; i<3; i++)
      *E0 +=  mov_xlist[3*n+i] * mov_xlist[3*n+i]  
            + ref_xlist[3*n+i] * ref_xlist[3*n+i];
    
    /*
     * correlation matrix R:   
     *   R[i,j) = sum(over n): y(n,i) * x(n,j)  
     *   where x(n) and y(n) are two vector sets   
     */
    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
        R[i][j] += mov_xlist[3*n+i] * ref_xlist[3*n+j];
    }
  }
  *E0 *= 0.5;
  }



#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }
/*   
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
//__device__ int jacobi3(float a[3][3], float d[3], float v[3][3], int* n_rot)
 

__global__ void jacobiTest(float a[3][3], float d[3], float v[3][3])
{
  int count, k, i, j;
  float tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

  /*Initialize v to the identity matrix.*/
  for (i=0; i<3; i++) 
  { 
    for (j=0; j<3; j++) 
      v[i][j] = 0.0;
    v[i][i] = 1.0;
  }

  /* Initialize b and d to the diagonal of a */
  for (i=0; i<3; i++) 
    b[i] = d[i] = a[i][i];

  /* z will accumulate terms */
  for (i=0; i<3; i++) 
    z[i] = 0.0; 
  
  //*n_rot = 0;

  /* 50 tries */
  for (count=0; count<50; count++)     
  {

    /* sum off-diagonal elements */
    sum = 0.0;
    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++)
         sum += fabsf(a[i][j]);
    }

    /* if converged to machine underflow */
    if (sum == 0.0) 
      return;

    /* on 1st three sweeps... */
    if (count < 3) 
      tresh = sum * 0.2 / 9.0;    
    else       
      tresh = 0.0;      

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        g = 100.0 * fabsf(a[i][j]);

        /*  after four sweeps, skip the rotation if
         *   the off-diagonal element is small 
         */
        if ( count > 3  &&  fabsf(d[i])+g == fabsf(d[i])
              &&  fabsf(d[j])+g == fabsf(d[j]) ) 
        {
          a[i][j] = 0.0;
        } 
        else if (fabsf(a[i][j]) > tresh) 
        {
          h = d[j] - d[i];
          
          if (fabsf(h)+g == fabsf(h))
          {
            t = a[i][j] / h;
          }
          else 
          {
            theta = 0.5 * h / (a[i][j]);
            t = 1.0 / ( fabsf(theta) +
                        (float)sqrtf(1.0 + theta*theta) );
            if (theta < 0.0) 
              t = -t;
          }
          
          c = 1.0 / (float) sqrtf(1 + t*t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

         /* for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)*/

          //++(*n_rot);
        }
      }
    }

    for (i=0; i<3; i++) 
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }

  printf("Too many iterations in jacobi3\n");
  //return (0);
}  



/* 
 * diagonalize_symmetric 
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */

__global__ void diagonalize_Test(float matrix[3][3], float vec[3][3], float eigenval[3])
//__device__ int diagonalize_symmetric(float matrix[3][3], float eigenval[3])
{
  //int n_rot, i, j, k;
  int i, j, k;
  //float vec[3][3];
  float val; 
  
  //if (!jacobi3(matrix, eigenval, vec, &n_rot))
 /* if (!jacobi3(matrix, eigenval, vec))
  {
    //printf("convergence failed\n");
    return ;
  }*/

  /* sort solutions by eigenval */
  for (i=0; i<3; i++) 
  {
    k = i;
    val = eigenval[i];
    
    for (j=i+1; j<3; j++)
      if (eigenval[j] >= val)
      { 
        k = j;
        val = eigenval[k];
      }
       
    if (k != i) 
    {
      eigenval[k] = eigenval[i];
      eigenval[i] = val;
      for (j=0; j<3; j++) 
      {
        val = vec[j][i];
        vec[j][i] = vec[j][k];
        vec[j][k] = val;
      }
    }
  }

  /* transpose such that first index refers to solution index */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      vec[i][j] = vec[j][i];

  //return (1);
}



/*
 * calculate_rotation_matrix() 
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
__global__ void calculate_rotation_matrix(float R[3][3],
                              float& E0,
                              float* residual)
{
  int i, j, k;
  float Rt[3][3], RtR[3][3];
  //float RtR[3][3];
  //float eigenval[3];
  float left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
  float v[3];
  float sigma;

  /* build Rt, transpose of R  */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Rt[i][j] = R[j][i];

  /* make symmetric RtR = Rt X R */
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++)
    {
      RtR[i][j] = 0.0;
      for (k = 0; k<3; k++)
        RtR[i][j] += Rt[k][i] * R[j][k];
        //RtR[i][j] += R[i][k] * R[j][k];
    }

  //if (!diagonalize_symmetric(RtR, right_eigenvec, eigenval))
    //return(0);

  // right_eigenvec's should be an orthogonal system but could be left
  // or right-handed. Let's force into right-handed system.
   
 // cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

  // From the Kabsch algorithm, the eigenvec's of RtR
  // are identical to the right_eigenvec's of R.
  // This means that left_eigenvec = R x right_eigenvec 
   
  /*for (i=0; i<3; i++) 
    for (j=0; j<3; j++) 
      left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);*/
      //left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &R[0][j]);

  //for (i=0; i<3; i++) 
   // normalize(&left_eigenvec[i][0]);

   
   // Force left_eigenvec[2] to be orthogonal to the other vectors.
   // First check if the rotational matrices generated from the 
   // orthogonal eigenvectors are in a right-handed or left-handed
   // co-ordinate system - given by sigma. Sigma is needed to
   // resolve this ambiguity in calculating the RMSD.
   
  /*cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
  if (dot(v, &left_eigenvec[2][0]) < 0.0)
    sigma = -1.0;
  else 
    sigma = 1.0;*/
  for (i=0; i<3; i++)
    left_eigenvec[2][i] = v[i];

    
  *residual = E0 - (float) sqrt(fabs(eigenval[0])) 
                 - (float) sqrt(fabs(eigenval[1]))
                 - sigma * (float) sqrt(fabs(eigenval[2]));

  //return (1);
}

