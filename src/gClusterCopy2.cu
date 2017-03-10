#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
	//including wordom.h generates compiler warnings
	//we do not have to worry about this, since this only happens because nvcc generates c++ object files and the c++ compiler is not completely satisfied with some stuff in wordom.h
#include "wordom.h"
#include "fileio.h"
#include "tools.h"
#include "analysis.h"
#include "cluster.h"

__device__ void normalize(float a[3])
{
  float  b;

  b = sqrtf((float)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] /= b;
  a[1] /= b;
  a[2] /= b;
}



__device__ float dot(float a[3], float b[3])
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}



__device__ static void cross(float a[3], float b[3], float c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}



/*
 * setup_rotation() 
 *
 *      given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
__device__ void setup_rotation(const float *ref_xlist, const float *mov_xlist, int& n_list, float R[3][3], float* E0, volatile int& i, volatile int&j, volatile int& n)
{
  //int i, j, n;

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
__device__ int jacobi3(float a[3][3], float d[3], float v[3][3], volatile int& k, volatile int& i, volatile int& j)
{
  //int count, k, i, j;
  //float tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];
  //float t, sum, s, h, g, c, b[3], z[3];
  float b[3], z[3];

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
  int count;
  for (count=0; count<50; count++)     
  {
	
	float tresh;
    /* sum off-diagonal elements */
    {
    float sum = 0.0;
    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++)
         sum += fabsf(a[i][j]);
    }

    /* if converged to machine underflow */
    if (sum == 0.0) 
      return(1);

    /* on 1st three sweeps... */
    if (count < 3) 
      tresh = sum * 0.2 / 9.0;    
    else       
      tresh = 0.0;
      
  } 

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        float g = 100.0 * fabsf(a[i][j]);

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
          float h = d[j] - d[i];
          float t;
          
          if (fabsf(h)+g == fabsf(h))
          {
            t = a[i][j] / h;
          }
          else 
          {
            float theta = 0.5 * h / (a[i][j]);
            t = 1.0 / ( fabsf(theta) +
                        (float)sqrtf(1.0 + theta*theta) );
            if (theta < 0.0) 
              t = -t;
          }
          
          float c = 1.0 / (float) sqrtf(1 + t*t);
          float s = t * c;
          float tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

          for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)

         // ++(*n_rot);
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
  return (0);
}  



/* 
 * diagonalize_symmetric 
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
__device__ int diagonalize_symmetric(float matrix[3][3], float vec[3][3], float eigenval[3], volatile int& i, volatile int& j, volatile int& k)
{
  //int n_rot, i, j, k;
  //float vec[3][3];
  //float val; 
  
  if (!jacobi3(matrix, eigenval, vec, i,j,k)) 
  {
    //printf("convergence failed\n");
    return (0);
  }

  /* sort solutions by eigenval */
  for (i=0; i<3; i++) 
  {
	  float val;
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

  return (1);
}



/*
 * calculate_rotation_matrix() 
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
//__device__ int calculate_rotation_matrix(float R[3][3], float E0, float* residual, volatile int& i, volatile int& j, volatile int& k)
__device__ int calculate_rotation_matrix(float R[3][3], float& residual, volatile int& i, volatile int& j, volatile int& k)
{
  //int i, j, k;
  float Rt[3][3];
  float right_eigenvec[3][3], eigenval[3];
  //float v[3];
  //float sigma;

  /* build Rt, transpose of R  */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Rt[i][j] = R[j][i];

  /* make symmetric RtR = Rt X R */
  
  {
  float RtR[3][3];
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++)
    {
      RtR[i][j] = 0.0;
      for (k = 0; k<3; k++)
        RtR[i][j] += Rt[k][i] * R[j][k];
    }

  if (!diagonalize_symmetric(RtR, right_eigenvec, eigenval,i,j,k))
    return(0);
  }

  /* right_eigenvec's should be an orthogonal system but could be left
   * or right-handed. Let's force into right-handed system.
   */
  cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

  /* From the Kabsch algorithm, the eigenvec's of RtR
   * are identical to the right_eigenvec's of R.
   * This means that left_eigenvec = R x right_eigenvec 
   */
   
   float left_eigenvec[3][3];
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++) 
      left_eigenvec[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);

  for (i=0; i<3; i++) 
    normalize(&left_eigenvec[i][0]);

  /* 
   * Force left_eigenvec[2] to be orthogonal to the other vectors.
   * First check if the rotational matrices generated from the 
   * orthogonal eigenvectors are in a right-handed or left-handed
   * co-ordinate system - given by sigma. Sigma is needed to
   * resolve this ambiguity in calculating the RMSD.
   */
  //cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
  cross(&right_eigenvec[0][0], &left_eigenvec[0][0], &left_eigenvec[1][0]);
  
  //if (dot(v, &left_eigenvec[2][0]) < 0.0)
  float sigma;
  if (dot(&right_eigenvec[0][0], &left_eigenvec[2][0]) < 0.0)
    sigma = -1.0;
  else 
    sigma = 1.0;
    
  //for (i=0; i<3; i++)
    //left_eigenvec[2][i] = v[i]; 

    
  //*residual = E0 - (float) sqrtf(fabsf(eigenval[0])) - (float) sqrtf(fabsf(eigenval[1]))- sigma * (float) sqrtf(fabsf(eigenval[2]));
  residual = residual - (float) sqrtf(fabsf(eigenval[0])) - (float) sqrtf(fabsf(eigenval[1]))- sigma * (float) sqrtf(fabsf(eigenval[2]));
                 

  return (1);
}



/*__device__ void calculate_rotation_rmsd(const float* ref_xlist,
                             const float* mov_xlist, 
                             int n_list,
                             float* rmsd)
{
  float Eo, residual;
  float R[3][3];
  
  setup_rotation(ref_xlist, mov_xlist, n_list, R, &Eo);
                 
  if(calculate_rotation_matrix(R, Eo, &residual)) {
  
		residual = fabsf(residual); // avoids the awkward case of -0.0 
		*rmsd = sqrtf( fabsf((float) (residual)*2.0/((float)n_list)) ); 
  
	}               
                 
}*/
 
 

__device__ float RmsdCalc_nosup(const float *refcoor, const float *movcoor, int nato)
{
  // compute rmsd and return value
  int          ii;
  float        rmsd,di;
  
  rmsd=0;
  di=0;
  
  for ( ii=0; ii<3*nato; ii++ ) {
		di= refcoor[ii]-movcoor[ii];
		rmsd += di*di;
             
	}
  
  rmsd /= nato;
  rmsd = sqrt ( rmsd );
  
  return rmsd;
}


/*__device__ float RmsdCalc(const float *refcoor, const float *movcoor, int nato, int super )
{
  float rmsd=-1;
  
  if( super == 0 )
    return RmsdCalc_nosup( refcoor, movcoor, nato );
  
  calculate_rotation_rmsd(refcoor,movcoor, nato, &rmsd);
  
  return rmsd;
  
}*/

__global__ void gRmsdMaxSuper(int nato, const int nframes, const float* gclust_coords, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance) {
	
	
	// nato = number of atoms, nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated rmsd because it will be needed in the post processing
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with the maxspeed flag frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {

		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			return;
		}
		
		//float rmsd=-1;
	
		//calculate_rotation_rmsd(&gclust_coords[cluster*3*nato],&gclust_coords[index*3*nato], nato, &rmsd);
		
				//calculate_rotation_rmsd(&gclust_coords[cluster*3*nato],&gclust_coords[index*3*nato], nato, &rmsd);
		float rmsd;
		float R[3][3];
		
		//volatile int i,j,k,n;
		volatile int i,j,n;	
		setup_rotation(&gclust_coords[cluster*3*nato],&gclust_coords[index*3*nato], nato, R, &rmsd, i, j, n);
		
		//setup_rotation(&gclust_coords[cluster*3*nato],&gclust_coords[index*3*nato], nato, R, &E0);		
			
		//float rmsd=-1;
		//if(calculate_rotation_matrix(R, E0, &rmsd, i, j, n)) {
		if(calculate_rotation_matrix(R, rmsd, i, j, n)) {
		//if(calculate_rotation_matrix(R, E0, &rmsd)) {
  
			rmsd = fabsf(rmsd); // avoids the awkward case of -0.0 
			rmsd = sqrtf( fabsf((float) (rmsd)*2.0/((float)nato)) ); 
			
			if (rmsd<cutoff){
				frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
				distance[index] = rmsd;
			return;
			}
  
		}	
			
		/*if (rmsd<cutoff && rmsd != -1){
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}*/
		
	}
	frameapp_write[index] = frameapp_read[index];
			  
			
}

/*__global__ void gRmsdMax(const int nato, const int nframes, float* gclust_coordsX, float* gclust_coordsY, float* gclust_coordsZ, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const int super) {
	
	
	// nato = number of atoms, nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated rmsd because it will be needed in the post processing
	
	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//if(index <100) printf("This works %d\n",index);
	
	//with the maxspeed flag frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
		
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			return;
		}
		
		float *g3_index[3];
		float *g3_cluster[3];
		
		g3_index[0]=&gclust_coordsX[index*nato];
		g3_index[1]=&gclust_coordsY[index*nato];
		g3_index[2]=&gclust_coordsZ[index*nato];
		
		g3_cluster[0]=&gclust_coordsX[cluster*nato];
		g3_cluster[1]=&gclust_coordsY[cluster*nato];
		g3_cluster[2]=&gclust_coordsZ[cluster*nato];
		
		float rmsd;//=0.; , di;
		
		rmsd = RmsdCalc(g3_cluster,g3_index,nato,super);
		
		fast_rmsd(float ref_xlist[][3], float mov_xlist[][3], int n_list, float* rmsd)
		
		//int ii;
		
	    
	    // compute rmsd and return value
			  
		//for ( ii=0; ii<nato; ii++ ) {
			//rmsd += ( (gclust_coordsX[index*nato+ii]-gclust_coordsX[cluster*nato+ii])*(gclust_coordsX[index*nato+ii]-gclust_coordsX[cluster*nato+ii]) + 
             //(gclust_coordsY[index*nato+ii]-gclust_coordsY[cluster*nato+ii])*(gclust_coordsY[index*nato+ii]-gclust_coordsY[cluster*nato+ii]) + 
             //(gclust_coordsZ[index*nato+ii]-gclust_coordsZ[cluster*nato+ii])*(gclust_coordsZ[index*nato+ii]-gclust_coordsZ[cluster*nato+ii]) );
			//}
			  
		//rmsd /= nato;
		//rmsd = sqrt ( rmsd );
	
		if (rmsd<cutoff){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}
		
	}
	frameapp_write[index] = frameapp_read[index];
			  
			
}*/

__global__ void gRmsdClosest(const int nato, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance) {
	
	// nato = number of atoms; nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms for comparison
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//frames can only be assigned to clusters with a lower frame number
	if(cluster > index) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
	
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			return;
		}		
		float rmsd=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<3*nato;jj++){
			di=(gclust_coords[index*3*nato+jj]-gclust_coords[cluster*3*nato+jj]);
			rmsd+=di*di; 
		}			
		rmsd /= nato;
		rmsd = sqrt ( rmsd );
	
		//at the beginning distance is set to the cutoff, by always comparing the drms to the current value of distance instead of only the cutoff we can reassign the frame if we find a closer cluster
		if (rmsd < distance[index]){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}
	}
	
	frameapp_write[index] = frameapp_read[index];
}


__global__ void gRmsdClustersMax(const int nato, const int nframes, const float* gclust_coords, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const int nclusters, const int clustercenter) {
	

	// nato = number of atoms; nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = number of the cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	// distance stores the calculated drms because it will be needed in the post processing;
	// nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x; 	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with maxspeed frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}		

	float rmsd=0. , di;
	int ii;
	
	    // compute rmsd and return value			  

	for ( ii=0; ii<3*nato; ii++ ) {
			di=gclust_coords[(index+nclusters)*3*nato+ii]-gclust_coords[cluster*3*nato+ii];
			rmsd+=di*di;
			/*di=gclust_coords[nframes*nato+index*nato+ii]-gclust_coords[nframes*nato+cluster*nato+ii];
			rmsd+=di*di;
			di=gclust_coords[2*nframes*nato+index*nato+ii]-gclust_coords[2*nframes*nato+cluster*nato+ii];
			rmsd+=di*di;*/
		}
			  
		rmsd /= nato;
		rmsd = sqrt ( rmsd );

		if (rmsd<cutoff){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}
	
	frameapp_write[index] = frameapp_read[index];
}

__global__ void gRmsdFramesMax(const int nato, const int nframes, const float* gclust_coords, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const int framesFinished, const int nclusters, int* newClusters, int* clusterCenters) {
	
	
	// nato = number of atoms; nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms because it will be needed in the post processing;
	// framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	// newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with maxspeed frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
		
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = framesFinished + cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			clusterCenters[*newClusters] = framesFinished + cluster + 1;
			(*newClusters)++;
			return;
		}
		
		float rmsd=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<3*nato;jj++){
			di=gclust_coords[(index+nclusters)*3*nato+jj]-gclust_coords[(cluster+nclusters)*3*nato+jj];
			rmsd+=di*di; 
		}
					
		rmsd /= nato;
		rmsd = sqrt ( rmsd );
	
		if (rmsd<cutoff){		
			frameapp_write[index] = framesFinished + cluster + 1;
			distance[index] = rmsd;
			return;
		}

	}
	frameapp_write[index] = frameapp_read[index];
}


__global__ void gRmsdClustersClosest(const int nato, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int nclusters, const int clustercenter) {
	
	// nato = number of atoms; nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = number of the cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	// distance stores the calculated for comparison;
	// nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	float rmsd=0. , di;
	int jj;
	
	//calculate the rmsd of cluster and index
	for (jj=0;jj<3*nato;jj++){
		di=(gclust_coords[(index+nclusters)*3*nato+jj]-gclust_coords[cluster*3*nato+jj]);
		rmsd+=di*di; 
	}
		
	rmsd /= nato;
	rmsd = sqrt ( rmsd );

	//at the beginning distance is set to cutoff, by always comparing the drms to the current value of distance instead of the cutoff we can reassign the frame if we find a closer cluster
	if (rmsd < distance[index]){		
		frameapp_write[index] = clustercenter;
		distance[index] = rmsd;
		return;
	}
	frameapp_write[index] = frameapp_read[index];
	
}

__global__ void gRmsdFramesClosest(const int nato, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int framesFinished, const int nclusters, int* newClusters, int* clusterCenters) {
	
	// nato = number of atoms; nframes = number of frames used for clustering;
	// gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms for comparison;
	// framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	// newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters	
	
	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//frames can only be assigned to clusters with a lower frame number
	if(cluster > index) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
	
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = framesFinished + cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			clusterCenters[*newClusters] = framesFinished + cluster + 1;
			(*newClusters)++;
			return;
		}		
		float rmsd=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<3*nato;jj++){
			di=(gclust_coords[(index+nclusters)*3*nato+jj]-gclust_coords[(cluster+nclusters)*3*nato+jj]);
			rmsd+=di*di; 
		}			
		rmsd /= nato;
		rmsd = sqrt ( rmsd );
	
		//at the beginning distance is set to cutoff, by always comparing the drms to the current value of distance instead of the cutoff we can reassign the frame if we find a closer cluster
		if (rmsd < distance[index]){	
			frameapp_write[index] = framesFinished + cluster + 1;
			distance[index] = rmsd;
			return;
		}	
	}
	frameapp_write[index] = frameapp_read[index];
	
}

// a kernel for the maxspeed flag
__global__ void gDrmsMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact) {
	
	// msize = size of each distance matrix, nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms because it will be needed in the post processing
	
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with the maxspeed flag frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
		
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			return;
		}
		
		
		float drms=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<msize;jj++){
			di=(gclust_dmtx[index*msize+jj]-gclust_dmtx[cluster*msize+jj]);
			drms+=di*di; 
		}		
		drms = sqrtf(drms);
	    drms *= nointrasegm_corr_fact/sqrtf((float)msize); //Renormalize the distance properly
	
		if (drms<cutoff){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = drms;
			return;
		}
		
	}
	frameapp_write[index] = frameapp_read[index];
}

// a kernel for the lfull flag
__global__ void gDrmsClosest(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact) {
	
	// msize = size of each distance matrix; nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms for comparison

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//frames can only be assigned to clusters with a lower frame number
	if(cluster > index) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
	
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			return;
		}		
		float drms=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<msize;jj++){
			di=(gclust_dmtx[index*msize+jj]-gclust_dmtx[cluster*msize+jj]);
			drms+=di*di; 
		}			
		drms = sqrtf(drms);
	    drms *= nointrasegm_corr_fact/sqrtf((float)msize); //Renormalize the distance properly
	
		//at the beginning distance is set to the cutoff, by always comparing the drms to the current value of distance instead of only the cutoff we can reassign the frame if we find a closer cluster
		if (drms < distance[index]){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = drms;
			return;
		}
	}
	
	frameapp_write[index] = frameapp_read[index];
}

// a kernel for the maxspeed flag for calculation with limited memory; this is for comparing the frames of the chunk to the previously found clusters
__global__ void gDrmsClustersMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int nclusters, const int clustercenter) {
	
	// msize = size of each distance matrix; nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = number of the cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	// distance stores the calculated drms because it will be needed in the post processing;
	// nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x; 	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with maxspeed frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}		
	float drms=0. , di;
	int jj;
	
	//calculate the drms of cluster and index
	for (jj=0;jj<msize;jj++){
		di=(gclust_dmtx[(index+nclusters)*msize+jj]-gclust_dmtx[cluster*msize+jj]);
		drms+=di*di; 
	}			
	drms = sqrtf(drms);
    drms *= nointrasegm_corr_fact/sqrtf(msize); //Renormalize the distance properly

	if (drms<cutoff){		
		frameapp_write[index] = clustercenter;
		distance[index] = drms;
		return;
	}
	
	frameapp_write[index] = frameapp_read[index];
}

// a kernel for the lfull flag for calculation with limited memory; this is for comparing the frames of the chunk to the previously found clusters
__global__ void gDrmsClustersClosest(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int nclusters, const int clustercenter) {
	
	// msize = size of each distance matrix; nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = number of the cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	// distance stores the calculated for comparison;
	// nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}		
	float drms=0. , di;
	int jj;
	
	//calculate the drms of cluster and index
	for (jj=0;jj<msize;jj++){
		di=(gclust_dmtx[(index+nclusters)*msize+jj]-gclust_dmtx[cluster*msize+jj]);
		drms+=di*di; 
	}			
	drms = sqrtf(drms);
    drms *= nointrasegm_corr_fact/sqrtf(msize); //Renormalize the distance properly

	//at the beginning distance is set to cutoff, by always comparing the drms to the current value of distance instead of the cutoff we can reassign the frame if we find a closer cluster
	if (drms < distance[index]){		
		frameapp_write[index] = clustercenter;
		distance[index] = drms;
		return;
	}
	frameapp_write[index] = frameapp_read[index];
}

// a kernel for the maxspeed flag for calculation with limited memory; this is for comparing the frames of the chunk among themselves
__global__ void gDrmsFramesMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int framesFinished, const int nclusters, int* newClusters, int* clusterCenters ) {
	
	// msize = size of each distance matrix; nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms because it will be needed in the post processing;
	// framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	// newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//with maxspeed frames are assigned to the first cluster -> do not check frame again if it was already assigned
	if(frameapp_read[index] != -1) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
		
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = framesFinished + cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			clusterCenters[*newClusters] = framesFinished + cluster + 1;
			(*newClusters)++;
			return;
		}		
		float drms=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<msize;jj++){
			di=(gclust_dmtx[(index+nclusters)*msize+jj]-gclust_dmtx[(cluster+nclusters)*msize+jj]);
			drms+=di*di; 
		}			
		drms = sqrtf(drms);
	    drms *= nointrasegm_corr_fact/sqrtf((float)msize); //Renormalize the distance properly
	
		if (drms<cutoff){		
			frameapp_write[index] = framesFinished + cluster + 1;
			distance[index] = drms;
			return;
		}
	}
	frameapp_write[index] = frameapp_read[index];
}

// a kernel for the lfull flag for calculation with limited memory; this is for comparing the frames of the chunk among themselves
__global__ void gDrmsFramesClosest(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int framesFinished, const int nclusters, int* newClusters, int* clusterCenters ) {
	
	// msize = size of each distance matrix; nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated drms for comparison;
	// framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	// newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters	
	
	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	//frames can only be assigned to clusters with a lower frame number
	if(cluster > index) {
		frameapp_write[index] = frameapp_read[index];
		return;
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {
	
		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = framesFinished + cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0;
			clusterCenters[*newClusters] = framesFinished + cluster + 1;
			(*newClusters)++;
			return;
		}		
		float drms=0. , di;
		int jj;
		
		//calculate the drms of cluster and index
		for (jj=0;jj<msize;jj++){
			di=(gclust_dmtx[(index+nclusters)*msize+jj]-gclust_dmtx[(cluster+nclusters)*msize+jj]);
			drms+=di*di; 
		}			
		drms = sqrtf(drms);
	    drms *= nointrasegm_corr_fact/sqrtf((float)msize); //Renormalize the distance properly
	
		//at the beginning distance is set to cutoff, by always comparing the drms to the current value of distance instead of the cutoff we can reassign the frame if we find a closer cluster
		if (drms < distance[index]){	
			frameapp_write[index] = framesFinished + cluster + 1;
			distance[index] = drms;
			return;
		}	
	}
	frameapp_write[index] = frameapp_read[index];
}

__global__ void shiftToCenter(float* gclust_coords, const int nato, const int nframes) {
  
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	float cms[3];
  
	int ii,jj;

  // calculate the centre of mass 
	for (ii=0; ii<3; ii++)
		cms[ii] = 0.0;
  
	for (ii=0; ii<nato; ii++) 
		for (jj=0; jj<3; jj++)
			cms[jj] += gclust_coords[3*nato*index+3*ii+jj];
			
    
	for (ii=0; ii<3; ii++)
		cms[ii] /= nato;


  // shift mov_xlist and ref_xlist to centre of mass
	for (ii=0; ii<nato; ii++) 
		for (jj=0; jj<3; jj++) 
			gclust_coords[3*nato*index+3*ii+jj] -= cms[jj];

}

//A macro for handling CUDA errors
void errorHandler  (cudaError_t error, int line){
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    fprintf(stderr,"CUDA error: %s in line number %d\n", cudaGetErrorString(error),line);
    exit(-1);
  }
}

// the CUDA compiler generates C++ object files, thus the main procedure has to be an extern "C" for usage in wordom
extern "C" int gClusterRmsd (struct inp_Cluster *inp_cluster,float *distance) {
	
	
	int ii;
    float cutoff = inp_cluster->threshold;
   	int nato = inp_cluster->nato;
	int totframe = inp_cluster->totframe;
	//float **gclust_coords = inp_cluster->gclust_coords;
	float *gclust_coords = inp_cluster->gclust_coords;
    int *frameapp = inp_cluster->frameapp;
    int super = inp_cluster->super;
    int step = inp_cluster->step;
    int frames = totframe/step+(totframe%step == 0 ? 0 : 1); //the number of frames that have to be analysed 

	size_t coords_size = 3*nato*sizeof(float); //memory size for coords of a single frame in one dimension
	size_t memsize= frames * coords_size; //memory size for the array of coords of all frames
	size_t cmemsize= frames * sizeof(int); //memory size for the frameapp array
	size_t dmemsize= frames * sizeof(float); //memory size for the distance array
	size_t totalmemsize = memsize + 2*cmemsize + dmemsize;
	
	//float *devPtr_gclust_coordsX;
	//float *devPtr_gclust_coordsY;
	//float *devPtr_gclust_coordsZ;
	float *devPtr_gclust_coords;
	float *devPtr_distance;
	int *devPtr_frameapp1;
	int *devPtr_frameapp2;
			
	int deviceCount; // number of devices, i.e. gpus
	int device;
	int threadsPerBlock;
	struct cudaDeviceProp properties;		
	cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
	if (cudaResultCode != cudaSuccess)
		deviceCount = 0;
		
	fprintf(stderr,"Starting GPU calculation, devicecount : %d\n", deviceCount);
		
	// machines with no GPUs can still report one emulation device 	
	for (device = 0; device < deviceCount; ++device) {
		cudaGetDeviceProperties(&properties, device);
		if (properties.major != 9999) // 9999 means emulation only
			if (device==0){
				fprintf(stderr,"multiProcessorCount %d\n",properties.multiProcessorCount);
				fprintf(stderr,"maxThreadsPerMultiProcessor %d\n",properties.maxThreadsPerMultiProcessor);
				
				if(properties.major == 2)
					threadsPerBlock = 448;
				else
					threadsPerBlock = 320;
			}
	}
	
	fprintf(stderr,"threads per block: %d\n",threadsPerBlock);
	
	size_t freemem;
	size_t total;
		
	//cuda API functions always return some type of error, but if no error occured, this error is just a cudaSuccess
	//errorHandler terminates program in case there was no cudaSuccess reported
	errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
	
		
	fprintf(stderr,"Available memory on device: %u\n Total memory necessary on device for calculation: %u\n",freemem,totalmemsize);
		

	/*if(freemem < totalmemsize) {
		
				//DEBUG fprintf(stderr,"Available memory on device: %u\n Total memory necessary on device for calculation: %u\n .. splitting up calculation\n",freemem,totalmemsize);
		
		int framesRemaining = frames; //the number of frames that still have to be analysed
		int framesFinished = 0;
		int nclusters = 0; //the number of clusters already found
		int newClusters = 0;
		int *cluster = (int*)calloc(frames,sizeof(int)); //clustercenters of existing clusters
		int *devPtr_cluster; //stores the clustercenter of new found clusters
		int *devPtr_newClusters; //stores the number of new found clusters
		int nframes;	

		while(framesRemaining > 0) {
				
			size_t clust_coords_mem = nclusters * coords_size; //additional memory for the clusters' distance matrices
			errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
			
			//number of frames that fit into memory; 2MB of the total memory reported to freemem have to remain free, allocations fail otherwise (value found by trial and error)
			nframes = (freemem -2000000 - clust_coords_mem - sizeof(int))/(coords_size+3*sizeof(int)+sizeof(float));
			//DEBUG fprintf(stderr,"Free memory: %u, Frames remaining: %d, Frames fitting into memory: %d, Number of clusters: %d\n",freemem,framesRemaining,nframes,nclusters);
			
			//nframes is either the number of frames that fit into gpu memory, or the number of remaining frames
			if(nframes >= framesRemaining) {
				nframes = framesRemaining;
			} else {
				
				//if the number of clusters gets too high the calculation has to be stopped
				if(nclusters > nframes) {
					fprintf(stderr,"Number of clusters has exceeded number of frames that fit on GPU memory, calculation is getting too slow!\n Quitting calculation... Please choose a greater cutoff!\n");
					exit(-1);
				}
			}		
			
			//recalculate the memory sizes
			size_t memsize= nframes * coords_size;
			size_t cmemsize= nframes * sizeof(int);
			size_t dmemsize= nframes * sizeof(float);	
				
			//allocating memory on the GPU
			errorHandler(cudaMalloc((void**)&devPtr_gclust_coords,memsize+clust_coords_mem),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_distance,dmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_frameapp1,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_frameapp2,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_cluster,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_newClusters,sizeof(int)),__LINE__);
					
			//if there were already clusters found, copy their distance matrices
			if(clust_coords_mem>0){
				//because of the overhead of a single copy instruction we prefer to copy one large data packet over lots of small ones, we use a temporary array for this
				float *clusters_coords;
				clusters_coords=(float *)malloc(clust_coords_mem);
				
				for(ii = 0; ii < nclusters; ii++)
					memcpy(clusters_coords + ii*3*nato,gclust_coords + cluster[ii]*3*nato,coords_size);
				
				errorHandler(cudaMemcpy(devPtr_gclust_coords,clusters_coords,clust_coords_mem,cudaMemcpyHostToDevice),__LINE__);
				free(clusters_coords);
			}		
			
			//copy distance matrices of the frames to gpu, they are copied right after the distance matrices of the clusters
			errorHandler(cudaMemcpy(devPtr_gclust_coords + nclusters*3*nato, gclust_coords + (framesFinished + 1)*3*nato, memsize, cudaMemcpyHostToDevice),__LINE__);	
				
			//Sets all indices to -1
			errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
			
			//in order to find the closest cluster we set the distances to the cutoff for the start
			if(!inp_cluster->maxspeed){
				for(ii=framesFinished; ii<=nframes+framesFinished;ii++) {
					distance[ii] = cutoff;
				}
				errorHandler(cudaMemcpy(devPtr_distance, distance+framesFinished+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
			}	
				
			//set number of new clusters to 0
			errorHandler(cudaMemset((void*)devPtr_newClusters,0,sizeof(int)),__LINE__);
					
			int blocks = nframes/threadsPerBlock +1; //in total we want 1 thread for each frame
						
			if(inp_cluster->maxspeed) {			
				//compare the frames of the chunk to the previously found clusters first
				for(ii = 0; ii < nclusters; ii++) {		
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gRmsdClustersMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nclusters, cluster[ii]);
					else gRmsdClustersMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
									
				}
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other
				for(ii=0; ii < nframes; ii++){	
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gRmsdFramesMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords,ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gRmsdFramesMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Calculating Stage %% %f\r",(double)(framesFinished + ii)/frames*100.0);//just a progress bar
					
				}
				
			} else {	
				//compare the frames of the chunk to the previously found clusters first						
				for(ii = 0; ii < nclusters; ii++) {		
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gRmsdClustersClosest<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nclusters, cluster[ii]);
					else gRmsdClustersClosest<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
					
				}
				
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other	
				for(ii=0; ii < nframes; ii++){
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gRmsdFramesClosest<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gRmsdFramesClosest<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Calculating Stage %% %f\r",(double)(framesFinished + ii)/frames*100.0);//just a progress bar
					
				}
				
			}
			
			printf("\n");	
			//DEBUG fprintf(stderr,"Copying to Host ..\n");
			
			//copy back to host, by adding framesFinished/nclusters to the pointers we make sure not to overwrite the results from previous runs
			if((nclusters + ii + 1)%2) errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(distance+framesFinished + 1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(&newClusters, devPtr_newClusters, sizeof(int), cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(cluster+nclusters, devPtr_cluster,newClusters*sizeof(int),cudaMemcpyDeviceToHost),__LINE__);
									
			//update number of clusters and processed frames
			nclusters += newClusters;
			framesFinished += nframes;
			framesRemaining -= nframes;	
			
			//free all the GPU memory
			errorHandler(cudaFree(devPtr_gclust_coords),__LINE__);
			errorHandler(cudaFree(devPtr_distance),__LINE__);
			errorHandler(cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler(cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler(cudaFree(devPtr_cluster),__LINE__);
			errorHandler(cudaFree(devPtr_newClusters),__LINE__);			
		}	
		
			
		return 0;
	} else {*/
			
		//allocate gpu memory
		errorHandler(cudaMalloc((void**)&devPtr_gclust_coords, memsize),__LINE__);
		//errorHandler(cudaMalloc((void**)&devPtr_gclust_coordsX, memsize/3),__LINE__);
		//errorHandler(cudaMalloc((void**)&devPtr_gclust_coordsY, memsize/3),__LINE__);
		//errorHandler(cudaMalloc((void**)&devPtr_gclust_coordsZ, memsize/3),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_distance, dmemsize),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_frameapp1, cmemsize),__LINE__); 
		errorHandler(cudaMalloc((void**)&devPtr_frameapp2, cmemsize),__LINE__); 
						
		//copy coords to gpu
		errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords + 3*nato, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords+frames*nato, gclust_coords[1] + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords+2*frames*nato, gclust_coords[2] + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords + 3*nato, memsize, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coordsX, gclust_coords[0] + nato, memsize/3, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coordsY, gclust_coords[1] + nato, memsize/3, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coordsZ, gclust_coords[2] + nato, memsize/3, cudaMemcpyHostToDevice),__LINE__);
		
		//in order to find the closest cluster we set the distances to the cutoff for the start
		if(!inp_cluster->maxspeed){
			for(ii=0; ii<=frames;ii++) {
					distance[ii] = cutoff;
				}
			errorHandler(cudaMemcpy(devPtr_distance, distance+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
		}
			
		//Set all indices to -1
		errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
				
	
		int blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame
		
		if(super) {
			
			shiftToCenter<<<blocks, threadsPerBlock>>>(devPtr_gclust_coords, nato, frames);
			errorHandler( cudaPeekAtLastError(),__LINE__);
		
					
			if(inp_cluster->maxspeed) {
				for(ii=0;ii< frames;ii++){
							
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					//if((ii+1)%2) gRmsdMax<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coordsX, devPtr_gclust_coordsY, devPtr_gclust_coordsZ, ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, super);
					//else gRmsdMax<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coordsX, devPtr_gclust_coordsY, devPtr_gclust_coordsZ, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, super);
					
					if((ii+1)%2) gRmsdMaxSuper<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance);
					else gRmsdMaxSuper<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance);
							
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
							
					}
						
				} else {
					for(ii=0;ii< frames;ii++){
							
						//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
						//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
						//if((ii+1)%2) gRmsdClosest<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance);
						//else gRmsdClosest<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance);
							
						errorHandler( cudaPeekAtLastError(),__LINE__);
						fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
							
					}
					
				}
			
		}	
			printf("\n");
				
			//DEBUG fprintf(stderr,"Copying results to Host ..\n");
			
			if((ii+1)%2) errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
				
			errorHandler( cudaMemcpy(distance+1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
				
			//free GPU memory
			errorHandler( cudaFree(devPtr_gclust_coords),__LINE__);
			//errorHandler( cudaFree(devPtr_gclust_coordsX),__LINE__);
			//errorHandler( cudaFree(devPtr_gclust_coordsY),__LINE__);
			//errorHandler( cudaFree(devPtr_gclust_coordsZ),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler( cudaFree(devPtr_distance),__LINE__);		
			return 0;
		
	//}
		
}
	

// the CUDA compiler generates C++ object files, thus the main procedure has to be an extern "C" for usage in wordom
extern "C" int gClusterDrms (struct inp_Cluster *inp_cluster,float *distance)
{
    int ii;
    float cutoff = inp_cluster->threshold;
   	int msize = inp_cluster->msize;
	int totframe = inp_cluster->totframe;
	float *gclust_dmtx = inp_cluster->gclust_dmtx;
    float nointrasegm_corr_fact = 1.0;
    int *frameapp = inp_cluster->frameapp;
    int step = inp_cluster->step;
    int frames = totframe/step+(totframe%step == 0 ? 0 : 1); //the number of frames that have to be analysed 
     
    //change correction factor if correction should be applied
	if( inp_cluster->nointrasegm != 0)
		nointrasegm_corr_fact = inp_cluster->nointrasegm_corr_fact;

	size_t dmtx_size = msize*sizeof(float); //memory size for a single distance matrix
	size_t memsize= frames * dmtx_size; //memory size for the array of distance matrices
	size_t cmemsize= frames * sizeof(int); //memory size for the frameapp array
	size_t dmemsize= frames * sizeof(float); //memory size for the distance array
	size_t totalmemsize = memsize + cmemsize + dmemsize;
	
	float *devPtr_gclust_dmtx;
	float *devPtr_distance;
	int *devPtr_frameapp1;
	int *devPtr_frameapp2;
			
	int deviceCount; // number of devices, i.e. gpus
	int device;
	int threadsPerBlock;
	struct cudaDeviceProp properties;		
	cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
	if (cudaResultCode != cudaSuccess)
		deviceCount = 0;
		
	fprintf(stderr,"Starting GPU calculation, devicecount : %d\n", deviceCount);
		
	// machines with no GPUs can still report one emulation device 	
	for (device = 0; device < deviceCount; ++device) {
		cudaGetDeviceProperties(&properties, device);
		if (properties.major != 9999) // 9999 means emulation only
			if (device==0){
				fprintf(stderr,"multiProcessorCount %d\n",properties.multiProcessorCount);
				fprintf(stderr,"maxThreadsPerMultiProcessor %d\n",properties.maxThreadsPerMultiProcessor);
				
				if(properties.major == 2)
					threadsPerBlock = 192;
				else
					threadsPerBlock = 256;
			}
	}
	size_t freemem;
	size_t total;
		
	//cuda API functions always return some type of error, but if no error occured, this error is just a cudaSuccess
	//errorHandler terminates program in case there was no cudaSuccess reported
	errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
	
	//check if there is enough gpu memory for the job and split up the calculation if not
	if(freemem < totalmemsize) {		
		//DEBUG fprintf(stderr,"Available memory on device: %u\n Total memory necessary on device for calculation: %u\n .. splitting up calculation\n",freemem,totalmemsize);
		
		int framesRemaining = frames; //the number of frames that still have to be analysed
		int framesFinished = 0;
		int nclusters = 0; //the number of clusters already found
		int newClusters = 0;
		int *cluster = (int*)calloc(frames,sizeof(int)); //clustercenters of existing clusters
		int *devPtr_cluster; //stores the clustercenter of new found clusters
		int *devPtr_newClusters; //stores the number of new found clusters
		int nframes;	

		while(framesRemaining > 0) {
				
			size_t clust_dmtx_mem = nclusters * dmtx_size; //additional memory for the clusters' distance matrices
			errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
			
			//number of frames that fit into memory; 2MB of the total memory reported to freemem have to remain free, allocations fail otherwise (value found by trial and error)
			nframes = (freemem -2000000 - clust_dmtx_mem - sizeof(int))/(dmtx_size+3*sizeof(int)+sizeof(float));
			//DEBUG fprintf(stderr,"Free memory: %u, Frames remaining: %d, Frames fitting into memory: %d, Number of clusters: %d\n",freemem,framesRemaining,nframes,nclusters);
			
			//nframes is either the number of frames that fit into gpu memory, or the number of remaining frames
			if(nframes >= framesRemaining) {
				nframes = framesRemaining;
			} else {
				
				//if the number of clusters gets too high the calculation has to be stopped
				if(nclusters > nframes) {
					fprintf(stderr,"Number of clusters has exceeded number of frames that fit on GPU memory, calculation is getting too slow!\n Quitting calculation... Please choose a greater cutoff!\n");
					exit(-1);
				}
			}		
			
			//recalculate the memory sizes
			size_t memsize= nframes * dmtx_size;
			size_t cmemsize= nframes * sizeof(int);
			size_t dmemsize= nframes * sizeof(float);	
				
			//allocating memory on the GPU
			errorHandler(cudaMalloc((void**)&devPtr_gclust_dmtx,memsize+clust_dmtx_mem),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_distance,dmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_frameapp1,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_frameapp2,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_cluster,cmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_newClusters,sizeof(int)),__LINE__);
					
			//if there were already clusters found, copy their distance matrices
			if(clust_dmtx_mem>0){
				//because of the overhead of a single copy instruction we prefer to copy one large data packet over lots of small ones, we use a temporary array for this
				float *clusters_dmtx;
				clusters_dmtx=(float *)malloc(clust_dmtx_mem);
				
				for(ii = 0; ii < nclusters; ii++)
					memcpy(clusters_dmtx + ii*msize,gclust_dmtx + cluster[ii]*msize,dmtx_size);
				
				errorHandler(cudaMemcpy(devPtr_gclust_dmtx,clusters_dmtx,clust_dmtx_mem,cudaMemcpyHostToDevice),__LINE__);
				free(clusters_dmtx);
			}		
			
			//copy distance matrices of the frames to gpu, they are copied right after the distance matrices of the clusters
			errorHandler(cudaMemcpy(devPtr_gclust_dmtx + nclusters*msize, gclust_dmtx + (framesFinished + 1)*msize, memsize, cudaMemcpyHostToDevice),__LINE__);	
				
			//Sets all indices to -1
			errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
			
			//in order to find the closest cluster we set the distances to the cutoff for the start
			if(!inp_cluster->maxspeed){
				for(ii=framesFinished; ii<=nframes+framesFinished;ii++) {
					distance[ii] = cutoff;
				}
				errorHandler(cudaMemcpy(devPtr_distance, distance+framesFinished+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
			}	
				
			//set number of new clusters to 0
			errorHandler(cudaMemset((void*)devPtr_newClusters,0,sizeof(int)),__LINE__);
					
			int blocks = nframes/threadsPerBlock +1; //in total we want 1 thread for each frame
						
			if(inp_cluster->maxspeed) {			
				//compare the frames of the chunk to the previously found clusters first
				for(ii = 0; ii < nclusters; ii++) {		
							
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gDrmsClustersMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					else gDrmsClustersMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
									
				}
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other
				for(ii=0; ii < nframes; ii++){	
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gDrmsFramesMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx,ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gDrmsFramesMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Calculating Stage %% %f\r",(double)(framesFinished + ii)/frames*100.0);//just a progress bar
					
				}
				
			} else {	
				//compare the frames of the chunk to the previously found clusters first						
				for(ii = 0; ii < nclusters; ii++) {		
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gDrmsClustersClosest<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					else gDrmsClustersClosest<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
					
				}
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other	
				for(ii=0; ii < nframes; ii++){
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gDrmsFramesClosest<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gDrmsFramesClosest<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Calculating Stage %% %f\r",(double)(framesFinished + ii)/frames*100.0);//just a progress bar
					
				}
			}
			
			printf("\n");	
			//DEBUG fprintf(stderr,"Copying to Host ..\n");
			
			//copy back to host, by adding framesFinished/nclusters to the pointers we make sure not to overwrite the results from previous runs
			if((nclusters + ii + 1)%2) errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(distance+framesFinished + 1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(&newClusters, devPtr_newClusters, sizeof(int), cudaMemcpyDeviceToHost),__LINE__);
			errorHandler( cudaMemcpy(cluster+nclusters, devPtr_cluster,newClusters*sizeof(int),cudaMemcpyDeviceToHost),__LINE__);
									
			//update number of clusters and processed frames
			nclusters += newClusters;
			framesFinished += nframes;
			framesRemaining -= nframes;	
			
			//free all the GPU memory
			errorHandler(cudaFree(devPtr_gclust_dmtx),__LINE__);
			errorHandler(cudaFree(devPtr_distance),__LINE__);
			errorHandler(cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler(cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler(cudaFree(devPtr_cluster),__LINE__);
			errorHandler(cudaFree(devPtr_newClusters),__LINE__);			
		}	
		
			
		return 0;		
		
		} else {
			
			//allocate gpu memory
			errorHandler(cudaMalloc((void**)&devPtr_gclust_dmtx, memsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_distance, dmemsize),__LINE__);
			errorHandler(cudaMalloc((void**)&devPtr_frameapp1, cmemsize),__LINE__); 
			errorHandler(cudaMalloc((void**)&devPtr_frameapp2, cmemsize),__LINE__); 
					
			//copy distance matrices to gpu
			errorHandler(cudaMemcpy(devPtr_gclust_dmtx, gclust_dmtx + msize, memsize, cudaMemcpyHostToDevice),__LINE__);
		
			//Set all indices to -1
			errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
			
			//in order to find the closest cluster we set the distances to the cutoff for the start
			if(!inp_cluster->maxspeed){
				for(ii=0; ii<=frames;ii++) {
					distance[ii] = cutoff;
				}
				errorHandler(cudaMemcpy(devPtr_distance, distance+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
			}

			int blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame
				
			if(inp_cluster->maxspeed) {
				for(ii=0;ii< frames;ii++){
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gDrmsMax<<<blocks, threadsPerBlock>>>(msize, frames, devPtr_gclust_dmtx,ii, cutoff, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact);
					else gDrmsMax<<<blocks, threadsPerBlock>>>(msize, frames, devPtr_gclust_dmtx,ii, cutoff, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
					
				}
				
			} else {
				for(ii=0;ii< frames;ii++){
					
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gDrmsClosest<<<blocks, threadsPerBlock>>>(msize, frames, devPtr_gclust_dmtx, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact);
					else gDrmsClosest<<<blocks, threadsPerBlock>>>(msize, frames, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
					
				}	
			}	
			printf("\n");
			
			//DEBUG fprintf(stderr,"Copying results to Host ..\n");
		
			if((ii+1)%2) errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
			errorHandler( cudaMemcpy(distance+1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
			//free GPU memory
			errorHandler( cudaFree(devPtr_gclust_dmtx),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler( cudaFree(devPtr_distance),__LINE__);		
			return 0;
	   }
	}
