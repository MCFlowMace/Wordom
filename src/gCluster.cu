// ------------------------------------------------------------------
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ------------------------------------------------------------------


#include <math.h>
#include <stdio.h> 
#include <stdlib.h>

#define GCLUSTER_INCLUDES

#include "fileio.h"
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
 * given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
__device__ void setup_rotation(const float *ref_xlist, const float *mov_xlist, int& n_list, float R[3][3], float& E0, volatile int& i, volatile int&j, volatile int& n)
{
  // initialize
  for (i=0; i<3; i++)
    for (j=0; j<3; j++) 
      R[i][j] = 0.0f;
  E0 = 0.0f;

  for (n=0; n<n_list; n++) 
  {
    /* 
     * E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) 
     */
    for (i=0; i<3; i++)
      E0 +=  mov_xlist[3*n+i] * mov_xlist[3*n+i]  
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
  E0 *= 0.5f;
  }



#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }
/*   
 * jacobi3
 *
 * computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
__device__ int jacobi3(float a[3][3], float d[3], float v[3][3], volatile int& k, volatile int& i, volatile int& j)
{
  float b[3], z[3];

  //Initialize v to the identity matrix.
  for (i=0; i<3; i++) 
  { 
    for (j=0; j<3; j++) 
      v[i][j] = 0.0f;
    v[i][i] = 1.0f;
  }

  // Initialize b and d to the diagonal of a
  for (i=0; i<3; i++) 
    b[i] = d[i] = a[i][i];

  // z will accumulate terms
  for (i=0; i<3; i++) 
    z[i] = 0.0f; 

  // 50 tries
  int count;
  for (count=0; count<50; count++)     
  {	
	float tresh;
    // sum off-diagonal elements
    {
	    float sum = 0.0f;
	    for (i=0; i<2; i++) 
	      for (j=i+1; j<3; j++)
	         sum += fabsf(a[i][j]);
	
	    // if converged to machine underflow
	    if (sum == 0.0f) 
	      return(1);
	
	    // on 1st three sweeps...
	    if (count < 3) 
	      tresh = sum * 0.2f / 9.0f;    
	    else       
	      tresh = 0.0f;   
	}

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        float g = 100.0f * fabsf(a[i][j]);

        /*  after four sweeps, skip the rotation if
         *   the off-diagonal element is small 
         */
       if ( count > 3  &&  fabsf(d[i])+g == fabsf(d[i]) &&  fabsf(d[j])+g == fabsf(d[j]) ) 
        {
          a[i][j] = 0.0f;
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
            float theta = 0.5f * h / (a[i][j]);
            t = 1.0f / ( fabsf(theta) +
                        (float)sqrtf(1.0f + theta*theta) );
            if (theta < 0.0f) 
              t = -t;
          }
          
          float c = 1.0f / (float) sqrtf(1.0f + t*t);
          float s = t * c;
          float tau = s / (1.0f + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0f;

          for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)

        }
      }
    }

    for (i=0; i<3; i++) 
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0f;
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
  
	if (!jacobi3(matrix, eigenval, vec, i,j,k))
		return (0);
	
	// sort solutions by eigenval
	for (i=0; i<3; i++) 
	{
		k = i;
		matrix[0][0] = eigenval[i];
		
		for (j=i+1; j<3; j++)
			if (eigenval[j] >= matrix[0][0])
			{ 
			k = j;
			matrix[0][0] = eigenval[k];
			}
		   
		if (k != i) 
		{
			eigenval[k] = eigenval[i];
			eigenval[i] = matrix[0][0];
			for (j=0; j<3; j++) 
			{
				matrix[0][0] = vec[j][i];
				vec[j][i] = vec[j][k];
				vec[j][k] = matrix[0][0];
			}
		}
	}
	
	// transpose such that first index refers to solution index
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			vec[i][j] = vec[j][i];
	
	return (1);
}



/*
 * calculate_rotation() 
 *
 * calculates the residual from the R matrix and E0:
 * to reduce the number of used variables E0 has to be passed as the value of residual
 */

__device__ int calculate_rotation(float R[3][3], float& residual, volatile int& i, volatile int& j, volatile int& k)
{
	// build Rt, transpose of R 
	float Rt[3][3];
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			Rt[i][j] = R[j][i];

	// make symmetric R = Rt X R
	float right_eigenvec[3][3];
		for (i=0; i<3; i++) 
			for (j=0; j<3; j++)
			{
				right_eigenvec[i][j] = 0.0f;
				for (k = 0; k<3; k++)
					right_eigenvec[i][j] += Rt[k][i] * R[j][k];
		    }
    
	for (i=0; i<3; i++) 
		for (j=0; j<3; j++)
			R[i][j] = right_eigenvec[i][j];

	float eigenval[3];
	if (!diagonalize_symmetric(R, right_eigenvec, eigenval,i,j,k))
	    return(0);
  

	/* right_eigenvec's should be an orthogonal system but could be left
	* or right-handed. Let's force into right-handed system.
	*/
	cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);
	
	/* From the Kabsch algorithm, the eigenvec's of RtR
	* are identical to the right_eigenvec's of R.
	* This means that left_eigenvec = R x right_eigenvec 
	*/	

	for (i=0; i<3; i++) 
		for (j=0; j<3; j++) 
			R[i][j] = dot(&right_eigenvec[i][0], &Rt[j][0]);
      	
	for (i=0; i<3; i++) 
		normalize(&R[i][0]);

	/* 
	* Force left_eigenvec[2] to be orthogonal to the other vectors.
	* First check if the rotational matrices generated from the 
	* orthogonal eigenvectors are in a right-handed or left-handed
	* co-ordinate system - given by sigma. Sigma is needed to
	* resolve this ambiguity in calculating the RMSD.
	*/

	cross(&right_eigenvec[0][0], &R[0][0], &R[1][0]);
	
	float sigma;
	if (dot(&right_eigenvec[0][0], &R[2][0]) < 0.0)
		sigma = -1.0f;
	
	else 
		sigma = 1.0f;
	
	residual = residual - (float) sqrtf(fabsf(eigenval[0])) - (float) sqrtf(fabsf(eigenval[1]))- sigma * (float) sqrtf(fabsf(eigenval[2]));	
	return (1);
}

__global__ void rmsd_kernel(int nato, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int mode) {
	
	/* nato = number of atoms, nframes = number of frames used for clustering;
	 * gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated rmsd because it will be needed in the post processing
	 * mode is used to switch between maxspeed=1 and closest=0
	 */
	
	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	if(mode) {
		//with the maxspeed flag frames are assigned to the first cluster -> do not check frame again if it was already assigned
		if(frameapp_read[index] != -1) {
			frameapp_write[index] = frameapp_read[index];
			return;
		}
	} else {
		//frames can only be assigned to clusters with a lower frame number
		if(cluster > index) {
			frameapp_write[index] = frameapp_read[index];
			return;
		}	
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {

		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0f;
			return;
		}
		
		int ii;
		float rmsd,di;
  
		rmsd=0.0f;
		di=0.0f;
  
		for ( ii=0; ii<3*nato; ii++ ) {
			di= gclust_coords[cluster*3*nato+ii]-gclust_coords[index*3*nato+ii];
			rmsd += di*di;            
		}
  
		rmsd /= nato;
		rmsd = sqrtf ( rmsd );
	
		if (rmsd<distance[index]){
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}	
	}
	
	frameapp_write[index] = frameapp_read[index];
}

__global__ void rmsd_super_kernel(int nato, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int mode) {
	
	
	/* nato = number of atoms, nframes = number of frames used for clustering;
	 * gclust_coords = array of all coordinates of the frames; cluster = center of a possible cluster; frameapp stores the cluster center 
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated rmsd because it will be needed in the post processing
	 * mode is used to switch between maxspeed=1 and closest=0 
	 */
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	if(mode) {
		//with the maxspeed flag frames are assigned to the first cluster -> do not check frame again if it was already assigned
		if(frameapp_read[index] != -1) {
			frameapp_write[index] = frameapp_read[index];
			return;
		}
	} else {
		//frames can only be assigned to clusters with a lower frame number
		if(cluster > index) {
			frameapp_write[index] = frameapp_read[index];
			return;
		}	
	}
	
	//if cluster has already been assigned to another cluster then it cannot be a new cluster
	if(frameapp_read[cluster] == -1) {

		//spare the distance calculation if index and cluster are the same, by wordoms conventions the distance is set to -1 instead of 0 in this case
		if(cluster == index) {
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = -1.0f;
			return;
		}
		
		float rmsd;
		float R[3][3];
		volatile int i,j,n;
		
		setup_rotation(&gclust_coords[cluster*3*nato],&gclust_coords[index*3*nato], nato, R, rmsd, i, j, n);
		
		if(calculate_rotation(R, rmsd, i, j, n)) {
  
			rmsd = fabsf(rmsd); // avoids the awkward case of -0.0 
			rmsd = sqrtf( fabsf((float) (rmsd)*2.0f/((float)nato)) );
			
			if (rmsd<distance[index]){
				frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
				distance[index] = rmsd;
				return;
			}
		}		
	}
	
	frameapp_write[index] = frameapp_read[index];			
}

// a kernel for the maxspeed flag
__global__ void gDrmsMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const float cutoff, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact) {
	
	/* msize = size of each distance matrix, nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated drms because it will be needed in the post processing 
	 */
	
	

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
	
	/* msize = size of each distance matrix; nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated drms for comparison 
	 */

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
	
	/* msize = size of each distance matrix; nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = number of the cluster; frameapp stores the cluster center 
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	 * distance stores the calculated drms because it will be needed in the post processing;
	 * nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster 
	 */

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
	
	/* msize = size of each distance matrix; nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = number of the cluster; frameapp stores the cluster center 
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays; 
	 * distance stores the calculated for comparison;
	 * nclusters passes the number of already found clusters; clustercenter passes the center of the current cluster 
	 */

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
	
	/* msize = size of each distance matrix; nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated drms because it will be needed in the post processing;
	 * framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	 * newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters 
	 */
	

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
	
	/* msize = size of each distance matrix; nframes = number of frames used for clustering;
	 * gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	 * of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	 * distance stores the calculated drms for comparison;
	 * framesFinished passes the number of processed frames in a previous chunk; nclusters passes the number of already found clusters;
	 * newClusters stores the number of new found clusters in this chunk; clusterCenters stores the centers of these clusters
	 */
	
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

//shifts the center of mass to the origin
__global__ void shiftToCenter(float* gclust_coords, const int nato, const int nframes) {
  
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
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

//wrapper function for handling CUDA errors
void errorHandler  (cudaError_t error, int line){
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    fprintf(stderr,"CUDA error: %s in line number %d\n", cudaGetErrorString(error),line);
    exit(-1);
  }
}

extern "C" int find_GPUs() {
	
	int deviceCount, device, realGPUs;
	struct cudaDeviceProp prop;
	cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
	if (cudaResultCode != cudaSuccess)
		deviceCount = 0;
		
	// look for the number of real GPUs not counting CUDA emulation devices
	for (device = 0; device < deviceCount; ++device) {
		errorHandler(cudaGetDeviceProperties(&prop, device),__LINE__);
		if (prop.major != 9999) { // 9999 is an emulation device
			realGPUs++;
			fprintf(stderr,"Device %d:\n",device);
			fprintf(stderr,"%s @%dMHz\n",prop.name,prop.clockRate/1000);
			fprintf(stderr,"Total memory: %dMBytes\n",prop.totalGlobalMem/(1000*1000));
			fprintf(stderr,"\n");
		}
	}
	
	return realGPUs;
}

// the CUDA compiler generates C++ object files, thus the main procedure has to be an extern "C" for usage in wordom
extern "C" int gClusterRmsd (struct inp_Cluster *inp_cluster,float *distance) {
	
	int ii;
	float cutoff = inp_cluster->threshold;
	int nato = inp_cluster->nato;
	int totframe = inp_cluster->totframe;
	float *gclust_coords = inp_cluster->gclust_coords;
	int *frameapp = inp_cluster->frameapp;
	int super = inp_cluster->super;
	int step = inp_cluster->step;
	int device = inp_cluster->device;
	int frames = totframe/step+(totframe%step == 0 ? 0 : 1); //the number of frames that have to be analysed 
	
	size_t coords_size = 3*nato*sizeof(float); //memory size for coords of a single frame in one dimension
	size_t memsize= frames * coords_size; //memory size for the array of coords of all frames
	size_t cmemsize= frames * sizeof(int); //memory size for the frameapp array
	size_t dmemsize= frames * sizeof(float); //memory size for the distance array
	size_t totalmemsize = memsize + 2*cmemsize + dmemsize;
	
	float *devPtr_gclust_coords;
	float *devPtr_distance;
	int *devPtr_frameapp1;
	int *devPtr_frameapp2;
			
	int threadsPerBlock;
	int blocks;
	struct cudaDeviceProp properties;
	errorHandler(cudaSetDevice(device),__LINE__);
	errorHandler(cudaGetDeviceProperties(&properties, device),__LINE__);
		
	fprintf(stderr,"Starting GPU calculation\n");
	
	if(properties.kernelExecTimeoutEnabled)
		fprintf(stderr,"WARNING! The GPU you are using was set up with a run time limit for kernels. Therefore it is likely that the calculation of large proteins fails!\n");
		
	// machines with no GPUs can still report one emulation device 	
	/*for (device = 0; device < deviceCount; ++device) {
		cudaGetDeviceProperties(&properties, device);
		if (properties.major != 9999) // 9999 means emulation only
			if (device==0){
				fprintf(stderr,"multiProcessorCount %d\n",properties.multiProcessorCount);
				fprintf(stderr,"maxThreadsPerMultiProcessor %d\n",properties.maxThreadsPerMultiProcessor);
			}
	}*/
	
	//fprintf(stderr,"threads per block: %d\n",threadsPerBlock);
	
	size_t freemem;
	size_t total;
		
	//cuda API functions always return some type of error, but if no error occured, this error is just a cudaSuccess
	//errorHandler terminates program in case there was no cudaSuccess reported
	errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
	if(freemem < totalmemsize) {
		fprintf(stderr,"Available graphics memory: %8.3f MBytes\nRequired memory for calculation: %8.3f MBytes\n",(float)freemem/1000000,(float)totalmemsize/1000000);
		fprintf(stderr,"Terminating calculation. Maybe use the --STEP option?\n");
		exit(-1);
	}

	//allocate gpu memory
	errorHandler(cudaMalloc((void**)&devPtr_gclust_coords, memsize),__LINE__);
	errorHandler(cudaMalloc((void**)&devPtr_distance, dmemsize),__LINE__);
	errorHandler(cudaMalloc((void**)&devPtr_frameapp1, cmemsize),__LINE__);
	errorHandler(cudaMalloc((void**)&devPtr_frameapp2, cmemsize),__LINE__);
	
	//copy coords to gpu
	errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords + 3*nato, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);

	//distances are set to the cutoff for the start
	for(ii=0; ii<=frames;ii++) {
			distance[ii] = cutoff;
		}
	errorHandler(cudaMemcpy(devPtr_distance, distance+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);

	//Set all indices to -1
	errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);

	threadsPerBlock = 256; //tuned for 100% occupancy with CUDA occupancy calculator
	blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame

	if(super) {

		//shit the center of mass to the origin first for all frames
		shiftToCenter<<<blocks, threadsPerBlock>>>(devPtr_gclust_coords, nato, frames);
		errorHandler( cudaPeekAtLastError(),__LINE__);
	
		//change block size for the main kernel function for 100% occupancy again
		if(properties.major == 2)
			threadsPerBlock = 448;
		else
			threadsPerBlock = 256;
			
		blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame
				
		for(ii=0;ii< frames;ii++){
						
			//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
			//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
				
			if((ii+1)%2) rmsd_super_kernel<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, inp_cluster->maxspeed);
			else rmsd_super_kernel<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, inp_cluster->maxspeed);
						
			errorHandler( cudaPeekAtLastError(),__LINE__);
			fprintf(stderr,"Stage %% %f\r",(float)ii/frames*100.0);//just a progress bar
						
		}
		
		
	} else {
		
		for(ii=0;ii< frames;ii++){
						
			//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
			//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
				
			if((ii+1)%2) rmsd_kernel<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, inp_cluster->maxspeed);
			else rmsd_kernel<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, inp_cluster->maxspeed);
						
			errorHandler( cudaPeekAtLastError(),__LINE__);
			fprintf(stderr,"Stage %% %f\r",(float)ii/frames*100.0);//just a progress bar
						
		}
	}	
		printf("\n");
			
		//DEBUG fprintf(stderr,"Copying results to Host ..\n");
		
		//make sure to copy the correct frameapp array back
		if((ii+1)%2) errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
		else errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
		errorHandler( cudaMemcpy(distance+1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
		//free GPU memory
		errorHandler( cudaFree(devPtr_gclust_coords),__LINE__);
		errorHandler( cudaFree(devPtr_frameapp1),__LINE__);
		errorHandler( cudaFree(devPtr_frameapp2),__LINE__);
		errorHandler( cudaFree(devPtr_distance),__LINE__);		
		return 0;		
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
    int device = inp_cluster->device;
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
				
	// machines with no GPUs can still report one emulation device 	
	/*for (device = 0; device < deviceCount; ++device) {
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
	}*/
	
	struct cudaDeviceProp properties;
	errorHandler(cudaSetDevice(device),__LINE__);
	errorHandler(cudaGetDeviceProperties(&properties, device),__LINE__);
		
	fprintf(stderr,"Starting GPU calculation\n");
	
	if(properties.kernelExecTimeoutEnabled)
		fprintf(stderr,"WARNING! The GPU you are using was set up with a run time limit for kernels. Therefore it is likely that the calculation of large proteins fails!\n");
		
	int threadsPerBlock;
	if(properties.major == 2)
		threadsPerBlock = 192;
	else
		threadsPerBlock = 256;
	
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
			nframes = (freemem -10000000 - clust_dmtx_mem - sizeof(int))/(dmtx_size+3*sizeof(int)+sizeof(float));
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
			if((nclusters + ii + 1)%2) errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+framesFinished + 1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
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
		
			if((ii+1)%2) errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp1, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			else errorHandler( cudaMemcpy(frameapp+1, devPtr_frameapp2, cmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
			errorHandler( cudaMemcpy(distance+1, devPtr_distance, dmemsize, cudaMemcpyDeviceToHost),__LINE__);
			
			//free GPU memory
			errorHandler( cudaFree(devPtr_gclust_dmtx),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler( cudaFree(devPtr_distance),__LINE__);		
			return 0;
	   }
	}
