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

__global__ void gRmsdMax(const int msize, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance) {
	
	
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
		
		float rmsd=0. , di;
		int ii;
		
	    
	    // compute rmsd and return value
			  
		for ( ii=0; ii<msize; ii++ ) {
				di=gclust_coords[index*msize+ii]-gclust_coords[cluster*msize+ii];
				rmsd+=di*di;
				/*di=gclust_coords[nframes*nato+index*nato+ii]-gclust_coords[nframes*nato+cluster*nato+ii];
				rmsd+=di*di;
				di=gclust_coords[2*nframes*nato+index*nato+ii]-gclust_coords[2*nframes*nato+cluster*nato+ii];
				rmsd+=di*di;*/
			}
			  
		rmsd /= nato;
		rmsd = sqrt ( rmsd );
	
		if (rmsd<distance[index]){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}
		
	}
	frameapp_write[index] = frameapp_read[index];
			  
			
}

__global__ void gRmsdClosest(const int msize, const int nframes, const float* gclust_coords, const int cluster, const float* frameapp_read, float* frameapp_write, float* distance) {
	
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
		for (jj=0;jj<msize;jj++){
			di=(gclust_coords[index*msize+jj]-gclust_coords[cluster*msize+jj]);
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


__global__ void gRmsdClustersMax(const int msize, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int nclusters, const int clustercenter) {
	

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

	for ( ii=0; ii<msize; ii++ ) {
			di=gclust_coords[(index+nclusters)*msize+ii]-gclust_coords[cluster*msize+ii];
			rmsd+=di*di;
			/*di=gclust_coords[nframes*nato+index*nato+ii]-gclust_coords[nframes*nato+cluster*nato+ii];
			rmsd+=di*di;
			di=gclust_coords[2*nframes*nato+index*nato+ii]-gclust_coords[2*nframes*nato+cluster*nato+ii];
			rmsd+=di*di;*/
		}
			  
		rmsd /= nato;
		rmsd = sqrt ( rmsd );

		if (rmsd<distance[index]){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = rmsd;
			return;
		}
	
	frameapp_write[index] = frameapp_read[index];
}

__global__ void gRmsdFramesMax(const int msize, const int nframes, const float* gclust_coords, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const int framesFinished, const int nclusters, int* newClusters, int* cluster) {
	
	
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
		for (jj=0;jj<msize;jj++){
			di=gclust_coords[(index+nclusters)*msize+jj]-gclust_coords[(cluster+nclusters)*msize+jj]);
			rmsd+=di*di; 
		}
					
		rmsd /= nato;
		rmsd = sqrt ( rmsd );
	
		if (rmsd<distance[index]){		
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
		drms+=di*di; 
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
			di=(gclust_coords[(index+nclusters)*3*nato+jj]-gclust_dmtx[(cluster+nclusters)*3*nato+jj]);
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
__global__ void gDrms(const int max, const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact) {
	
	// max determines to run either the maxspeed (=1) or closest cluster search (=0)
	// msize = size of each distance matrix, nframes = number of frames used for clustering;
	// gclust_dmtx = array of ALL distance matrices; cluster = center of a possible cluster; frameapp stores the cluster center 
	// of each frame. There is one frameapp array for reading only and one for writing only, because otherwise this could go wrong for very long arrays;
	// distance stores the calculated distance for comparison and later usage
	

	// threadIdx.x is a built-in variable provided by CUDA at runtime, this gives each of the threads a unique index that corresponds to a frame number 
	int index = blockIdx.x * blockDim.x + threadIdx.x;	
	// if there are more threads than frames then stop these
	if (index>=nframes){
		return;
	}
	
	if(max) {
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
		if (drms<distance[index]){		
			frameapp_write[index] = cluster + 1; //+1 because in wordom the frames are counted starting with 1
			distance[index] = drms;
			return;
		}
		
	}
	frameapp_write[index] = frameapp_read[index];

}

// a kernel for the maxspeed flag for calculation with limited memory; this is for comparing the frames of the chunk to the previously found clusters
__global__ void gDrmsClustersMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int nclusters, const int clustercenter) {
	
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

	if (drms<distance[index]){		
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
__global__ void gDrmsFramesMax(const int msize, const int nframes, const float* gclust_dmtx, const int cluster, const int* frameapp_read, int* frameapp_write, float* distance, const float nointrasegm_corr_fact, const int framesFinished, const int nclusters, int* newClusters, int* clusterCenters ) {
	
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
	
		if (drms<distance[index]){		
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
    int step = inp_cluster->step;
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
	
		
	fprintf(stderr,"Available memory on device: %u\n Total memory necessary on device for calculation: %u\n",freemem,totalmemsize);
		

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
				
			size_t clust_coords_mem = nclusters * coords_size; //additional memory for the clusters' distance matrices
			errorHandler(cudaMemGetInfo(&freemem, &total),__LINE__);
			
			//number of frames that fit into memory; 2MB of the total memory reported to freemem have to remain free, allocations fail otherwise (value found by trial and error)
			nframes = (freemem -2000000 - clust_dmtx_mem - sizeof(int))/(coords_size+3*sizeof(int)+sizeof(float));
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
				free(clusters_dmtx);
			}		
			
			//copy distance matrices of the frames to gpu, they are copied right after the distance matrices of the clusters
			errorHandler(cudaMemcpy(devPtr_gclust_coords + nclusters*3*nato, gclust_coords + (framesFinished + 1)*3*nato, memsize, cudaMemcpyHostToDevice),__LINE__);	
				
			//Sets all indices to -1
			errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
			
			//setting the distances to the cutoff for comparison on the gpu
			for(ii=framesFinished; ii<=nframes+framesFinished;ii++) {
				distance[ii] = cutoff;
			}
			errorHandler(cudaMemcpy(devPtr_distance, distance+framesFinished+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
	
				
			//set number of new clusters to 0
			errorHandler(cudaMemset((void*)devPtr_newClusters,0,sizeof(int)),__LINE__);
					
			int blocks = nframes/threadsPerBlock +1; //in total we want 1 thread for each frame
						
			if(inp_cluster->maxspeed) {			
				//compare the frames of the chunk to the previously found clusters first
				for(ii = 0; ii < nclusters; ii++) {		
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gRmsdClustersMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nclusters, cluster[ii]);
					else gRmsdClustersMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
									
				}
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other
				for(ii=0; ii < nframes; ii++){	
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gRmsdFramesMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords,ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gRmsdFramesMax<<<blocks, threadsPerBlock>>>(nato, nframes, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
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
	} else {
			
		//allocate gpu memory
		errorHandler(cudaMalloc((void**)&devPtr_gclust_coords, memsize),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_distance, dmemsize),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_frameapp1, cmemsize),__LINE__); 
		errorHandler(cudaMalloc((void**)&devPtr_frameapp2, cmemsize),__LINE__); 
						
		//copy coords to gpu
		//errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords[0] + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords+frames*nato, gclust_coords[1] + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		//errorHandler(cudaMemcpy(devPtr_gclust_coords+2*frames*nato, gclust_coords[2] + coords_size, frames*coords_size, cudaMemcpyHostToDevice),__LINE__);
		errorHandler(cudaMemcpy(devPtr_gclust_coords, gclust_coords + 3*nato, memsize, cudaMemcpyHostToDevice),__LINE__);
		
		//setting the distances to the cutoff for comparison on the gpu
		for(ii=0; ii<=frames;ii++) {
				distance[ii] = cutoff;
			}
		errorHandler(cudaMemcpy(devPtr_distance, distance+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);

			
		//Set all indices to -1
		errorHandler(cudaMemset((void*)devPtr_frameapp1,-1,cmemsize),__LINE__);
				
	
		int blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame
		
					
		if(inp_cluster->maxspeed) {
			for(ii=0;ii< frames;ii++){
						
				//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
				//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
				if((ii+1)%2) gRmsdMax<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords,ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance);
				else gRmsdMax<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords,ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance);
						
				errorHandler( cudaPeekAtLastError(),__LINE__);
				fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
						
				}
					
			} else {
				for(ii=0;ii< frames;ii++){
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gRmsdClosest<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance);
					else gRmsdClosest<<<blocks, threadsPerBlock>>>(nato, frames, devPtr_gclust_coords, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance);
						
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
			errorHandler( cudaFree(devPtr_gclust_coords),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp1),__LINE__);
			errorHandler( cudaFree(devPtr_frameapp2),__LINE__);
			errorHandler( cudaFree(devPtr_distance),__LINE__);		
			return 0;
		
	}
		
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
			
			//setting the distances to the cutoff for comparison on the gpu
			for(ii=framesFinished; ii<=nframes+framesFinished;ii++) {
				distance[ii] = cutoff;
			}
			errorHandler(cudaMemcpy(devPtr_distance, distance+framesFinished+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);
				
			//set number of new clusters to 0
			errorHandler(cudaMemset((void*)devPtr_newClusters,0,sizeof(int)),__LINE__);
					
			int blocks = nframes/threadsPerBlock +1; //in total we want 1 thread for each frame
						
			if(inp_cluster->maxspeed) {			
				//compare the frames of the chunk to the previously found clusters first
				for(ii = 0; ii < nclusters; ii++) {		
							
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((ii+1)%2) gDrmsClustersMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					else gDrmsClustersMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, nclusters, cluster[ii]);
					
					errorHandler( cudaPeekAtLastError(),__LINE__);
					fprintf(stderr,"Comparing to previous clusters %% %f\r",(double)ii/nclusters*100.0);//just a progress bar
									
				}
				
				fprintf(stderr,"\n");
				
				//then check the remaining frames against each other
				for(ii=0; ii < nframes; ii++){	
						
					//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
					//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
					if((nclusters + ii + 1)%2) gDrmsFramesMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx,ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					else gDrmsFramesMax<<<blocks, threadsPerBlock>>>(msize, nframes, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact, framesFinished, nclusters, devPtr_newClusters, devPtr_cluster);
					
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
			
			//setting the distances to the cutoff for comparison on the gpu
			for(ii=0; ii<=frames;ii++) {
				distance[ii] = cutoff;
			}
			errorHandler(cudaMemcpy(devPtr_distance, distance+1, dmemsize, cudaMemcpyHostToDevice),__LINE__);

			int blocks = frames/threadsPerBlock +1; //in total we want 1 thread for each frame
				

			for(ii=0;ii< frames;ii++){
					
				//the kernel ensures that frameapp_read has been written to frameapp_write entirely after one iteration
				//to prevent wasting time on copying frameapp_write back to frameapp_read the kernel simply gets launched with both interchanged in every second iteration
				if((ii+1)%2) gDrms<<<blocks, threadsPerBlock>>>(inp_cluster->maxspeed,msize, frames, devPtr_gclust_dmtx, ii, devPtr_frameapp1, devPtr_frameapp2, devPtr_distance, nointrasegm_corr_fact);
				else gDrms<<<blocks, threadsPerBlock>>>(inp_cluster->maxspeed,msize, frames, devPtr_gclust_dmtx, ii, devPtr_frameapp2, devPtr_frameapp1, devPtr_distance, nointrasegm_corr_fact);
					
				errorHandler( cudaPeekAtLastError(),__LINE__);
				fprintf(stderr,"Stage %% %f\r",(double)ii/frames*100.0);//just a progress bar
					
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
