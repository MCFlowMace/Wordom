// ------------------------------------------------------------------
// Copyright (C) 2003  University of Zurich
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
/*! \file cluster.h
 \brief Clustering module
 
 Headers for all clustering procedures
*/
#ifndef CLUSTER
#define CLUSTER

// ------------------------------------------------------------------
// ** Clustering data
typedef struct cluster_
{
   float        threshold;
   int          ncluster;
   int         *cluspop;
   int         *clusstart;
   int         *cluscenter;
   int        **clusstrs;
   float       *rmsdspread;
   float       *cluscenterrmsd;
} _cluster;
// ------------------------------------------------------------------
// ** Clustering data - single cluster
typedef struct scluster_
{
   int          cluspop;
   int          clusstart;
   int          cluscenter;
   int         *clusstrs;
   float        rmsdspread;
   float        cluscenterrmsd;
   float       *distmtx;
   float        framedist;       // distance clustercenter-frame, when computing
} _scluster;
// ------------------------------------------------------------------
// ** linked list of cluster members
typedef struct _confs
{
  int             conf_index;
  struct _confs  *next;
} confs;
// ------------------------------------------------------------------
// ** Clustering data - single cluster
typedef struct clusters
{
  int              center;          // Frame of cluster center
  int              nelements;       // Number of elements
  int             *clusconfs;       // conformations in cluster (unused)
  float           *distance;        // Distance Matrix of cluster (x DRMS)
  float          **cluscoors;       // coordinates of cluster center  
  struct clusters *next;            // Pointer to next element of linked list
  struct clusters *last;
  confs           *ll_confs;
  confs           *ll_confs_tail;
  int             *clus_confs;
  // neighbors list added for QT post-optimization
  int                   n_clus_neig;
  struct ll_clus_neig  *llneig;       // linked list of neigbor clusters
  struct ll_clus_neig  *llneiglast;   // pointer to last element of llneig
} _clusters;
// ------------------------------------------------------------------
struct ll_clus_neig
{
  int                     cl_number;
  float                   distance;
  struct ll_clus_neig    *next;
};
// ------------------------------------------------------------------
struct ll_conf_neig
{
  int                   conf_number;
  float                 distance;
  struct ll_conf_neig  *next;
};
// ------------------------------------------------------------------

struct inp_Cluster
{
  char           title[64];
  char           inmol[256];
  char           intrj[256];
  char           header[512];
  FILE          *output;
  FILE			*gOutput;
  _cluster       cluster1;      // structure for clustering output
  int            nclusters;
  _scluster     *cluster;      // structure for each cluster coming out of the output
  Selection      sele;
  Molecule      *mol;
  int            nato;          // equal to sele.nselatm
  int            step;
  float          threshold;
  float			 threshold_bh;	//threshold for the Barnes-Hut approximation in tsne
  float			 perplexity;	//for probability distribution in tsne
  int			 max_iter;		//max iterations for tsne
  int			 dimension;		//reduced dimension for tsne
  int            method;
  int            distance;
  int            totframe;
//int            skip;
  int            begin;
  // know whether to operate from memory or temp file for coor
  short          himem;
  FILE          *tmpfile;
  short          threaded;
  int            nthreads;
  // if multiple clustering or distance-matrix required as output, save it:
  short          outmatrix;
  short          outmatrixonly;
  short          inmatrix;
  char           matrixfilename[256];
  char           tmpfilename[256];
  CoorSet       *tmpfilecoor;
  
  // if dist == rsmd, it can be minimized or not
  int            super;
  
  // in method 3 we can go faster by taking the first cluster < cutoff
  int            maxspeed;
  // even faster yet more accurate by comparing last frames' clusters first
  int            markov;
  
  // if method requires coors to be stored, allocate mem
  float        **xcoor;
  float        **ycoor;
  float        **zcoor;
//float       ***coors;
  float        **matrix;
  int            nframe;
  
  // if fitting to a selection and then rmsd-calc to another is required
  short          fit;
  float        **fitref;
  float        **fitmov;
  Selection       fitsele;

  //if method (3) requires immediate handling of coors, allocate mem
  float        **refcoor;
  //if method (3) would not store which cluster a frame belongs to, take care
  int           *frameapp;
  int           *cluster_to_check;
  
  int            pbcflag;
  float         *pbcbox;
  
  _clusters     *head;
  _clusters     *tail;
  _clusters    **clusters_array;
  _clusters    **clusterlist;
//struct clusters       *head;
  
  // if distance == drms, we need a distmtx for the ii^{th} frame
  int            msize;         // memory size of the distance matrix
  float         *distmtx;       // distance matrix
  float         *dist_mtx;
  int            dist_mtx_size;
  short          nointrasegm;
  float          nointrasegm_corr_fact;
  float         *nointrasegm_msk;
  short          cmapd;
  float          cmapd_cutoff;
  // if using drms and method 3, use a set of coordinates to pass to DistMtxCalc
  float        **tmpcoor;
  // if distance == rmsd
  struct ll_clus_neig  *tmp_llneig;
  struct ll_conf_neig  **llconfneig;
  int                   *conf_nneig;
  // for Gcluster array of distance matrices or array of all coords
  float 		*gclust_dmtx;
  float			*gclust_coords;
  int			 device; //number of device to run on for multi GPU systems
};

// ------------------------------------------------------------------
struct  thrFillMat
{
   int            thread_rank;
   int            thread_size;
   float        **matrix;
   struct inp_Cluster *inp_cluster;
};
// ------------------------------------------------------------------
// ** Cluster Assignment Module structure
typedef struct drmsdata
  {
   int            msize;         // memory size of the distance matrix
   float         *distmtx;       // distance matrix
   float         *dist_mtx;
   int            dist_mtx_size;
  }     _drmsdata;
// ------------------------------------------------------------------
struct inp_CAssign
{
  char           title[64];
  char           cluslistname[64];
  FILE          *output;
  FILE          *cluslist;
  Selection       sele;
  int            nato;          // equal to sele.nselatm
  int            step;
  int            distance;
  int            totframe;
  int            skip;
  int            nframe;
  float          cutoff;
  
  int            nclusters;
  int           *clustercenters;
  
  int            super;
  int            pbcflag;
  float         *pbcbox;
  
  // if distance == rmsd, array of centers coors
  float        **xcoor;
  float        **ycoor;
  float        **zcoor;
  struct coords
  {
   float       **coor;
  }    *refcoors;
  float       ***refcoor;
  float        **movcoor;
  // if distance == drms, array of centers distmtx
  _drmsdata     *drmsdata;
  // and a distmatrix for the moving frame (not really moving)
  float         *movdist_mtx;
};
// ------------------------------------------------------------------

#ifndef GCLUSTER_INCLUDES
// ------------------------------------------------------------------
int Read_iCluster(char **input, int inp_index, struct inp_Cluster * , char * , Molecule * , CoorSetData *coorsetdata, int nframe );
//----------------------------------------------------------
int Compute_Cluster ( struct inp_Cluster *, struct sopt *, CoorSet *, char *, int   );
// ------------------------------------------------------------------
int Post_Cluster (struct inp_Cluster * );
// ------------------------------------------------------------------
int Read_iCAssign ( char **input, int inp_index, struct inp_CAssign *, char *, Molecule *, CoorSetData *, int nframe );
int Compute_CAssign ( struct inp_CAssign *, struct sopt *, CoorSet *, char *, int  );
// ------------------------------------------------------------------
void ReadCluslist ( struct inp_CAssign * );
// ------------------------------------------------------------------
void FillDistMatrix( struct inp_Cluster *);
// ------------------------------------------------------------------
void FillDistMatrix_threaded( struct inp_Cluster *);
// -------
void * thr_FillRMSDMatrix( void * pFMdata );
// -------
void * thr_FillDRMSMatrix( void * pFMdata );
// ------------------------------------------------------------------
int ClusterQT ( float **, int ,  _cluster * );
// ------------------------------------------------------------------
void ClusterHiero ( float **, int ,  _cluster * );
// ------------------------------------------------------------------
void MatrixCluster ( float **, int , int , _cluster * );
//------------------------------------------------------------------------------
int CheckClusters( float **, _cluster * );
// ------------------------------------------------------------------
int FileClustering( char **ppcInput, int iInputLineNum );
// ------------------------------------------------------------------
#endif

#endif
