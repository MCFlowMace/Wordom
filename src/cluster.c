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
/*! \file cluster.c
 \brief Clustering module
 
 Source code for all clustering procedures
*/

#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "analysis.h"
#include "cluster.h"
#include "gCluster.h"

#define C_LINKING
#include "tsne.h"

//---------------------------------------------------------------------
// Clustering module - work in progress...
// ------------------------------------------------------------------
//void Read_iRmsd_c1 ( FILE *inpfile, struct inp_rms_c1 *inp_rms_c1, FILE *outfile, Molecule *molecule )
int Read_iCluster ( char **input, int inp_index, struct inp_Cluster *inp_cluster, char *printout, Molecule *molecule, CoorSetData *coorsetdata,int nframe )
{
   int         ii, jj, kk;
   char        buffer[1024];
   char        title[64];
   char        distance[64];
   char        method[64];
   int         step=1;          // skip step for frames in trajectory
   int         nsframe;         // frames to be clustered = nframe/step
   float       cutoff;
   int         gotit;
   
   extern short int  no_frame_par;
   no_frame_par = 1;
   
   inp_cluster->method = 0;
   inp_cluster->distance = 0;
   inp_cluster->threshold = 0.0;
   inp_cluster->maxspeed = 0;
   inp_cluster->super = 1;
   inp_cluster->step = 1;
   inp_cluster->himem = 1;
   inp_cluster->threaded = 0;
   inp_cluster->outmatrix = 0;
   inp_cluster->outmatrixonly = 0;
   inp_cluster->inmatrix = 0;
   inp_cluster->markov = 1;
   inp_cluster->cmapd = 0;
   inp_cluster->cmapd_cutoff = 0;
   inp_cluster->nointrasegm = 0;
   inp_cluster->perplexity = 0.0;
   inp_cluster->max_iter = 0;
   inp_cluster->dimension = 0;
   inp_cluster->threshold_bh = -1.0;
   inp_cluster->rand_seed = -1;
   
   
   strcpy( inp_cluster->inmol, molecule->rawmol.name );
   strcpy( inp_cluster->intrj, coorsetdata->datasetname );
   inp_cluster->mol=molecule;
   // scan the input file options
   memset (  title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
//  fgets(buffer, 2048, inpfile);
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_cluster->sele.selestring);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--FIT", 5))
    {
      inp_cluster->fit = 1;
      sscanf(buffer, "--FIT%[^\n]%*c ", inp_cluster->fitsele.selestring);
      fprintf( stderr, "Sorry, --FIT not supported for clustering yet\n");
      exit(0);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CUTOFF", 8))
    {
      sscanf (buffer, "--CUTOFF %f", &cutoff);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
      sscanf( buffer, "--TITLE %s", title);
      gotit = 1;
    } 
    else if ( !strncmp(buffer, "--DISTANCE", 10))
    {
      sscanf( buffer, "--DISTANCE %s", distance);
      if( !strncmp(distance, "rmsd", 4 )) inp_cluster->distance = 1;
      if( !strncmp(distance, "drms", 4 )) inp_cluster->distance = 2;
      if( !strncmp(distance, "tani", 4 )) inp_cluster->distance = 3;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CMAPD", 7))
    {
      inp_cluster->cmapd=1;
      sscanf( buffer, "--CMAPD %f", &(inp_cluster->cmapd_cutoff));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--SEED", 6))
    {
      sscanf( buffer, "--SEED %d", &(inp_cluster->rand_seed));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--PERPLEXITY", 12))
    {
      sscanf( buffer, "--PERPLEXITY %lf", &(inp_cluster->perplexity));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--THETABH", 9))
    {
      sscanf( buffer, "--THETABH %lf", &(inp_cluster->threshold_bh));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--DIM", 5))
    {
      sscanf( buffer, "--DIM %d", &(inp_cluster->dimension));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--IMAX", 6))
    {
      sscanf( buffer, "--IMAX %d", &(inp_cluster->max_iter));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOINTRASEGM", 13))
    {
      inp_cluster->nointrasegm=1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--METHOD", 8))
    {
      sscanf( buffer, "--METHOD %s", method);
      if( !strncmp(method, "hiero" , 5 )) inp_cluster->method = 1;
      if( !strncmp(method,   "qt"  , 2 )) inp_cluster->method = 2;
      if( !strncmp(method, "leader", 6 )) inp_cluster->method = 3;
      if( !strncmp(method, "fastqtlike", 10 )) inp_cluster->method = 4;
      if( !strncmp(method, "tsne", 4 )) inp_cluster->method = 5;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--STEP", 6))
    {
      sscanf( buffer, "--STEP %d", &step);
      sscanf( buffer, "--STEP %d", &inp_cluster->step);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--LOMEM", 7))
    {
      //fprintf( stderr, "Sorry, --LOMEM option is not fully implemented yet\n");
      //exit(0);
      inp_cluster->himem = 0;
      sscanf( buffer, "--LOMEM %s", inp_cluster->tmpfilename);
      inp_cluster->tmpfile = fopen( inp_cluster->tmpfilename, "wb" );
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOSUPER", 9))
    {
      inp_cluster->super=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MAXSPEED", 10))
    {
      inp_cluster->maxspeed=1;
      inp_cluster->markov=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--LFULL", 7))
    {
      inp_cluster->maxspeed=0;
      inp_cluster->markov=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--BACKW", 7) || !strncmp(buffer, "--MARKOV", 8))
    {
      inp_cluster->markov=1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NT", 4))
    {
      #ifndef THREADED
       fprintf( stderr, "this binary not compiled with threading\n");
       exit(0);
      #endif
      inp_cluster->threaded=1;
      sscanf(buffer, "--NT %d", &inp_cluster->nthreads);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--OUTMATRIX", 11))
    {
      if( inp_cluster->inmatrix )
      {
        fprintf( stderr, "use --INMATRIX _or_ --OUTMATRIX please\n");
        exit(0);
      }
      inp_cluster->outmatrix = 1;
      sscanf( buffer, "--OUTMATRIX %s", inp_cluster->matrixfilename);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MATRIXONLY", 12))
    {
      inp_cluster->outmatrixonly = 1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--INMATRIX", 10))
    {
      if( inp_cluster->outmatrix )
      {
        fprintf( stderr, "use --INMATRIX _or_ --OUTMATRIX please\n");
        exit(0);
      }
      inp_cluster->inmatrix = 1;
      sscanf( buffer, "--INMATRIX %s", inp_cluster->matrixfilename);
      gotit = 1;
    }
    if( gotit==0 )
    {
      fprintf( stderr, "Could not understand option: %s\n", buffer);
      exit(5);
    }
    inp_index++;
   }
   
   if( inp_cluster->method == 0 || (inp_cluster->distance == 0 && inp_cluster->inmatrix == 0) )
   {
     fprintf( stderr, "!!!: both method (%d) and distance (%d) required in cluster module\n", inp_cluster->method, inp_cluster->distance );
     exit(12);
   }
   
   if( inp_cluster->outmatrixonly == 1 && inp_cluster->outmatrix == 0 )
   {
     fprintf( stderr, "--MATRIXONLY requires --OUTMATRIX\n");
     exit(12);
   }
   
   if( inp_cluster->maxspeed && inp_cluster->markov )
   {
     fprintf( stderr, "non-markov procedure already implies maxspeed treatment\n");
     exit(12);
   }

   if( inp_cluster->cmapd == 1 && inp_cluster->distance != 2 )
   {
     fprintf( stderr, "CMAPD option is only compatible with drms distance \n");
     exit(12);
   }
   
   if( inp_cluster->nointrasegm == 1 && inp_cluster->distance != 2 )
   {
     fprintf( stderr, "NOINTRASEGM option is only compatible with drms distance \n");
     exit(12);
   }
   
   if( inp_cluster->method == 5)
   {
	   if(inp_cluster->distance != 2 && inp_cluster->distance != 1)
	   {
		   fprintf( stderr, "tsne only works with drms or rmsd\n");
		   exit(0);
	   }
	   
	   if( inp_cluster->max_iter <= 0)
	   {
		   fprintf( stderr, "Please supply a valid number of maximum iterations (use --IMAX)\n");
		   exit(0);
	   }

	   if( inp_cluster->perplexity <= 0)
	   {
		   fprintf( stderr, "Please supply a valid perplexity (use --PERPLEXITY)\n");
		   exit(0);
	   }
	   fprintf( stderr, "%f\n",inp_cluster->threshold_bh);
	   if( inp_cluster->threshold_bh < 0)
	   {
		   fprintf( stderr, "Please supply a valid Barnes-Hut threshold (use --THETABH)\n");
		   
		   exit(0);
	   }
	   
	   if( inp_cluster->dimension <= 0)
	   {
		   fprintf( stderr, "Please supply a valid dimension (use --DIM)\n");
		   exit(0);
	   }
	   
	   
   } else {
	   
		if(inp_cluster->max_iter != 0 || inp_cluster->threshold_bh != -1.0 || inp_cluster->perplexity != 0.0 || inp_cluster->dimension != 0)
		{
			fprintf( stderr, "IMAX, THETABH, DIM and PERPLEXITY options are only for T-SNE clustering!\n");
			exit(12);
		}

		if(cutoff == 0) 
		{
			fprintf(stderr,"Please supply a cutoff (use --CUTOFF)\n");
			exit(12);
		}
	   
   }
   
   if( title[0]=='\0' )
   {
     fprintf( stderr, "Please supply a title (use --TITLE)\n");
     exit(12);
   }
   
   GetSele ( inp_cluster->sele.selestring, &inp_cluster->sele, molecule);
   if( inp_cluster->sele.nselatm < 3 && inp_cluster->distance == 1 )
   {
     fprintf( stderr, "Cannot use rmsd distance with less than 3 atoms selected\n");
     exit(0);
   }
   else if ( inp_cluster->sele.nselatm < 2 )
   {
     fprintf( stderr, "Please provide at least 2 atoms (more if possible) as a selection\n");
     exit(0);
   }

   fprintf(stderr,"# Clustering using method %d , distance %d (cmapd_cutoff %f), cutoff %8.3f on %d selected atoms (%s) nointrasegm %d cmapd %d\n",inp_cluster->method,inp_cluster->distance,inp_cluster->cmapd_cutoff,cutoff,inp_cluster->sele.nselatm, inp_cluster->sele.selestring,inp_cluster->nointrasegm, inp_cluster->cmapd);
   sprintf( inp_cluster->header,"# Clustering using method %d , distance %d, cutoff %8.3f on %d selected atoms (%s)\n",inp_cluster->method,inp_cluster->distance,cutoff,inp_cluster->sele.nselatm, inp_cluster->sele.selestring);
   if( inp_cluster->method == 3 && inp_cluster->markov )
     fprintf(stderr, "# Using backward cluster search (non-markovian data set, default)\n");
   inp_cluster->nato = inp_cluster->sele.nselatm;
   nsframe = (int) nframe/step;
   fprintf( stdout, "DEBUG - numbers %d %d %d\n", nsframe, nframe, step); 
   
   inp_cluster->totframe = nsframe;
   inp_cluster->threshold = cutoff;
   inp_cluster->step = step;
   inp_cluster->nframe = 0;
   
   if( inp_cluster->method == 1 || inp_cluster->method == 2 )
   {
     fprintf( stdout, "DEBUG - alloc for %d frames\n", nsframe); fflush(stdout);
     // there actually is no particular reason why we should have 
     // consequent memory areas to store x-values for different frames
     // so calloc has been modified
     inp_cluster->xcoor = (float **)calloc( nsframe+1, sizeof(float *));
     for ( ii =0; ii<nsframe+1; ii++)
      inp_cluster->xcoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
     
     inp_cluster->ycoor = (float **)calloc( nsframe+1, sizeof(float *));
     for ( ii =0; ii<nsframe+1; ii++)
      inp_cluster->ycoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
     
     inp_cluster->zcoor = (float **)calloc( nsframe+1, sizeof(float *));
     for ( ii =0; ii<nsframe+1; ii++)
      inp_cluster->zcoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
     
     // as before, requiring consecutive memory is probably just a hassle
     // and makes optimizing memory size more difficult
     // so we're back to "normal" (albeit triangular) matrix allocation
     inp_cluster->matrix    = (float **)calloc( nsframe+1, sizeof(float *));
     for( ii=0; ii<nsframe+1; ii++)
       inp_cluster->matrix[ii] = (float *)calloc( ii, sizeof(float));
   }
   
   if( inp_cluster->method == 3 )
   {
    inp_cluster->totframe = nframe;
    //inp_cluster->frameapp = (int *)calloc( inp_cluster->totframe+1 , sizeof( int)); // old, not taking the step into account and therefore taking more memory
    //inp_cluster->frameapp = (int*) malloc((inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int));
    //memset(inp_cluster->frameapp,-1,(inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int)); //Memsetting to -1 for debugging purposes
    inp_cluster->frameapp = (int *) calloc((inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1), sizeof(int)); // +1 is needed because we will use the real number of frame, not the index
    
    inp_cluster->refcoor = calloc( 3, sizeof(float *));
    inp_cluster->refcoor[0] = calloc ( inp_cluster->nato*3, sizeof ( float  ));
    for ( ii =0; ii<3; ii++)
      inp_cluster->refcoor[ii] = inp_cluster->refcoor[0] + inp_cluster->nato*ii;
    /*
    inp_cluster->refcoor    = calloc(   inp_cluster->nato, sizeof(float *));
    inp_cluster->refcoor[0] = calloc( 3*inp_cluster->nato, sizeof(float  ));
    for ( ii =0; ii<inp_cluster->nato; ii++)
     inp_cluster->refcoor[ii] = inp_cluster->refcoor[0] + 3*ii;
    */
    inp_cluster->nclusters = 0;
    inp_cluster->clusterlist = calloc( inp_cluster->totframe+1, sizeof(_clusters *));
    
    if( inp_cluster->markov )
      inp_cluster->cluster_to_check = (int *)calloc( inp_cluster->totframe+1 , sizeof( int));
     
    
    if( inp_cluster->distance == 2 )
    {
      inp_cluster->msize = (int)(inp_cluster->nato*(inp_cluster->nato-1)/2);
      inp_cluster->dist_mtx = calloc ( inp_cluster->msize , sizeof (float));
      if(  inp_cluster->nointrasegm == 1 )
      {
        inp_cluster->nointrasegm_msk=calloc(inp_cluster->msize , sizeof (float));
        kk=0;
        inp_cluster->nointrasegm_corr_fact=0;
        for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
        {
          for ( jj=0; jj<ii; jj++)
          {
            if(!strcmp(inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[ii]-1], inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[jj]-1]))
              inp_cluster->nointrasegm_msk[kk] = 0.0;
            else
            {
              inp_cluster->nointrasegm_msk[kk] = 1.0; 
              inp_cluster->nointrasegm_corr_fact+=1.0;
            } 
            kk++;
          }
        }
         
        inp_cluster->nointrasegm_corr_fact = sqrt( kk/inp_cluster->nointrasegm_corr_fact );
      }
           

      inp_cluster->dist_mtx_size = inp_cluster->msize ;
      inp_cluster->tmpcoor = (float **) walloc( 3, sizeof(float *));
      inp_cluster->tmpcoor[0] = (float *) walloc( 3*inp_cluster->nato, sizeof(float ));
      for( ii=0; ii<3; ii++)
        inp_cluster->tmpcoor[ii] = inp_cluster->tmpcoor[0] + ii*inp_cluster->nato;
      for( ii=0; ii<3*inp_cluster->nato; ii++)
        inp_cluster->tmpcoor[0][ii] = 0.0;
    }

   }
  
  if( inp_cluster->method == 4 )
  {
    // new, optimized, gp method, using leader-like to find neighbours
    // and then apply further methods.
    // linked lists in leader-like, coupled to reallocated p-array
    
    inp_cluster->xcoor = (float **)calloc( nsframe+1, sizeof(float *));
    for ( ii =0; ii<nsframe+1; ii++)
     inp_cluster->xcoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
    
    inp_cluster->ycoor = (float **)calloc( nsframe+1, sizeof(float *));
    for ( ii =0; ii<nsframe+1; ii++)
     inp_cluster->ycoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
    
    inp_cluster->zcoor = (float **)calloc( nsframe+1, sizeof(float *));
    for ( ii =0; ii<nsframe+1; ii++)
     inp_cluster->zcoor[ii] = (float *)calloc( inp_cluster->sele.nselatm, sizeof ( float  ));
    
    inp_cluster->totframe = nframe;    
    //inp_cluster->frameapp = (int *)calloc( inp_cluster->totframe+1 , sizeof( int)); // old, not taking the step into account and therefore taking more memory
    //inp_cluster->frameapp = (int*) malloc((inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int));
    //memset(inp_cluster->frameapp,-1,(inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int)); //Memsetting to -1 for debugging purposes
    inp_cluster->frameapp = (int *) calloc((inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1), sizeof(int)); // +1 is needed because we will use the real number of frame, not the index
    
    inp_cluster->refcoor = calloc( 3, sizeof(float *));
    inp_cluster->refcoor[0] = calloc ( inp_cluster->nato*3, sizeof ( float  ));
    for ( ii =0; ii<3; ii++)
      inp_cluster->refcoor[ii] = inp_cluster->refcoor[0] + inp_cluster->nato*ii;
    
    inp_cluster->nclusters = 0;
    //inp_cluster->clusterlist = calloc( inp_cluster->totframe+1, sizeof(_clusters *));
    inp_cluster->head = NULL;
    inp_cluster->tail = NULL;
    
    if( inp_cluster->markov )
      inp_cluster->cluster_to_check = (int *)calloc( inp_cluster->totframe+1 , sizeof( int));
    
    if( inp_cluster->distance == 1 )
    {
     
    }
    
    if( inp_cluster->distance == 2 )
    {
      fprintf( stderr, "Sorry, fastqtlike only works with RMSD\n" );
      exit(0);
      inp_cluster->msize = (int)(inp_cluster->nato*(inp_cluster->nato-1)/2);
      inp_cluster->dist_mtx = calloc ( inp_cluster->msize , sizeof (float));
      if(  inp_cluster->nointrasegm == 1 )
      {
        // build mask to ignore intra-segment distances
        inp_cluster->nointrasegm_msk=calloc(inp_cluster->msize , sizeof (float));
        kk=0;
        inp_cluster->nointrasegm_corr_fact=0;
        for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
        {
          for ( jj=0; jj<ii; jj++)
          {
            if(!strcmp(inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[ii]-1], inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[jj]-1]))
              inp_cluster->nointrasegm_msk[kk] = 0.0;
            else
            {
              inp_cluster->nointrasegm_msk[kk] = 1.0; 
              inp_cluster->nointrasegm_corr_fact+=1.0;
            } 
            kk++;
          }
        }
        // correction factor for drms obtained with "masked" distances
        inp_cluster->nointrasegm_corr_fact = sqrt( kk/inp_cluster->nointrasegm_corr_fact );
      }
           
      inp_cluster->dist_mtx_size = inp_cluster->msize ;
      inp_cluster->tmpcoor = (float **) walloc( 3, sizeof(float *));
      inp_cluster->tmpcoor[0] = (float *) walloc( 3*inp_cluster->nato, sizeof(float ));
      for( ii=0; ii<3; ii++)
        inp_cluster->tmpcoor[ii] = inp_cluster->tmpcoor[0] + ii*inp_cluster->nato;
      for( ii=0; ii<3*inp_cluster->nato; ii++)
        inp_cluster->tmpcoor[0][ii] = 0.0;
    }
    
    
    // end clustering method 5 (new1) init
  }
  
  if(inp_cluster->method == 5)
  {
	  inp_cluster->totframe = nframe;
    
    /*inp_cluster->refcoor = calloc( 3, sizeof(float *));
    inp_cluster->refcoor[0] = calloc ( inp_cluster->nato*3, sizeof ( float  ));
    for ( ii =0; ii<3; ii++)
      inp_cluster->refcoor[ii] = inp_cluster->refcoor[0] + inp_cluster->nato*ii;*/
     
     
	if(inp_cluster->distance ==1) {
		        
			inp_cluster->tsne_coords = calloc ( (size_t)inp_cluster->nato *3*(size_t)(inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1), sizeof (double));
	         
			if(!inp_cluster->tsne_coords) {
				fprintf(stderr,"Not enough memory!\n");
				exit(-1);
			}
	}
    
    if( inp_cluster->distance == 2 )
    {
      inp_cluster->msize = (int)(inp_cluster->nato*(inp_cluster->nato-1)/2);
      inp_cluster->dist_mtx = calloc ( inp_cluster->msize , sizeof (float));
      
      if(  inp_cluster->nointrasegm == 1 )
      {
        inp_cluster->nointrasegm_msk=calloc(inp_cluster->msize , sizeof (float));
        kk=0;
        inp_cluster->nointrasegm_corr_fact=0;
        for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
        {
          for ( jj=0; jj<ii; jj++)
          {
            if(!strcmp(inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[ii]-1], inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[jj]-1]))
              inp_cluster->nointrasegm_msk[kk] = 0.0;
            else
            {
              inp_cluster->nointrasegm_msk[kk] = 1.0; 
              inp_cluster->nointrasegm_corr_fact+=1.0;
            } 
            kk++;
          }
        }
         
        inp_cluster->nointrasegm_corr_fact = sqrt( kk/inp_cluster->nointrasegm_corr_fact );
      }
      
      inp_cluster->dist_mtx_size = inp_cluster->msize ;
      inp_cluster->tsne_coords = calloc ( (size_t)inp_cluster->msize * (size_t)(inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1), sizeof (double));
      //inp_cluster->tsne_coords = malloc ( (size_t)inp_cluster->msize * (size_t)(inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1)* sizeof (double));
      
      //for(ii=0; ii < inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1; ii++)
		//inp_cluster->tsne_coords[ii]=-1.0;
       
      if(!inp_cluster->tsne_coords) {
		  fprintf(stderr,"Not enough memory!\n");
		  exit(-1);
	  }

      inp_cluster->tmpcoor = (float **) walloc( 3, sizeof(float *));
      inp_cluster->tmpcoor[0] = (float *) walloc( 3*inp_cluster->nato, sizeof(float ));
      for( ii=0; ii<3; ii++)
        inp_cluster->tmpcoor[ii] = inp_cluster->tmpcoor[0] + ii*inp_cluster->nato;
      for( ii=0; ii<3*inp_cluster->nato; ii++)
        inp_cluster->tmpcoor[0][ii] = 0.0;
    }
  }
  
   if( inp_cluster->himem == 0 )
   {
     inp_cluster->tmpfilecoor = InitCoor( inp_cluster->sele.nselatm );
   }
   
   inp_cluster->output = O_File( title, "w");
   if(inp_cluster->method == 3)
   {
     sprintf( printout, " %10s cluster# distance ", title);
     return 30;
   }
   else
     return 0;
}

#ifdef CUDA
int Read_iGcluster ( char **input, int inp_index, struct inp_Cluster *inp_cluster, char *printout, Molecule *molecule, CoorSetData *coorsetdata,int nframe )
{
   int         ii, jj, kk;
   char        buffer[1024];
   char        title[64];
   char        gTitle[68];
   char        distance[64];
   char        method[64];
   int         step=1;          // skip step for frames in trajectory
   int         nsframe;         // frames to be clustered = nframe/step
   float       cutoff;
   int         gotit;

   
   int realGPUs = find_GPUs();
   
   if(realGPUs >0 ) {
		
		fprintf(stderr,"Number of GPUs: %d\n",realGPUs);
	   
	} else {
		fprintf(stderr,"No CUDA GPUs found on the system. Cannot run gCluster!\n");
		exit(-1);
	}

   
   extern short int  no_frame_par;
   no_frame_par = 1;
   
   inp_cluster->method = 0;
   inp_cluster->distance = 0;
   inp_cluster->threshold = 0.0;
   inp_cluster->maxspeed = 1;
   inp_cluster->super = 1;
   inp_cluster->step = 1;
   inp_cluster->himem = 1;
   inp_cluster->threaded = 0;
   inp_cluster->outmatrix = 0;
   inp_cluster->outmatrixonly = 0;
   inp_cluster->inmatrix = 0;
   inp_cluster->markov = 0;
   inp_cluster->cmapd = 0;
   inp_cluster->cmapd_cutoff = 0;
   inp_cluster->nointrasegm = 0;
   inp_cluster->device = 0;
   
   
   strcpy( inp_cluster->inmol, molecule->rawmol.name );
   strcpy( inp_cluster->intrj, coorsetdata->datasetname );
   inp_cluster->mol=molecule;
   // scan the input file options
   memset (  title, '\0', sizeof(title));
   memset (  gTitle, '\0', sizeof(gTitle));
   memset ( buffer, '\0', sizeof(buffer));
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
//  fgets(buffer, 2048, inpfile);
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--SELE", 6))
    {
      sscanf(buffer, "--SELE %[^\n]%*c ", inp_cluster->sele.selestring);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--DEVICE", 6))
    {
      sscanf(buffer, "--DEVICE %u", &(inp_cluster->device));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--FIT", 5))
    {
      inp_cluster->fit = 1;
      sscanf(buffer, "--FIT%[^\n]%*c ", inp_cluster->fitsele.selestring);
      fprintf( stderr, "Sorry, --FIT not supported for clustering yet\n");
      exit(0);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CUTOFF", 8))
    {
      sscanf (buffer, "--CUTOFF %f", &cutoff);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
		sscanf( buffer, "--TITLE %s", title);
		snprintf(gTitle, sizeof gTitle, "%s.out", title);
		gotit = 1;
    } 
    else if ( !strncmp(buffer, "--DISTANCE", 10))
    {
      sscanf( buffer, "--DISTANCE %s", distance);
      if( !strncmp(distance, "rmsd", 4 )) inp_cluster->distance = 1;
      if( !strncmp(distance, "drms", 4 )) inp_cluster->distance = 2;
      if( !strncmp(distance, "tani", 4 )) inp_cluster->distance = 3;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CMAPD", 7))
    {
      inp_cluster->cmapd=1;
      sscanf( buffer, "--CMAPD %f", &(inp_cluster->cmapd_cutoff));
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOINTRASEGM", 13))
    {
      inp_cluster->nointrasegm=1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--METHOD", 8))
    {
      sscanf( buffer, "--METHOD %s", method);
      if( !strncmp(method, "hiero" , 5 )) inp_cluster->method = 1;
      if( !strncmp(method,   "qt"  , 2 )) inp_cluster->method = 2;
      if( !strncmp(method, "leader", 6 )) inp_cluster->method = 3;
      if( !strncmp(method, "fastqtlike", 10 )) inp_cluster->method = 4;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--STEP", 6))
    {
      sscanf( buffer, "--STEP %d", &step);
      sscanf( buffer, "--STEP %d", &inp_cluster->step);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--LOMEM", 7))
    {
      fprintf( stderr, "Sorry, --LOMEM option is not fully implemented yet for GPU clustering\n");
      exit(0);
      //inp_cluster->himem = 0;
      //sscanf( buffer, "--LOMEM %s", inp_cluster->tmpfilename);
      //inp_cluster->tmpfile = fopen( inp_cluster->tmpfilename, "wb" );
      //gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOSUPER", 9))
    {
      inp_cluster->super=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MAXSPEED", 10))
    {
      inp_cluster->maxspeed=1;
      inp_cluster->markov=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--LFULL", 7))
    {
      inp_cluster->maxspeed=0;
      inp_cluster->markov=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--BACKW", 7) || !strncmp(buffer, "--MARKOV", 8))
    {
      inp_cluster->markov=1;
      inp_cluster->maxspeed=0;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NT", 4))
    {
      #ifndef THREADED
       fprintf( stderr, "this binary not compiled with threading\n");
       exit(0);
      #endif
      inp_cluster->threaded=1;
      sscanf(buffer, "--NT %d", &inp_cluster->nthreads);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--OUTMATRIX", 11))
    {
      if( inp_cluster->inmatrix )
      {
        fprintf( stderr, "use --INMATRIX _or_ --OUTMATRIX please\n");
        exit(0);
      }
      inp_cluster->outmatrix = 1;
      sscanf( buffer, "--OUTMATRIX %s", inp_cluster->matrixfilename);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MATRIXONLY", 12))
    {
      inp_cluster->outmatrixonly = 1;
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--INMATRIX", 10))
    {
      if( inp_cluster->outmatrix )
      {
        fprintf( stderr, "use --INMATRIX _or_ --OUTMATRIX please\n");
        exit(0);
      }
      inp_cluster->inmatrix = 1;
      sscanf( buffer, "--INMATRIX %s", inp_cluster->matrixfilename);
      gotit = 1;
    }
    if( gotit==0 )
    {
      fprintf( stderr, "Could not understand option: %s\n", buffer);
      exit(5);
    }
    inp_index++;
   }
   
   if((inp_cluster->device+1)>realGPUs || inp_cluster->device<0) {
	   fprintf(stderr,"Please supply the number of an existing device (--DEVICE)\n");
	   exit(12);
   }
   
   if(cutoff == 0) 
   {
	   fprintf(stderr,"Please supply a cutoff (use --CUTOFF)\n");
	   exit(12);
   }
   
   if(inp_cluster->markov) 
   {
	   fprintf(stderr, "Sorry, non-markov procedure not yet implemented for GPU clustering\n");
	   exit(12);
   }
   
   if( inp_cluster->method == 0 || (inp_cluster->distance == 0 && inp_cluster->inmatrix == 0) )
   {
     fprintf( stderr, "!!!: both method (%d) and distance (%d) required in cluster module\n", inp_cluster->method, inp_cluster->distance );
     exit(12);
   }
   
   if( inp_cluster->outmatrixonly == 1 && inp_cluster->outmatrix == 0 )
   {
     fprintf( stderr, "--MATRIXONLY requires --OUTMATRIX\n");
     exit(12);
   }
   
   if( inp_cluster->maxspeed && inp_cluster->markov )
   {
     fprintf( stderr, "non-markov procedure already implies maxspeed treatment\n");
     exit(12);
   }

   if( inp_cluster->cmapd == 1 && inp_cluster->distance != 2 )
   {
     fprintf( stderr, "CMAPD option is only compatible with drms distance \n");
     exit(12);
   }

   if( inp_cluster->nointrasegm == 1 && inp_cluster->distance != 2 )
   {
     fprintf( stderr, "NOINTRASEGM option is only compatible with drms distance \n");
     exit(12);
   }
   
   if( title[0]=='\0' )
   {
     fprintf( stderr, "Please supply a title (use --TITLE)\n");
     exit(12);
   }
   
   if( inp_cluster->method != 3)
   {
		fprintf(stderr, "Sorry, GPU clustering only supports method leader so far\n");
		exit(12);
   }
   
   if( inp_cluster->distance == 3) 
   {
		fprintf(stderr, "Sorry, GPU clustering does not support distance tani so far\n");
		exit(12);
   }
   
   /*if(inp_cluster->super && inp_cluster->distance ==1) {
	   fprintf(stderr, "Sorry, GPU clustering with rmsd only works without superposition so far\n");
	   exit(12);
   }*/
   
   GetSele ( inp_cluster->sele.selestring, &inp_cluster->sele, molecule);
   if( inp_cluster->sele.nselatm < 3 && inp_cluster->distance == 1 )
   {
     fprintf( stderr, "Cannot use rmsd distance with less than 3 atoms selected\n");
     exit(0);
   }
   else if ( inp_cluster->sele.nselatm < 2 )
   {
     fprintf( stderr, "Please provide at least 2 atoms (more if possible) as a selection\n");
     exit(0);
   }

	   fprintf(stderr,"# GPU Clustering using method %d , distance %d (cmapd_cutoff %f), cutoff %8.3f on %d selected atoms (%s) nointrasegm %d cmapd %d on device %d\n",inp_cluster->method,inp_cluster->distance,inp_cluster->cmapd_cutoff,cutoff,inp_cluster->sele.nselatm, inp_cluster->sele.selestring,inp_cluster->nointrasegm, inp_cluster->cmapd, inp_cluster->device);
	   sprintf( inp_cluster->header,"# Clustering using method %d , distance %d, cutoff %8.3f on %d selected atoms (%s)\n",inp_cluster->method,inp_cluster->distance,cutoff,inp_cluster->sele.nselatm, inp_cluster->sele.selestring);
	   if( inp_cluster->method == 3 && inp_cluster->maxspeed )
	     fprintf(stderr, "# Using maxspeed cluster search (default for GPU clustering)\n");
	   inp_cluster->nato = inp_cluster->sele.nselatm;
	   nsframe = (int) nframe/step;
	   fprintf( stdout, "DEBUG - numbers %d %d %d\n", nsframe, nframe, step); 
	   
	   inp_cluster->totframe = nsframe;
	   inp_cluster->threshold = cutoff;
	   inp_cluster->step = step;
	   inp_cluster->nframe = 0;
   
   
		//always uses method 3, any other chosen method already exited at this point
	    inp_cluster->totframe = nframe;
	    //inp_cluster->frameapp = (int *)calloc( (inp_cluster->totframe/inp_cluster->step+1), sizeof( int)); // +1 is needed because we will use the real number of frame, not the index (maybe)
	    inp_cluster->frameapp = (int*) malloc((inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int));
	    memset(inp_cluster->frameapp,-1,(inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1)+1)*sizeof(int));
	    
	    inp_cluster->nclusters = 0;
	    inp_cluster->clusterlist = calloc( inp_cluster->totframe+1, sizeof(_clusters *));
	    
	    if(inp_cluster->distance ==1) {
			
		    //refcoor[0] are x coords of all atoms of reference frame, refcoor[1] y coords ...
		    /*inp_cluster->refcoor = calloc( 3, sizeof(float *));
		    inp_cluster->refcoor[0] = calloc ( inp_cluster->nato*3, sizeof ( float  ));
		    for ( ii =0; ii<3; ii++)
				inp_cluster->refcoor[ii] = inp_cluster->refcoor[0] + inp_cluster->nato*ii;*/     
		    
			inp_cluster->gclust_coords = calloc ( (size_t)inp_cluster->nato *3*(size_t)(inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1), sizeof (float));
	         
			if(!inp_cluster->gclust_coords) {
				fprintf(stderr,"Not enough memory!\n");
				exit(-1);
			}
	    }    
	    
	    if( inp_cluster->distance == 2 )
	    {
		      inp_cluster->msize = (int)(inp_cluster->nato*(inp_cluster->nato-1)/2);
		      inp_cluster->dist_mtx = calloc ( inp_cluster->msize , sizeof (float));
		      if(  inp_cluster->nointrasegm == 1 )
		      {
		        inp_cluster->nointrasegm_msk=calloc(inp_cluster->msize , sizeof (float));
		        kk=0;
		        inp_cluster->nointrasegm_corr_fact=0;
		        for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
		        {
		          for ( jj=0; jj<ii; jj++)
		          {
		            if(!strcmp(inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[ii]-1], inp_cluster->mol->rawmol.segId[inp_cluster->sele.selatm[jj]-1]))
		              inp_cluster->nointrasegm_msk[kk] = 0.0;
		            else
		            {
		              inp_cluster->nointrasegm_msk[kk] = 1.0; 
		              inp_cluster->nointrasegm_corr_fact+=1.0;
		            } 
		            kk++;
		          }
		        }
		         
		        inp_cluster->nointrasegm_corr_fact = sqrt( kk/inp_cluster->nointrasegm_corr_fact );
		      }
		         
		
		      inp_cluster->dist_mtx_size = inp_cluster->msize ;
		      inp_cluster->gclust_coords = calloc ( (size_t)inp_cluster->msize * (size_t)(inp_cluster->totframe/inp_cluster->step+ (inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1) +1), sizeof (float));
		       
		      //gclust_dmtx can get very big and thus the allocation can fail, if there is not enough memory
		      if(!inp_cluster->gclust_coords) {
				  fprintf(stderr,"Not enough memory!\n");
				  exit(-1);
			  }
			  
		      inp_cluster->tmpcoor = (float **) walloc( 3, sizeof(float *));
		      inp_cluster->tmpcoor[0] = (float *) walloc( 3*inp_cluster->nato, sizeof(float ));
		      
		      for( ii=0; ii<3; ii++)
		        inp_cluster->tmpcoor[ii] = inp_cluster->tmpcoor[0] + ii*inp_cluster->nato;
		      for( ii=0; ii<3*inp_cluster->nato; ii++)
		        inp_cluster->tmpcoor[0][ii] = 0.0;
	    }
	      
	   
	   inp_cluster->output = O_File( title, "w");
	   inp_cluster->gOutput = O_File( gTitle, "w");
	
	   return 30;

}
#endif


// ------------------------------------------------------------------
int ClusterLeaderRmsd( struct inp_Cluster *inp_cluster, CoorSet *trj_crd, int frame, float *finalrmsd )
{
   int          ii, index;
   float        rmsd, minrmsd;
   float      **framecoor;
   int          checked_clusters;
   int          lframe;
   
   lframe = (frame-1)/inp_cluster->step +1;
   //lframe = (frame)/inp_cluster->step; probably wrong
   //lframe = frame;
   if( inp_cluster->nclusters == 0 )
   { // at the beginning, get first frame as first cluster
     inp_cluster->clusterlist[0] = calloc( 1, sizeof(_clusters));
     inp_cluster->clusterlist[0][0].center = frame;
     inp_cluster->clusterlist[0][0].nelements = 1;
     /*
     inp_cluster->clusterlist[0][0].cluscoors    = calloc(   inp_cluster->nato, sizeof(float *));
     inp_cluster->clusterlist[0][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
     for ( ii =0; ii<inp_cluster->nato; ii++)
       inp_cluster->clusterlist[0][0].cluscoors[ii] = inp_cluster->clusterlist[0][0].cluscoors[0] + 3*ii;
     for( ii=0; ii<inp_cluster->nato; ii++ )
     {
       inp_cluster->clusterlist[0][0].cluscoors[ii][0] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[ii][1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[ii][2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
     }
     */
     inp_cluster->clusterlist[0][0].cluscoors    = calloc(   3, sizeof(float *));
     inp_cluster->clusterlist[0][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
     for ( ii =0; ii<3; ii++)
       inp_cluster->clusterlist[0][0].cluscoors[ii] = inp_cluster->clusterlist[0][0].cluscoors[0] + ii*inp_cluster->nato;
     for( ii=0; ii<inp_cluster->nato; ii++ )
     {
       inp_cluster->clusterlist[0][0].cluscoors[0][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[1][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[2][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
     }
     inp_cluster->nclusters = 1;
     inp_cluster->frameapp[lframe] = 0;
     finalrmsd[0] = -1.0;
     return inp_cluster->clusterlist[0][0].center ;
   }
   
   // if we really have to deal with the frame, let's copy the coordinates
   framecoor = inp_cluster->refcoor;
   for( ii=0; ii<inp_cluster->nato; ii++ )
   {
     framecoor[0][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     framecoor[1][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     framecoor[2][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   if( inp_cluster->markov)     // comparing with previous frames' cluster first
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
       inp_cluster->cluster_to_check[ii] = 1;
     
     checked_clusters = 0;
     
     for( ii=frame-1; ii>=0; ii--)
     {
       if( checked_clusters == inp_cluster->nclusters )
         break; // if we have already checked all possible clusters, get out of loop
       index = inp_cluster->frameapp[ii];
       if( inp_cluster->cluster_to_check[index] == 0 )
         continue;      // if this cluster is already checked, continue tp next (previous) frame
       
       rmsd = RmsdCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato , inp_cluster->super );
       //printf( "DEBUG - rmsd: %8.3f\n", rmsd);
       if( rmsd<inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finalrmsd[0] = rmsd;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
       else
       {
         checked_clusters ++;
         inp_cluster->cluster_to_check[index] = 0;
       }
     }
   }
   else if( inp_cluster->maxspeed)      // no markov process, but first close cluster is taken
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       rmsd = RmsdCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato , inp_cluster->super );
       //printf( "DEBUG - rmsd: %8.3f\n", rmsd);
       if( rmsd<inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finalrmsd[0] = rmsd;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
     }
   }
   else         // default is to compare with all existing clusters 
   {            // (order is then irrelevant) and assign to closest
     minrmsd = inp_cluster->threshold;
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       rmsd = RmsdCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato , inp_cluster->super );
       //printf( "DEBUG - rmsd: %8.3f\n", rmsd);
       if( rmsd<minrmsd )
       {
         minrmsd = rmsd;
         inp_cluster->frameapp[lframe] = index;
       }
     }
     if( minrmsd < inp_cluster->threshold )
     {
       inp_cluster->clusterlist[inp_cluster->frameapp[lframe]]->nelements += 1;
       finalrmsd[0] = minrmsd;
       return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
     }
   }
   
   // if no cluster is close, create a new one
   index = inp_cluster->nclusters;
   inp_cluster->clusterlist[index] = calloc( 1, sizeof(_clusters));
   inp_cluster->clusterlist[index][0].center = frame;
   inp_cluster->clusterlist[index][0].nelements = 1;
   /*
   inp_cluster->clusterlist[index][0].cluscoors    = calloc(   inp_cluster->nato, sizeof(float *));
   inp_cluster->clusterlist[index][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
   for ( ii =0; ii<inp_cluster->nato; ii++)
     inp_cluster->clusterlist[index][0].cluscoors[ii] = inp_cluster->clusterlist[index][0].cluscoors[0] + 3*ii;
   for( ii=0; ii<inp_cluster->nato; ii++ )
   {
     inp_cluster->clusterlist[index][0].cluscoors[ii][0] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[ii][1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[ii][2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   */
   inp_cluster->clusterlist[index][0].cluscoors    = calloc(   3, sizeof(float *));
   inp_cluster->clusterlist[index][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
   for ( ii =0; ii<3; ii++)
     inp_cluster->clusterlist[index][0].cluscoors[ii] = inp_cluster->clusterlist[index][0].cluscoors[0] + inp_cluster->nato*ii;
   for( ii=0; ii<inp_cluster->nato; ii++ )
   {
     inp_cluster->clusterlist[index][0].cluscoors[0][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[1][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[2][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   for( ii=0; ii<3*inp_cluster->nato; ii++ )
     inp_cluster->clusterlist[index][0].cluscoors[0][ii] = framecoor[0][ii];
   inp_cluster->nclusters += 1;
   inp_cluster->frameapp[lframe] = index;
   
   return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
}

// ------------------------------------------------------------------
//this method only stores the coordinates of the current frame in an array the clustering takes place later
void saveCoords( struct inp_Cluster *inp_cluster, CoorSet *trj_crd, int frame )
{
   int          ii, lframe;

   lframe = (frame-1)/inp_cluster->step+1;
   
   // trj_crd contains coords of current frame
   
   if( inp_cluster->method==5) {
	   
	   for( ii=0; ii<inp_cluster->nato; ii++ )
	   {
	     inp_cluster->tsne_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
	     inp_cluster->tsne_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3+1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
	     inp_cluster->tsne_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3+2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
	   }
	} else {
		
	   for( ii=0; ii<inp_cluster->nato; ii++ )
	   {
	     inp_cluster->gclust_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
	     inp_cluster->gclust_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3+1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
	     inp_cluster->gclust_coords[(size_t)lframe*(size_t)inp_cluster->nato*3+(size_t)ii*3+2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
	   }
   }
   
}


// ------------------------------------------------------------------
int ClusterLeaderDrms( struct inp_Cluster *inp_cluster, CoorSet *trj_crd, int frame, float *finaldrms )
{
   int          ii, index;
   float        drms, mindrms;
   int          checked_clusters;
   int          lframe;
   
   //lframe = frame/inp_cluster->step; probably wrong 
   lframe = (frame-1)/inp_cluster->step+1;
   
   

/* GetCoor( trj_crd, inp_cluster->tmpcoor[0], inp_cluster->tmpcoor[1], inp_cluster->tmpcoor[2], inp_cluster->sele );*/
   for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
   {
     inp_cluster->tmpcoor[0][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->tmpcoor[1][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->tmpcoor[2][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   DistMtxCalc( inp_cluster->pbcflag , inp_cluster->pbcbox, inp_cluster->nato , inp_cluster->tmpcoor, inp_cluster->dist_mtx );
   if( inp_cluster->nointrasegm == 1)
     for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
       inp_cluster->dist_mtx[ii] *= inp_cluster->nointrasegm_msk[ii];
   
   if( inp_cluster->cmapd == 1)         //Compute the smoothed contact map instead of the distance matrix 
     for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
       inp_cluster->dist_mtx[ii] = fret( inp_cluster->dist_mtx[ii] , inp_cluster->cmapd_cutoff );
   
   if( inp_cluster->nclusters == 0 )
   { 
     // at the beginning, get first frame as first cluster
     inp_cluster->clusterlist[0] = calloc( 1, sizeof(_clusters));
     inp_cluster->clusterlist[0][0].center = frame;
     inp_cluster->clusterlist[0][0].nelements = 1;
     inp_cluster->clusterlist[0][0].distance = calloc ( inp_cluster->msize , sizeof (float));
     inp_cluster->dist_mtx_size = inp_cluster->msize ;
     for( ii=0; ii<inp_cluster->dist_mtx_size; ii++ )
       inp_cluster->clusterlist[0][0].distance[ii] = inp_cluster->dist_mtx[ii];

     inp_cluster->nclusters = 1;
     inp_cluster->frameapp[lframe] = 0;
     finaldrms[0] = -1.0;
     
     return inp_cluster->clusterlist[0][0].center ;
   }
   
   if( inp_cluster->markov)     // comparing with previous frames' cluster first
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
       inp_cluster->cluster_to_check[ii] = 1;
     
     checked_clusters = 0;
     
     for( ii=frame-1; ii>=0; ii--)
     {
       if( checked_clusters == inp_cluster->nclusters )
         break; // if we have already checked all possible clusters, get out of loop
       index = inp_cluster->frameapp[ii];
       if( inp_cluster->cluster_to_check[index] == 0 )
         continue;      // if this cluster is already checked, continue tp next (previous) frame
       
       Drms( inp_cluster->dist_mtx_size, inp_cluster->clusterlist[index][0].distance, inp_cluster->dist_mtx, &drms );
       if( inp_cluster->nointrasegm == 1 ) //Renormalize the distance properly
         drms *= inp_cluster->nointrasegm_corr_fact;

       if( drms<inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finaldrms[0] = drms;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
       else
       {
         checked_clusters ++;
         inp_cluster->cluster_to_check[index] = 0;
       }
     }
   }
   else if( inp_cluster->maxspeed)      // no markov process, but first close cluster is taken
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       Drms( inp_cluster->dist_mtx_size, inp_cluster->clusterlist[ii][0].distance, inp_cluster->dist_mtx, &drms );
 
       if( inp_cluster->nointrasegm == 1 ) //Renormalize the distance properly
         drms *= inp_cluster->nointrasegm_corr_fact;

       if( drms<inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finaldrms[0] = drms;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
     }
   }
   else         // default is to compare with all existing clusters 
   {            // (order is then irrelevant) and assign to closest
     mindrms = inp_cluster->threshold;
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       Drms( inp_cluster->dist_mtx_size, inp_cluster->clusterlist[ii][0].distance, inp_cluster->dist_mtx, &drms );

       if( inp_cluster->nointrasegm == 1 ) //Renormalize the distance properly
         drms *= inp_cluster->nointrasegm_corr_fact;

       if( drms<mindrms )
       {
         mindrms = drms;
         inp_cluster->frameapp[lframe] = index;
       }
     }
     if( mindrms < inp_cluster->threshold )
     {
       inp_cluster->clusterlist[inp_cluster->frameapp[lframe]]->nelements += 1;
       finaldrms[0] = mindrms;
       return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
     }
   } 
   
   // if no cluster is close, create a new one
   index = inp_cluster->nclusters;
   inp_cluster->nclusters += 1;
   inp_cluster->frameapp[lframe] = index;
   inp_cluster->clusterlist[index] = calloc( 1, sizeof(_clusters));
   inp_cluster->clusterlist[index][0].center = frame;
   inp_cluster->clusterlist[index][0].nelements = 1;
   inp_cluster->clusterlist[index][0].distance = calloc ( inp_cluster->msize , sizeof (float));
   inp_cluster->dist_mtx_size = inp_cluster->msize ;
   for( ii=0; ii<inp_cluster->dist_mtx_size; ii++ )
     inp_cluster->clusterlist[index][0].distance[ii] = inp_cluster->dist_mtx[ii];
     
   
   return inp_cluster->clusterlist[index][0].center ;

}

//for the GPU clustering this method has been changed to reading only the distance matrix of a frame and storing it in an array
//the real clustering takes place later in the Post_Gcluster part
void saveDistMtx( struct inp_Cluster *inp_cluster, CoorSet *trj_crd, int frame )
{
   int          ii,lframe;
   
   lframe = (frame-1)/inp_cluster->step+1;

   for ( ii=0; ii<inp_cluster->sele.nselatm; ii++)
   {
     inp_cluster->tmpcoor[0][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->tmpcoor[1][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->tmpcoor[2][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   
   DistMtxCalc( inp_cluster->pbcflag , inp_cluster->pbcbox, inp_cluster->nato , inp_cluster->tmpcoor, inp_cluster->dist_mtx );
   if( inp_cluster->nointrasegm == 1)
     for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
       inp_cluster->dist_mtx[ii] *= inp_cluster->nointrasegm_msk[ii];
   
   //Compute the smoothed contact map instead of the distance matrix 
   if( inp_cluster->cmapd == 1)
     for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
       inp_cluster->dist_mtx[ii] = fret( inp_cluster->dist_mtx[ii] , inp_cluster->cmapd_cutoff );
       
	if( inp_cluster->method==5) {
	
		for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
			inp_cluster->tsne_coords[(size_t)lframe*(size_t)inp_cluster->dist_mtx_size+(size_t)ii]=inp_cluster->dist_mtx[ii];
	} else {
		
		for ( ii=0; ii<inp_cluster->dist_mtx_size; ii++)
			inp_cluster->gclust_coords[(size_t)lframe*(size_t)inp_cluster->dist_mtx_size+(size_t)ii]=inp_cluster->dist_mtx[ii];
	}

}

// ------------------------------------------------------------------
float TaniCalc ( float **cluscoors , float **framecoor , int nato )
{
  int     ii;
  float   sum1=0, sum2=0, sum3=0, sim=0;
  
  for( ii=0; ii<nato; ii++)
  {
    sum1 += cluscoors[ii][0]*framecoor[ii][0];
    sum2 += cluscoors[ii][0]*cluscoors[ii][0];
    sum3 += framecoor[ii][0]*framecoor[ii][0];
  }
  
  sim = sum1/(sum2+sum3-sum1);
  
  return sim;
}
// ------------------------------------------------------------------
int ClusterLeaderTanimoto( struct inp_Cluster *inp_cluster, CoorSet *trj_crd, int frame, float *finaltani )
{
   int          ii, index;
   float        tani, maxtani;
   float      **framecoor;
   int          checked_clusters;
   int          lframe;
   
//   lframe = (frame-1)/inp_cluster->step;
   lframe = (frame)/inp_cluster->step;
   //lframe = frame;
   if( inp_cluster->nclusters == 0 )
   { // at the beginning, get first frame as first cluster
     inp_cluster->clusterlist[0] = calloc( 1, sizeof(_clusters));
     inp_cluster->clusterlist[0][0].center = frame;
     inp_cluster->clusterlist[0][0].nelements = 1;
     inp_cluster->clusterlist[0][0].cluscoors    = calloc(   inp_cluster->nato, sizeof(float *));
     inp_cluster->clusterlist[0][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
     for ( ii =0; ii<inp_cluster->nato; ii++)
       inp_cluster->clusterlist[0][0].cluscoors[ii] = inp_cluster->clusterlist[0][0].cluscoors[0] + 3*ii;
     for( ii=0; ii<inp_cluster->nato; ii++ )
     {
       inp_cluster->clusterlist[0][0].cluscoors[ii][0] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[ii][1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
       inp_cluster->clusterlist[0][0].cluscoors[ii][2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
     }
     inp_cluster->nclusters = 1;
     inp_cluster->frameapp[lframe] = 0;
     finaltani[0] = -1.0;
     return inp_cluster->clusterlist[0][0].center ;
   }
   
   // if we really have to deal with the frame, let's copy the coordinates
   framecoor = inp_cluster->refcoor;
   for( ii=0; ii<inp_cluster->nato; ii++ )
   {
     framecoor[ii][0] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     framecoor[ii][1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     framecoor[ii][2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   if( inp_cluster->markov)     // comparing with previous frames' cluster first
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
       inp_cluster->cluster_to_check[ii] = 1;
     
     checked_clusters = 0;
     
     for( ii=frame-1; ii>=0; ii--)
     {
       if( checked_clusters == inp_cluster->nclusters )
         break; // if we have already checked all possible clusters, get out of loop
       index = inp_cluster->frameapp[ii];
       if( inp_cluster->cluster_to_check[index] == 0 )
         continue;      // if this cluster is already checked, continue to next (ie previous) frame
       
       tani = TaniCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato );
       if( tani>inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finaltani[0] = tani;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
       else
       {
         checked_clusters ++;
         inp_cluster->cluster_to_check[index] = 0;
       }
     }
   }
   else if( inp_cluster->maxspeed)      // no markov process, but first close cluster is taken
   {
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       tani = TaniCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato );
       if( tani>inp_cluster->threshold )
       {
         inp_cluster->frameapp[lframe] = index;
         inp_cluster->clusterlist[index]->nelements += 1;
         finaltani[0] = tani;
         return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
       }
     }
   }
   else         // default is to compare with all existing clusters 
   {            // (order is then irrelevant) and assign to closest
     maxtani = inp_cluster->threshold;
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       index = ii;
       tani = TaniCalc ( inp_cluster->clusterlist[index]->cluscoors , framecoor , inp_cluster->nato );
       if( tani>maxtani )
       {
         maxtani = tani;
         inp_cluster->frameapp[lframe] = index;
       }
     }
     if( maxtani > inp_cluster->threshold )
     {
       inp_cluster->clusterlist[inp_cluster->frameapp[lframe]]->nelements += 1;
       finaltani[0] = maxtani;
       return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
     }
   } 
   
   // if no cluster is close, create a new one
   index = inp_cluster->nclusters;
   inp_cluster->clusterlist[index] = calloc( 1, sizeof(_clusters));
   inp_cluster->clusterlist[index][0].center = frame;
   inp_cluster->clusterlist[index][0].nelements = 1;
   inp_cluster->clusterlist[index][0].cluscoors    = calloc(   inp_cluster->nato, sizeof(float *));
   inp_cluster->clusterlist[index][0].cluscoors[0] = calloc( 3*inp_cluster->nato, sizeof ( float  ));
   for ( ii =0; ii<inp_cluster->nato; ii++)
     inp_cluster->clusterlist[index][0].cluscoors[ii] = inp_cluster->clusterlist[index][0].cluscoors[0] + 3*ii;
   for( ii=0; ii<inp_cluster->nato; ii++ )
   {
     inp_cluster->clusterlist[index][0].cluscoors[ii][0] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[ii][1] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
     inp_cluster->clusterlist[index][0].cluscoors[ii][2] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
   }
   for( ii=0; ii<3*inp_cluster->nato; ii++ )
     inp_cluster->clusterlist[index][0].cluscoors[0][ii] = framecoor[0][ii];
   inp_cluster->nclusters += 1;
   inp_cluster->frameapp[lframe] = index;
   
   return inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center;
}
// ------------------------------------------------------------------
int Compute_Cluster ( struct inp_Cluster *inp_cluster, struct sopt *OPT, CoorSet *trj_crd, char *outprint, int frame )
{
   int     ii;
   int     lframe;
   float   distance;
   int     clusterapp;
   
   if( fmod( frame-1, inp_cluster->step) )
   { // if you are to skip some frames (--STEP option), do so
     if ( inp_cluster->method == 1 || inp_cluster->method == 2 )
       return 0;
     else
       return 21;
   }
   if ( inp_cluster->method == 1 || inp_cluster->method == 2 )
   { 
     lframe = (frame-1)/inp_cluster->step; 
     //fprintf( stdout, "DEBUG - lframe: %d\n", lframe); fflush(stdout);
     // choosen method requires all-vs-all frame comparison: store coordinates for post-processing
     if( inp_cluster->himem )
     {
       for( ii=0; ii < inp_cluster->sele.nselatm; ii++ )
       {
         inp_cluster->xcoor[lframe][ii] = trj_crd->xcoor[inp_cluster->sele.selatm[ii]-1];
         inp_cluster->ycoor[lframe][ii] = trj_crd->ycoor[inp_cluster->sele.selatm[ii]-1];
         inp_cluster->zcoor[lframe][ii] = trj_crd->zcoor[inp_cluster->sele.selatm[ii]-1];
       }
     }
     else
     {
       GetSeleCoor( trj_crd, inp_cluster->tmpfilecoor, &inp_cluster->sele);
       fwrite( inp_cluster->tmpfilecoor->cords, sizeof(float), 3*inp_cluster->sele.nselatm, inp_cluster->tmpfile );
     }
     
     inp_cluster->nframe ++;
     return 0;
   }
 
   if ( inp_cluster->method == 3 )
   {
     lframe = (frame-1)/inp_cluster->step+1; 
     //lframe = frame/inp_cluster->step; probably wrong 
     if( inp_cluster->distance == 1 )
     {
       clusterapp = -1;
       distance   = -1;
       clusterapp = ClusterLeaderRmsd( inp_cluster, trj_crd, frame, &distance);
       sprintf( outprint, " %10d %8d %8.3f ", clusterapp, inp_cluster->frameapp[lframe]+1, distance);
       return 25;
     }
     if( inp_cluster->distance == 2 )
     {
       clusterapp = -1;
       distance   = -1;
       clusterapp = ClusterLeaderDrms( inp_cluster, trj_crd, frame, &distance);
       sprintf( outprint, " %10d %8d %8.3f ", clusterapp, inp_cluster->frameapp[lframe]+1, distance);
       return 25;
     }
     if( inp_cluster->distance == 3 )
     {
       clusterapp = -1;
       distance   = -1;
       clusterapp = ClusterLeaderTanimoto( inp_cluster, trj_crd, frame, &distance);
       sprintf( outprint, " %10d %8d %8.3f ", clusterapp, inp_cluster->frameapp[lframe]+1, distance);
       return 25;
     }
   }
   
	if ( inp_cluster->method == 5 )
	{
		
		if( inp_cluster->distance == 1 )
		{
			//only saves coordinates of current frame for later processing
			saveCoords( inp_cluster, trj_crd, frame);
			return 25;
	    }
	    if( inp_cluster->distance == 2 )
	    {
			//only saves distance matrix of current frame for later processing
			saveDistMtx( inp_cluster, trj_crd, frame);
			return 25;
	    }
	    
	    if( inp_cluster->distance == 3 )
	    {
			//not supported
	    }
	}
   
   return 0;
}

#ifdef CUDA
int Compute_Gcluster ( struct inp_Cluster *inp_cluster, struct sopt *OPT, CoorSet *trj_crd, char *outprint, int frame )
{
   int     ii, lframe;
   
   if( fmod( frame-1, inp_cluster->step) )
   { // if you are to skip some frames (--STEP option), do so
     if ( inp_cluster->method == 1 || inp_cluster->method == 2 )
       return 0;
     else
       return 21;
   }
   
   if ( inp_cluster->method == 1 || inp_cluster->method == 2 )
   { 
     //Not supported
   }
 
   if ( inp_cluster->method == 3 )
   {
     lframe = (frame-1)/inp_cluster->step+1;
     if( inp_cluster->distance == 1 )
     {
		//only saves coordinates of current frame in inp_cluster->gclust_coords, the clustering takes place in Post_Gcluster
       saveCoords( inp_cluster, trj_crd, frame);
       return 25;
     }
     if( inp_cluster->distance == 2 )
     {
		 //only saves distance matrix of current frame in inp_cluster->gclust_coords, the clustering takes place in Post_Gcluster
       saveDistMtx( inp_cluster, trj_crd, frame);
       return 25;
     }
     
     if( inp_cluster->distance == 3 )
     {
		//Not supported
     }
   }
   
   return 0;
}
#endif

// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ------------------------------------------------------------------
int llist_cmp(const void *AA, const void *BB)
{
  _clusters **aa = (_clusters **)AA;
  _clusters **bb = (_clusters **)BB;

  return (int)(bb[0]->nelements - aa[0]->nelements);
}
// ------------------------------------------------------------------
int MkDistList( struct inp_Cluster *inp_cluster )
{
  int          ii, jj, kk;
  int          nframe=0, nato=0;
  float      **refcoor=NULL, **movcoor=NULL;
  float        rmsd;
  int          thisconfindex ;
  _clusters    *thiscluster;
  _clusters    *thatcluster;
  struct ll_conf_neig *local_llconfneig = NULL;
  struct ll_clus_neig *ll_neig_cursor = NULL;
  
  int           tcounter;
  int          *n_conf_neig;
  
  confs        *cursor_on_confs;
  struct ll_conf_neig *cursor_conf_neig = NULL;

  nframe = inp_cluster->nframe;
  nato = inp_cluster->sele.nselatm;
  
  refcoor = calloc( 3, sizeof(float *));
  refcoor[0] = calloc ( nato*3, sizeof ( float  ));
  for ( ii =0; ii<3; ii++)
    refcoor[ii] = refcoor[0] + nato*ii;
  
  movcoor = (float **)calloc( 3, sizeof(float *));
  movcoor[0] = (float *)calloc ( nato*3, sizeof ( float  ));
  for ( ii =0; ii<3; ii++)
    movcoor[ii] = movcoor[0] + nato*ii;
  
  n_conf_neig = (int *)calloc( nframe, sizeof( int));
  
  inp_cluster->llconfneig = (struct ll_conf_neig **) calloc( nframe, sizeof( struct ll_conf_neig *));
  for( ii=0; ii<nframe; ii++ )
    inp_cluster->llconfneig[ii] = NULL;
  inp_cluster->conf_nneig = (int *) calloc( nframe, sizeof( int ));
  
  //ccenterindex = (int *)calloc(inp_cluster->nclusters, sizeof( int ));
  
  // build list of elements for each cluster
  for( ii=0; ii<inp_cluster->nclusters; ii++ )
  {
    inp_cluster->clusters_array[ii]->clus_confs = (int *)calloc( inp_cluster->clusters_array[ii]->nelements, sizeof(int));
    cursor_on_confs = inp_cluster->clusters_array[ii]->ll_confs;
    for( jj=0; jj<inp_cluster->clusters_array[ii]->nelements; jj++ )
    {
      if( cursor_on_confs == NULL )
      {
        fprintf( stderr, "Error while traveling on confs linked list for cluster #%d\n", ii );
        exit(0);
      }
      inp_cluster->clusters_array[ii]->clus_confs[jj] = cursor_on_confs->conf_index;
      cursor_on_confs = cursor_on_confs->next;
    }
  }
  
  for( ii=0; ii<nframe; ii++ )
  {
    for( kk=0; kk<nato; kk++ )
    {
      refcoor[0][kk] = inp_cluster->xcoor[ii][kk];
      refcoor[1][kk] = inp_cluster->ycoor[ii][kk];
      refcoor[2][kk] = inp_cluster->zcoor[ii][kk];
    }
    thiscluster = inp_cluster->clusters_array[inp_cluster->frameapp[ii]];
    // compare with frames from its native cluster
    for( jj=0; jj<thiscluster->nelements; jj++ )
    {
      thisconfindex = thiscluster->clus_confs[jj]; //-1;
      if( ii == thisconfindex )
        rmsd = 0.0;
      else
      {
        for( kk=0; kk<nato; kk++ )
        {
          movcoor[0][kk] = inp_cluster->xcoor[thisconfindex][kk];
          movcoor[1][kk] = inp_cluster->ycoor[thisconfindex][kk];
          movcoor[2][kk] = inp_cluster->zcoor[thisconfindex][kk];
        }
        rmsd = RmsdCalc( refcoor, movcoor, nato, inp_cluster->super );
      }
      
      if( rmsd<inp_cluster->threshold )
      {
        n_conf_neig[ii] += 1;
        
        if( inp_cluster->llconfneig[ii] == NULL )
        {
          inp_cluster->llconfneig[ii] = (struct ll_conf_neig *) calloc( 1, sizeof( struct ll_conf_neig ));
          local_llconfneig = inp_cluster->llconfneig[ii];
        }
        else
        {
          local_llconfneig->next = (struct ll_conf_neig *) calloc( 1, sizeof( struct ll_conf_neig ));
          local_llconfneig = local_llconfneig->next;
        }
        local_llconfneig->conf_number = thisconfindex;
        local_llconfneig->distance = rmsd;
        local_llconfneig->next = NULL;
      }
    }
    // then the nearby clusters
    ll_neig_cursor = thiscluster->llneig;
    if( ll_neig_cursor == NULL )
      continue;
    tcounter = 0;
    if( ll_neig_cursor->cl_number == -1 ) // weird check for weird cases ...
      if( ll_neig_cursor->next != NULL )
        ll_neig_cursor = ll_neig_cursor->next;
        
    while( ll_neig_cursor != NULL )
    {
      thatcluster = inp_cluster->clusters_array[ll_neig_cursor->cl_number];
      for( jj=0; jj<thatcluster->nelements; jj++ )
      {
        thisconfindex = thatcluster->clus_confs[jj];//-1;
        for( kk=0; kk<nato; kk++ )
        {
          movcoor[0][kk] = inp_cluster->xcoor[thisconfindex][kk];
          movcoor[1][kk] = inp_cluster->ycoor[thisconfindex][kk];
          movcoor[2][kk] = inp_cluster->zcoor[thisconfindex][kk];
        }
        rmsd = RmsdCalc( refcoor, movcoor, nato, inp_cluster->super );
        
        if( rmsd<inp_cluster->threshold )
        {
          n_conf_neig[ii] += 1;
          
          local_llconfneig->next = (struct ll_conf_neig *) calloc( 1, sizeof( struct ll_conf_neig ));
          local_llconfneig = local_llconfneig->next;
          local_llconfneig->conf_number = thisconfindex;
          local_llconfneig->distance = rmsd;
          local_llconfneig->next = NULL;
        }
      }
      tcounter ++;
      ll_neig_cursor = ll_neig_cursor->next;
    }
  }
  
  for( ii=0; ii<inp_cluster->nframe; ii++ )
  {
    inp_cluster->conf_nneig[ii] = n_conf_neig[ii];
    // also check neighbors' lists
    cursor_conf_neig = inp_cluster->llconfneig[ii];
    jj = 0;
    while( cursor_conf_neig != NULL )
    {
      //fprintf( stdout, "DEBUG neigh #%d: %d\n", jj+1, cursor_conf_neig->conf_number );
      cursor_conf_neig = cursor_conf_neig->next;
      jj++;
    }
  }
  
  return 0;
}
// ------------------------------------------------------------------
int DistListCluster( struct inp_Cluster *inp_cluster )
{
  int       ii;
  int       nn, clusstarter;
  int      *unclustered;
  int       max_neighbors;
  int      *local_conf_nneig;
  
  struct ll_conf_neig  *cursor_conf_neig = NULL;
  _cluster             *cluster1;
  
  
  nn = inp_cluster->nframe;
  cluster1 = (_cluster *) calloc( 1, sizeof(_cluster ));
  cluster1->ncluster = 0;
  cluster1->cluspop = (int *)calloc( nn, sizeof(int));
  for(ii=0;ii<nn;ii++)
   cluster1->cluspop[ii] = 0;
  cluster1->clusstart = (int *)calloc( nn, sizeof(int));
  cluster1->clusstrs = (int **)calloc( nn, sizeof(int *));
  
  local_conf_nneig = (int *)calloc( nn, sizeof(int));
  for( ii=0; ii<nn; ii++ )
    local_conf_nneig[ii] = inp_cluster->conf_nneig[ii];
    
  unclustered = (int *)calloc( nn, sizeof(int));
  for( ii=0; ii<nn; ii++)
   unclustered[ii] = 1;
  
  // loop until it finds any worthy (pop>1) cluster
  while ( 1 )
  {
    //for( ii=0; ii<inp_cluster->nframe; ii++ )
      //fprintf( stdout, "DEBUG frame %d has %d neighbors\n", ii, local_conf_nneig[ii]); fflush(stdout);
    // find biggest neighbors cluster
    max_neighbors = 0;
    for( ii=0; ii<inp_cluster->nframe; ii++ )
    {
      if( local_conf_nneig[ii] > max_neighbors )
      {
        max_neighbors = local_conf_nneig[ii];
        clusstarter = ii;
      }
    }
    // keep going only if biggest clusters' not a single structure
    if( max_neighbors < 2 )
      break;
    // pass all infos about this cluster to the cluster1 structure
    cluster1->ncluster ++;
    cluster1->cluspop[cluster1->ncluster] = max_neighbors;
    cluster1->clusstart[cluster1->ncluster] = clusstarter;
    cluster1->clusstrs[cluster1->ncluster] = calloc(max_neighbors, sizeof(int));
    
    ii = 0;
    cursor_conf_neig = inp_cluster->llconfneig[clusstarter];
    while( cursor_conf_neig != NULL )
    {
      if( unclustered[cursor_conf_neig->conf_number] )
      {
        cluster1->clusstrs[cluster1->ncluster][ii] = cursor_conf_neig->conf_number;
        unclustered[cursor_conf_neig->conf_number] = 0;
        ii++;
      }
      cursor_conf_neig = cursor_conf_neig->next;
    }
    
    unclustered[clusstarter] = 0;
    local_conf_nneig[clusstarter] = 0;
    
    // update neighbors counts ignoring already-clustered confs
    // first, set neig_count=0 for frames included in cluster
    for( ii=0; ii<max_neighbors; ii++ )
      local_conf_nneig[cluster1->clusstrs[cluster1->ncluster][ii]] = 0;
    // then, re-count neigs for each frame
    for( ii=0; ii<inp_cluster->nframe; ii++ )
    {
      if( local_conf_nneig[ii] == 0 )
        continue;
      cursor_conf_neig = inp_cluster->llconfneig[ii];
      local_conf_nneig[ii] = 0;
      while( cursor_conf_neig != NULL )
      {
        if( unclustered[cursor_conf_neig->conf_number] )
          local_conf_nneig[ii] ++;
        cursor_conf_neig = cursor_conf_neig->next;
      }
    }
  }
  
  fprintf( stdout, "DEBUG Listing clusters ...\n");
  cluster1->cluspop[0] = 0;
  for( ii=0; ii<inp_cluster->nframe; ii++ )
    if( unclustered[ii] == 1 )
      cluster1->cluspop[0] ++;
  cluster1->clusstrs[0] = calloc( cluster1->cluspop[0], sizeof(int));
  cluster1->clusstrs[0][0] = -1;
  for( ii=0; ii<cluster1->ncluster+1 ; ii++ )
  {
    fprintf( stdout, "DEBUG Cluster # %4d ; structures: %4d ; origin: %4d\n", ii, cluster1->cluspop[ii], cluster1->clusstart[ii] );
  }
  
  return 0;
}


// ------------------------------------------------------------------
int Post_Cluster ( struct inp_Cluster *inp_cluster )
{
	
   int      ii, jj;
   int      index;
   int      nframe;
   //clock_t  start, end;
   FILE    *matrixfile;
   int      matrix_nframe[1];
   
   
   
	//for(ii = 0; ii <= (inp_cluster->totframe)/inp_cluster->step +1;ii++)
	  //fprintf(stdout,"DEBUG frame %d -> cluster nr %d -> center %d -> elements %d\n",ii,inp_cluster->frameapp[ii],inp_cluster->frameapp[ii] >= 0 ? inp_cluster->clusterlist[inp_cluster->frameapp[ii]][0].center : -1,inp_cluster->frameapp[ii] >= 0 ? inp_cluster->clusterlist[inp_cluster->frameapp[ii]][0].nelements : -1);
   
   //start = clock();
   if(inp_cluster->himem == 0 )
   {
     fclose(inp_cluster->tmpfile);
     inp_cluster->tmpfile = fopen( inp_cluster->tmpfilename, "rb");
   }
     
   nframe = inp_cluster->nframe;
   
   if( inp_cluster->method == 1 || inp_cluster->method == 2 )
   {
     
     if( inp_cluster->inmatrix )
     {
       fprintf( stdout,"Reading distance matrix from %s\n", inp_cluster->matrixfilename);
       matrixfile = O_File( inp_cluster->matrixfilename, "r");
       fread( matrix_nframe, sizeof(int), 1, matrixfile );
       if( matrix_nframe[0] != nframe )
       {
         fprintf( stderr, "Warning: provided nframe (%d) differs from matrix nframe (%d)\n", nframe, matrix_nframe[0]);
         fprintf( stderr, "Correcting ...\n");
         for( ii=0; ii<nframe; ii++)
           free(inp_cluster->matrix[ii]);
         free(inp_cluster->matrix);
         nframe = inp_cluster->nframe = matrix_nframe[0];
         inp_cluster->matrix = calloc( inp_cluster->nframe, sizeof(float *));
         for( ii=0; ii<inp_cluster->nframe; ii++)
           inp_cluster->matrix[ii] = calloc( ii, sizeof(float));
       }
       for( ii=0; ii<nframe; ii++)
         fread( inp_cluster->matrix[ii], sizeof(float), ii, matrixfile );
       fprintf( stdout,"Done reading %s\n", inp_cluster->matrixfilename); fflush(stdout);
       fclose( matrixfile );
     }
     else
     {
       if( inp_cluster->threaded == 0 )
         FillDistMatrix( inp_cluster );
       #ifdef THREADED
       else
         FillDistMatrix_threaded( inp_cluster );
       #endif
       
       for( ii=0; ii<nframe; ii++)
         free(inp_cluster->xcoor[ii]);
       free(inp_cluster->xcoor);
       for( ii=0; ii<nframe; ii++)
         free(inp_cluster->ycoor[ii]);
       free(inp_cluster->ycoor);
       for( ii=0; ii<nframe; ii++)
         free(inp_cluster->zcoor[ii]);
       free(inp_cluster->zcoor);
     }
     
     if( inp_cluster->outmatrix )
     {
       fprintf( stdout,"Writing distance matrix to %s\n", inp_cluster->matrixfilename);
       matrixfile = O_File( inp_cluster->matrixfilename, "w");
       matrix_nframe[0] = nframe;
       fwrite( matrix_nframe, sizeof(int), 1, matrixfile );
       for( ii=0; ii<nframe; ii++)
         fwrite( inp_cluster->matrix[ii], sizeof(float), ii, matrixfile);
       fclose( matrixfile );
       fprintf( stdout,"Done writing %s\n", inp_cluster->matrixfilename); fflush(stdout);
       if( inp_cluster->outmatrixonly )
       {
         fprintf( stdout, "Only distance matrix required: task completed\n");
         return 0;
       }
     }
     
     inp_cluster->cluster1.threshold = inp_cluster->threshold;
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     //start = end;
     MatrixCluster( inp_cluster->matrix, inp_cluster->nframe, inp_cluster->method, &inp_cluster->cluster1 );
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     //start = end;
     
     // bugged ? TODO
     CheckClusters( inp_cluster->matrix, &inp_cluster->cluster1 );
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
     //start = end;
     
     if( inp_cluster->inmatrix == 0 )
     {
       //fprintf( inp_cluster->output, "# Clustering run on %s and %s (%d frames)\n", inp_cluster->inmol, inp_cluster->intrj, inp_cluster->nframe );
       fprintf( inp_cluster->output, "%s", inp_cluster->header);
     }
     else
     {
       //fprintf( inp_cluster->output, "# Clustering run on %s (%d frames)\n", inp_cluster->matrixfilename, inp_cluster->nframe );
       fprintf( inp_cluster->output, "%s", inp_cluster->header);
     }
     fprintf( inp_cluster->output, "# Nclusters (method %d): %d (cluster #0 => unclustered conformations)\n", inp_cluster->method, inp_cluster->cluster1.ncluster);
     
     for( ii=0; ii<inp_cluster->cluster1.ncluster; ii++)
     {
//     fprintf(inp_cluster->output, "Cluster# %4d ; structures: %d ; center: %5d\n", ii, inp_cluster->cluster1.cluspop[ii], inp_cluster->cluster1.cluscenter[ii]+1);
       fprintf(inp_cluster->output, "Cluster# %4d ; structures: %d ; center: %5d\n", ii, inp_cluster->cluster1.cluspop[ii], (inp_cluster->cluster1.cluscenter[ii]*inp_cluster->step)+1);
       
       for( jj=0; jj<inp_cluster->cluster1.cluspop[ii]; jj++)
        fprintf(inp_cluster->output, "%d\n", inp_cluster->cluster1.clusstrs[ii][jj]*inp_cluster->step +1);
     }
     
     //end = clock();
     //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
   }
   else if (inp_cluster->method == 3)
   {
     //fprintf( inp_cluster->output, "# Clustering run on %s and %s (%d frames)\n", inp_cluster->inmol, inp_cluster->intrj, inp_cluster->nframe );
     fprintf( inp_cluster->output, "%s", inp_cluster->header);
     inp_cluster->cluster = (_scluster *)calloc(inp_cluster->nclusters+1, sizeof(_scluster));
     //ccenterindex = (int *)calloc(inp_cluster->nclusters, sizeof( int ));
	 int frames = inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1);
   
     for( ii=0; ii<inp_cluster->nclusters; ii++ )
     {
       inp_cluster->cluster[ii].cluspop = inp_cluster->clusterlist[ii][0].nelements;
       inp_cluster->cluster[ii].clusstart = inp_cluster->cluster[ii].cluscenter = inp_cluster->clusterlist[ii][0].center;
       inp_cluster->cluster[ii].clusstrs = calloc( inp_cluster->cluster[ii].cluspop+1, sizeof(int));
       index = 0;  		
       for( jj=0; jj<frames; jj++ )
       {
         if( inp_cluster->frameapp[jj+1] == ii )
         {
           //inp_cluster->cluster[ii].clusstrs[index] = (jj + 1) * inp_cluster->step; //seems to be wrong
           inp_cluster->cluster[ii].clusstrs[index] = jj * inp_cluster->step + 1;
           index++;
         }
       }
     }

     fprintf(inp_cluster->output, "# Nclusters (method %d): %d \n", inp_cluster->method, inp_cluster->nclusters);
     fprintf(inp_cluster->output, "Cluster#    0 ; structures: 0 ; header: leader-like has no isolated structures\n" );
     
     for( ii=0; ii<inp_cluster->nclusters; ii++)
     {
       fprintf(inp_cluster->output, "Cluster# %4d ; structures: %d ; header: %5d\n", ii+1, inp_cluster->cluster[ii].cluspop, inp_cluster->cluster[ii].cluscenter );
       for( jj=0; jj<inp_cluster->cluster[ii].cluspop; jj++)
         fprintf(inp_cluster->output, "%5d\n", inp_cluster->cluster[ii].clusstrs[jj]);
     }
   }
   else if (inp_cluster->method == 4)
   {
     // DEBUG: calling new distance optimization function
     fprintf( stderr, "DEBUG About to start mkdistlist (%d)\n", inp_cluster->nframe ); fflush(stderr);
     MkDistList( inp_cluster );
     fprintf( stderr, "DEBUG About to start distlistcluster (%d)\n", inp_cluster->nframe ); fflush(stderr);
     DistListCluster( inp_cluster );
   }
   else if (inp_cluster->method == 5)
   {
     
     
    int frames = inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1); //number of datapoints
	//int msize = inp_cluster->msize; //original dimension
	//int no_dims = inp_cluster->dimension; //output dimension
	//int max_iter = inp_cluster->max_iter; 
	//double perplexity = inp_cluster->perplexity;
	//double theta = inp_cluster->threshold_bh;
	     
	//fprintf(stderr,"Reached post cluster\n");
	// Now fire up the SNE implementation
	
	//fprintf(stderr, "postcluster_start\n");
	double* Y = (double*) malloc(frames * inp_cluster->dimension * sizeof(double));
    if(Y == NULL) { fprintf(stderr, "Memory allocation failed!\n"); exit(1); }
    //fprintf(stderr, "postcluster_start2\n");
    
    if(inp_cluster->distance ==1) {
		
		if(inp_cluster->super) {
			setup_tsne(inp_cluster->tsne_coords+3*inp_cluster->nato, frames, 3*inp_cluster->nato, Y, inp_cluster->dimension, inp_cluster->perplexity, inp_cluster->threshold_bh, inp_cluster->rand_seed, inp_cluster->max_iter, true, true);
		} else {
			//fprintf(stderr, "Nosuper!\n");
			setup_tsne(inp_cluster->tsne_coords+3*inp_cluster->nato, frames, 3*inp_cluster->nato, Y, inp_cluster->dimension, inp_cluster->perplexity, inp_cluster->threshold_bh, inp_cluster->rand_seed, inp_cluster->max_iter, false, true);
		}
    } else {
		
		setup_tsne(inp_cluster->tsne_coords+inp_cluster->msize, frames, inp_cluster->msize, Y, inp_cluster->dimension, inp_cluster->perplexity, inp_cluster->threshold_bh, inp_cluster->rand_seed, inp_cluster->max_iter, false, false);
    }
    
    //output
    fprintf(inp_cluster->output, " Result T-Sne, #frames: %d, dimension: %d, perplexity: %f, theta: %f, iterations: %d\n",frames, inp_cluster->dimension, inp_cluster->perplexity, inp_cluster->threshold_bh, inp_cluster->max_iter);
    fprintf(inp_cluster->output, " %6s %12s\n","frame","coords");
    
    for (ii=0; ii<frames; ii++) {
		fprintf(inp_cluster->output, " %10d ",ii);
		for (jj=0; jj<inp_cluster->dimension; jj++) {
			fprintf(inp_cluster->output, "%15.5f ",Y[jj+inp_cluster->dimension*ii]);
		}
		fprintf(inp_cluster->output, "\n");
	}

    // Clean up the memory
	free(Y); Y = NULL;
    free(inp_cluster->tsne_coords); inp_cluster->tsne_coords = NULL;
    
   }

   return 0;
}

//----------------------
#ifdef CUDA
int Post_Gcluster(struct inp_Cluster *inp_cluster,FILE * oA_f)
{
		float *distance;
		int frames = inp_cluster->totframe/inp_cluster->step+(inp_cluster->totframe%inp_cluster->step == 0 ? 0 : 1);
		distance = (float*)calloc((frames + 1),sizeof(float));
		
		//DEBUG fprintf(oA_f,"Results gcluster:\n");
		//DEBUG fprintf( oA_f, " %6s %11s %5s %6s\n","frame","center","cluster#","distance");
	    
	    //calculate the distances and fill the frameapp array using CUDA

	    if(inp_cluster->distance == 1) gClusterRmsd(inp_cluster,distance);
		if(inp_cluster->distance == 2) gClusterDrms(inp_cluster,distance);
		    
		int ii,jj,index; 	
		
		//for(ii = 0; ii <= frames;ii++)
		//	fprintf(stdout,"DEBUG frame %d  -> center %d -> distance %f\n",ii,inp_cluster->frameapp[ii],distance[ii]);
	
		fprintf(stderr,"GPU calculation finished. Postprocessing ..\n");
	
		//generating the actual clusters from the frameapp array
		for(ii = 1; ii <= frames; ii++) {
			
			//after gpu calculation frameapp contains not the cluster numbers but the clustercenters
			//so we get clustercenter from frameapp and correct it according to the stepsize
			int clustercenter = (inp_cluster->frameapp[ii] - 1)*inp_cluster->step +1;
			
			//we check the value of frameapp of all frames, everytime we find a clustercenter that matches the current frame we found a new cluster
			if(inp_cluster->frameapp[ii] == ii) {
				
				index = inp_cluster->nclusters;
				inp_cluster->nclusters += 1;
				inp_cluster->clusterlist[index] = (_clusters*) calloc( 1, sizeof(_clusters));
				inp_cluster->clusterlist[index][0].center = clustercenter;
				//inp_cluster->clusterlist[index][0].distance = (float*) calloc (inp_cluster->msize , sizeof (float));
				inp_cluster->clusterlist[index][0].nelements = 1;
				inp_cluster->frameapp[ii] = index; //change frameapp from clustercenter to the number of the newly created cluster
									
			} else {		
				//if ii is no new cluster check the existing clusters and add it to the one that matches its frameapp value
				for(jj = 0; jj<inp_cluster->nclusters;jj++) {		
					if(inp_cluster->clusterlist[jj][0].center == clustercenter) {
						inp_cluster->clusterlist[jj][0].nelements += 1;
						inp_cluster->frameapp[ii] = jj;
					}
				}
			}
		}	
			
		//generating the output file
		fprintf(stderr,"Generating outputfiles ..\n");
		fprintf(inp_cluster->gOutput,"Results gcluster - numbers %u %u %u\n",inp_cluster->totframe/inp_cluster->step, inp_cluster->totframe, inp_cluster->step);
		fprintf(inp_cluster->gOutput, " %6s %11s %5s %6s\n","frame","center","cluster#","distance");
		for(ii = 0; ii < inp_cluster->totframe; ii++) {
			fprintf(inp_cluster->gOutput,"%7d ", ii+1);
			int lframe = ii/inp_cluster->step + 1;
			if(ii%inp_cluster->step == 0)
				fprintf(inp_cluster->gOutput, " %10d %8d %8.3f ", inp_cluster->clusterlist[inp_cluster->frameapp[lframe]][0].center, inp_cluster->frameapp[lframe]+1, distance[lframe]);
			fprintf(inp_cluster->gOutput,"\n");
			fflush(inp_cluster->gOutput);
		}
		
		fclose(inp_cluster->gOutput);
		free(distance);
		distance = NULL;
		Post_Cluster(inp_cluster);
		return 0;	
}
#endif

//---------------------------------------------------------------------
int FileClustering( char **ppcInput, int iInputLineNum )
{
   int                ii;
   char               buffer[1024];
   char               title[64];
   char               method[64];
   struct inp_Cluster inp_cluster;
   short              gotit;
   short              gottitle;
   
//inp_cluster = calloc( 1, sizeof(struct inp_Cluster));
   inp_cluster.method = 0;
   inp_cluster.himem = 1;
   inp_cluster.threaded = 0;
   inp_cluster.inmatrix = 0;
   inp_cluster.nframe = 0;
   inp_cluster.totframe = 0;
   inp_cluster.step = 1;
   gottitle = 0;
   
       
   // scan the input file options
   memset (  title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
//  fgets(buffer, 2048, inpfile);
    sprintf( buffer, "%s", ppcInput[iInputLineNum]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
     sscanf( buffer, "--TITLE %s", title);
     gottitle = 1;
     gotit = 1;
    } 
    else if ( !strncmp(buffer, "--CUTOFF", 8))
    {
     sscanf (buffer, "--CUTOFF %f", &inp_cluster.threshold);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--NFRAME", 8))
    {
     sscanf (buffer, "--NFRAME %d", &inp_cluster.nframe);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--METHOD", 8))
    {
     sscanf( buffer, "--METHOD %s", method);
     if( !strncmp(method, "hiero" , 5 )) inp_cluster.method = 1;
     else if( !strncmp(method,   "qt"  , 2 )) inp_cluster.method = 2;
     else
     {
       fprintf(stderr, "Unrecognized method label (leader not allowed for file clustering)\n");
       exit(0);
     }
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--INMATRIX", 10))
    {
      inp_cluster.inmatrix = 1;
      sscanf( buffer, "--INMATRIX %s", inp_cluster.matrixfilename);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--NT", 4))
    {
      #ifndef THREADED
       fprintf( stderr, "this binary not compiled with threading\n");
       exit(0);
      #endif
      inp_cluster.threaded=1;
      sscanf(buffer, "--NT %d", &inp_cluster.nthreads);
      gotit = 1;
    }
    if( gotit==0 )
    {
     fprintf( stderr, "Could not understand option: %s\n", buffer);
     exit(5);
    }
    iInputLineNum++;
   }
   
   if( !inp_cluster.inmatrix )
   {
     fprintf( stderr, "You need to specify a distance matrix\n(output from previous clustering run with --OUTMATRIX option)\n");
     exit(0);
   }
   if( inp_cluster.method == 1 )
   {
    fprintf( stderr, "Please run hiero clustering through -iA path with --INMATRIX option\n");
    fprintf( stderr, "complete with mol and trj files (may give -end 1 to speed trj reading)\n");
    fprintf( stderr, "ongoing work on this bug\n");
    //exit(12);
   }
   
   if( inp_cluster.nframe == 0 )
   {
     fprintf( stderr, "You need to specify a frame number (and be sure it is right for the given distance matrix\n");
     exit(0);
   }
   if( gottitle == 0 )
   {
     fprintf( stderr, "You didn't provide a title: a generic \"temp_cluster.out\" will be used\n");
     sprintf( title, "temp_cluster.out");
   }
   
// fprintf( stdout, "DEBUG nframe:%d; threshold: %8.3f\n", inp_cluster->nframe, inp_cluster->threshold);fflush(stdout);
   inp_cluster.matrix = calloc( inp_cluster.nframe, sizeof(float *));
   for( ii=0; ii<inp_cluster.nframe; ii++)
     inp_cluster.matrix[ii] = calloc( ii, sizeof(float));
   
   inp_cluster.totframe = inp_cluster.nframe;
   inp_cluster.output = O_File( title, "w");
   
   Post_Cluster( &inp_cluster );
   
//free(inp_cluster);
   return 0;
}
//---------------------------------------------------------------------
// Clustering module - work in progress...
// ------------------------------------------------------------------

int Read_iCAssign ( char **input, int inp_index, struct inp_CAssign *inp_cassign, char *printout, Molecule *molecule, CoorSetData *coorsetdata, int nframe )
{
   int          ii, jj;
   char         buffer[1024];
   char         title[64];
   char         distance[64];
   int          step=1;          // skip step for frames in trajectory
   int          nsframe;         // frames to be clustered = nframe/step
   float        cutoff;
   int          gotit;
   float      **coor;
   CoorSet      crd;
   
   extern short int  no_frame_par;
   no_frame_par = 1;
   
   sprintf(inp_cassign->sele.selestring, "/*/*/*");
   // scan the input file options
   
   memset (  title, '\0', sizeof(title));
   memset ( buffer, '\0', sizeof(buffer));
   
   inp_cassign->distance = 0;
   inp_cassign->cutoff = 0;
   inp_cassign->super = 1;
   CalloCoor( &crd, molecule->nato );
   
   while( strncmp (buffer, "END", 3))
   {
    gotit = 0;
//  fgets(buffer, 2048, inpfile);
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--SELE", 6))
    {
     sscanf(buffer, "--SELE %[^\n]%*c ", inp_cassign->sele.selestring);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--TITLE", 7))
    {
     sscanf( buffer, "--TITLE %s", title);
     gotit = 1;
    } 
    else if ( !strncmp(buffer, "--FILE", 6))
    {
     sscanf( buffer, "--FILE %s", inp_cassign->cluslistname);
     gotit = 1;
    } 
    else if ( !strncmp(buffer, "--CUTOFF", 8))
    {
     sscanf (buffer, "--CUTOFF %f", &cutoff);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--DISTANCE", 10))
    {
     sscanf( buffer, "--DISTANCE %s", distance);
     if( !strncmp(distance, "rmsd", 4 )) inp_cassign->distance = 1;
     if( !strncmp(distance, "drms", 4 )) inp_cassign->distance = 2;
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--STEP", 6))
    {
     sscanf( buffer, "--STEP %d", &step);
     gotit = 1;
    }
    else if ( !strncmp(buffer, "--NOSUPER", 9))
    {
     inp_cassign->super=0;
     gotit = 1;
    }
    if( gotit==0 )
    {
     fprintf( stderr, "Could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_index++;
   }
   
   if( inp_cassign->distance == 0 )
   {
    fprintf( stderr, "distance is required in cluster assign module\n");
    exit(12);
   }
   if( cutoff == 0 )
   {
    fprintf( stderr, "cutoff is required in cluster assign module\n");
    exit(12);
   }
   
   GetSele ( inp_cassign->sele.selestring, &inp_cassign->sele, molecule);

   inp_cassign->nato = inp_cassign->sele.nselatm;
   nsframe = nframe/step;
   
   inp_cassign->totframe = nsframe;
   inp_cassign->step = step;
   inp_cassign->nframe = 0;
   inp_cassign->cutoff = cutoff;
   
   inp_cassign->output = O_File( title, "w");
   
   ReadCluslist( inp_cassign );
   
   inp_cassign->movcoor = wrd_fmatAlloc( 3, inp_cassign->sele.nselatm );
   
   if( inp_cassign->distance == 1 ) // store coordinates of centers
   {
//   inp_cassign->refcoors = calloc( inp_cassign->nclusters+1, sizeof( struct coords));
     // let's allocate for more than the existing clusters, so new ones can be defined
     inp_cassign->refcoors = calloc( inp_cassign->totframe, sizeof( struct coords));
     for( ii=1; ii<inp_cassign->nclusters+1; ii++ )
     {
       inp_cassign->refcoors[ii].coor = wrd_fmatAlloc( 3, inp_cassign->sele.nselatm );
     }
     
     
     for( ii=1; ii<inp_cassign->nclusters; ii++ )
     {
       GetThisCoorSet( coorsetdata, inp_cassign->clustercenters[ii], &crd);
       for( jj=0; jj<inp_cassign->sele.nselatm; jj++ )
       {
         inp_cassign->refcoors[ii].coor[0][jj] = crd.xcoor[inp_cassign->sele.selatm[jj]-1];
         inp_cassign->refcoors[ii].coor[1][jj] = crd.ycoor[inp_cassign->sele.selatm[jj]-1];
         inp_cassign->refcoors[ii].coor[2][jj] = crd.zcoor[inp_cassign->sele.selatm[jj]-1];
       }
     }
     
   }
   
   if( inp_cassign->distance == 2 ) // store distmtx of centers
   {
//  inp_cassign->drmsdata = malloc( (inp_cassign->nclusters+1)*sizeof(struct drmsdata) );
    // agina, le's allocate space for all possible clusters
    inp_cassign->drmsdata = malloc( (inp_cassign->totframe)*sizeof(struct drmsdata) );
    for( ii=0; ii<(inp_cassign->nclusters+1); ii++ )
    {
     inp_cassign->drmsdata[ii].msize = (int)(inp_cassign->nato*(inp_cassign->nato-1)/2);
     inp_cassign->drmsdata[ii].dist_mtx = calloc ( inp_cassign->drmsdata[ii].msize , sizeof (float));
     inp_cassign->drmsdata[ii].dist_mtx_size = inp_cassign->drmsdata[ii].msize ;
    }
    
    coor = calloc( inp_cassign->sele.nselatm, sizeof(float *));
    coor[0] = calloc ( inp_cassign->sele.nselatm*3, sizeof ( float  ));
    for ( ii =0; ii<inp_cassign->sele.nselatm; ii++)
     coor[ii] = coor[0] + 3*ii;
    
    inp_cassign->movdist_mtx = calloc ( inp_cassign->drmsdata[0].msize , sizeof (float));
    
    for( ii=1; ii<inp_cassign->nclusters; ii++ )
    {
     GetThisCoorSet( coorsetdata, inp_cassign->clustercenters[ii], &crd);
     for( jj=0; jj<inp_cassign->sele.nselatm; jj++ )
     {
      coor[jj][0] = crd.xcoor[inp_cassign->sele.selatm[jj]-1];
      coor[jj][1] = crd.ycoor[inp_cassign->sele.selatm[jj]-1];
      coor[jj][2] = crd.zcoor[inp_cassign->sele.selatm[jj]-1];
     }
     DistMtxCalc ( inp_cassign->pbcflag, inp_cassign->pbcbox, inp_cassign->nato , coor , inp_cassign->drmsdata[ii].dist_mtx ) ;
     printf("cluster #%d (%d) - %8.3f%8.3f%8.3f\n", ii, inp_cassign->clustercenters[ii], inp_cassign->drmsdata[ii].dist_mtx[0], inp_cassign->drmsdata[ii].dist_mtx[1], inp_cassign->drmsdata[ii].dist_mtx[2]);
    }
   }

   sprintf( printout, " %13s ", title);
   
   if( coorsetdata->islistoftrajs )
   {
     CloseTrj( &coorsetdata->traj );
     OpenTrj( coorsetdata->trj_list.trjname[0], &coorsetdata->traj, "r" );
     ReadTrjHeader ( &coorsetdata->traj );
   }
   
   return 18;
}
// ------------------------------------------------------------------
void ReadCluslist ( struct inp_CAssign *inp_cassign )
{
   int          ii, jj;
   FILE        *inpfile;
   char         buffer[128];
   char         nullbuffer[1024];
      
   inpfile = O_File( inp_cassign->cluslistname, "r");
// fscanf(inpfile, "# Nclusters (method %*d): %d %*[^\n}%*c", &inp_cassign->nclusters);
   // grezzly remove first two lines (comments)
   if(fgets( nullbuffer, 1024, inpfile )==NULL){
      fprintf(stderr, "Error! Premature end of the stream reached!\n");
   }
   //fgets( nullbuffer, 1024, inpfile );
   fscanf(inpfile, "%*23c %d %*[^\n]%*c", &inp_cassign->nclusters);
// inp_cassign->nclusters--;    // written ncluster include #0 (isolated frames)
// printf("Number of clusters found: %d\n", inp_cassign->nclusters);

   inp_cassign->clustercenters = malloc( (inp_cassign->nclusters+1)*sizeof( int) );
   
   ii=0;
   while( ii<inp_cassign->nclusters )
   {
    if(fgets(buffer, 128, inpfile)==NULL){
      fprintf(stderr, "Error! Premature end of the stream reached!\n");
    }
    if( !strncmp(buffer, "Cluster#", 8))
    {
     sscanf( buffer, "Cluster# %d ; structures: %*d ; center: %d", &jj, &inp_cassign->clustercenters[ii]);
     
     if( ii != jj )
     {
      fprintf( stderr, "Cluster number unclear for cluster# %d!\n", jj );
      exit(92);
     }
     ii++;
    }
   }
   
   fclose(inpfile);
   
   return;
}
// ------------------------------------------------------------------
int Compute_CAssign ( struct inp_CAssign *inp_cassign, struct sopt *OPT, CoorSet *dcd_crd, char *outprint, int frame )
{
   int          ii, jj;
   int          closestcluster;
   float        nato, min, rmsd, drms;
   
   if( fmod( frame, inp_cassign->step) )
    return 0;
   
   //lframe = (frame-1)/inp_cassign->step;
     
   for( ii=0; ii< inp_cassign->sele.nselatm; ii++ )
   {
    inp_cassign->movcoor[0][ii] = dcd_crd->xcoor[inp_cassign->sele.selatm[ii]-1];
    inp_cassign->movcoor[1][ii] = dcd_crd->ycoor[inp_cassign->sele.selatm[ii]-1];
    inp_cassign->movcoor[2][ii] = dcd_crd->zcoor[inp_cassign->sele.selatm[ii]-1];
   }
   
   if( inp_cassign->distance == 1 )
   {
    nato = inp_cassign->sele.nselatm;
    min=inp_cassign->cutoff;
    closestcluster = 0;
    for( ii=1; ii<inp_cassign->nclusters; ii++ )
    {
     rmsd = RmsdCalc( inp_cassign->refcoors[ii].coor, inp_cassign->movcoor, nato, inp_cassign->super );
     if( rmsd < min ) // difference using <= is insignificant
     {
      min = rmsd;
      closestcluster = ii;
     }
    }
    // if no suitable cluster is found, create a new one
    if( closestcluster == 0 )
    {
      min = 0.0;
      inp_cassign->refcoors[inp_cassign->nclusters].coor = wrd_fmatAlloc( 3, inp_cassign->sele.nselatm );
      for( jj=0; jj<inp_cassign->sele.nselatm*3; jj++ )
        inp_cassign->refcoors[inp_cassign->nclusters].coor[0][jj] = inp_cassign->movcoor[0][jj];
      closestcluster = inp_cassign->nclusters;
      inp_cassign->nclusters++;
    }
    
    sprintf( outprint, " %8d (%8.3f) ", closestcluster, min ); 
   }

   if( inp_cassign->distance == 2 )
   {
    DistMtxCalc ( inp_cassign->pbcflag, inp_cassign->pbcbox, inp_cassign->nato , inp_cassign->movcoor , inp_cassign->movdist_mtx ) ;
    min=inp_cassign->cutoff;
    closestcluster = 0;
//  printf("Frame #%d\n", frame);
    for( ii=1; ii<inp_cassign->nclusters+1; ii++ )
    {
     Drms ( inp_cassign->drmsdata[ii].msize , inp_cassign->drmsdata[ii].dist_mtx, inp_cassign->movdist_mtx, &drms ) ;
//   printf("cluster #%d - %8.3f\n", ii, drms);
     if( drms < min ) // difference using <= is insignificant
     {
      min = drms;
      closestcluster = ii;
     }
    } 
    // if no suitable cluster is found, create a new one
    if( closestcluster == 0 )
    {
      inp_cassign->nclusters++;
      min = 0.0;
      closestcluster = inp_cassign->nclusters;
      inp_cassign->drmsdata[inp_cassign->nclusters].msize = (int)(inp_cassign->nato*(inp_cassign->nato-1)/2);
      inp_cassign->drmsdata[inp_cassign->nclusters].dist_mtx = calloc ( inp_cassign->drmsdata[inp_cassign->nclusters].msize , sizeof (float));
      inp_cassign->drmsdata[inp_cassign->nclusters].dist_mtx_size = inp_cassign->drmsdata[inp_cassign->nclusters].msize ;
      DistMtxCalc ( inp_cassign->pbcflag, inp_cassign->pbcbox, inp_cassign->nato , inp_cassign->movcoor , inp_cassign->drmsdata[inp_cassign->nclusters].dist_mtx ) ;
    }
    sprintf( outprint, " %8d (%5.3f) ", closestcluster, min ); 
   }
   
   return 18;
}
// ------------------------------------------------------------------
//void FillDistMatrix_new(float **matrix, int nato, int nframe, float ***coors)
void FillDistMatrix( struct inp_Cluster *inp_cluster)
{
   
   int          ii, jj, kk;
   int          nframe=0, nato=0;
   float      **refcoor=NULL, **movcoor=NULL;
   
   nframe = inp_cluster->nframe;
   nato = inp_cluster->sele.nselatm;
   
   refcoor = calloc( 3, sizeof(float *));
   refcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<3; ii++)
     refcoor[ii] = refcoor[0] + nato*ii;
   
   movcoor = calloc( 3, sizeof(float *));
   movcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<3; ii++)
     movcoor[ii] = movcoor[0] + nato*ii;
   
   switch ( inp_cluster->distance )
   {
     case 1 :
       // NOTE: during matrix fillup RmsdCalc will alter frame position if super is selected.
       // as superimpose is run every time the rmsd is computed, this should not matter

       for( ii=0; ii<nframe; ii++ )
       {  
         for( kk=0; kk<nato; kk++ )
         {
           refcoor[0][kk] = inp_cluster->xcoor[ii][kk];
           refcoor[1][kk] = inp_cluster->ycoor[ii][kk];
           refcoor[2][kk] = inp_cluster->zcoor[ii][kk];
         }
         for( jj=0; jj<ii; jj++ )
         {
           for( kk=0; kk<nato; kk++ )
           {
             movcoor[0][kk] = inp_cluster->xcoor[jj][kk];
             movcoor[1][kk] = inp_cluster->ycoor[jj][kk];
             movcoor[2][kk] = inp_cluster->zcoor[jj][kk];
           }
           inp_cluster->matrix[ii][jj] = RmsdCalc( refcoor, movcoor, nato, inp_cluster->super );
         }
       }
       break;
     case 2 :
       for( ii=0; ii<nframe; ii++ )
       {
         for( kk=0; kk<nato; kk++ )
         {
           refcoor[0][kk] = inp_cluster->xcoor[ii][kk];
           refcoor[1][kk] = inp_cluster->ycoor[ii][kk];
           refcoor[2][kk] = inp_cluster->zcoor[ii][kk];
         }
//       Calc_dist_mtx ( OPT, n , inp_c->selatm , dcd_crd , dist_mtx ) ;
         for( jj=0; jj<ii; jj++ )
         {
           for( kk=0; kk<nato; kk++ )
           {
             movcoor[0][kk] = inp_cluster->xcoor[jj][kk];
             movcoor[1][kk] = inp_cluster->ycoor[jj][kk];
             movcoor[2][kk] = inp_cluster->zcoor[jj][kk];
           }
           inp_cluster->matrix[ii][jj] = DrmsCalc( refcoor, movcoor, nato, inp_cluster->pbcflag, inp_cluster->pbcbox );
         }
       }
       break;
     default :
       fprintf( stderr, "Unclear distance selection for clustering");
       exit(0);
       break;
   }
   
   free(refcoor[0]);
   free(refcoor);
   free(movcoor[0]);
   free(movcoor);
   
   return;
}
// ------------------------------------------------------------------
//void FillDistMatrix( struct inp_Cluster *inp_cluster)
/*void FillDistMatrix_old( struct inp_Cluster *inp_cluster)
{
   
   int          ii, jj, kk;
   int          nframe, nato;
   float      **refcoor, **movcoor;
   int          nfloats;
   off_t        offset;
   
   nframe = inp_cluster->nframe;
   nato = inp_cluster->sele.nselatm;
   
   switch ( inp_cluster->distance )
   {
     case 1 :
       // NOTE: during matrix fillup RmsdCalc will alter frame position if super is selected.
       // as superimpose is run every time the rmsd is computed, this should not matter
       if( inp_cluster->himem == 0 )
       {
         refcoor = calloc( 3, sizeof(float *));
         refcoor[0] = calloc ( 3*nato, sizeof ( float  ));
         for ( ii =0; ii<3; ii++)
           refcoor[ii] = refcoor[0] + nato*ii;
         
         movcoor = calloc( 3, sizeof(float *));
         movcoor[0] = calloc ( 3*nato, sizeof ( float  ));
         for ( ii =0; ii<3; ii++)
           movcoor[ii] = movcoor[0] + nato*ii;
         
         nfloats = 3*nato;
         for( ii=0; ii<nframe; ii++ )
         { 
           offset = ii*nfloats*sizeof( float); 
           fseeko64( inp_cluster->tmpfile, offset, SEEK_SET );
           fread( refcoor[0], sizeof( float ), nfloats, inp_cluster->tmpfile );
           for( jj=0; jj<ii; jj++ )
           {
             offset = jj*nfloats*sizeof( float); 
             fseeko64(inp_cluster->tmpfile, offset, SEEK_SET );
             fread( movcoor[0], sizeof( float ), nfloats, inp_cluster->tmpfile );
             inp_cluster->matrix[ii][jj] = RmsdCalcKabsch3n( refcoor, movcoor, nato, inp_cluster->super );
           }
         }
       }
       else if( inp_cluster->himem == 1 )
       {
         refcoor = calloc( nato, sizeof(float *));
         refcoor[0] = calloc ( nato*3, sizeof ( float  ));
         for ( ii =0; ii<nato; ii++)
           refcoor[ii] = refcoor[0] + 3*ii;
         
         movcoor = calloc( nato, sizeof(float *));
         movcoor[0] = calloc ( nato*3, sizeof ( float  ));
         for ( ii =0; ii<nato; ii++)
           movcoor[ii] = movcoor[0] + 3*ii;
         
         for( ii=0; ii<nframe; ii++ )
         {  
           for( kk=0; kk<nato; kk++ )
           {
             refcoor[kk][0] = inp_cluster->xcoor[ii][kk];
             refcoor[kk][1] = inp_cluster->ycoor[ii][kk];
             refcoor[kk][2] = inp_cluster->zcoor[ii][kk];
           }
           for( jj=0; jj<ii; jj++ )
           {
             for( kk=0; kk<nato; kk++ )
             {
               movcoor[kk][0] = inp_cluster->xcoor[jj][kk];
               movcoor[kk][1] = inp_cluster->ycoor[jj][kk];
               movcoor[kk][2] = inp_cluster->zcoor[jj][kk];
             }
             inp_cluster->matrix[ii][jj] = RmsdCalc( refcoor, movcoor, nato, inp_cluster->super );
           }
         }
       }
       break;
     case 2 :
       refcoor = calloc( nato, sizeof(float *));
       refcoor[0] = calloc ( 3*nato, sizeof ( float  ));
       for ( ii =0; ii<nato; ii++)
         refcoor[ii] = refcoor[0] + 3*ii;
       
       movcoor = calloc( nato, sizeof(float *));
       movcoor[0] = calloc ( 3*nato, sizeof ( float  ));
       for ( ii =0; ii<nato; ii++)
         movcoor[ii] = movcoor[0] + 3*ii;
       for( ii=0; ii<nframe; ii++ )
       {
//       Calc_dist_mtx ( OPT, n , inp_c->selatm , dcd_crd , dist_mtx ) ;
         for( kk=0; kk<nato; kk++ )
         {
           refcoor[kk][0] = inp_cluster->xcoor[ii][kk];
           refcoor[kk][1] = inp_cluster->ycoor[ii][kk];
           refcoor[kk][2] = inp_cluster->zcoor[ii][kk];
         }
         for( jj=0; jj<ii; jj++ )
         {
           for( kk=0; kk<nato; kk++ )
           {
             movcoor[kk][0] = inp_cluster->xcoor[jj][kk];
             movcoor[kk][1] = inp_cluster->ycoor[jj][kk];
             movcoor[kk][2] = inp_cluster->zcoor[jj][kk];
           }
           inp_cluster->matrix[ii][jj] = DrmsCalc( refcoor, movcoor, nato, inp_cluster->pbcflag, inp_cluster->pbcbox );
         }
       }
       break;
     default :
       fprintf( stderr, "Unclear distance selection for clustering");
       exit(0);
       break;
   }
   
   free(refcoor[0]);
   free(refcoor);
   free(movcoor[0]);
   free(movcoor);
   
// fprintf( stdout, "Distance Matrix\n");
// for( ii=0; ii<nframe; ii++)
// {
//   for( jj=0; jj<ii; jj++)
//     fprintf( stdout, "%8.3f  ", inp_cluster->matrix[ii][jj]);
//   fprintf( stdout, "\n");
// }
// fflush(stdout);
   
   return;
}*/
// ------------------------------------------------------------------
#ifdef THREADED
void FillDistMatrix_threaded( struct inp_Cluster *inp_cluster )
{
   int    TT;
   struct thrFillMat    *FMdata;
   pthread_t            *pThreads;
   int                  *piExitStatus;
   
   pThreads = (pthread_t *) calloc ( inp_cluster->nthreads, sizeof(pthread_t));
   piExitStatus = (int *) calloc ( inp_cluster->nthreads, sizeof(int));
   FMdata = calloc( inp_cluster->nthreads, sizeof(struct threading));
   //printf( "DEBUG: thr called\n"); fflush(stdout);
   
   for( TT=0; TT<inp_cluster->nthreads; TT++)
   {
     FMdata[TT].inp_cluster = inp_cluster;
     FMdata[TT].thread_size = inp_cluster->nthreads;
     FMdata[TT].thread_rank = TT;
     if( inp_cluster->distance == 1 )
     {
       //printf( "DEBUG: going with thr_FillRMSDMatrix\n"); fflush(stdout);
       piExitStatus[TT] = pthread_create( &pThreads[TT], NULL, thr_FillRMSDMatrix, (void *) &FMdata[TT] );
     }
     else if( inp_cluster->distance == 2 )
       piExitStatus[TT] = pthread_create( &pThreads[TT], NULL, thr_FillDRMSMatrix, (void *) &FMdata[TT] );
   }
   //printf( "DEBUG: thr launched\n"); fflush(stdout);
     
   for(TT=0; TT<inp_cluster->nthreads; TT++)
     pthread_join(pThreads[TT], NULL);
   
   return;
}
// -----
void * thr_FillRMSDMatrix( void * pFMdata )
{
   struct thrFillMat *FMdata;
   int          ii, jj, kk;
   int          nframe, nato, t_rank, t_size;
   float      **refcoor, **movcoor;
   struct inp_Cluster *inp_cluster;
   
   //printf( "DEBUG: thr started\n"); fflush(stdout);
   FMdata = (struct thrFillMat *)pFMdata;
   inp_cluster = FMdata->inp_cluster;
   
   nframe = inp_cluster->nframe;
   nato = inp_cluster->sele.nselatm;
   t_rank = FMdata->thread_rank;
   t_size = FMdata->thread_size;
   
   //printf( "DEBUG: thr # %d/%d started\n", t_rank, t_size); fflush(stdout);
   
   refcoor = calloc( 3, sizeof(float *));
   refcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<3; ii++)
     refcoor[ii] = refcoor[0] + nato*ii;
   
   movcoor = calloc( 3, sizeof(float *));
   movcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<3; ii++)
     movcoor[ii] = movcoor[0] + nato*ii;
   
   
   for( ii=t_rank; ii<nframe; ii+=t_size )
   {  
     for( kk=0; kk<nato; kk++ )
     {
       refcoor[0][kk] = inp_cluster->xcoor[ii][kk];
       refcoor[1][kk] = inp_cluster->ycoor[ii][kk];
       refcoor[2][kk] = inp_cluster->zcoor[ii][kk];
     }
     for( jj=0; jj<ii; jj++ )
     {
       for( kk=0; kk<nato; kk++ )
       {
         movcoor[0][kk] = inp_cluster->xcoor[jj][kk];
         movcoor[1][kk] = inp_cluster->ycoor[jj][kk];
         movcoor[2][kk] = inp_cluster->zcoor[jj][kk];
       }
       inp_cluster->matrix[ii][jj] = RmsdCalc( refcoor, movcoor, nato, inp_cluster->super );
     }
   }
  
   return NULL;
}
// -----
void * thr_FillDRMSMatrix( void * pFMdata )
{
   struct thrFillMat *FMdata;
   int          ii, jj, kk;
   int          nframe, nato, t_rank, t_size;
   float      **refcoor, **movcoor;
   struct inp_Cluster *inp_cluster;
   
   FMdata = (struct thrFillMat *)pFMdata;
   inp_cluster = FMdata->inp_cluster;
   
   nframe = inp_cluster->nframe;
   nato = inp_cluster->sele.nselatm;
   t_rank = FMdata->thread_rank;
   t_size = FMdata->thread_size;
   
   refcoor = calloc( nato, sizeof(float *));
   refcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<nato; ii++)
     refcoor[ii] = refcoor[0] + 3*ii;
   
   movcoor = calloc( nato, sizeof(float *));
   movcoor[0] = calloc ( nato*3, sizeof ( float  ));
   for ( ii =0; ii<nato; ii++)
     movcoor[ii] = movcoor[0] + 3*ii;
   
   for( ii=t_rank; ii<nframe; ii+=t_size )
   {
     for( jj=0; jj<ii; jj++ )
     {
       for( kk=0; kk<nato; kk++ )
       {
         refcoor[kk][0] = inp_cluster->xcoor[ii][kk];
         refcoor[kk][1] = inp_cluster->ycoor[ii][kk];
         refcoor[kk][2] = inp_cluster->zcoor[ii][kk];
         movcoor[kk][0] = inp_cluster->xcoor[jj][kk];
         movcoor[kk][1] = inp_cluster->ycoor[jj][kk];
         movcoor[kk][2] = inp_cluster->zcoor[jj][kk];
         inp_cluster->matrix[ii][jj] = DrmsCalc( refcoor, movcoor, nato, inp_cluster->pbcflag, inp_cluster->pbcbox );
       }
     }
   }
   
   return NULL;
}
#endif
// ------------------------------------------------------------------
void MatrixCluster ( float **matrix, int msize, int method, _cluster *cluster1 )
{
   //int  ii, jj, kk;
   
   switch ( method )
   {
    case  1 :
     ClusterHiero( matrix, msize, cluster1);
     break;
    case  2 :
     ClusterQT( matrix, msize, cluster1);
     break;
    default :
     fprintf( stderr, "Unclear algorithm selection for clustering");
     exit(0);
     break;
   }
   
   return;
}
// -----
void ClusterHiero ( float **datamatrix, int nn, _cluster *cluster1 )
{
   int     ii, jj, /*jjtop,*/ kk, ll, /*mm,*/ hh, c1, c2, *sols, /*ncluster,*/ *list1 /**list2*/;
   float   min=0, max=0, cmax=0;
   float **matrix, threshold=0;
   //short **checkmatrix;
   int **checkmatrix;

   threshold = cluster1->threshold;
   // to spare memory, clustering is destructive
   // if multiple clustering is needed, use the temp file
   matrix = datamatrix;
   
   checkmatrix    = (int **) calloc( nn, sizeof(int *));
   if( checkmatrix == NULL )
   {
     fprintf( stderr, "Error while allocating memory (c0)\n");
     exit(0);
   }
   for( ii=0; ii<nn; ii++)
   {
     checkmatrix[ii] = (int *) calloc( ii, sizeof(int));
     if( checkmatrix[ii] == NULL )
     {
       fprintf( stderr, "Error while allocating memory (c1)\n");
       exit(0);
     }
   }
     
// for(ii=0; ii<nn; ii++ )
//  for(jj=0; jj<ii; jj++)
//   checkmatrix[ii][jj]  = (short) datamatrix[ii][jj];
   
   for(ii=0; ii<nn; ii++ )
     for(jj=0; jj<ii; jj++)
       checkmatrix[ii][jj]  = 1;
   
   sols = calloc( nn, sizeof(int));
   for(ii=0;ii<nn;ii++)
     sols[ii]=0;
   
   max=-1;
   /* find maximum distance */
   for(ii=0; ii<nn; ii++ )
     for(jj=0; jj<ii; jj++)
       if(matrix[ii][jj] > max)
         max = matrix[ii][jj];

   ii = jj = kk = ll = 0;
   // kk will be the number of clusters ?
   while(min<threshold && min<max)
   {
     min = max+1;
     for(ii=0; ii<nn-1; ii++ )
       for(jj=0; jj<ii; jj++)
         if(matrix[ii][jj] < min && checkmatrix[ii][jj]!=-1)
         {
           /* find minimum distance and involved structure pair */
           min = matrix[ii][jj];
           c1 = ii;
           c2 = jj;
         }

     if(min<threshold)
     {
      if( sols[c1]==0 && sols[c2]==0 )
      {
       /* if none of the structures is in a cluster, establish a new one */
       kk++;
       sols[c1] = kk;
       sols[c2] = kk;
       checkmatrix[c1][c2] = -1;
       // now update all other distances 
       for( ii=0; ii<nn; ii++)
       {
         // find max distance from frame ii to either c1 or c2
         if( c1>ii ) {
           cmax=matrix[c1][ii];
         } else {
           cmax = matrix[ii][c1];
         }
         if( c2>ii ) {
           if(matrix[c2][ii]>cmax) cmax = matrix[c2][ii];
         } else {
           if(matrix[ii][c2]>cmax) cmax = matrix[ii][c2];
         }

         // set all distances ii-c1/c2 to the max
         if( c1>ii ) {
           if(matrix[c1][ii]>0) matrix[c1][ii] = cmax;
         } else {
           if(matrix[ii][c1]>0) matrix[ii][c1] = cmax;
         }
         if( c2>ii ) {
           if(matrix[ii][c2]>0) matrix[ii][c2] = cmax;
         } else {
           if(matrix[c2][ii]>0) matrix[c2][ii] = cmax;
         }
       }
      }
      else if( (sols[c1]!=0 || sols[c2]!=0) && !(sols[c1]!=0 && sols[c2]!=0))
      {
       /* if one of structures belongs to a cluster, both now belong to that */
       if(sols[c1]!=0) {
         sols[c2]=sols[c1];
       } else if(sols[c2]!=0) {
         sols[c1]=sols[c2];
       }
       for( ii=0; ii<nn; ii++)
       {
         if( c1>ii ) {
           cmax=matrix[c1][ii];
         } else {
           cmax = matrix[ii][c1];
         }
         if( c2>ii ) {
           if(matrix[c2][ii]>cmax) cmax = matrix[c2][ii];
         } else {
           if(matrix[ii][c2]>cmax) cmax = matrix[ii][c2];
         }

//       cmax=matrix[c1][ii];
//       if(matrix[ii][c1]>cmax) cmax = matrix[ii][c1];
//       if(matrix[c2][ii]>cmax) cmax = matrix[c2][ii];
//       if(matrix[ii][c2]>cmax) cmax = matrix[ii][c2];
         for(jj=0; jj<nn; jj++)
         {
           if(sols[jj]==sols[c1])
           {
             if( ii>jj && matrix[ii][jj]>0) matrix[ii][jj] = cmax;
             if( ii<jj && matrix[jj][ii]>0) matrix[jj][ii] = cmax;
           }
         }
       }
       
       checkmatrix[c1][c2]=-1;
      }
      else if(sols[c1]!=0 && sols[c2]!=0 && sols[c1]!=sols[c2])
      {
        // if both already belong to some cluster, combine the clusters 
        // assign elements of cluster sols[c2] to cluster sols[c1]
        hh = sols[c2];
        for(ii=0; ii<nn; ii++)
          if(sols[ii]==hh)
            sols[ii]=sols[c1];
        // set every distance to the maximum
        for(ii=0; ii<nn-1; ii++)
        {
//       cmax=matrix[c1][ii];
//       if(matrix[ii][c1]>cmax) cmax = matrix[ii][c1];
//       if(matrix[c2][ii]>cmax) cmax = matrix[c2][ii];
//       if(matrix[ii][c2]>cmax) cmax = matrix[ii][c2];
          if( c1>ii ) {
            cmax=matrix[c1][ii];
          } else {
            cmax = matrix[ii][c1];
          }
          if( c2>ii ) {
            if(matrix[c2][ii]>cmax) cmax = matrix[c2][ii];
          } else {
            if(matrix[ii][c2]>cmax) cmax = matrix[ii][c2];
          }
         for(jj=0; jj<nn; jj++)
         {
           if(sols[jj]==sols[c1])
           {
             if( ii>jj && matrix[ii][jj]>0) matrix[ii][jj] = cmax;
             if( ii<jj && matrix[jj][ii]>0) matrix[jj][ii] = cmax;
           }
         }
        }
        checkmatrix[c1][c2]=-1;
        // COMPLETE THIS!
      }
      else if(sols[c1]!=0 && sols[c2]!=0 && sols[c1]==sols[c2])
      {
        checkmatrix[c1][c2]=-1;
      }
     }
     ll++;
   }
   
   
   // now fill up the cluster1 structure
   list1 = calloc( kk+1, sizeof(int));
   //list2 = calloc( kk+1, sizeof(int));
   cluster1->ncluster=0;
   //mm=0;
   
   list1[0] = 0;
   cluster1->ncluster = 1;
   for( ii=1; ii<=kk; ii++)
    for( jj=0; jj<nn; jj++ )
     if(sols[jj]==ii)
     {
      list1[cluster1->ncluster] = ii;
      cluster1->ncluster++;
      break;
     }

   
   cluster1->cluspop = calloc( cluster1->ncluster, sizeof(int));
   for( ii=0; ii<cluster1->ncluster; ii++)
    cluster1->cluspop[ii] = 0;
   
   for( ii=0; ii<cluster1->ncluster; ii++)
    for( jj=0; jj<nn; jj++ )
     if( sols[jj]==list1[ii] )
      cluster1->cluspop[ii]++;
   
   cluster1->clusstrs = calloc( cluster1->ncluster, sizeof( int *));
   for ( ii=0; ii < cluster1->ncluster ; ii++)
    cluster1->clusstrs[ii] = calloc( cluster1->cluspop[ii], sizeof(int));
   
   for( ii=0; ii<cluster1->ncluster; ii++)
   {
    hh=0;
    for( jj=0; jj<nn; jj++ )
     if(sols[jj]==list1[ii])
     {
      cluster1->clusstrs[ii][hh]=jj;
      hh++;
     }
   }

// printing (as a check) unclustered frames   
// for( ii=0; ii<cluster1->cluspop[0]; ii++)
//  printf("%d\n", cluster1->clusstrs[0][ii]+1);
   for( ii=0; ii<nn; ii++ )
     free(checkmatrix[ii]);
   free( checkmatrix );
   
   return;
}
// -----
int ClusterQT ( float **datamatrix, int nn, _cluster *cluster1 )
{
   int     ii, jj, kk, /*keepongoing,*/ clusstarter;
   //int     min, max, cmax, isolated;
   int    *unclustered;
   int     uu;
   float **matrix, threshold;//, checksum;
   short **checkmatrix;
   int     neighbors_here, *neighbors_list;
   int     max_neighbors, *max_neighbors_list;
   
   threshold = cluster1->threshold;

   // to spare memory, clustering is destructive
   // if multiple clustering is needed, use the temp file
   matrix = datamatrix;
   
   checkmatrix    = calloc( nn, sizeof(short *));
   for( ii=0; ii<nn; ii++)
     checkmatrix[ii] = calloc( ii, sizeof(short));
   
   for(ii=0; ii<nn; ii++ )
     for(jj=0; jj<ii; jj++)
       checkmatrix[ii][jj]  = 1;
   
   neighbors_list = calloc( nn, sizeof(int));
   max_neighbors_list = calloc( nn, sizeof(int));
   
   cluster1->ncluster = 0;
   cluster1->cluspop = calloc( nn, sizeof(int));
   for(ii=0;ii<nn;ii++)
    cluster1->cluspop[ii]=0;
   cluster1->clusstrs = calloc( nn, sizeof(int *));

   unclustered = calloc( nn, sizeof(int));
   for( ii=0; ii<nn; ii++)
    unclustered[ii] = 1;
   uu=nn;
   
   //ll=1;
   // loops until it finds any worthy cluster
   while ( 1 )
   {
     // reset neighbors count
     for( ii=0; ii<nn; ii++ )
       neighbors_list[ii] = 0;
     // count neighbors (dist<threshold) for each conformation and remember max
     max_neighbors = 0;
     for( ii=0; ii<nn; ii++ )
     {
       neighbors_list[0]=ii;
       neighbors_here = 1;      // frame itself is at least its own neighbor
       for( jj=0; jj<ii; jj++)
       {
         if(matrix[ii][jj]<threshold && checkmatrix[ii][jj]!=0)
         {
           neighbors_list[neighbors_here]=jj;
           neighbors_here++;
         }
       }
       for( jj=ii+1; jj<nn; jj++)
       {
         if(matrix[jj][ii]<threshold && checkmatrix[jj][ii]!=0)
         {
           neighbors_list[neighbors_here]=jj;
           neighbors_here++;
         }
       }
       if( max_neighbors < neighbors_here )
       {
         clusstarter = ii;
         max_neighbors = neighbors_here;
         for( jj=0; jj<max_neighbors; jj++)
           max_neighbors_list[jj] = neighbors_list[jj];
       }
     }
     // if bigger neighbors'list is with 0 member only, let's stop here
     if( max_neighbors<2 ) break;
     
     // pass to cluster1 structure infos about this cluster
     cluster1->ncluster ++;
     cluster1->cluspop[cluster1->ncluster] = max_neighbors;
     cluster1->clusstrs[cluster1->ncluster] = calloc(max_neighbors, sizeof(int));
     unclustered[clusstarter] = 0;
     for( ii=0; ii<max_neighbors; ii++ )
     {
      cluster1->clusstrs[cluster1->ncluster][ii] = max_neighbors_list[ii];
      unclustered[max_neighbors_list[ii]] = 0;
      uu--;
     }
     // erase confs from checkmatrix to proceed to next (smaller) cluster
     for( ii=0; ii<max_neighbors; ii++ )
     {
       for( jj=0; jj<max_neighbors_list[ii]; jj++ )
         checkmatrix[max_neighbors_list[ii]][jj] = 0;
       for( jj=max_neighbors_list[ii]+1; jj<nn; jj++ )
         checkmatrix[jj][max_neighbors_list[ii]] = 0;
       max_neighbors_list[ii] = 0;
     }
   }
   
   cluster1->ncluster ++;
   
   // take care of the "unclustered" cluster (#0)
   cluster1->cluspop[0] = 0;
   for( ii=0; ii<nn; ii++ )
     if( unclustered[ii] )
       cluster1->cluspop[0] ++;
   cluster1->clusstrs[0] = calloc( cluster1->cluspop[0], sizeof(int));
   kk = 0;
   for( ii=0; ii<nn; ii++ )
     if( unclustered[ii] )
     {
       cluster1->clusstrs[0][kk] = ii;
       kk++;
     }
   
   return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int CheckClusters( float **datamatrix, _cluster *cluster1 )
{
   int     ii, jj, kk;
   int     center;
   float  *rmsd, *centerrmsd, min, ndata;
   
   // bugging call for small nframes && method==qt
   cluster1->rmsdspread = calloc( cluster1->ncluster, sizeof(float));
   cluster1->cluscenter = calloc( cluster1->ncluster, sizeof( int ));
   for( ii=0; ii<cluster1->ncluster; ii++ )
    cluster1->cluscenter[ii] = -1;
   
   rmsd = calloc( cluster1->ncluster, sizeof(float));
   for( ii=0; ii<cluster1->ncluster; ii++ )
   rmsd[ii] = 0;
   
   centerrmsd = calloc( 1, sizeof(float));
   
   for( ii=0; ii<(cluster1->ncluster); ii++ )
   {
    free(centerrmsd);
    centerrmsd = calloc( cluster1->cluspop[ii], sizeof(float));
    for( jj=0; jj<cluster1->cluspop[ii]; jj++ )
     centerrmsd[jj] = 0;
    
    for( jj=0; jj<cluster1->cluspop[ii]; jj++ )
    {
     for( kk=jj; kk<cluster1->cluspop[ii]; kk++)
     {
      rmsd[ii] += datamatrix[cluster1->clusstrs[ii][jj]][cluster1->clusstrs[ii][kk]];
      rmsd[ii] += datamatrix[cluster1->clusstrs[ii][kk]][cluster1->clusstrs[ii][jj]];
     }
     for( kk=00; kk<cluster1->cluspop[ii]; kk++)
     {
      centerrmsd[jj] += datamatrix[cluster1->clusstrs[ii][jj]][cluster1->clusstrs[ii][kk]];
      centerrmsd[jj] += datamatrix[cluster1->clusstrs[ii][kk]][cluster1->clusstrs[ii][jj]];
     }
    }
    
    ndata = (float) ((((cluster1->cluspop[ii]*cluster1->cluspop[ii])-cluster1->cluspop[ii])/2));
    cluster1->rmsdspread[ii] = rmsd[ii]/ndata;
    
    min = rmsd[ii];
    center = -1;
    for( jj=0; jj<cluster1->cluspop[ii]; jj++ )
    {
     if( min>=centerrmsd[jj] )
     {
      min = centerrmsd[jj];
      center = cluster1->clusstrs[ii][jj];
     }
    }
    cluster1->cluscenter[ii] = center;
   }
   
   return 0;
}
//------------------------------------------------------------------------------

