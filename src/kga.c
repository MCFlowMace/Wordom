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
/*! \file kga.c
 \brief Kinetic Grouping Analysis module
 
 Source code for the Kinetic Grouping Analysis (KGA) and 
 Free Energy Profile (FEP) module.
*/
#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "kga.h"
//
int LogBin ( char **input, int inp_index, FILE *output )
{
   int    ii, jj;
   int    lmax = 17;
   int    lmin = 0;
   int    nbin;
   int    snapshot;
   float  lD;
   char   buffer[64];
   int    nlines;
   int    count;
   float *bins;
   float  raw;
   FILE  *data;
   
   struct inp_lb        inp_lb;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--BPD", 5))
     {
       sscanf(buffer, "--BPD %d ", &inp_lb.BpD);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CLUSFILE", 10))
     {
       sscanf(buffer, "--CLUSFILE %[^\n]%*c ", inp_lb.clusfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &inp_lb.target);
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   nbin= inp_lb.BpD * (lmax-lmin);
   lD = 1.0/inp_lb.BpD;
   bins = calloc ( 1+(lmax-lmin)*inp_lb.BpD, sizeof(float));
   
   data = O_File( inp_lb.clusfile, "r" );
   nlines = -1;
   while ( !feof(data))
   {
    nlines++;
    if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
   }
   rewind(data);
   
   count = 0;
   for(ii=0; ii<nlines; ii++)
   {
    fscanf( data, "%d", &snapshot);
    count ++;
    // added by Settanni as debug
    if(snapshot == 0 ) {count=0; continue;}
    // debug ends here
    if(snapshot == inp_lb.target) count = 0;
    if (count!=0)
     raw = log( count )/log(10);
    else 
     raw = log( 1 )/log(10);
    jj  = ( (raw-lmin)/lD );
    bins [ jj ]++;
   }
   
   for ( ii=0; ii<=nbin; ii++ )
    if ( bins [ ii ] > 0 )
     printf( "%8.3f %8.3f\n", pow( 10, lmin+ii*lD), -log(bins[ ii ]/nlines));
   
   fclose (data);
   exit(0);
   return 0;
}                
// ------------------------------------------------------------------
int Ka ( char **input, int inp_index, FILE *output )
{
   int    ii, jj, kk;
   char   buffer[64];
   int    nlines;
   int   *cluster;
   int   *target;
   int    listlength=0;
   
   float  pcomm;
   int   *start;
   int  **counter;
   int  **nn;
   int  **success;
   
   FILE  *data;
   
   int    index;
   int    change;
   int   *col;
   float  cutoff = 0.5;
   int    minnn = 99;
   int   *node1, *node2;
   
   int    allnodes;		// all the nodes (lots of them!)
   int    label;		// color to which the snapshot is converting
   int   *farben;		// array of colors of all 
   int  **nodeinfo;
   int    nbasins;
   int   *basinlist;
   int   *basinlist_1;
   int    max;
   int    maxcolor;
   int    population;
   
   struct inp_ka        inp_ka;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TCOMM", 7))
     {
       sscanf(buffer, "--TCOMM %d ", &inp_ka.tcomm );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CLUSFILE", 10))
     {
       sscanf(buffer, "--CLUSFILE %[^\n]%*c ", inp_ka.clusfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NNODES", 8))
     {
       sscanf(buffer, "--NNODES %d ", &inp_ka.nnodes );
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   printf ("tcomm: %6d\n",inp_ka.tcomm);
   
   data = O_File( inp_ka.clusfile, "r" );
   nlines = 0;
   while ( !feof(data))
   {
    nlines++;
    if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
    {
       fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
   }
   nlines--;

   rewind(data);
   cluster = calloc( nlines, sizeof(int));
   for(ii=0; ii<nlines; ii++)
    fscanf( data, "%d", &cluster[ii]);
   
   target = calloc( inp_ka.nnodes, sizeof(int));
   for( ii=0; ii<inp_ka.nnodes; ii++)
    target[ii] = ii+1;
   
   listlength = inp_ka.nnodes;
   
   start = calloc( listlength, sizeof(int));
   
   counter = calloc( listlength, sizeof(int *));
   for( ii=0; ii<listlength; ii++)
    counter[ii] = calloc( listlength, sizeof(int));
   
   nn = calloc( listlength, sizeof(int *));
   for( ii=0; ii<listlength; ii++)
    nn[ii] = calloc( listlength, sizeof(int));
   
   success = calloc( listlength, sizeof(int *));
   for( ii=0; ii<listlength; ii++)
    success[ii] = calloc( listlength, sizeof(int));
   
   for( ii=nlines-1; ii>-1; ii-- )
   {
    for( jj=0; jj<listlength; jj++ )
    {
     for( kk=0; kk<listlength; kk++ )
      counter[jj][kk]++;
     if(cluster[ii] == target[jj])
     {
      start[jj] = 1;
      for( kk=0; kk<listlength; kk++ )
      {
       if(start[kk] == 1)
       {
        nn[jj][kk]++;
	if( counter[jj][kk] <= inp_ka.tcomm )
	 success[jj][kk]++;
       }
       counter[kk][jj] = 0;
      }
     }
    }
    if (cluster[ii]==0)
    {
      for( jj=0; jj<listlength; jj++ )
      {
      	 for( kk=0; kk<listlength; kk++ )
      	{
      	   counter[kk][jj] += inp_ka.tcomm;
      	}
      }
    }
   }
   
   for( ii=0; ii<listlength; ii++)
    for( jj=0; jj<listlength; jj++)
     if( ii!=jj && nn[ii][jj]>0)
     {
      pcomm = (float)success[ii][jj]/(float)nn[ii][jj];
      printf("%8d %8d %8.5f %8d\n", target[ii], target[jj], pcomm, nn[ii][jj]);
     }
   
// now the painting (coloring!) of  part
   
   col = calloc( inp_ka.nnodes+1, sizeof(int));
   for( ii=1; ii<listlength+1; ii++ )
    col[ii] = ii;
   
   node1 = calloc( inp_ka.nnodes*inp_ka.nnodes, sizeof(int));
   node2 = calloc( inp_ka.nnodes*inp_ka.nnodes, sizeof(int));
   
   index = 0;
   for( ii=0; ii<listlength; ii++ )
    for( jj=0; jj<listlength; jj++ )
    {
     pcomm = (float)success[ii][jj]/(float)nn[ii][jj];
     if( pcomm>=cutoff && nn[ii][jj]>minnn )
     {
      node1[index] = target[ii];
      node2[index] = target[jj];
      index++;
     }
    }
   
   change = 1;
   while( change == 1 )
   {
    change=0;
    for( ii=0; ii<index; ii++ )
    {
     if( col[node1[ii]]>col[node2[ii]])
     {
      change = 1;
      col[node1[ii]]=col[node2[ii]];
     }
     else if( col[node1[ii]]<col[node2[ii]])
     {
      change = 1;
      col[node2[ii]]=col[node1[ii]];
     }
    }
   }
   
   fflush(stdout);
   
// all-node assignement to main basins (found in step above)
   
 // get the number of the basins (colors) and set up 2 maps
   nbasins = 0;
   basinlist = calloc(inp_ka.nnodes+1, sizeof(int));
   for( ii=1; ii<inp_ka.nnodes+1; ii++ )
    if( col[ii] == ii )
    {
     nbasins ++;
     basinlist[ii] = nbasins;
    }
   
   basinlist_1 = calloc( nbasins+1, sizeof(int));
   for( ii=1; ii<inp_ka.nnodes+1; ii++ )
    if( basinlist[ii] != 0 )
     basinlist_1[basinlist[ii]] = ii;
   
   
   allnodes = 0;
   for( ii=0; ii<nlines; ii++ )
   {
    if( cluster[ii] > allnodes )
     allnodes = cluster[ii];
   }
   
   nodeinfo = calloc( allnodes+1, sizeof(int *));
   for( ii=0; ii<allnodes+1; ii++ )
    nodeinfo[ii] = calloc( nbasins+1, sizeof(int));
   
   farben = calloc( allnodes+1, sizeof(int)); 
   for( ii=1; ii<allnodes+1; ii++ )
   {
    if( ii <= inp_ka.nnodes )
     farben[ii] = col[ii];
    else
     farben[ii] = 0;
   } 

   label = 0;
   index = 0;
   for( ii=nlines-1; ii>-1; ii-- )
   {
    
    if( farben[cluster[ii]] == 0 && label !=0 )
    {
     index++;
     if( index > inp_ka.tcomm || cluster[ii]==0) // CHANGED!! || cluster[ii]==0
      label = 0;
    }
    
    if( farben[cluster[ii]] != 0 )
    {
     index = 0;
     label = farben[cluster[ii]];
    }
    
    nodeinfo[cluster[ii]][basinlist[label]]++;
    
   }
   
   for( ii=0; ii<allnodes+1; ii++ )
   {
    population = 0;
    for( jj=0; jj<nbasins+1; jj++ )
     population+=nodeinfo[ii][jj];
    max = 0;
    for( jj=0; jj<nbasins+1; jj++ )
     if( nodeinfo[ii][jj]>max )
     {
      max=nodeinfo[ii][jj];
      maxcolor = basinlist_1[jj];
     }
    if( (float)max/(float)population >= 0.5 )
     farben[ii] = maxcolor;
    else if( (float)max/(float)population < 0.5 )
     farben[ii] = 0;
   }
   
   printf("+++++++++++++\n");
   for( ii=0; ii<allnodes+1; ii++ )
    printf("%6d %8d\n", ii, farben[ii]);
   
   
   
   exit(0);
}                
// ------------------------------------------------------------------
int Basin ( char **input, int inp_index, FILE *output ) 
{
   int     ii;
   char    buffer[64];
   int     max=0;
   int     nnodes=0; // total number of nodes
   int     nlines=0; // number of lines
   int     counter=0; 
   int     start=0;
           
   float   pcomm;
   int     *cluster; // of length 'nlines' 
   int     *success; // of length 'nnodes'
   int     *nn; // of length 'nnodes'
   
   FILE    *data;
   
   struct inp_basin     inp_basin;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TCOMM", 7))
     {
       sscanf(buffer, "--TCOMM %d ", &inp_basin.tcomm );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CLUSFILE", 10))
     {
       sscanf(buffer, "--CLUSFILE %[^\n]%*c ", inp_basin.clusfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &inp_basin.target );
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   data = O_File ( inp_basin.clusfile, "r" );
   nlines=0;
   while ( !feof(data))
   {
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   //nlines--;
   
   rewind(data);   
   
   cluster = calloc( nlines, sizeof(int));
   max=0;
   for(ii=0; ii<nlines; ii++)
   {
     fscanf( data, "%d", &cluster[ii]); // fill the array of the timeseries of clusters (as in earlier functions)
     if (cluster[ii]>max)
       max=cluster[ii];
   }
   
   nnodes=max; // the maximal rank is the number of nodes
   
   fclose (data);
   
   nn = calloc( nnodes+1, sizeof(int));
   success = calloc( nnodes+1, sizeof(int));
   
   for (ii=nlines-1;ii>=0;ii--)
   {
           if (cluster[ii] == inp_basin.target )
           {
             start=1;
             counter=0;
           }
           counter++;
           if (start==1)
           {
             nn[cluster[ii]]++;
                            
             if (counter <= inp_basin.tcomm)
             {
               success[cluster[ii]]++;
             }
           }
   }
   
   // print pcommit of all nodes to the target node, MAYBE WE CAN SKIP THAT?
   for (ii=1;ii<=nnodes;ii++)
   {
     if (nn[ii]>0)
     {
       pcomm = (float)success[ii]/(float)nn[ii];
       printf("%8d %8d %8.5f %8d\n", ii, inp_basin.target, pcomm, nn[ii]);
     }
   }
   
   printf("=============\n");
   // print only nodes that belong to the investigated basin, i.e., pcomm>=0.5
   
   printf("Nodes in basin %8d :\n", inp_basin.target);
   
   for (ii=1;ii<=nnodes;ii++)
   {
     if (nn[ii]>0)
     {
       pcomm = (float)success[ii]/(float)nn[ii];
       if (pcomm>=0.5)
         printf("%8d\n", ii);
     }
   }

   return 0; 
}
// ------------------------------------------------------------------
int MFPT ( char **input_text, int inp_index, FILE *output )
{
   int 	ii, jj, kk;
   char   	buffer[512];
   float 	max=0.0;
   int 	nnodes=0; // total number of nodes
   int 	nlines=0; // number of lines
   
   
   double   sum;
   double   cutoff, cut;
   double   totweight=0.0;
   
   int	*cluster; // of length 'nlines'
   
   int	*xx;
   int	*yy;
   
   double	*weight; // population of nodes
   double 	*solution;
   double	*solution_tmp;	
   
   FILE  	*data;
   
   int 	    nit     = 50000;  // (default) Number of iterations in the solvation of equations. User can change this with the -nit option (see below) 
   float    temp    = 300;    // (default) Temperature. User can change this with the -temp option (see below) 
   int 	    nonsymm = 1;
   int      target  = 1;
	 char     clusfile[512];
   
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   
   memset( buffer, '\0', 512);
   memset( clusfile, '\0', 512);
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input_text[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TEMP", 6))
     {
       sscanf(buffer, "--TEMP %f ", &temp );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CLUSFILE", 10))
     {
       sscanf(buffer, "--CLUSFILE %[^\n]%*c ", clusfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &target );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NIT", 5))
     {
       sscanf(buffer, "--NIT %d ", &nit );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NONSYMM", 9))
     {
       nonsymm = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SYMM", 6))
     {
       nonsymm = 0;
       gotit = 1;
     }
     if( gotit==0 )
     {
       fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_index++;
   }
   
   fprintf( stdout, "# clusfile: %s\n", clusfile);
   fprintf( stdout, "# target: %d\n", target);
   fprintf( stdout, "# temp: %8.3f\n", temp);
   fprintf( stdout, "# nit: %d\n", nit);

   data = O_File ( clusfile, "r" );

   while ( !feof(data))
   { 
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   nlines--;

   rewind(data);
   cluster = calloc( nlines, sizeof(int));

   for(ii=0; ii<nlines; ii++)
   {
     fscanf( data, "%d", &cluster[ii]); // fill the array of the timeseries of clusters (as in earlier functions)
     
     if (cluster[ii]>nnodes)
       nnodes=cluster[ii];
   }

   fclose (data);

   xx = calloc(2*nlines, sizeof(int));
   yy = calloc(2*nlines, sizeof(int));
   weight = calloc( nnodes+1, sizeof(double));

   weight[cluster[0]]+=0.5;
   for(ii=1; ii<nlines-1; ii++)
   {
           if (cluster[ii+1]==0 || cluster[ii-1]==0)
           {
              weight[cluster[ii]]+=0.5;
              totweight+=0.5;
           }
           
           else if( cluster[ii]!=0 )
           {
              weight[cluster[ii]]+=1.0;
              totweight+=1.0;
           }
   }

   weight[0]=0.0;


   for (ii=1;ii<nlines;ii++)
   {
     
     xx[2*ii]=cluster[ii-1];
     yy[2*ii]=cluster[ii];

     xx[2*ii+1]=cluster[ii];
     yy[2*ii+1]=cluster[ii-1];
     
     if ( nonsymm )
     {     
        xx[2*ii+1]=0;
        yy[2*ii+1]=0;
     }
                   
   }



   solution = calloc( nnodes+1, sizeof(double));
   solution_tmp = calloc( nnodes+1, sizeof(double));
           

   for (ii=0;ii<nit;ii++) 
   {
           
           for (jj=0;jj<nnodes+1;jj++)
              solution_tmp[jj]=1.0;

           for (kk=2;kk<2*nlines;kk++)
           {      
                  if (xx[kk]==0)
                     solution_tmp[xx[kk]]=0.0;
                  
                  
                  if (yy[kk]!=target && weight[xx[kk]]!=0 && weight[yy[kk]]!=0 && xx[kk]!=0 && yy[kk]!=0)
                  {
                     
                     solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]])*solution[yy[kk]];
                     
                     if (nonsymm==1)
                     {
                       solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]])*solution[yy[kk]];
                     }
                  }

           }
           solution_tmp[target]=0;
           
           max=0.0;
           
           for (kk=1;kk<nnodes+1;kk++){
              if (solution_tmp[kk]-solution[kk]>max)
              {
                 max=solution_tmp[kk]-solution[kk];
              }
              solution[kk]=solution_tmp[kk];
           }
           
           
           if (max==0.0)
           {
              ii=nit; // again: 1000 replaced by "inp_mfpt.nit"
           }               
   }



   for (ii=1;ii<nnodes+1;ii++)
   {
      cutoff=solution[ii];
      cut=0.0;
      sum=0;
      for (jj=2;jj<2*nlines-1;jj++)
      {
        if (xx[jj]!=0 && yy[jj]!=0)
        {
          if ((solution[xx[jj]]>=cutoff &&  solution[yy[jj]]<cutoff) || (solution[yy[jj]]>= cutoff &&  solution[xx[jj]]<cutoff))
          {
            cut+=0.5;
          }
        }
      }
      for (jj=1;jj<nnodes+1;jj++)
      {
        if (solution[jj]<cutoff)
        {
          sum+=(double)weight[jj];
        }
      }
           
      fprintf(stdout, "%8.4f %10.4f %10.4f %6d\n",  sum/totweight, -temp*0.002*log((double)cut/totweight), solution[ii],  ii); 
   }
           
    


   return 0;
}
// ------------------------------------------------------------------
int MFPTnet ( char **input, int inp_index, FILE *output )
{
   int     ii, jj, kk;
   int     aa, bb;
   float   cc;
   char    buffer[64];
   float   max=0.0;
   int     nnodes=0; // total number of nodes
   int     nlines=0; // number of lines


   double   sum;
   double   totweight=0.0;
   double   cutoff, cut;
           
   int     *xx;
   int     *yy;

   double  *weight; // population of nodes
   double   *linkweight; // weight of links
   double   *solution;
   double   *solution_tmp;  

   FILE    *data;

   int     nit=50000; // (default) Number of iterations in the solvation of equations. User can change this with the -nit option (see below) 
   float   temp = 300;// (default) Temperature. User can change this with the -temp option (see below) 
   int     nonsymm=1;


   struct inp_mfptnet   inp_mfptnet;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TEMP", 6))
     {
       sscanf(buffer, "--TEMP %f ", &inp_mfptnet.temp );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--LINKFILE", 10))
     {
       sscanf(buffer, "--LINKFILE %[^\n]%*c ", inp_mfptnet.linkfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &inp_mfptnet.target );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NIT", 5))
     {
       sscanf(buffer, "--NIT %d ", &inp_mfptnet.nit );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NONSYMM", 9))
     {
       inp_mfptnet.nonsymm = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SYMM", 6))
     {
       inp_mfptnet.nonsymm = 0;
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   if(inp_mfptnet.nit)
    nit = inp_mfptnet.nit;
   if(inp_mfptnet.temp)
    temp = inp_mfptnet.temp;
   if (inp_mfptnet.nonsymm)
    nonsymm=1;
    
   
   data = O_File ( inp_mfptnet.linkfile, "r" );

   while ( !feof(data))
   { 
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   nlines--;

   rewind(data);

   xx = calloc(2*nlines, sizeof(int));
   yy = calloc(2*nlines, sizeof(int));
   linkweight = calloc(2*nlines, sizeof(double));

   for(ii=0; ii<nlines; ii++)
   {
     fscanf( data, "%d %d %8f", &aa, &bb, &cc); // fill the array of the links
     xx[2*ii]=aa;
     yy[2*ii]=bb;
     linkweight[2*ii]=cc/2;
     
     if (nonsymm==0)
     {
        xx[2*ii+1]=bb;
        yy[2*ii+1]=aa;
        linkweight[2*ii+1]=cc/2;
     }
     if (nonsymm==1)
     {
        xx[2*ii+1]=aa;
        yy[2*ii+1]=bb;
        linkweight[2*ii+1]=cc/2;
     }
     
     
     if (aa>nnodes)
       nnodes=aa;
     if (bb>nnodes)
       nnodes=bb;
       
      totweight+=cc;
   }

   fclose (data);

   weight = calloc(nnodes+1, sizeof(double));

   for(ii=0; ii<2*nlines; ii++)
   {
      weight[xx[ii]]+=linkweight[ii];
   }

   solution = calloc( nnodes+1, sizeof(double));
   solution_tmp = calloc( nnodes+1, sizeof(double));

   for (ii=0;ii<nit;ii++) 
   {
           
           for (jj=0;jj<nnodes+1;jj++)
              solution_tmp[jj]=1.0;

           for (kk=0;kk<2*nlines;kk++)
           {      
                  
                  if (yy[kk]!=inp_mfptnet.target && weight[xx[kk]]!=0 && weight[yy[kk]]!=0)
                  {
                     solution_tmp[xx[kk]]+=linkweight[kk]/(weight[xx[kk]])*solution[yy[kk]];
                  }

           }
           solution_tmp[inp_mfptnet.target]=0;
           
           
           max=0.0;
           
           for (kk=1;kk<nnodes+1;kk++){
              if (solution_tmp[kk]-solution[kk]>max)
              {
                 max=solution_tmp[kk]-solution[kk];
              }
              solution[kk]=solution_tmp[kk];
           }
           
           
           if (max==0.0)
           {
              ii=nit; // again: 1000 replaced by "inp_mfptnet.nit"
           }               
   }

   for (ii=1;ii<nnodes+1;ii++)
   {
      cutoff=solution[ii];
      cut=0.0;
      sum=0;
      for (jj=0;jj<2*nlines-1;jj++)
      {
           
         if ((solution[xx[jj]]>=cutoff &&  solution[yy[jj]]<cutoff) || (solution[yy[jj]]>= cutoff &&  solution[xx[jj]]<cutoff))
         {
            cut+=linkweight[jj];
         }
      }
      for (jj=1;jj<nnodes+1;jj++)
      {
         if (solution[jj]<cutoff)
         {
            sum+=weight[jj];
         }
      }
      if (weight[ii]>0)
      {
         fprintf(stdout, "%8.4f %10.4f %10.4f %6d\n",  sum/totweight, -temp*0.002*log((double)cut/totweight), solution[ii],  ii); 
      }
   }
           
   return 0;
}
// ------------------------------------------------------------------
int PFoldf ( char **input, int inp_index, FILE *output )
{
   int     ii, jj, kk;
   char    buffer[64];
   float   max=0.0;
   int     nnodes=0; // total number of nodes
   int     nlines=0; // number of lines


   double   sum;
   double   cutoff, cut;
   double   totweight=0.0;
   float   lambda=0.0001; // (default)

   int     *cluster; // of length 'nlines'

   int     *xx;
   int     *yy;
           
   double   *weight; // population of nodes
   double   *solution;
   double   *solution_tmp;  

   FILE    *data;

   int     nit=50000; // (default) Number of iterations in the solvation of equations. User can change this with the -nit option (see below) 
   float   temp = 300;// (default) Temperature. User can change this with the -temp option (see below) 
   int     target2=0;// (default)
   int     nonsymm=1;// (default)

   struct inp_pfoldf    inp_pfoldf;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TEMP", 6))
     {
       sscanf(buffer, "--TEMP %f ", &inp_pfoldf.temp );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--LAMBDA", 8))
     {
       sscanf(buffer, "--LAMBDA %f ", &inp_pfoldf.lambda );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CLUSFILE", 10))
     {
       sscanf(buffer, "--CLUSFILE %[^\n]%*c ", inp_pfoldf.clusfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET2", 9))
     {
       sscanf(buffer, "--TARGET2 %d ", &inp_pfoldf.target2 );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &inp_pfoldf.target );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NIT", 5))
     {
       sscanf(buffer, "--NIT %d ", &inp_pfoldf.nit );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NONSYMM", 9))
     {
       inp_pfoldf.nonsymm = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SYMM", 6))
     {
       inp_pfoldf.nonsymm = 0;
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   if(inp_pfoldf.nit)
    nit = inp_pfoldf.nit;
   if(inp_pfoldf.temp)
    temp = inp_pfoldf.temp;
   if (inp_pfoldf.nonsymm)
    nonsymm=1;
   
   lambda=inp_pfoldf.lambda; // default set in inp_pfoldf.lambda is 0.0001
   if( inp_pfoldf.target2 != 0 )
    lambda = 0.0;
   
   data = O_File ( inp_pfoldf.clusfile, "r" );

   while ( !feof(data))
   { 
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   nlines--;

   rewind(data);
   cluster = calloc( nlines, sizeof(int));

   for(ii=0; ii<nlines; ii++)
   {
     fscanf( data, "%d", &cluster[ii]); // fill the array of the timeseries of clusters (as in earlier functions)
     
     if (cluster[ii]>nnodes)
       nnodes=cluster[ii];
   }

   fclose (data);

   if ( inp_pfoldf.target2 == 0 )
    target2 = nnodes+1;
   else target2 = inp_pfoldf.target2;
    

   xx = calloc(4*nlines, sizeof(int));
   yy = calloc(4*nlines, sizeof(int));
   weight = calloc( nnodes+2, sizeof(double));

// weight[cluster[0]]+=0.5;
   for(ii=0; ii<nlines-1; ii++)
   {
           if (nonsymm==1){
              if (cluster[ii+1]!=0 && cluster[ii]!=0)
              {
                 weight[cluster[ii]]+=1.0;
                 totweight+=1.0;
              }
           }
           if (nonsymm==0){
           
              if (cluster[ii+1]==0 || cluster[ii-1]==0)
              {
                 weight[cluster[ii]]+=0.5;
                 totweight+=0.5;
              }
           
              else if( cluster[ii]!=0 )
              {
                 weight[cluster[ii]]+=1.0;
                 totweight+=1.0;
              }
           }
   }

   weight[0]=0.0;
   weight[nnodes+1]=totweight*lambda;


   for (ii=1;ii<nlines;ii++)
   {
     
     xx[4*ii]=cluster[ii-1];
     yy[4*ii]=cluster[ii];

     xx[4*ii+1]=cluster[ii];
     yy[4*ii+1]=cluster[ii-1];
      
     xx[4*ii+2]=cluster[ii];
     yy[4*ii+2]=nnodes+1;
     
     xx[4*ii+3]=nnodes+1;
     yy[4*ii+3]=cluster[ii];

     if (nonsymm==1)
     {     
       xx[4*ii+1]=0;
       yy[4*ii+1]=0;
     
     }
   }


   solution = calloc( nnodes+2, sizeof(double));
   solution_tmp = calloc( nnodes+2, sizeof(double));
           
   solution[inp_pfoldf.target]=1.0; //inp_pfoldf.target1 instead of inp_pfoldf.target

   for (ii=0;ii<nit;ii++) 
   {
           
           for (jj=0;jj<nnodes+2;jj++)
           {
              solution_tmp[jj]=0.0;
           }


           for (kk=4;kk<4*nlines;kk++)
           {
                  
                  if (xx[kk]==0)
                     solution_tmp[xx[kk]]=0.0;
                  
                  if (weight[xx[kk]]!=0 && weight[yy[kk]]!=0 && xx[kk]!=0 && yy[kk]!=0)
                  {
                     if (xx[kk]!=nnodes+1 && yy[kk]!=nnodes+1)
                        //solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]])*solution[yy[kk]];
                        solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]]*(1.0+lambda))*solution[yy[kk]];
                     /*
                     if (xx[kk]==nnodes+1)
                        solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]])*solution[yy[kk]];
                     if (yy[kk]==nnodes+1)
                        solution_tmp[xx[kk]]+=0.5*lambda/(weight[xx[kk]])*solution[yy[kk]];
                     */
                     if (nonsymm==1)
                     {
                        if (xx[kk]!=nnodes+1 && yy[kk]!=nnodes+1)
                           //solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]])*solution[yy[kk]];
                           solution_tmp[xx[kk]]+=0.5/(weight[xx[kk]]*(1.0+lambda))*solution[yy[kk]];
                     }
                  }
                   
           }
           solution_tmp[inp_pfoldf.target]=1.0; //inp_pfoldf.target1 instead of inp_pfoldf.target
           
           solution_tmp[target2]=0; //inp_pfoldf.target2 instead of nnodes+1
           
           
           max=0.0;
           
           for (kk=1;kk<nnodes+2;kk++){
              if (solution_tmp[kk]-solution[kk]>max)
              {
                 max=solution_tmp[kk]-solution[kk];
              }
              solution[kk]=solution_tmp[kk];
           }
           
           
           if (max==0.0)
           {
              ii=nit; // again: 1000 replaced by "inp_pfoldf.nit"
           }               
   }



   for (ii=1;ii<nnodes+1;ii++)
   {
      cutoff=solution[ii];
      cut=0.0;
      sum=0;
      for (jj=4;jj<4*nlines-1;jj++)
      {
           if (xx[jj]!=0 && yy[jj]!=0 && xx[jj]!=nnodes+1 && yy[jj]!=nnodes+1)
           {
              if ((solution[xx[jj]]>cutoff &&  solution[yy[jj]]<=cutoff) || (solution[yy[jj]]> cutoff &&  solution[xx[jj]]<=cutoff))
              {
                 cut+=0.5;
              }
           }
      }
      for (jj=1;jj<nnodes+1;jj++)
      {
         if (solution[jj]>cutoff)
         {
            sum+=(double)weight[jj]/(double)(1+lambda); //The final plot is made on the original weights, not those with extra node
         }
      }
           
      printf("%8.4lf %10.4lf %16.10lf %6d\n",  sum/(totweight), -temp*0.002*log((double)cut/(totweight)), solution[ii], ii);  
   }
           
    


   return 0;
}
// ------------------------------------------------------------------
int PFoldfnet ( char **input, int inp_index, FILE *output )
{
   int     ii, jj, kk;
   int     aa, bb;
   float   cc;
   char    buffer[64];
   float   max=0.0;
   int     nnodes=0; // total number of nodes
   int     nlines=0; // number of lines

   double   sum=0.0;
   double   cutoff, cut;
   double   totweight=0.0;
   float   lambda=0.0001; // (default)
           
   int     *xx;
   int     *yy;
           
   double   *weight; // population of nodes
   double   *linkweight; // population of nodes
   double   *solution;
   double   *solution_tmp;  

   FILE    *data;

   int     nit=50000; // (default) Number of iterations in the solvation of equations. User can change this with the -nit option (see below) 
   float   temp = 300;// (default) Temperature. User can change this with the -temp option (see below) 
   int     target2=0;// (default)
   int     nonsymm=1;// (default)

   struct inp_pfoldfnet      inp_pfoldfnet;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TEMP", 6))
     {
       sscanf(buffer, "--TEMP %f ", &inp_pfoldfnet.temp );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--LAMBDA", 8))
     {
       sscanf(buffer, "--LAMBDA %f ", &inp_pfoldfnet.lambda );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--LINKFILE", 10))
     {
       sscanf(buffer, "--LINKFILE %[^\n]%*c ", inp_pfoldfnet.linkfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET2", 9))
     {
       sscanf(buffer, "--TARGET2 %d ", &inp_pfoldfnet.target2 );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--TARGET", 8))
     {
       sscanf(buffer, "--TARGET %d ", &inp_pfoldfnet.target );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NIT", 5))
     {
       sscanf(buffer, "--NIT %d ", &inp_pfoldfnet.nit );
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NONSYMM", 9))
     {
       inp_pfoldfnet.nonsymm = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SYMM", 6))
     {
       inp_pfoldfnet.nonsymm = 0;
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   if(inp_pfoldfnet.nit)
    nit = inp_pfoldfnet.nit;
   if(inp_pfoldfnet.temp)
    temp = inp_pfoldfnet.temp;
   if (inp_pfoldfnet.nonsymm)
    nonsymm=1;
   
   lambda=inp_pfoldfnet.lambda; // default set in inp_pfoldfnet.lambda is 0.0001
   if( inp_pfoldfnet.target2 != 0 )
    lambda = 0.0;
   
   data = O_File ( inp_pfoldfnet.linkfile, "r" );

   while ( !feof(data))
   { 
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   nlines--;

   rewind(data);

   xx = calloc(4*nlines, sizeof(int));
   yy = calloc(4*nlines, sizeof(int));
   linkweight = calloc(4*nlines, sizeof(double));


   for(ii=0; ii<nlines; ii++)
   {
     fscanf( data, "%d %d %8f", &aa, &bb, &cc); // fill the array of the links
     xx[2*ii]=aa;
     yy[2*ii]=bb;
     linkweight[2*ii]=cc/2;
     
     
     if (nonsymm==0)
     {
        xx[2*ii+1]=bb;
        yy[2*ii+1]=aa;
        linkweight[2*ii+1]=cc/2;
     }
     if (nonsymm==1)
     {
        xx[2*ii+1]=aa;
        yy[2*ii+1]=bb;
        linkweight[2*ii+1]=cc/2;
     }
     
     if (aa>nnodes)
        nnodes=aa;
     if (bb>nnodes)
        nnodes=bb;
       
     totweight+=cc;
             
   }

   fclose (data);



   weight = calloc(nnodes+2, sizeof(double));

   for(ii=0; ii<2*nlines; ii++)
   {
      weight[xx[ii]]+=linkweight[ii];
   }
   weight[nnodes+1]=totweight*lambda;


   for (ii=1;ii<nnodes+1;ii++){
           xx[2*nlines+2*ii-1]=ii;
           yy[2*nlines+2*ii-1]=nnodes+1;
           linkweight[2*nlines+2*ii-1]=lambda*weight[ii];
           
           
           xx[2*nlines+2*ii]=nnodes+1;
           yy[2*nlines+2*ii]=ii;
           linkweight[2*nlines+2*ii]=lambda*weight[ii];
           
   }


   sum=0.0;
   for(ii=1; ii<nnodes+1; ii++)
   {
      weight[ii]+=lambda*weight[ii];
   }


   if ( inp_pfoldfnet.target2 == 0 )
    target2 = nnodes+1;
   else target2 = inp_pfoldfnet.target2;



   solution = calloc( nnodes+2, sizeof(double));
   solution_tmp = calloc( nnodes+2, sizeof(double));
           
   solution[inp_pfoldfnet.target]=1.0; //inp_pfoldfnet.target1 instead of inp_pfoldfnet.target

   for (ii=0;ii<nit;ii++) 
   {
           
           for (jj=0;jj<nnodes+2;jj++)
              solution_tmp[jj]=0.0;

           for (kk=0;kk<4*nlines;kk++)
           {      
                  
                  if (weight[xx[kk]]!=0 && weight[yy[kk]]!=0)
                  {
                     solution_tmp[xx[kk]]+=linkweight[kk]/(weight[xx[kk]])*solution[yy[kk]];
                  }

           }
           solution_tmp[inp_pfoldfnet.target]=1.0;
           solution_tmp[target2]=0.0;
           
           
           max=0.0;
           
           for (kk=1;kk<nnodes+2;kk++){
              if (solution_tmp[kk]-solution[kk]>max)
              {
                 max=solution_tmp[kk]-solution[kk];
              }
              solution[kk]=solution_tmp[kk];
           }
           
           if (max==0.0)
           {
              ii=nit; // again: 1000 replaced by "inp_pfoldfnet.nit"
           }               
   }




   for (ii=1;ii<nnodes+1;ii++)
   {
      cutoff=solution[ii];
      cut=0.0;
      sum=0;
      for (jj=0;jj<4*nlines-1;jj++)
      {
           if (xx[jj]!=nnodes+1 && yy[jj]!=nnodes+1)
           {
              if ((solution[xx[jj]]>cutoff &&  solution[yy[jj]]<=cutoff) || (solution[yy[jj]]> cutoff &&  solution[xx[jj]]<=cutoff))
              {
                 cut+=linkweight[jj];
              }
           }
      }
      for (jj=1;jj<nnodes+1;jj++)
      {
         if (solution[jj]>cutoff)
         {
            sum+=(double)weight[jj]/(double)(1+lambda); //The final plot is made on the original weights, not those with extra node
         }
      }
           
//    printf("%8.4f %10.4f %6d\n",  sum/(totweight), -temp*0.002*log((float)cut/(totweight)), ii); 
      printf("%8.4lf %10.4lf %16.10lf %6d\n",  sum/(totweight), -temp*0.002*log((double)cut/(totweight)), solution[ii], ii); 
   }
           
    


   return 0;
}
// ------------------------------------------------------------------
int Equil ( char **input, int inp_index, FILE *output )
{
   int aa, bb;
   int ii, jj, ll;
   int nnodes=0;
   int nlines=0;
   double cc;
   double totweight;

   char    buffer[64];

   int *xx;
   int *yy;

   double max=0.0;
   double soltot=0.0;

   double *weight;
   double *linkweight;
   double *p_ij;
   double *solution;
   double *solution_tmp;


   int    nit=100000;
   
   struct inp_equil      inp_equil;
   int                  gotit;
   
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--LINKFILE", 10))
     {
       sscanf(buffer, "--LINKFILE %[^\n]%*c ", inp_equil.linkfile);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NIT", 5))
     {
       sscanf(buffer, "--NIT %d ", &inp_equil.nit );
       gotit = 1;
     }
     if( gotit==0 )
     {
      fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
      exit(5);
     }
     inp_index++;
   }
   
   if(inp_equil.nit)
    nit = inp_equil.nit;

   FILE    *data;

   data = O_File ( inp_equil.linkfile, "r" );

   while ( !feof(data))
   { 
     nlines++;
     if(fgets(buffer, 64, data)==NULL && ( !feof(data) || ferror(data) ))
     {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
     }
   }
   nlines--;

   rewind(data);

   xx = calloc(nlines, sizeof(int));
   yy = calloc(nlines, sizeof(int));
   linkweight = calloc(nlines, sizeof(double));


   for(ii=0; ii<nlines; ii++)
   {
      fscanf( data, "%d %d %12lf", &aa, &bb, &cc);
      xx[ii]=aa;
      yy[ii]=bb;
      linkweight[ii]=cc;
      totweight+=linkweight[ii];
      if (xx[ii]>nnodes)
      {
           nnodes=aa;
      }
      if (yy[ii]>nnodes){
           nnodes=bb;
      }
      
   }

   fclose (data);



   weight = calloc(nnodes+1, sizeof(double));
   solution = calloc(nnodes+1, sizeof(double));
   solution_tmp = calloc(nnodes+1, sizeof(double));
   p_ij = calloc(nlines, sizeof(double));

   for (ii=0;ii<nlines;ii++)
   {
      weight[xx[ii]]+=linkweight[ii];      
   }

   for (ii=0;ii<nlines;ii++)
   {
      p_ij[ii]=linkweight[ii]/weight[xx[ii]];
   }

   for (ii=0;ii<nnodes+1;ii++)
   {
           solution[ii]=weight[ii]/totweight; 
   }


   for (ll=1;ll<=nit;ll++)
   {
           
           for (ii=0;ii<nnodes+1;ii++)
           {
              solution_tmp[ii]=0;
           }
           soltot=0.0;

           for (jj=0;jj<nlines;jj++)
           {
                   
              solution_tmp[yy[jj]]+=p_ij[jj]*solution[xx[jj]];
              soltot+=p_ij[jj]*solution[xx[jj]];   
           }
                                   
           
           max=0.0;
           
           for (ii=0;ii<nnodes+1;ii++)
           {
              if (solution_tmp[ii]-solution[ii]>max)
              {
                 max=solution_tmp[ii]-solution[ii];
              }
              solution[ii]=solution_tmp[ii];
           }
           
           
           if (max==0.0)
           {
              ll=nit+1;
           }

   }
// for (ii=0;ii<nlines;ii++)
// {
//    solution[ii]=(double)solution[ii]/(double)soltot;
// }

   printf("%8f\n",max);

   for (ii=1;ii<=nnodes;ii++)
   {
      printf("NODE  %6d %16.8lf %14.4lf\n", ii, solution[ii], solution[ii]*totweight);
   } 
                                                           
   for (ii=0;ii<nlines;ii++)
   {
      printf("LINK  %6d %8d %14.2lf\n", xx[ii], yy[ii], (p_ij[ii]*solution[xx[ii]]*totweight));    
   }

           
   return 0; 
   
}
// =======================================================
intlist * readSymbTrj( char *filename )
{
  int  ii, nlines;
  FILE *inputfile;
  char  buffer[64];
  intlist  *strj;
  
  inputfile = fopen( filename, "r");
  nlines = -1;
  while ( !feof(inputfile))
  {
   nlines++;
   if(fgets(buffer, 64, inputfile)==NULL && ( !feof(inputfile) || ferror(inputfile) ))
   {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
   }
  }
  rewind(inputfile);
   
  strj = (intlist *)calloc( 1, sizeof(intlist));
  strj->data = (int *)calloc( nlines, sizeof(int));
  strj->ndata = nlines;
  for( ii=0; ii<nlines; ii++)
    fscanf( inputfile, "%d", &strj->data[ii]);
  
  return strj;
}
//--------------------------------------------------------------------==
graph * sTrj2Graph( intlist *strj )
{
  int  ii;
  int  nnodes=0;
  graph *grapho;
  
  for( ii=0; ii<strj->ndata; ii++)
    if( strj->data[ii] > nnodes )
      nnodes = strj->data[ii];
  
  grapho = (graph *)calloc( 1, sizeof(graph));
  grapho->nodes = NULL;
  grapho->leenklist = NULL;
  
  for( ii=0; ii<strj->ndata-1; ii++)
    if( strj->data[ii]!= 0 && strj->data[ii+1] != 0 )
      addtograph1( grapho, strj->data[ii], strj->data[ii+1] );
  
  if( nnodes != grapho->nnodes )
  {
    printf( "something's amiss\n");
    exit(0);
  }
  
  return grapho;
}
//--------------------------------------------------------------------==
int imin( int num1, int num2 )
{
  if( num1>num2 )
    return num2;
  else
    return num1;
}
//--------------------------------------------------------------------==
int whereinlist( intlist *list, int number )
{
  int   ii;
  
  for( ii=0; ii<list->ndata; ii++)
    if( list->data[ii] == number )
      return ii;
  return -1;
}
//--------------------------------------------------------------------==
int addtolist( intlist *list, int number )
{
  list->ndata ++;
  //list->data = realloc( list->data, list->ndata * sizeof(int) );
  list->data[list->ndata-1] = number;
  return 0;
}
//--------------------------------------------------------------------==
int popfromlist( intlist *list )
{
  list->ndata --;
  return list->data[list->ndata];
}
//--------------------------------------------------------------------==
void printlist( intlist *ilist )
{
  int   ii;
  
  if( ilist == NULL)
  {
    printf(" NULL found!\n");
    exit(0);
  }
  
  printf(">>>{%d}:   [", ilist->ndata);
  for(ii=0; ii<ilist->ndata; ii++)
   printf("%d, ", ilist->data[ii] );
  printf("\b\b]\n");
  
  return;
}
//--------------------------------------------------------------------==
int addtodlist( intdlist *dlist, int key, int value )
{
  int     ii;
  //int   * tmp;
  
  for( ii=0; ii<dlist->ndata; ii++)
  {
    if( dlist->keys[ii] == key )
    {
      dlist->values[ii] = value;
      return 0;
    }
  }
  
  dlist->ndata ++;
  //tmp = realloc( dlist->keys, dlist->ndata );
  //if( !tmp )
  //{
    //printf("Erorr jdhfvgf\n");
    //exit(0);
  //}
  //dlist->keys = tmp;
  //tmp = NULL;
  //tmp = realloc( dlist->values, dlist->ndata );
  //if( !tmp )
  //{
    //printf("Erorr jdhfvgf2\n");
    //exit(0);
  //}
  //dlist->values = tmp;
  
  dlist->keys[dlist->ndata-1] = key;
  dlist->values[dlist->ndata-1] = value;
  
  return 0;
}
//--------------------------------------------------------------------==
int getfromdlist( intdlist *dlist, int key )
{
  int     ii, value = -1;
  
  for( ii=0; ii<dlist->ndata; ii++ )
  {
    if( dlist->keys[ii] == key )
    {
      value = dlist->values[ii];
      break;
    }
  }
  
  return value;
}
//--------------------------------------------------------------------==
void printdlist( intdlist *dlist )
{
  int   ii;
  
  printf("{");
  for(ii=0; ii<dlist->ndata; ii++)
   printf("%d: %d, ", dlist->keys[ii], dlist->values[ii]);
  printf("\b\b}\n");
  
  return;
}
//--------------------------------------------------------------------==
int whereingraph( graph *graph1, int number )
{
  int   ii;
  
  for( ii=0; ii<graph1->nnodes; ii++)
    if( graph1->nodes[ii] == number )
      return ii;
  return -1;
}
//--------------------------------------------------------------------==
int addtograph1( graph *lgraph, int key, int value) // this adds to the current value of the 'value'
{
  int     ii;
  
  for ( ii=0; ii<lgraph->nnodes; ii++)
  {
    if( lgraph->nodes[ii] == key )
    {
      if( whereinlist( &lgraph->leenklist[ii], value) != -1 )
        return 0;
      lgraph->leenklist[ii].ndata +=1;
      lgraph->leenklist[ii].data = realloc( lgraph->leenklist[ii].data, lgraph->leenklist[ii].ndata * sizeof(int) );
      lgraph->leenklist[ii].data[lgraph->leenklist[ii].ndata-1] = value;
      return 0;
    }
  }

  lgraph->nnodes ++;
  lgraph->nodes = realloc( lgraph->nodes, lgraph->nnodes*sizeof(int) );
  lgraph->nodes[lgraph->nnodes-1] = key;
  lgraph->leenklist = realloc( lgraph->leenklist, lgraph->nnodes*sizeof(intlist) );
  lgraph->leenklist[lgraph->nnodes-1].ndata = 1;
  lgraph->leenklist[lgraph->nnodes-1].data = (int *)calloc( 1, sizeof(int));
  lgraph->leenklist[lgraph->nnodes-1].data[0] = value;
  
  return 0;
}
//--------------------------------------------------------------------==
int addtograph2( graph *lgraph, int key, int value) // this changes the current value of the 'value'
{
  int       ii;
  int     * itmp;
  intlist * tmp;
  
  for ( ii=0; ii<lgraph->nnodes; ii++)
  {
    if( lgraph->nodes[ii] == key )
    {
      lgraph->leenklist[ii].data[0] = value;
      return 0;
    }
  }
  
  lgraph->nnodes ++;
  itmp = realloc( lgraph->nodes, lgraph->nnodes*sizeof(int) );
  if( itmp == NULL) 
  {
    printf( "error!\n");
    exit(0);
  }
  lgraph->nodes = itmp;
  
  lgraph->nodes[lgraph->nnodes-1] = key;
  
  tmp = (intlist *)realloc( lgraph->leenklist, (lgraph->nnodes+1)*sizeof(int) );
  if( tmp != NULL)
    lgraph->leenklist = tmp;
  else
    printf( "error 3\n");
  lgraph->leenklist[lgraph->nnodes-1].ndata = 1;
  lgraph->leenklist[lgraph->nnodes-1].data = (int *)calloc( 1, sizeof(int));
  lgraph->leenklist[lgraph->nnodes-1].data[0] = value;
  
  return 0;
}
//--------------------------------------------------------------------==
intlist * getgraphval( graph *lgraph, int key )
{
  int       ii;
  
  for ( ii=0; ii<lgraph->nnodes; ii++)
  {
    if( lgraph->nodes[ii] == key )
    {
      return &lgraph->leenklist[ii];
    }
  }
  return NULL;
}
//----------------------------------------------------------------------
void printgraph( graph *grapho )
{
  int ii, jj;
  
  for( ii=0; ii<grapho->nnodes; ii++)
  {
    printf(" %d:\n", grapho->nodes[ii]);
    printf("\t");
    for( jj=0; jj<grapho->leenklist[ii].ndata; jj++)
      printf(" %d; ", grapho->leenklist[ii].data[jj]);
    printf("\n");
  }
  
  return;
}
//----------------------------------------------------------------------
graph * Tarjan_ite( graph *grapho )
{
  int                 ii, jj, hh, ss, min;
  graph *             components;
  graph               todo;
  intdlist            number;   // list of int couples, used as dictionary (int key, int value)
  intdlist            low;      // list of int couples, used as dictionary (int key, int value)
  intlist             explored, stack;
  
  // set up some of the structures
  //number.ndata = 0;
  //number.keys = NULL;
  //number.values = NULL;
  //low.ndata = 0;
  //low.keys = NULL;
  //low.values = NULL;
  number.ndata = 0;
  number.keys = (int *)calloc( grapho->nnodes, sizeof(int));
  number.values = (int *)calloc( grapho->nnodes, sizeof(int));
  low.ndata = 0;
  low.keys = (int *)calloc( grapho->nnodes, sizeof(int));
  low.values = (int *)calloc( grapho->nnodes, sizeof(int));
  explored.ndata = 0;
  explored.data = (int *)calloc( grapho->nnodes, sizeof(int));
  stack.ndata = 0;
  stack.data = (int *)calloc( grapho->nnodes, sizeof(int));
  
  components = calloc( 1, sizeof( graph));

  // copy grapho to local copy, ie todo
  todo.nnodes = grapho->nnodes;
  todo.nodes = calloc( todo.nnodes, sizeof(int));
  todo.leenklist = (intlist *)calloc( todo.nnodes, sizeof( intlist));
  for( ii=0; ii<todo.nnodes; ii++)
  {
    todo.nodes[ii] = grapho->nodes[ii];
    todo.leenklist[ii].ndata = grapho->leenklist[ii].ndata;
    todo.leenklist[ii].data = (int *)calloc( todo.leenklist[ii].ndata, sizeof(int));
    for( jj=0; jj<grapho->leenklist[ii].ndata; jj++ )
      todo.leenklist[ii].data[jj] = grapho->leenklist[ii].data[jj];
  }
  
  for( ii=0; ii<todo.nnodes; ii++)
  {
    if( getfromdlist( &number, todo.nodes[ii]) == -1 )
    {
      addtolist( &explored, todo.nodes[ii] );
      addtodlist( &low, todo.nodes[ii], number.ndata);
      addtodlist( &number, todo.nodes[ii], number.ndata);
      addtolist( &stack, todo.nodes[ii]);
      while( explored.ndata != 0 )
      {
        hh = explored.data[explored.ndata-1];
        if( getgraphval( &todo, hh)->ndata != 0 )
        {
          ss = popfromlist( getgraphval( &todo, hh )  );
          if( getfromdlist( &number, ss) == -1 )
          {
            addtodlist( &low, ss, number.ndata);
            addtodlist( &number, ss, number.ndata);
            addtolist( &stack, ss );
            addtolist( &explored, ss);
            fflush(stdout);
          }
          else if( getfromdlist( &number, ss) < getfromdlist( &number, hh) && whereinlist( &stack, ss)!= -1 )
          {
            min = imin( getfromdlist( &number, ss ), getfromdlist( &low, hh));
            addtodlist( &low, hh, min);
          }
        }
        else
        {
          popfromlist( &explored );
          if( getfromdlist( &low, hh ) == getfromdlist( &number, hh ) )
          {
            while( stack.ndata != 0 && getfromdlist( &number, stack.data[stack.ndata-1] ) >= getfromdlist( &number, hh ) )
              addtograph1( components, hh, popfromlist( &stack) );
          }
          if( explored.ndata != 0 )
          {
            min = imin( getfromdlist( &low, hh), getfromdlist( &low, explored.data[explored.ndata-1]) );
            addtodlist( &low, explored.data[explored.ndata-1], min);
          }
        }
      }
    }
  }
  
  return components;
}
// ---------------------------------------------------------------------
int ELEC( char **input, int inp_index )
{
  int           ii, jj;
  int           gotit;
  char          buffer[64];
  char          filename[512], outfilename[512];    
  intlist      *sTrj;
  graph        *test1, *out;
  
  while( strncmp (buffer, "END", 3))
  {
    gotit = 0;
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--STRJ", 6))
    {
      sscanf(buffer, "--STRJ %s ", filename );
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--OUT", 6))
    {
      sscanf(buffer, "--STRJ %s ", outfilename );
      gotit = 1;
    }
    if( gotit==0 )
    {
     fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_index++;
  }
  
  sTrj = readSymbTrj( filename );
  test1 = sTrj2Graph( sTrj );
  
  out = Tarjan_ite( test1 );
  
  printf("strongly connected components (%d clusters)\n", out->nnodes);
  
  for( ii=0; ii<out->nnodes; ii++ )
  {
    for( jj=0; jj<out->leenklist[ii].ndata; jj++)
      printf( "%d ", out->leenklist[ii].data[jj] );
    printf("\n");
  }
  
  exit(0);

  
}

//======================================================================
