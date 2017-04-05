// ------------------------------------------------------------------
// Copyright (C) 2003  University of Zurich and Univ. of Modena & Reggio Emilia
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
/*! \file tools.c
 \brief mathematical and misc tools
 
 Mathematical and other basic functions used by several files and
 functions in the wordom program. When adding a function that might be
 useful to many, place it here
*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "fileio.h"
#include "tools.h"
#include "qcprot.h"

/*#define ROTATE(a,i,j,k,l) g=a[i][j];\
                          h=a[k][l];\
                          a[i][j]=g-s*(h+g*tau);\
	                  a[k][l]=h+s*(g-h*tau);
*/
extern short int verbose;

//---------------------------------------------------------------------
/*int compare_ints( const void* a, const void* b )
{
  int* arg1 = (int*) a;
  int* arg2 = (int*) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}*/
//---------------------------------------------------------------------
void DiagMatrix ( _eigen *eigen)
{
   int ii, jj;
   // SSYEV variables
   char JOBZ = 'V';
   char UPLO = 'U';

   static int No = 0;
   static int LDA = 0;
   static int INFO = 0;
   static float *W = NULL;
   static float *A = NULL;
   static float *WORK = NULL;
   static float LWORK = 0;
   
// static matrix;
   //fprintf(stderr, "got into diagmatrix (size %d)\n", eigen->size); fflush(stderr);
   
  #ifdef LAPACK
  
   if ( W == NULL )
   {
    No    = eigen->size;
    LDA   = eigen->size;
    LWORK = LDA*No;
    WORK  = calloc(LWORK,sizeof(float));
    W     = calloc(No,sizeof(float));
    A = (float *)walloc ( eigen->size*eigen->size , sizeof(float));
   }
   else if ( No != eigen->size )
   {
    free( WORK ) ;
    free( W ) ;
    free( A ) ;
    No    = eigen->size;
    LDA   = eigen->size;
    LWORK = LDA*No;
    WORK  = calloc(LWORK,sizeof(float));
    W     = calloc(No,sizeof(float));
    A = (float *)walloc ( eigen->size*eigen->size , sizeof(float));
   }
   
   for( ii=0; ii<eigen->size; ii++)
     for ( jj=0; jj<eigen->size; jj++)
       A[ii*eigen->size+jj] = eigen->inpmat[jj][ii];
   
   //fprintf(stderr, "launching ssyev\n"); fflush(stderr);
   ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);
   //fprintf(stderr, "done ssyev\n"); fflush(stderr);
   
   // copy eigenvalues and eigenvectors to output variables
   // ssyev lists both in increasing order (biggest last)
   // we want biggest first (seems more logical)
   eigen->eigval = (float *) calloc( eigen->size, sizeof(float));
   for( ii=0; ii<eigen->size; ii++)
    eigen->eigval[ii] = W[eigen->size-(ii+1)];
   
   eigen->eigvec =  (float **) calloc( eigen->size, sizeof(float *));
   eigen->eigvec[0] = (float *) calloc( eigen->size*eigen->size, sizeof(float));
   for( ii=0; ii<eigen->size; ii++)
    eigen->eigvec[ii] = eigen->eigvec[0] + ii*eigen->size;
   
   for( ii=0; ii<eigen->size; ii++)
    for( jj=0; jj<eigen->size; jj++)
     eigen->eigvec[ii][jj] = A[eigen->size*eigen->size - ((ii+1) * eigen->size) + jj];
   
  #else
   fprintf( stderr, "Sorry, all modules requiring matrix diagonalization with lapack libraries\n have been disabled in this executable\n");
   exit(0);
  #endif
   
   return;
}
// ------------------------------------------------------------------
void superimpose ( int npoints, float **referenceSet, float **tomoveSet, float rotMatrix[3][3], float transVec[3] )
{
   int ii, jj, kk;
   float centerRef[3], centerMov[3];
   float **refSet, **movSet;
   float R[3][3], RtR[3][3], eigVec[3][3];
   float a[3][3], b[3][3], norm;

    // SSYEV variables
   char JOBZ = 'V';
   char UPLO = 'U';

   int No;
   int LDA;
   int INFO;

   float *W;
   float *A;
   float *WORK;
   float LWORK;
   
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
       RtR[ii][jj]    = 0;

   
    // Copy referenceSet,toMoveSet on refSet,movSet
   refSet    = calloc (   npoints, sizeof ( float *));
   refSet[0] = calloc ( 3*npoints, sizeof ( float  ));
   for ( ii=0; ii<npoints; ii++ )
    refSet[ii] = refSet[0] + 3*ii;
   for ( ii=0; ii<npoints; ii++)
    for ( jj=0; jj<3; jj++)
     refSet[ii][jj] = 0.0;
   
   movSet    = calloc (   npoints, sizeof ( float *));
   movSet[0] = calloc ( 3*npoints, sizeof ( float  ));
   for ( ii=0; ii<npoints; ii++ )
    movSet[ii] = movSet[0] + 3*ii;
   
   for ( ii=0; ii<npoints; ii++)
    for ( jj=0; jj<3; jj++)
    {
     refSet[ii][jj] = referenceSet[ii][jj];
     movSet[ii][jj] =    tomoveSet[ii][jj];
    }
   
   
    // Translate the vectors with respect to their centroid
   for ( ii=0; ii<3; ii++)
    {
     centerRef[ii] = 0;
     centerMov[ii] = 0;
     for ( jj=0; jj<npoints; jj++)
      {
       centerRef[ii] += refSet[jj][ii];
       centerMov[ii] += movSet[jj][ii];
      }
     centerRef[ii] /= npoints;
     centerMov[ii] /= npoints;
    }
   
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<npoints; jj++)
     {
      refSet[jj][ii] -= centerRef[ii];
      movSet[jj][ii] -= centerMov[ii];
     }
   
    // Compute the matrix R
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     {
      R[ii][jj]=0.0;
      for ( kk=0; kk<npoints; kk++)
        R[ii][jj] += refSet[kk][ii]*movSet[kk][jj];
     }
   //printf("======R matrix========\n");
   //printf("%8.3f %8.3f %8.3f\n", R[0][0], R[1][0], R[2][0]);
   //printf("%8.3f %8.3f %8.3f\n", R[0][1], R[1][1], R[2][1]);
   //printf("%8.3f %8.3f %8.3f\n", R[0][2], R[1][2], R[2][2]);
   //printf("--------------\n");
   
    // Compute the matrix RtR ( (R transposed)*R )
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     for ( kk=0; kk<3; kk++)
      RtR[ii][jj] += R[kk][ii]*R[kk][jj];
    
   //printf("======inpmatrix========\n");
   //printf("%8.3f %8.3f %8.3f\n", RtR[0][0], RtR[1][0], RtR[2][0]);
   //printf("%8.3f %8.3f %8.3f\n", RtR[0][1], RtR[1][1], RtR[2][1]);
   //printf("%8.3f %8.3f %8.3f\n", RtR[0][2], RtR[1][2], RtR[2][2]);
   //printf("--------------\n");
    // Compute the eigenvalues and eigenvectors of the matrix RtR
    // ssyev_ from lapack
   No    = 3;
   LDA   = 3;
   LWORK = 9;                           // LDA*No;
   WORK  = calloc(LWORK,sizeof(float));
   W     = calloc(No,sizeof(float));
   INFO  = 0;

   A = calloc ( 3*3, sizeof (float) );
   for( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
       A[ii*3+jj] = RtR[jj][ii];

   #ifdef LAPACK
    ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);
   #else
    fprintf( stderr, "Sorry, all modules requiring matrix diagonalization with lapack libraries\n have been disabled in this executable\n");
    exit(0);
   #endif
   for ( ii=0; ii<3; ii++)
    {
     //eigVal[ii] = W[ii];
     for ( jj=0; jj<3; jj++)
      eigVec[ii][jj] = A[ii+jj*3];
    }
   
   //printf("======eigvec matrix========\n");
   //printf("%8.3f %8.3f %8.3f\n", eigVec[0][0], eigVec[1][0], eigVec[2][0]);
   //printf("%8.3f %8.3f %8.3f\n", eigVec[0][1], eigVec[1][1], eigVec[2][1]);
   //printf("%8.3f %8.3f %8.3f\n", eigVec[0][2], eigVec[1][2], eigVec[2][2]);
   //printf("--------------\n");
   
   // Determine 'a' 
   // ssyev_ gives eigenvecs in (reverse) order of eigenval, no need to sort
   // a1
   for ( ii=0; ii<3; ii++)
     a[ii][0]=eigVec[ii][2];
   // a2
   for ( ii=0; ii<3; ii++)
     a[ii][1]=eigVec[ii][1];
   
   // a3=a1^a2 to be sure to have a right-handed system
   a[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
   a[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
   a[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

   // Determine 'b'
   // b1 and b2
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<2; jj++)
       b[ii][jj]=0.0;
   for ( ii=0; ii<2; ii++)
     for ( jj=0; jj<3; jj++)
       for ( kk=0; kk<3; kk++)
         b[jj][ii] += R[jj][kk]*a[kk][ii];
   // normalize b1 and b2
   for ( ii=0; ii<2; ii++)
   {
    norm = sqrt (b[0][ii]*b[0][ii] + b[1][ii]*b[1][ii] + b[2][ii]*b[2][ii]);
    b[0][ii] /= norm;
    b[1][ii] /= norm;
    b[2][ii] /= norm;
   }
   //normalizeVector(&(b[0][ii]),&(b[1][ii]),&(b[2][ii]));
   // b3=b1^b2
   b[0][2]=b[1][0]*b[2][1]-b[2][0]*b[1][1];
   b[1][2]=b[2][0]*b[0][1]-b[0][0]*b[2][1];
   b[2][2]=b[0][0]*b[1][1]-b[1][0]*b[0][1];


    // Compute the orthogonal matrix rotationMatrix
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     {
      rotMatrix[ii][jj]=0.0;
      
      for ( kk=0; kk<3; kk++)
       rotMatrix[ii][jj] += b[ii][kk]*a[jj][kk];
     }
   
   
    // Compute the translation vector "translationVector"
   for ( ii=0; ii<3; ii++)
    {
     transVec[ii] = centerRef[ii];
     for ( jj=0; jj<3; jj++)
      transVec[ii] -= rotMatrix[ii][jj]*centerMov[jj];
    }
   
   free (refSet[0]);
   free (refSet);
   free (movSet[0]);
   free (movSet);
   free (A);
   free (W);
   free (WORK);
   
   return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
/*void superimpose_old ( int npoints, float **referenceSet, float **tomoveSet, float rotMatrix[3][3], float transVec[3] )
{
   int ii, jj, kk;
   float centerRef[3], centerMov[3];
   float **refSet, **movSet;
   float R[3][3], RtR[3][3], eigVec[3][3];
   float a[3][3], b[3][3], norm;

    // SSYEV variables
   char JOBZ = 'V';
   char UPLO = 'U';

   int No;
   int LDA;
   int INFO;

   float *W;
   float *A;
   float *WORK;
   float LWORK;
   
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
       RtR[ii][jj]    = 0;

   
    // Copy referenceSet,toMoveSet on refSet,movSet
   refSet    = calloc (   npoints, sizeof ( float *));
   refSet[0] = calloc ( 3*npoints, sizeof ( float  ));
   for ( ii=0; ii<npoints; ii++ )
    refSet[ii] = refSet[0] + 3*ii;
   for ( ii=0; ii<npoints; ii++)
    for ( jj=0; jj<3; jj++)
     refSet[ii][jj] = 0.0;
   
   movSet    = calloc (   npoints, sizeof ( float *));
   movSet[0] = calloc ( 3*npoints, sizeof ( float  ));
   for ( ii=0; ii<npoints; ii++ )
    movSet[ii] = movSet[0] + 3*ii;
   
   for ( ii=0; ii<npoints; ii++)
    for ( jj=0; jj<3; jj++)
    {
     refSet[ii][jj] = referenceSet[ii][jj];
     movSet[ii][jj] =    tomoveSet[ii][jj];
    }
   
   
    // Translate the vectors with respect to their centroid
   for ( ii=0; ii<3; ii++)
    {
     centerRef[ii] = 0;
     centerMov[ii] = 0;
     for ( jj=0; jj<npoints; jj++)
      {
       centerRef[ii] += refSet[jj][ii];
       centerMov[ii] += movSet[jj][ii];
      }
     centerRef[ii] /= npoints;
     centerMov[ii] /= npoints;
    }
   
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<npoints; jj++)
     {
      refSet[jj][ii] -= centerRef[ii];
      movSet[jj][ii] -= centerMov[ii];
     }
   
    // Compute the matrix R
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     {
      R[ii][jj]=0.0;
      for ( kk=0; kk<npoints; kk++)
        R[ii][jj] += refSet[kk][ii]*movSet[kk][jj];
     }
   
    // Compute the matrix RtR ( (R transposed)*R )
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     for ( kk=0; kk<3; kk++)
      RtR[ii][jj] += R[kk][ii]*R[kk][jj];
    
    // Compute the eigenvalues and eigenvectors of the matrix RtR
    // ssyev_ from lapack
   No    = 3;
   LDA   = 3;
   LWORK = 9;                           // LDA*No;
   WORK  = calloc(LWORK,sizeof(float));
   W     = calloc(No,sizeof(float));
   INFO  = 0;

   A = calloc ( 3*3, sizeof (float) );
   for( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
       A[ii*3+jj] = RtR[jj][ii];

   #ifdef LAPACK
    ssyev_(&JOBZ, &UPLO, &No, A, &LDA, W, WORK, &LWORK, &INFO);
   #else
    fprintf( stderr, "Sorry, all modules requiring matrix diagonalization with lapack libraries\n have been disabled in this executable\n");
    exit(0);
   #endif
   for ( ii=0; ii<3; ii++)
    {
     //eigVal[ii] = W[ii];
     for ( jj=0; jj<3; jj++)
      eigVec[ii][jj] = A[ii+jj*3];
    }
   
   // Determine 'a' 
   // ssyev_ gives eigenvecs in (reverse) order of eigenval, no need to sort
   // a1
   for ( ii=0; ii<3; ii++)
     a[ii][0]=eigVec[ii][2];
   // a2
   for ( ii=0; ii<3; ii++)
     a[ii][1]=eigVec[ii][1];
   
   // a3=a1^a2 to be sure to have a right-handed system
   a[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
   a[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
   a[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

   // Determine 'b'
   // b1 and b2
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<2; jj++)
       b[ii][jj]=0.0;
   for ( ii=0; ii<2; ii++)
     for ( jj=0; jj<3; jj++)
       for ( kk=0; kk<3; kk++)
         b[jj][ii] += R[jj][kk]*a[kk][ii];
   // normalize b1 and b2
   for ( ii=0; ii<2; ii++)
   {
    norm = sqrt (b[0][ii]*b[0][ii] + b[1][ii]*b[1][ii] + b[2][ii]*b[2][ii]);
    b[0][ii] /= norm;
    b[1][ii] /= norm;
    b[2][ii] /= norm;
   }
   //normalizeVector(&(b[0][ii]),&(b[1][ii]),&(b[2][ii]));
   // b3=b1^b2
   b[0][2]=b[1][0]*b[2][1]-b[2][0]*b[1][1];
   b[1][2]=b[2][0]*b[0][1]-b[0][0]*b[2][1];
   b[2][2]=b[0][0]*b[1][1]-b[1][0]*b[0][1];


    // Compute the orthogonal matrix rotationMatrix
   for ( ii=0; ii<3; ii++)
    for ( jj=0; jj<3; jj++)
     {
      rotMatrix[ii][jj]=0.0;
      
      for ( kk=0; kk<3; kk++)
       rotMatrix[ii][jj] += b[ii][kk]*a[jj][kk];
     }
   
   
    // Compute the translation vector "translationVector"
   for ( ii=0; ii<3; ii++)
    {
     transVec[ii] = centerRef[ii];
     for ( jj=0; jj<3; jj++)
      transVec[ii] -= rotMatrix[ii][jj]*centerMov[jj];
    }
   
   free (refSet[0]);
   free (refSet);
   free (movSet[0]);
   free (movSet);
   free (A);
   free (W);
   free (WORK);
   
   return;
}*/
// ------------------------------------------------------------------
void CalcRMSDRotTrans(int npoints, float **referenceSet, float **tomoveSet, float rotMatrix[3][3], float transVec[3])
{
   // same as superimpose, attempting to use coor[3][nato] format
   int      ii, jj, kk;
   float    centerRef[3], centerMov[3];
   float  **refSet, **movSet;
   _eigen   matrix;
   float    R[3][3];
   float    a[3][3], b[3][3], norm;
   
    // Copy referenceSet,toMoveSet on refSet,movSet
   refSet = wrd_fmatAlloc( 3, npoints );
   movSet = wrd_fmatAlloc( 3, npoints );

   for ( ii=0; ii<3*npoints; ii++)
   {
     refSet[0][ii] = referenceSet[0][ii];
     movSet[0][ii] =    tomoveSet[0][ii];
   }
   
    // Translate the vectors with respect to their centroid
   for ( ii=0; ii<3; ii++)
   {
     centerRef[ii] = 0;
     centerMov[ii] = 0;
     for ( jj=0; jj<npoints; jj++)
     {
      centerRef[ii] += refSet[ii][jj];
      centerMov[ii] += movSet[ii][jj];
     }
     centerRef[ii] /= npoints;
     centerMov[ii] /= npoints;
   }
   
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<npoints; jj++)
     {
       refSet[ii][jj] -= centerRef[ii];
       movSet[ii][jj] -= centerMov[ii];
     }
   
    // Compute matrix R
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
     {
       R[ii][jj]=0.0;
       for ( kk=0; kk<npoints; kk++)
         R[ii][jj] += refSet[ii][kk]*movSet[jj][kk];
     }
   
   // WARNING DEBUG WARNING - free matrix.inpmat + matrix.eigvec + matrix.eigval !!!
   matrix.inpmat = wrd_fmatAlloc( 3, 3 );
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
     {
       matrix.inpmat[ii][jj]  = 0.0;
       for ( kk=0; kk<3; kk++)
         matrix.inpmat[ii][jj] += R[kk][ii]*R[kk][jj];
     }
   
   matrix.size = 3;
   DiagMatrix( &matrix );
   //fprintf( stderr, "DEBUG - got here2\n"); fflush(stderr);
   
   //printf("======eigvec matrix========\n");
   //printf("%8.3f %8.3f %8.3f\n", matrix.eigvec[0][0], matrix.eigvec[1][0], matrix.eigvec[2][0]);
   //printf("%8.3f %8.3f %8.3f\n", matrix.eigvec[0][1], matrix.eigvec[1][1], matrix.eigvec[2][1]);
   //printf("%8.3f %8.3f %8.3f\n", matrix.eigvec[0][2], matrix.eigvec[1][2], matrix.eigvec[2][2]);
   //printf("--------------\n");
   
   // Determine 'a' 
   // a1
   for ( ii=0; ii<3; ii++)
     //a[ii][0]=matrix.eigvec[ii][2];
     a[ii][0]=matrix.eigvec[2][ii];
   // a2
   for ( ii=0; ii<3; ii++)
     //a[ii][1]=matrix.eigvec[ii][1];
     a[ii][1]=matrix.eigvec[1][ii];
   
   // a3=a1^a2 to be sure to have a right-handed system
   a[0][2]=a[1][0]*a[2][1]-a[2][0]*a[1][1];
   a[1][2]=a[2][0]*a[0][1]-a[0][0]*a[2][1];
   a[2][2]=a[0][0]*a[1][1]-a[1][0]*a[0][1];

   // Determine 'b'
   // b1 and b2
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<2; jj++)
       b[ii][jj]=0.0;
   for ( ii=0; ii<2; ii++)
     for ( jj=0; jj<3; jj++)
       for ( kk=0; kk<3; kk++)
         b[jj][ii] += R[jj][kk]*a[kk][ii];
   // normalize b1 and b2
   for ( ii=0; ii<2; ii++)
   {
     norm = sqrt (b[0][ii]*b[0][ii] + b[1][ii]*b[1][ii] + b[2][ii]*b[2][ii]);
     b[0][ii] /= norm;
     b[1][ii] /= norm;
     b[2][ii] /= norm;
   }
   //normalizeVector(&(b[0][ii]),&(b[1][ii]),&(b[2][ii]));
   // b3=b1^b2
   b[0][2]=b[1][0]*b[2][1]-b[2][0]*b[1][1];
   b[1][2]=b[2][0]*b[0][1]-b[0][0]*b[2][1];
   b[2][2]=b[0][0]*b[1][1]-b[1][0]*b[0][1];


    // Compute the orthogonal matrix rotationMatrix
   for ( ii=0; ii<3; ii++)
     for ( jj=0; jj<3; jj++)
     {
       rotMatrix[ii][jj]=0.0;
       
       for ( kk=0; kk<3; kk++)
         rotMatrix[ii][jj] += b[ii][kk]*a[jj][kk];
     }
   
   
    // Compute the translation vector
   for ( ii=0; ii<3; ii++)
   {
     transVec[ii] = centerRef[ii];
     for ( jj=0; jj<3; jj++)
       transVec[ii] -= rotMatrix[ii][jj]*centerMov[jj];
   }
   
   return;
}
// ------------------------------------------------------------------
int applyRotTrans( float **moving, float rotmatrix[3][3], float transvec[3], int nato )
{
  int     ii;
  float   aa, bb, cc;
  
  for ( ii=0; ii<nato; ii++)
   {
    
    aa = moving[0][ii]*rotmatrix[0][0] + moving[1][ii]*rotmatrix[0][1] + moving[2][ii]*rotmatrix[0][2] ;
    bb = moving[0][ii]*rotmatrix[1][0] + moving[1][ii]*rotmatrix[1][1] + moving[2][ii]*rotmatrix[1][2] ;
    cc = moving[0][ii]*rotmatrix[2][0] + moving[1][ii]*rotmatrix[2][1] + moving[2][ii]*rotmatrix[2][2] ;
    moving[0][ii] = aa + transvec[0]; //+ xcenter
    moving[1][ii] = bb + transvec[1]; //+ ycenter
    moving[2][ii] = cc + transvec[2]; //+ zcenter
    
   }
   
   return 0;
}

// ------------------------------------------------------------------
float RmsdCalc(float **refcoor, float **movcoor, int nato, int super)
{

  float rmsd=-1;
  float rotmat[9];
  
  if( super == 0 )
    return RmsdCalc_nosup( refcoor, movcoor, nato );
  
  rmsd = RmsdCalcQCP( refcoor, movcoor, nato, super);
  if( rmsd == -1 )
    return RmsdCalcKabsch3n( refcoor, movcoor, nato, super );
  
  //rmsd = CalcRMSDRotationalMatrix(refcoor, movcoor, nato, rotmat, NULL);
  
  return rmsd;
}
// ------------------------------------------------------------------
float RmsdCalcKabsch(float **refcoor, float **movcoor, int nato, int super )
{
   int          ii;
   float        rmsd;
   float        rotmatrix[3][3], transvec[3];
   float        aa, bb, cc;
   
   // if rmsd optimization is called, compute superposition and rototraslate movcoor
   if(super)
   {
     superimpose( nato, refcoor, movcoor, rotmatrix, transvec );
     
     for ( ii=0; ii<nato; ii++)
     {
       aa = movcoor[ii][0]*rotmatrix[0][0] + movcoor[ii][1]*rotmatrix[0][1] + movcoor[ii][2]*rotmatrix[0][2] ;
       bb = movcoor[ii][0]*rotmatrix[1][0] + movcoor[ii][1]*rotmatrix[1][1] + movcoor[ii][2]*rotmatrix[1][2] ;
       cc = movcoor[ii][0]*rotmatrix[2][0] + movcoor[ii][1]*rotmatrix[2][1] + movcoor[ii][2]*rotmatrix[2][2] ;
       movcoor[ii][0] = aa + transvec[0]; //+ xcenter
       movcoor[ii][1] = bb + transvec[1]; //+ ycenter
       movcoor[ii][2] = cc + transvec[2]; //+ zcenter
     }
   }
   
   // compute rmsd and return value
   rmsd=0;
   for ( ii=0; ii<nato; ii++ )
    rmsd += ( (refcoor[ii][0]-movcoor[ii][0])*(refcoor[ii][0]-movcoor[ii][0]) + 
              (refcoor[ii][1]-movcoor[ii][1])*(refcoor[ii][1]-movcoor[ii][1]) + 
              (refcoor[ii][2]-movcoor[ii][2])*(refcoor[ii][2]-movcoor[ii][2]) );
   
   rmsd /= nato;
   rmsd = sqrt ( rmsd );
   return rmsd;
}
// ------------------------------------------------------------------
float RmsdCalcKabsch3n(float **refcoor, float **movcoor, int nato, int super )
{
   float        rmsd;
   float        rotmatrix[3][3], transvec[3];
   
  if ( super )
  {
    CalcRMSDRotTrans ( nato, refcoor, movcoor, rotmatrix, transvec);
    applyRotTrans( movcoor, rotmatrix, transvec, nato);
  }
  
  rmsd = RmsdCalc_nosup( refcoor, movcoor, nato );
  //fprintf(stderr, "RMSD: %f\n", rmsd);
  return rmsd;
}

// ------------------------------------------------------------------
float RmsdCalcQCP(float **refcoor, float **movcoor, int nato, int super)
{
  float   rotmat[9];
  float   rmsd=0.0;
  int     ii;
       
  for( ii=0; ii<9; ii++)
    rotmat[ii] = 0.0;
  if( super )
    rmsd = CalcRMSDRotationalMatrix(refcoor, movcoor, nato, rotmat, NULL);
  else
    rmsd =  RmsdCalc_nosup( refcoor, movcoor, nato );
  
  return rmsd;
}
// ------------------------------------------------------------------
float RmsdCalc_preserve(float **refcoor, float **movcoor, int nato, int super )
{
   float **hereref, **heremov;
   float   rmsd;
   int     ii;
   
   hereref = calloc( 3, sizeof(float *));
   hereref[0] = calloc( 3*nato, sizeof(float  ));
   for ( ii =0; ii<3; ii++)
     hereref[ii] = hereref[0] + nato*ii;
   
   heremov = calloc( 3, sizeof(float *));
   heremov[0] = calloc( 3*nato, sizeof(float  ));
   for ( ii =0; ii<3; ii++)
     heremov[ii] = heremov[0] + nato*ii;
   
   for( ii=0; ii<3*nato; ii++ )
   {
     hereref[ii] = refcoor[ii];
     heremov[ii] = movcoor[ii];
   }
   
   rmsd = RmsdCalc( hereref, heremov, nato, super );
   
   free(hereref[0]);
   free(hereref);
   free(heremov[0]);
   free(heremov);
   
   return rmsd;
}
// ------------------------------------------------------------------
float RmsdCalc_nosup(float **refcoor, float **movcoor, int nato)
{
  // compute rmsd and return value
  int          ii;
  float        rmsd;
  
  rmsd=0;
  for ( ii=0; ii<nato; ii++ )
   rmsd += ( (refcoor[0][ii]-movcoor[0][ii])*(refcoor[0][ii]-movcoor[0][ii]) + 
             (refcoor[1][ii]-movcoor[1][ii])*(refcoor[1][ii]-movcoor[1][ii]) + 
             (refcoor[2][ii]-movcoor[2][ii])*(refcoor[2][ii]-movcoor[2][ii]) );
  
  rmsd /= nato;
  rmsd = sqrt ( rmsd );
  
  return rmsd;
}
// ------------------------------------------------------------------
// python-serving function
float MolRmsd( Molecule *mol1, Molecule *mol2, char *selestring )
{
  int         ii;
  int         nato1, nato2;
  Selection   sele1, sele2;
  float     **refcoor, **movcoor;
  float       rmsd;
  
  nato1 = GetSele ( selestring, &sele1, mol1 );
  nato2 = GetSele ( selestring, &sele2, mol2 );
  if( nato1 != nato2 )
    return -1;
  for( ii=0; ii<nato1; ii++ )
    if( sele1.selatm[ii] != sele2.selatm[ii] )
      return -1;
  
  refcoor = wrd_fmatAlloc( 3, nato1 );
  movcoor = wrd_fmatAlloc( 3, nato2 );
  GetSeleCoor2Vecs( &mol1->coor, refcoor[0], refcoor[1], refcoor[2], &sele1 );
  GetSeleCoor2Vecs( &mol2->coor, movcoor[0], movcoor[1], movcoor[2], &sele2 );
  
  rmsd = RmsdCalc(refcoor, movcoor, nato1, 1);
  wrd_fmatFree( refcoor );
  wrd_fmatFree( movcoor );
  
  return rmsd;
}
// ------------------------------------------------------------------
// Compute the distance matrix
void DistMtxCalc( int pbcflag , float *pbcbox, int nato , float **incoor, float * dist_mtx )
{
  int                   ii, jj, kk;
  static float        **coor;
  static int            coor_size;
  if(!coor)
  {
   coor = calloc (6, sizeof(float *));
   for (ii=0;ii<6;ii++)
    coor[ii] = calloc ( (int)(nato*(nato-1)/2) , sizeof(float));
   coor_size = nato ;
   for( ii=0; ii<6; ii++)
    for ( jj=0; jj<nato; jj++)
     coor[ii][jj] = 0;
  }
  else if( coor_size < nato )
  {
   for (ii=0;ii<6;ii++)
    free( coor[ii] );
   
   for (ii=0; ii<6; ii++)
    coor[ii] = calloc ( (int)(nato*(nato-1)/2) , sizeof(float));
   
   coor_size = nato ;
   for( ii=0; ii<6; ii++)
    for ( jj=0; jj<nato; jj++)
     coor[ii][jj] = 0;
  }
  
  kk=0;
  for( ii=0; ii<nato; ii++ )
  {
   for ( jj=0; jj<ii; jj++ )
   {
    coor[0][kk] = incoor[0][ii];   // Assign the coordinates of the relevant pairs of
    coor[1][kk] = incoor[1][ii];   //  atoms to the vectors for distance calculation
    coor[2][kk] = incoor[2][ii];
    coor[3][kk] = incoor[0][jj];
    coor[4][kk] = incoor[1][jj];
    coor[5][kk] = incoor[2][jj];
    kk++;
   }
  }
  
  // Compute the distance matrix as a vector having k = n * (n-1) / 2 components
  if( pbcflag )
    vDistancePBC( kk, dist_mtx, coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], pbcbox); 
  else
    vDistance( kk, dist_mtx, coor[0], coor[1], coor[2], coor[3], coor[4], coor[5]); 
  return;
}
//-------------------------------------------------------------------
// Compute DRMS
void Drms( int k , float  * distance , float * dist_mtx , float * drms ) 
{
    float *  tmp_mtx;
    int i=1;
    int j=1;
    float diff=-1.0;

    tmp_mtx = calloc( k , sizeof(float) );
    
    scopy( &k , dist_mtx  , &j , tmp_mtx , &i );
    saxpy( &k , &diff , distance , &i , tmp_mtx  , &j ) ;
    (*drms) = (snrm2( &k ,  tmp_mtx  , &j ) / sqrt((float)k));
    
    free(tmp_mtx);
    return;
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------
float DrmsCalc(float **refcoor, float **movcoor, int nato, int pbcflag, float *pbcbox )
{
   int                  msize;      // counter + distance matrix size
   static float        *dist_mtx1;              // distance matrix
   static float        *dist_mtx2;              // distance matrix
   static int           dist_mtx_size;          // distance matrix size
   float                drms;
   
   msize = (int)(nato*(nato-1)/2);
   if ( ! dist_mtx1 && ! dist_mtx2 )
   {
    dist_mtx1 = calloc ( msize , sizeof (float));
    dist_mtx2 = calloc ( msize , sizeof (float));
    dist_mtx_size = msize ;
   }
   else if ( dist_mtx_size < msize )
   {
    free( dist_mtx1 ) ;
    free( dist_mtx2 ) ;
    dist_mtx1 = calloc ( msize , sizeof (float));
    dist_mtx2 = calloc ( msize , sizeof (float));
    dist_mtx_size = msize ;
   }
   
   DistMtxCalc ( pbcflag, pbcbox, nato , refcoor , dist_mtx1 ) ;
   DistMtxCalc ( pbcflag, pbcbox, nato , movcoor , dist_mtx2 ) ;
   
   drms=0.0;
   Drms ( msize , dist_mtx1 , dist_mtx2 , &drms ) ;
   
   return(drms);
}
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// Dimension distvcts structure
void FillDistvcts ( _distvcts *vcts )
{
       vcts->dist_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->dist= (float *) (((unsigned long int) vcts->dist_t + 8 ) & (~0xf));
       vcts->disty_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->disty= (float *) (((unsigned long int) vcts->disty_t + 8 ) & (~0xf));
       vcts->distz_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->distz= (float *) (((unsigned long int) vcts->distz_t + 8 ) & (~0xf));

       vcts->x0_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->x0= (float *) (((unsigned long int) vcts->x0_t + 8 ) & (~0xf));
       vcts->y0_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->y0= (float *) (((unsigned long int) vcts->y0_t + 8 ) & (~0xf));
       vcts->z0_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->z0= (float *) (((unsigned long int) vcts->z0_t + 8 ) & (~0xf));

       vcts->x1_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->x1= (float *) (((unsigned long int) vcts->x1_t + 8 ) & (~0xf));
       vcts->y1_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->y1= (float *) (((unsigned long int) vcts->y1_t + 8 ) & (~0xf));
       vcts->z1_t = ( float *) calloc ( vcts->ndata+1, sizeof(float) );
       vcts->z1= (float *) (((unsigned long int) vcts->z1_t + 8 ) & (~0xf));
}
// ------------------------------------------------------------------
// Compute Distance
/*void vdistance ( _distvcts *vcts)
{
   int ii;
   
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->dist[ii]   = vcts->x0[ii];
      vcts->dist[ii]  -= vcts->x1[ii];
      vcts->dist[ii]  *= vcts->dist[ii];
      vcts->disty[ii]  = vcts->y0[ii];
     }
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->disty[ii] -= vcts->y1[ii];
      vcts->disty[ii] *= vcts->disty[ii];
      vcts->distz[ii]  = vcts->z0[ii];
      vcts->distz[ii] -= vcts->z1[ii];
      vcts->distz[ii] *= vcts->distz[ii];
      vcts->dist[ii]  += vcts->disty[ii];
      vcts->dist[ii]  += vcts->distz[ii];
      vcts->dist[ii]   = sqrt(vcts->dist[ii]);
    }

   return;
}*/
// ------------------------------------------------------------------
// Compute Distance with Periodic Boundary Conditions
void vdistancePBC ( float *PBCBOX, _distvcts *vcts)
{
   int ii;
   
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->dist[ii]   = vcts->x0[ii];
      vcts->dist[ii]  -= vcts->x1[ii];
     }
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->disty[ii]  = vcts->y0[ii];
      vcts->disty[ii] -= vcts->y1[ii];
     }
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->distz[ii]  = vcts->z0[ii];
      vcts->distz[ii] -= vcts->z1[ii];
     }
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->dist[ii]   -= PBCBOX[0] * rintf(vcts->dist[ii]/PBCBOX[0]);
      vcts->disty[ii]  -= PBCBOX[1] * rintf(vcts->disty[ii]/PBCBOX[1]);
      vcts->distz[ii]  -= PBCBOX[2] * rintf(vcts->distz[ii]/PBCBOX[2]);
     }
   for ( ii=0; ii<vcts->ndata; ii++)
    {
      vcts->disty[ii] *= vcts->disty[ii];
      vcts->distz[ii] *= vcts->distz[ii];
      vcts->dist[ii]  *= vcts->dist[ii];
      vcts->dist[ii]  += vcts->disty[ii];
      vcts->dist[ii]  += vcts->distz[ii];
      vcts->dist[ii]   = sqrt(vcts->dist[ii]);
    }

   return;
}
// ------------------------------------------------------------------
// Compute Distance
void vDistance ( int ndata, float *dist, float *x0, float *y0, float *z0, float *x1, float *y1, float *z1 )
{
   int ii;
   static float *disty, *distz;
   static int    coor_size;
   
   if ( !disty )
    {
     disty = (float *)walloc ( ndata, sizeof(float));
     distz = (float *)walloc ( ndata, sizeof(float));
     coor_size = ndata;
    }
   else if ( coor_size < ndata )
    {
     free ( disty );
     free ( distz );
     disty = (float *)walloc ( ndata, sizeof(float));
     distz = (float *)walloc ( ndata, sizeof(float));
     coor_size = ndata;
    }
   for ( ii=0; ii<ndata; ii++)
    {
      dist[ii]   = x0[ii];
      dist[ii]  -= x1[ii];
      dist[ii]  *= dist[ii];
      disty[ii]  = y0[ii];
     }
   for ( ii=0; ii<ndata; ii++)
    {
      disty[ii] -= y1[ii];
      disty[ii] *= disty[ii];
      distz[ii]  = z0[ii];
      distz[ii] -= z1[ii];
      distz[ii] *= distz[ii];
      dist[ii]  += disty[ii];
      dist[ii]  += distz[ii];
      dist[ii]   = sqrt(dist[ii]);
    }

   return;
}
// ------------------------------------------------------------------
// Compute Distance with Periodic Boundary Conditions
void vDistancePBC ( int ndata, float *dist, float *x0, float *y0, float *z0, float *x1, float *y1, float *z1 , float *PBCBOX)
{
   int ii;
   static float *disty, *distz;
   static int    coor_size;
   
   if ( !disty )
    {
     disty = (float *)walloc ( ndata, sizeof(float));
     distz = (float *)walloc ( ndata, sizeof(float));
     coor_size = ndata;
    }
   else if ( coor_size < ndata )
    {
     free ( disty );
     free ( distz );
     disty = (float *)walloc ( ndata, sizeof(float));
     distz = (float *)walloc ( ndata, sizeof(float));
     coor_size = ndata;
    }
   for ( ii=0; ii<ndata; ii++)
    {
      dist[ii]   = x0[ii];
      dist[ii]  -= x1[ii];
    }
   for ( ii=0; ii<ndata; ii++)
    {
      disty[ii]  = y0[ii];
      disty[ii] -= y1[ii];
    }
   for ( ii=0; ii<ndata; ii++)
    {
      distz[ii]  = z0[ii];
      distz[ii] -= z1[ii];
    }
   for ( ii=0; ii<ndata; ii++)
    {
      dist[ii]   -= PBCBOX[0] * rintf(dist[ii]/PBCBOX[0]);
      disty[ii]  -= PBCBOX[1] * rintf(disty[ii]/PBCBOX[1]);
      distz[ii]  -= PBCBOX[2] * rintf(distz[ii]/PBCBOX[2]);
     }
   for ( ii=0; ii<ndata; ii++)
    {
      dist[ii]  *= dist[ii];
      disty[ii] *= disty[ii];
      distz[ii] *= distz[ii];
      dist[ii]  += disty[ii];
      dist[ii]  += distz[ii];
      dist[ii]   = sqrt(dist[ii]);
    }

   return;
}
//-----------------------------------------------------------------------  
float PBC1( CoorSet *coorset1, int atom1, CoorSet *coorset2, int atom2)
{
   /* in case a "Orthogonal" (CUBIC) pbc is encountered (all angles == 90) */
   float        value, temp;
   
   value = 0;
   temp = coorset1->xcoor[atom1-1] - coorset2->xcoor[atom2-1];
   temp -= coorset1->pbc->a_size * rintf( temp/coorset1->pbc->a_size );
   value += temp*temp;
   temp = coorset1->ycoor[atom1-1] - coorset2->ycoor[atom2-1];
   temp -= coorset1->pbc->b_size * rintf( temp/coorset1->pbc->b_size );
   value += temp*temp;
   temp = coorset1->zcoor[atom1-1] - coorset2->zcoor[atom2-1];
   temp -= coorset1->pbc->c_size * rintf( temp/coorset1->pbc->c_size );
   value += temp*temp;
   
   value = sqrt(value);
   
   return value;
}
//-----------------------------------------------------------------------  
//float PBC2( CoorSet *coorset1, int atom1, CoorSet *coorset2, int atom2)
//{
   /* in case a Truncated octahedron (TO) pbc is encountered () */
/*   int          ii, jj;
   float        value, temp;
   return 0.0;
}*/
//-----------------------------------------------------------------------  
//float PBC3( CoorSet *coorset1, int atom1, CoorSet *coorset2, int atom2)
//{
   /* in case a Rhombic dodecahedron (RHDO) pbc is encountered () */
//   return 0.0;
//}
//-----------------------------------------------------------------------
void DistanceAxis( float *distances, float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc )
{
   if( pbc == NULL )
   {
     distances[0] = xcoor1-xcoor2;
     distances[1] = ycoor1-ycoor2;
     distances[2] = zcoor1-zcoor2;
     return;
   }
   else if( pbc != NULL )
   {
     if( pbc->angle1 == 0.0 && pbc->angle2 == 0.0 && pbc->angle3 == 0.0 )
     {
       distances[0] = (xcoor1 - xcoor2) - (pbc->a_size * rintf((xcoor1 - xcoor2)/pbc->a_size));
       distances[1] = (ycoor1 - ycoor2) - (pbc->b_size * rintf((ycoor1 - ycoor2)/pbc->b_size));
       distances[2] = (zcoor1 - zcoor2) - (pbc->c_size * rintf((zcoor1 - zcoor2)/pbc->c_size));
       return;
     }
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
   }
   exit(1);
}
//-----------------------------------------------------------------------
void DistanceAxes( float *distances, float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc )
{
   if( pbc == NULL )
   {
     distances[0] = xcoor1-xcoor2;
     distances[1] = ycoor1-ycoor2;
     distances[2] = zcoor1-zcoor2;
     return;
   }
   else if( pbc != NULL )
   {
     if( pbc->angle1 == 0.0 && pbc->angle2 == 0.0 && pbc->angle3 == 0.0 )
     {
       distances[0] = (xcoor1 - xcoor2) - (pbc->a_size * rintf((xcoor1 - xcoor2)/pbc->a_size));
       distances[1] = (ycoor1 - ycoor2) - (pbc->b_size * rintf((ycoor1 - ycoor2)/pbc->b_size));
       distances[2] = (zcoor1 - zcoor2) - (pbc->c_size * rintf((zcoor1 - zcoor2)/pbc->c_size));
       return;
     }
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
   }
   exit(1);
}
//-----------------------------------------------------------------------
void SumAxis( float *sums, float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc )
{
   if( pbc == NULL )
   {
     sums[0] = xcoor1+xcoor2;
     sums[1] = ycoor1+ycoor2;
     sums[2] = zcoor1+zcoor2;
     return;
   }
   else if( pbc != NULL )
   {
     if( pbc->angle1 == 0.0 && pbc->angle2 == 0.0 && pbc->angle3 == 0.0 )
     {
       sums[0] = (xcoor1 + xcoor2) - (pbc->a_size * rintf((xcoor1 + xcoor2)/pbc->a_size));
       sums[1] = (ycoor1 + ycoor2) - (pbc->b_size * rintf((ycoor1 + ycoor2)/pbc->b_size));
       sums[2] = (zcoor1 + zcoor2) - (pbc->c_size * rintf((zcoor1 + zcoor2)/pbc->c_size));
       return;
     }
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
   }
   exit(1);
}
//-----------------------------------------------------------------------
float DistanceCoor( float xcoor1, float ycoor1, float zcoor1, float xcoor2, float ycoor2, float zcoor2, Pbc *pbc )
{
   float        distance, temp;
   
   if( pbc == NULL )
   {
     return sqrt( (xcoor1-xcoor2) * (xcoor1-xcoor2) + \
                  (ycoor1-ycoor2) * (ycoor1-ycoor2) + \
                  (zcoor1-zcoor2) * (zcoor1-zcoor2));
  }
   else if( pbc != NULL )
   {
    if( pbc->angle1 == 0.0 && pbc->angle2 == 0.0 && pbc->angle3 == 0.0 )
    {
      distance = 0;
      temp = xcoor1 - xcoor2;
      temp -= pbc->a_size * rintf(temp/pbc->a_size);
      distance += temp*temp;
      temp = ycoor1 - ycoor2;
      temp -= pbc->b_size * rintf(temp/pbc->b_size);
      distance += temp*temp;
      temp = zcoor1 - zcoor2;
      temp -= pbc->c_size * rintf(temp/pbc->c_size);
      distance += temp*temp;
      distance = sqrt(distance);
      //fprintf(stdout, "%8.3f %8.3f %8.3f => %10.5f\n", pbc->a_size, pbc->b_size, pbc->c_size, distance); fflush(stdout);
      return distance;//DistPbc1( coorset1, piSele1[ii], coorset2, piSele2[ii]);
    }
    else
    {
      fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
      exit(0);
    }
   }

   return -1.0;
}
//-----------------------------------------------------------------------
float *DistanceVectCoor( CoorSet *coorset1, int *piSele1, CoorSet *coorset2, int *piSele2, int iNumOfAtomPairs )
{
   int          ii;
   int          nato;
   float       *distances;
   
   nato = iNumOfAtomPairs;
   distances = calloc( nato, sizeof(float));

   if( (coorset1->pbc_flag == 0 && coorset2->pbc_flag == 0) || (coorset1->pbc_flag == -1 || coorset2->pbc_flag == -1))
   {
     for( ii=0; ii<nato; ii++)
     {
       distances[ii] += (coorset1->xcoor[piSele1[ii]-1]-coorset2->xcoor[piSele2[ii]-1]) *
                        (coorset1->xcoor[piSele1[ii]-1]-coorset2->xcoor[piSele2[ii]-1]);
       distances[ii] += (coorset1->ycoor[piSele1[ii]-1]-coorset2->ycoor[piSele2[ii]-1]) *
                        (coorset1->ycoor[piSele1[ii]-1]-coorset2->ycoor[piSele2[ii]-1]);
       distances[ii] += (coorset1->zcoor[piSele1[ii]-1]-coorset2->zcoor[piSele2[ii]-1]) *
                        (coorset1->zcoor[piSele1[ii]-1]-coorset2->zcoor[piSele2[ii]-1]);
       distances[ii] = sqrt(distances[ii]);
     }
     
     return distances;
   }
   
   else if( coorset1->pbc_flag == 1 && coorset2->pbc_flag == 0 )
   {
     if( coorset1->pbc->angle1 == 0.0 && coorset1->pbc->angle2 == 0.0 && coorset1->pbc->angle3 == 0.0 )
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset1, piSele1[ii], coorset2, piSele2[ii]);
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }
   else if( coorset1->pbc_flag == 0 && coorset2->pbc_flag == 1 )
   {
     if( coorset2->pbc->angle1 == 0.0 && coorset2->pbc->angle2 == 0.0 && coorset2->pbc->angle3 == 0.0 )
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset2, piSele2[ii], coorset1, piSele1[ii]);
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }
   else if( coorset1->pbc_flag == 1 && coorset2->pbc_flag == 1 )
   {
     if( coorset1->pbc->angle1 == 0.0 && coorset1->pbc->angle2 == 0.0 && coorset1->pbc->angle3 == 0.0 )
     {
       if(coorset1->pbc->a_size!=coorset2->pbc->a_size || coorset1->pbc->b_size!=coorset2->pbc->b_size || coorset1->pbc->c_size!=coorset2->pbc->c_size ||
          coorset1->pbc->angle1!=coorset2->pbc->angle1 || coorset1->pbc->angle2!=coorset2->pbc->angle2 || coorset1->pbc->angle3!=coorset2->pbc->angle3 )
         {
           fprintf( stderr, "Error in DistanceCoor! pbc data not equal!\n");
           exit(0);
         }
       
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset1, piSele1[ii], coorset2, piSele2[ii]);
     }
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }

   return 0;
}
/* ------------------------------------------------------------------ */
float * DistanceSelCoor( CoorSet *coorset1, Selection *sele1, CoorSet *coorset2, Selection *sele2 )
{
   int          ii;
   int          nato;
   float       *distances;

   if( sele1->nselatm != sele2->nselatm )
   {
     fprintf( stderr, "Error in DistanceCoor! nselatm not equal!\n");
     exit(0);
   }
   //printf("DEBUG just a check...%d ; %d\n", coorset1->pbc_flag, coorset2->pbc_flag); fflush(stdout);
   nato = sele1->nselatm;
   distances = calloc( nato, sizeof(float));
   //printf("DEBUG allocating %d: %p\n", nato, distances);
   for( ii=0; ii<nato; ii++ )
     distances[ii] = 0;
   
   if( (coorset1->pbc_flag == 0 && coorset2->pbc_flag == 0) || (coorset1->pbc_flag == -1 || coorset2->pbc_flag == -1))
   {
     for( ii=0; ii<nato; ii++)
     {
       distances[ii] += (coorset1->xcoor[sele1->selatm[ii]-1]-coorset2->xcoor[sele2->selatm[ii]-1]) *
                        (coorset1->xcoor[sele1->selatm[ii]-1]-coorset2->xcoor[sele2->selatm[ii]-1]);
       distances[ii] += (coorset1->ycoor[sele1->selatm[ii]-1]-coorset2->ycoor[sele2->selatm[ii]-1]) *
                        (coorset1->ycoor[sele1->selatm[ii]-1]-coorset2->ycoor[sele2->selatm[ii]-1]);
       distances[ii] += (coorset1->zcoor[sele1->selatm[ii]-1]-coorset2->zcoor[sele2->selatm[ii]-1]) *
                        (coorset1->zcoor[sele1->selatm[ii]-1]-coorset2->zcoor[sele2->selatm[ii]-1]);
       distances[ii] = sqrt(distances[ii]);
     }
     return distances;
   }
   else if( coorset1->pbc_flag == 1 && coorset2->pbc_flag == 0 )
   {
     if( coorset1->pbc->angle1 == 0 && coorset1->pbc->angle2 == 0 && coorset1->pbc->angle3 == 0.0 )
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset1, sele1->selatm[ii], coorset2, sele2->selatm[ii]);
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }
   else if( coorset1->pbc_flag == 0 && coorset2->pbc_flag == 1 )
   {
     if( coorset2->pbc->angle1 == 0 && coorset2->pbc->angle2 == 0 && coorset2->pbc->angle3 == 0.0 )
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset2, sele2->selatm[ii], coorset1, sele1->selatm[ii]);
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }
   else if( coorset1->pbc_flag == 1 && coorset2->pbc_flag == 1 )
   {
     if( coorset1->pbc->angle1 == 0.0 && coorset1->pbc->angle2 == 0.0 && coorset1->pbc->angle3 == 0.0 )
     {
       if(coorset1->pbc->a_size!=coorset2->pbc->a_size || coorset1->pbc->b_size!=coorset2->pbc->b_size || coorset1->pbc->c_size!=coorset2->pbc->c_size ||
          coorset1->pbc->angle1!=coorset2->pbc->angle1 || coorset1->pbc->angle2!=coorset2->pbc->angle2 || coorset1->pbc->angle3!=coorset2->pbc->angle3 )
           fprintf( stderr, "Warning in DistanceCoor! pbc data not equal between sets (taking set #1)!\n");
       
       for( ii=0; ii<nato; ii++)
         distances[ii] = PBC1( coorset1, sele1->selatm[ii], coorset2, sele2->selatm[ii]);
     }
     else
     {
       fprintf( stderr, "Only CUBIC Periodic Boundary Conditions accepted, sorry (try -nopbc)\n");
       exit(0);
     }
     return distances;
   }
   
   fprintf( stderr, "DistanceSelCoor could not find a workable option\n");
   exit(0);
   return distances;
}

double 
AngleCalc( double *a, double *b, double *c)
{
    double dist12, dist23, dist13, angle_cos;
   
    dist12 =  (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) ;
    dist13 =  (a[0]-c[0])*(a[0]-c[0]) + (a[1]-c[1])*(a[1]-c[1]) + (a[2]-c[2])*(a[2]-c[2]) ;
    dist23 =  (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) + (c[2]-b[2])*(c[2]-b[2]) ;
    angle_cos = (dist12 + dist23 - dist13)/(2.0*sqrt(dist12*dist23));

    if( abs(angle_cos) < 1 ) {
        return acos(angle_cos)*180.00/M_PI;
    } else {
        return 0.0;
    }
}

double DiheCalc (double *a, double *b, double *c, double *d)
{
   int       jj;
   double    angle, angle_cos, triple;
   double    Rij[3], Rkj[3], Rkl[3];
   double    Rim[3], Rln[3];
   double    Rkj_sq;
   double    fac_m, fac_n = 0.0;
   double    scal_m, scal_n = 0.0;

   // vector Rij
   Rij[0] = b[0] - a[0];
   Rij[1] = b[1] - a[1];
   Rij[2] = b[2] - a[2];
 
   // vector Rkj
   Rkj[0] = b[0] - c[0];
   Rkj[1] = b[1] - c[1];
   Rkj[2] = b[2] - c[2];
   Rkj_sq = wrd_vect_dotprod(Rkj, Rkj, 3);

    //printf ("%f \n", Rkj_sq); 

   // vector Rkl
   Rkl[0] = d[0] - c[0];
   Rkl[1] = d[1] - c[1];
   Rkl[2] = d[2] - c[2];
   
   
   // scal_m = Rij . Rkj
   // scal_n = Rkl . Rkj
   scal_m = wrd_vect_dotprod(Rij, Rkj, 3);
   scal_n = wrd_vect_dotprod(Rkl, Rkj, 3);
 
 
   // fac_m = (Rij . Rkj) / Rkj_sq
   // fac_n = (Rkl . Rkj) / Rkj_sq
   fac_m = scal_m / Rkj_sq;
   fac_n = scal_n / Rkj_sq;
 
 
   for (jj=0; jj<3; jj++)
   { 
    Rim[jj] =   Rij[jj] - fac_m * Rkj[jj];
    Rln[jj] = - Rkl[jj] + fac_n * Rkj[jj];	    
   }

   wrd_vect_norm(Rim, 3);
   wrd_vect_norm(Rln, 3);

   angle_cos = wrd_vect_dotprod(Rim, Rln, 3);
   
   if( abs(angle_cos) < 1 )
     angle = acos( angle_cos );
   else
     angle = 0.0;
  
   triple = wrd_vect_tripleprod(Rij, Rkl, Rkj);

   return wsign(triple) * (angle * 180.0) / M_PI;
}

//----------------------------------------------------------------------  
float * GetNearestImage( CoorSet *coorset, int atom1, int atom2 )
{
   float       *coords;
   float        dx, dy, dz;
   
   coords = calloc( 3, sizeof(float));
   
   if( coorset->pbc_flag == 0 )
   {
     coords[0] = coorset->xcoor[atom2-1];
     coords[1] = coorset->ycoor[atom2-1];
     coords[2] = coorset->zcoor[atom2-1];
     return 0;
   }
   
   if( coorset->pbc_flag == 1 )
   {
     if( coorset->pbc->angle1 == 0.0 && coorset->pbc->angle2 == 0.0 && coorset->pbc->angle3 == 0.0 )
     {
       /* dealing with CUBIC box */
       dx = coorset->xcoor[atom1-1] - coorset->xcoor[atom2-1];
       dx -= floor(dx/coorset->pbc->a_size) * coorset->pbc->a_size;
       dy = coorset->ycoor[atom1-1] - coorset->ycoor[atom2-1];
       dy -= floor(dy/coorset->pbc->b_size) * coorset->pbc->b_size;
       dz = coorset->zcoor[atom1-1] - coorset->zcoor[atom2-1];
       dz -= floor(dz/coorset->pbc->c_size) * coorset->pbc->c_size;
     }
   }
   else
   {
     fprintf( stderr, "Sorry, only CUBIC PBC are usable right now.");
     exit(1);
   }
   
   return 0;
}
//----------------------------------------------------------------------  
char * parseSlash( char *datastring )
{
  int       ii;
  char     *newstring;
  int       datastrlen;
  int       nslashes=0;
  int       slashloc[16];
  int       space1=0, space2=0, startpoint;
  
  datastrlen = strlen( datastring );
  
  for( ii=0; ii<datastrlen; ii++ )
  {
    if( datastring[ii] == '/')
    {
      slashloc[nslashes] = ii;
      nslashes ++;
    }
  }
  if( nslashes > 4 )
  {
    fprintf( stderr, "Selection string may have at most 4 sections (4 /s)\n" );
    return NULL;
  }
  space1 = (4-nslashes)*2;
  for( ii=0; ii<nslashes-1; ii++ )
  {
    //if( slashloc[ii] = slashloc[ii+1]-1 )  /* TODO Check this!! */
    if( slashloc[ii] == (slashloc[ii+1]-1) )
      space2++;
  }
  space2 *= 2;
  
  newstring = calloc( datastrlen+space1+space2+1, sizeof(char));
  for( ii=0; ii<(4-nslashes); ii++ )
  {
    sprintf( &newstring[ii*2], "/*");
  }
  startpoint = ii*2;
  for( ii=0; ii<datastrlen; ii++ )
  {
    newstring[startpoint] = datastring[ii] ;
    startpoint++;
    if( datastring[ii] == '/' && ( datastring[ii+1] == '\0' || datastring[ii+1] == '/' || datastring[ii+1] == '[') )
    {
      newstring[startpoint] = '*';
      startpoint ++;
    }
  }
  
  return newstring;
}

//----------------------------------------------------------------------  
void scopy( int * k, float * x1, int * n, float * x2,int * m)
{
  int i;
  for(i=0;i<(*k);i++){
    x2[i]=x1[i];
  }
  return;
}
void saxpy( int * k, float *a, float *x1, int *n, float *x2, int * m)
{
  int i;
  for(i=0;i<(*k);i++){
    x2[i]+=(*a)*x1[i];
  }
  return;
}
float snrm2( int * k, float * x1, int * n)
{
  int i;
  float snrm=0.0;
  for(i=0;i<(*k);i++){
    snrm+=x1[i]*x1[i];
  }
  snrm=sqrtf(snrm);
  return (snrm);
}
// ------------------------------------------------------------------
void DoubleSort_int( int nn, int *indices, int *data )
{
   int    ii;
   int    done, tmp;
   
   done=0;
   while( done!=1 )
   {
    done = 1;
    for( ii=1; ii<nn-1; ii++ )
     if( indices[ii]<indices[ii+1] )
     {
      tmp = indices[ii];
      indices[ii] = indices[ii+1];
      indices[ii+1] = tmp;
      tmp = data[ii];
      data[ii] = data[ii+1];
      data[ii+1] = tmp;
      done = 0;
     }
   }
}
// ------------------------------------------------------------------
int checkSorted( struct intlist * listA )
{
   int  ii;
   
   for( ii=1; ii<listA->nframes; ii++)
     if( listA->nframe[ii]<listA->nframe[ii-1] )
       return 0;
   return 1;
}
// ------------------------------------------------------------------
int wrd_isnumber( char *nstring )
{
  int     isnumber=1;
  int     number = 0;
  int     ii;
  
  number = atoi(nstring);
  if(number == 0)
  {
    isnumber = 0;
    return isnumber;
  }
  for( ii=0; ii<strlen(nstring); ii++)
  {
    if(isalpha(nstring[ii]) || !isalnum(nstring[ii]))
    {
      isnumber = 0;
      return isnumber;
    }
  }
  
  return isnumber;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------
// Mathematical Functions
// ------------------------------------------------------------------
// Compute the scalar product (a,b)
float dot_prod(float * a, float * b)
{
    float scalar=0.0;

    scalar = a[0] * b[0] + 
	     a[1] * b[1] + 
	     a[2] * b[2]  ;

    return scalar;
}

/* Compute the scalar product (a,b) between two vectors of size n
    Double precision, generic version */
double
wrd_vect_dotprod(double *a, double *b, size_t n)
{
    double scalar=0.0;
    int ii;

    for (ii=0; ii<n; ii++) {
        scalar += a[ii] * b[ii];
    }

    return scalar;
}

// ------------------------------------------------------------------
// Compute the vector product (a,b)
void cross_prod(float * a, float * b, float * res)
{
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}
void wrd_vect_crossprod(double *a, double *b, double *res)
{
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}
// ------------------------------------------------------------------
// Compute the triple product -> a.(b X c)
/*float triple_prod(float * a, float * b, float * c)
{
    float triple=0.0;
    float vec[3];

    cross_prod(b, c, vec);
    triple = dot_prod(a, vec);

    return triple;
}*/
double wrd_vect_tripleprod(double *a, double *b, double *c)
{
    double triple=0.0;
    double vec[3];

    wrd_vect_crossprod(b, c, vec);
    triple = wrd_vect_dotprod(a, vec, 3);

    return triple;
}
// ------------------------------------------------------------------
// Normalize a vector (vec)
void w_norm(float * vec)
{
    int ii;
    float square=0.0;
    float scalar=0.0;

    square = dot_prod(vec, vec);
    scalar = sqrt(square);

    for(ii=0; ii<3; ii++)
	vec[ii]/=scalar;
}


/* Normalize a vector (vec) of size n - double precision, generic version */
void 
wrd_vect_norm(double *vec, size_t n)
{
    size_t ii;
    double norm;

    norm = wrd_vect_dotprod(vec, vec, n);
    norm = sqrt(norm);
    
    for(ii=0; ii<n; ii++) {
	    vec[ii] /= norm;
    }

    return;
}


// ------------------------------------------------------------------
// Determine the sign of a float variable
float wsign(float value)
{
    const float eps = 1.0E-9;

    if(value <= eps)
	return -1.0;
    else
	return 1.0;
}
//----------------------------------------------------------------------
void Res3ToRes1(char *cResCode3, char *cResCode1)
{
  if     (strncmp(cResCode3, "ALA", 3)==0) strcpy(cResCode1, "A\0");
  else if(strncmp(cResCode3, "ARG", 3)==0) strcpy(cResCode1, "R\0");
  else if(strncmp(cResCode3, "ASN", 3)==0) strcpy(cResCode1, "N\0");
  else if(strncmp(cResCode3, "ASP", 3)==0) strcpy(cResCode1, "D\0");
  else if(strncmp(cResCode3, "CYS", 3)==0) strcpy(cResCode1, "C\0");
  else if(strncmp(cResCode3, "CYN", 3)==0) strcpy(cResCode1, "C\0");
  else if(strncmp(cResCode3, "GLU", 3)==0) strcpy(cResCode1, "E\0");
  else if(strncmp(cResCode3, "GLN", 3)==0) strcpy(cResCode1, "Q\0");
  else if(strncmp(cResCode3, "GLY", 3)==0) strcpy(cResCode1, "G\0");
  else if(strncmp(cResCode3, "HIS", 3)==0) strcpy(cResCode1, "H\0");
  else if(strncmp(cResCode3, "HSC", 3)==0) strcpy(cResCode1, "H\0");
  else if(strncmp(cResCode3, "HSD", 3)==0) strcpy(cResCode1, "H\0");
  else if(strncmp(cResCode3, "HSP", 3)==0) strcpy(cResCode1, "H\0");
  else if(strncmp(cResCode3, "HID", 3)==0) strcpy(cResCode1, "H\0");
  else if(strncmp(cResCode3, "ILE", 3)==0) strcpy(cResCode1, "I\0");
  else if(strncmp(cResCode3, "LEU", 3)==0) strcpy(cResCode1, "L\0");
  else if(strncmp(cResCode3, "LYS", 3)==0) strcpy(cResCode1, "K\0");
  else if(strncmp(cResCode3, "LYP", 3)==0) strcpy(cResCode1, "K\0");
  else if(strncmp(cResCode3, "MET", 3)==0) strcpy(cResCode1, "M\0");
  else if(strncmp(cResCode3, "PHE", 3)==0) strcpy(cResCode1, "F\0");
  else if(strncmp(cResCode3, "PRO", 3)==0) strcpy(cResCode1, "P\0");
  else if(strncmp(cResCode3, "SER", 3)==0) strcpy(cResCode1, "S\0");
  else if(strncmp(cResCode3, "THR", 3)==0) strcpy(cResCode1, "T\0");
  else if(strncmp(cResCode3, "TRP", 3)==0) strcpy(cResCode1, "W\0");
  else if(strncmp(cResCode3, "TYR", 3)==0) strcpy(cResCode1, "Y\0");
  else if(strncmp(cResCode3, "VAL", 3)==0) strcpy(cResCode1, "V\0");
  else if(strncmp(cResCode3, "RET", 3)==0) strcpy(cResCode1, "X\0");
  else if(strncmp(cResCode3, "GDP", 3)==0) strcpy(cResCode1, "d\0");
  else if(strncmp(cResCode3, "GTP", 3)==0) strcpy(cResCode1, "t\0");
  else if(strncmp(cResCode3, "Mg+", 3)==0) strcpy(cResCode1, "m\0");
  else if(strncmp(cResCode3, "ZMA", 3)==0) strcpy(cResCode1, "Z\0");
  else                                     strcpy(cResCode1, "U\0");
    return;
}


//======================================================================
// -------------- GROMACS SECTION --------------------------------------
