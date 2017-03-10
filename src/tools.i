%module tools

%{
# include "tools.h"
# include "fileio.h"
# include "qcprot.h"
%}

%include tools.h

%inline %{

// Coordinate Orientation.
// Given a reference coordinate structure and a selection, 
// it reorients the full molecule with respect to the latter
float CoorOrientRms ( CoorSet * refcoor,
                      CoorSet * oricoor,
                      Molecule *molecule,
                      char * selestring, 
                      char * fitrule )
{
   int          ii, jj, kk;
   float        rmsd;
   float        rotmatrix[3][3], transvec[3];
   float        aa, bb, cc, xcenter, ycenter, zcenter;
   int      nele, nato;
   float ** refarray;
   float ** oriarray;
   Selection sele;
   char   * fitstring;
   int      super;

    if      ( strcmp( fitrule, "True") == 0 ) 
    {
        super     = 1;
        fitstring = selestring;
    }
    else if ( strcmp( fitrule, "False") == 0 ) 
    {
        super     = 0;
    }
    else 
    {
        super     = 1;
        fitstring = fitrule;
    }

   if ( super )
   {
       // This creates the selection where you want
       // to make the fit between the two structures, 
       // in general, this is independent from
       // the selection used for the calculation
       // of the RMSD (see below)
       GetSele ( fitstring, &sele, molecule );
       
       // number of elements in the selection
       nele = sele.nselatm;

       // create the bidimensional arrays to give to
       // the superimpose function
       refarray = calloc ( nele, sizeof(float *));
       oriarray = calloc ( nele, sizeof(float *));
       for ( ii=0; ii<nele; ii++ ) 
       {
           refarray[ii]=calloc( 3, sizeof(float));
           oriarray[ii]=calloc( 3, sizeof(float));

           jj=sele.selatm[ii]-1;

           refarray[ii][0] = refcoor->xcoor[jj];
           refarray[ii][1] = refcoor->ycoor[jj];
           refarray[ii][2] = refcoor->zcoor[jj];

           oriarray[ii][0] = oricoor->xcoor[jj];
           oriarray[ii][1] = oricoor->ycoor[jj];
           oriarray[ii][2] = oricoor->zcoor[jj];
       } 
       // Generate the rotation matrix and translation vector
       // to induce the desidered orientation
       superimpose( nele, refarray, oriarray, rotmatrix, transvec );

       // Apply the coordinate transformation
       nato = oricoor->nato;
       for ( ii=0; ii<nato; ii++)
       {
           aa = oricoor->xcoor[ii]*rotmatrix[0][0] + oricoor->ycoor[ii]*rotmatrix[0][1] + oricoor->zcoor[ii]*rotmatrix[0][2] ;
           bb = oricoor->xcoor[ii]*rotmatrix[1][0] + oricoor->ycoor[ii]*rotmatrix[1][1] + oricoor->zcoor[ii]*rotmatrix[1][2] ;
           cc = oricoor->xcoor[ii]*rotmatrix[2][0] + oricoor->ycoor[ii]*rotmatrix[2][1] + oricoor->zcoor[ii]*rotmatrix[2][2] ;
           oricoor->xcoor[ii] = aa + transvec[0]; //+ xcenter
           oricoor->ycoor[ii] = bb + transvec[1]; //+ ycenter
           oricoor->zcoor[ii] = cc + transvec[2]; //+ zcenter
       }
       free(sele.selatm);
   }

   // Now take the selection where you want to compute
   // the RMSD, in general, this is independent from
   // the selection used for fitting the molecules
   GetSele ( selestring, &sele, molecule );

   // number of elements in the selection
   nele = sele.nselatm;

   // Compute the RMSD on the given selection
   rmsd = 0.0;
   for ( ii=0; ii<nele; ii++ ) 
   {
       jj=sele.selatm[ii]-1;
       rmsd += ((oricoor->xcoor[ jj ]-refcoor->xcoor[ jj ])*(oricoor->xcoor[ jj ]-refcoor->xcoor[ jj ]) +
                (oricoor->ycoor[ jj ]-refcoor->ycoor[ jj ])*(oricoor->ycoor[ jj ]-refcoor->ycoor[ jj ]) +
                (oricoor->zcoor[ jj ]-refcoor->zcoor[ jj ])*(oricoor->zcoor[ jj ]-refcoor->zcoor[ jj ]));
   } 
   rmsd /= nele;
   rmsd = sqrt ( rmsd );

   free(refarray);
   free(oriarray);

   free(sele.selatm);

   return(rmsd);
}

float MolRmsd(Molecule *mol1, Molecule *mol2, char *selestring );

// -------------------------------[ OBSOLETE ]-------------------------------------------



//float tra[3];
//float rot[3][3];
//
//
//float ** ArrayCoor ( struct coor * strcoor , struct selestr_ sele ) {
//    int      ii,jj;
//    int      nato;
//    float ** coor;
//
//    printf ("%10f\n",strcoor->xcoor[sele.selatm[ii]-1]);
//    nato = sele.nselatm;
//    coor = calloc ( nato, sizeof(float *));
//    for ( ii=0; ii<nato; ii++ ) {
//        coor[ii]=calloc( 3, sizeof(float));
//        coor[ii][0] = strcoor->xcoor[sele.selatm[ii]-1];
//        coor[ii][1] = strcoor->ycoor[sele.selatm[ii]-1];
//        coor[ii][2] = strcoor->zcoor[sele.selatm[ii]-1];
//        printf ("%10f %10f %10f\n",coor[ii][0],coor[ii][1],coor[ii][2]);
//        printf( "%10d %10d\n",sele.selatm[ii], nato );
//    }
//    return coor;
//
//}

// Coordinate Orientation.
// Given a reference coordinate structure and a selection, 
// it reorients the full molecule with respect to the latter
//float CoorOrient ( struct coor * refcoor, struct coor * oricoor, struct selestr_ sele, int super  )
//{
//   int          ii, jj, kk;
//   float        rmsd;
//   float        rotmatrix[3][3], transvec[3];
//   float        aa, bb, cc, xcenter, ycenter, zcenter;
//   int      nele, nato;
//   float ** refarray;
//   float ** oriarray;
//
//   // number of elements in the selection
//   nele = sele.nselatm;
//
//   if ( super )
//   {
//       // create the bidimensional arrays to give to
//       // the superimpose function
//       refarray = calloc ( nele, sizeof(float *));
//       oriarray = calloc ( nele, sizeof(float *));
//       for ( ii=0; ii<nele; ii++ ) 
//       {
//           refarray[ii]=calloc( 3, sizeof(float));
//           oriarray[ii]=calloc( 3, sizeof(float));
//
//           jj=sele.selatm[ii]-1;
//
//           refarray[ii][0] = refcoor->xcoor[jj];
//           refarray[ii][1] = refcoor->ycoor[jj];
//           refarray[ii][2] = refcoor->zcoor[jj];
//
//           oriarray[ii][0] = oricoor->xcoor[jj];
//           oriarray[ii][1] = oricoor->ycoor[jj];
//           oriarray[ii][2] = oricoor->zcoor[jj];
//       } 
//       // Generate the rotation matrix and translation vector
//       // to induce the desidered orientation
//       superimpose( nele, refarray, oriarray, rotmatrix, transvec );
//       free(refarray);
//       free(oriarray);
//
//       // Apply the coordinate transformation
//       nato = oricoor->nato;
//       //printf("NATO %d\n",nato);
//       for ( ii=0; ii<nato; ii++)
//       {
//           aa = oricoor->xcoor[ii]*rotmatrix[0][0] + oricoor->ycoor[ii]*rotmatrix[0][1] + oricoor->zcoor[ii]*rotmatrix[0][2] ;
//           bb = oricoor->xcoor[ii]*rotmatrix[1][0] + oricoor->ycoor[ii]*rotmatrix[1][1] + oricoor->zcoor[ii]*rotmatrix[1][2] ;
//           cc = oricoor->xcoor[ii]*rotmatrix[2][0] + oricoor->ycoor[ii]*rotmatrix[2][1] + oricoor->zcoor[ii]*rotmatrix[2][2] ;
//           oricoor->xcoor[ii] = aa + transvec[0]; //+ xcenter
//           oricoor->ycoor[ii] = bb + transvec[1]; //+ ycenter
//           oricoor->zcoor[ii] = cc + transvec[2]; //+ zcenter
//       }
//   }
//
//   // Compute the RMSD on the given selection
//   rmsd = 0.0;
//   for ( ii=0; ii<nele; ii++ ) 
//   {
//       jj=sele.selatm[ii]-1;
//       //printf ( "> %10d %10f\n",jj,oricoor->xcoor[ jj ]);
//       //printf ( "> %10d %10f\n",jj,refcoor->xcoor[ jj ]);
//       rmsd += ((oricoor->xcoor[ jj ]-refcoor->xcoor[ jj ])*(oricoor->xcoor[ jj ]-refcoor->xcoor[ jj ]) +
//                (oricoor->ycoor[ jj ]-refcoor->ycoor[ jj ])*(oricoor->ycoor[ jj ]-refcoor->ycoor[ jj ]) +
//                (oricoor->zcoor[ jj ]-refcoor->zcoor[ jj ])*(oricoor->zcoor[ jj ]-refcoor->zcoor[ jj ]));
//   } 
//   rmsd /= nele;
//   rmsd = sqrt ( rmsd );
//   return(rmsd);
//}
//

%}
