// ------------------------------------------------------------------
// Copyright (C) 2012  University of Modena and Reggio Emilia and
//                     University of Zurich
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
/*! \file volumes.c
 \brief volumes/superposition module
 
 Source code for the Volumes/Superposition module.
*/
#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "volumes.h"
//
// --------------------------------------------------------------------
int linecount( FILE *inpfile )
{
  int     nlines = 0;
  char    textline[1024];
  long long    offset;
  
  offset = ftell( inpfile );
  rewind( inpfile );
  while ( !feof(inpfile))
  {
    if ( textline[0] != '#' && textline[0] != '\n' )
      ++nlines;
    if(fgets(textline, 1024, inpfile)==NULL && ( !feof(inpfile) || ferror(inpfile) ))
    {
      fprintf(stderr, "Warning! Premature end of file reached!\n");
    }
  }
  fseek( inpfile, offset, SEEK_SET );
  
  return nlines-1;
}
// -------
BoxDim * buildBox( moldata *moldata_s, moldata **moldata_l, int nligs, float step )
{
  int             ii, jj;
  BoxDim         *boxdim;
  float           xorg=1000.0, yorg=1000.0, zorg=1000.0;
  float           xmax=-1000.0, ymax=-1000.0, zmax=-1000.0;
  
  boxdim = (BoxDim *)calloc( 1, sizeof(BoxDim));
  for( ii=0; ii<moldata_s->natoms; ii ++)
  {
    if( xorg > moldata_s->coords[0][ii] )
      xorg = moldata_s->coords[0][ii];
    if( yorg > moldata_s->coords[1][ii] )
      yorg = moldata_s->coords[1][ii];
    if( zorg > moldata_s->coords[2][ii] )
      zorg = moldata_s->coords[2][ii];
    if( xmax < moldata_s->coords[0][ii] )
      xmax = moldata_s->coords[0][ii];
    if( ymax < moldata_s->coords[1][ii] )
      ymax = moldata_s->coords[1][ii];
    if( zmax < moldata_s->coords[2][ii] )
      zmax = moldata_s->coords[2][ii];
  }
  
  for( jj=0; jj<nligs; jj++ )
    for( ii=0; ii<moldata_l[jj]->natoms; ii ++)
    {
      if( xorg > moldata_l[jj]->coords[0][ii] )
        xorg = moldata_l[jj]->coords[0][ii];
      if( yorg > moldata_l[jj]->coords[1][ii] )
        yorg = moldata_l[jj]->coords[1][ii];
      if( zorg > moldata_l[jj]->coords[2][ii] )
        zorg = moldata_l[jj]->coords[2][ii];
      if( xmax < moldata_l[jj]->coords[0][ii] )
        xmax = moldata_l[jj]->coords[0][ii];
      if( ymax < moldata_l[jj]->coords[1][ii] )
        ymax = moldata_l[jj]->coords[1][ii];
      if( zmax < moldata_l[jj]->coords[2][ii] )
        zmax = moldata_l[jj]->coords[2][ii];
    }
  boxdim->xorg = xorg - 2.0;
  boxdim->yorg = yorg - 2.0;
  boxdim->zorg = zorg - 2.0;
  boxdim->xsize = xmax - xorg + 4.0;
  boxdim->ysize = ymax - yorg + 4.0;
  boxdim->zsize = zmax - zorg + 4.0;
  
  boxdim->step = step;
  
  boxdim->xbins = ceil( boxdim->xsize/boxdim->step );
  boxdim->ybins = ceil( boxdim->ysize/boxdim->step );
  boxdim->zbins = ceil( boxdim->zsize/boxdim->step );
  boxdim->totalbins = boxdim->xbins*boxdim->ybins*boxdim->zbins;
  
  fprintf( stdout, "# Box origin: %8.3f%8.3f%8.3f                   #\n", boxdim->xorg, boxdim->yorg, boxdim->zorg);
  fprintf( stdout, "# Box size:   %8.3f%8.3f%8.3f                   #\n", boxdim->xsize, boxdim->ysize, boxdim->zsize);
  fprintf( stdout, "# Box bins:   %8d%8d%8d                   #\n", boxdim->xbins, boxdim->ybins, boxdim->zbins);
  fprintf( stdout, "# Box step:   %8.3f                                   #\n", boxdim->step);
  
  return boxdim;
}
// -------
BoxDim * deriveBox( struct inp_vol *inpvol )
{
  BoxDim         *boxdim;
  
  boxdim = (BoxDim *)calloc( 1, sizeof(BoxDim));
  if( inpvol->centerboxflag )
  {
    boxdim->xorg = inpvol->box_xcenter - inpvol->box_xsize/2;
    boxdim->yorg = inpvol->box_ycenter - inpvol->box_ysize/2;
    boxdim->zorg = inpvol->box_zcenter - inpvol->box_zsize/2;
  }
  else if( inpvol->cornerboxflag )
  {
    boxdim->xorg = inpvol->box_xcorner;
    boxdim->yorg = inpvol->box_ycorner;
    boxdim->zorg = inpvol->box_zcorner;
  }
  boxdim->xsize = inpvol->box_xsize;
  boxdim->ysize = inpvol->box_ysize;
  boxdim->zsize = inpvol->box_zsize;
  
  boxdim->step = inpvol->step;
  
  boxdim->xbins = ceil( boxdim->xsize/boxdim->step );
  boxdim->ybins = ceil( boxdim->ysize/boxdim->step );
  boxdim->zbins = ceil( boxdim->zsize/boxdim->step );
  boxdim->totalbins = boxdim->xbins*boxdim->ybins*boxdim->zbins;
  
  fprintf( stdout, "# Box origin: %8.3f%8.3f%8.3f                   #\n", boxdim->xorg, boxdim->yorg, boxdim->zorg);
  fprintf( stdout, "# Box size:   %8.3f%8.3f%8.3f                   #\n", boxdim->xsize, boxdim->ysize, boxdim->zsize);
  fprintf( stdout, "# Box bins:   %8d%8d%8d                   #\n", boxdim->xbins, boxdim->ybins, boxdim->zbins);
  fprintf( stdout, "# Box step:   %8.3f                                   #\n", boxdim->step);
  
  return boxdim;
}
// -------
short *** mkGrid( BoxDim *boxdim )
{
  int       ii, jj;
  short  ***gridarray;
  
  // build grid
  // one big vector with appropriate indexes
  gridarray = (short ***)calloc( boxdim->xbins, sizeof( short **));
  gridarray[0] = (short **)calloc( boxdim->xbins*boxdim->ybins, sizeof(short *));
  for( ii=0; ii<boxdim->xbins; ii++)
    gridarray[ii] = gridarray[0] + ii*boxdim->ybins;
  gridarray[0][0] = (short *)calloc( boxdim->xbins*boxdim->ybins*boxdim->zbins, sizeof(short));
  for( ii=0; ii<boxdim->xbins; ii++ )
    for( jj=0; jj<boxdim->ybins; jj++ )
      gridarray[ii][jj] = gridarray[0][0] + ii*boxdim->ybins*boxdim->zbins + jj*boxdim->zbins;
  
  /*
  // regular tri-dimentional matrix
  gridarray = (short ***)calloc( boxdim->xbins, sizeof( short **));
  for( ii=0; ii<boxdim->xbins; ii++ )
  {
    gridarray[ii] = (short **)calloc( boxdim->ybins, sizeof( short *));
    for( jj=0; jj<boxdim->ybins ; jj++ )
      gridarray[ii][jj] = (short *)calloc( boxdim->zbins, sizeof( short ));
  }*/
  
  
  return(gridarray);
}
// -------
int scanGrid(moldata *moldata1, short ***gridarray, BoxDim *boxdim)
{
  int         ii, jj, kk, ll;
  float       aaa, bbb, ccc;
  float       xcoor, ycoor, zcoor, rad2;
  float       thisx, thisy, thisz;
  float       dist;
  float       xgridorg, ygridorg, zgridorg;
  float       xgridstart, ygridstart, zgridstart;
  float       xgridstop, ygridstop, zgridstop;
  float       step = 0.5;
  
  step = boxdim->step;
  xgridorg = boxdim->xorg;
  ygridorg = boxdim->yorg;
  zgridorg = boxdim->zorg;

  for( ii=0; ii<boxdim->xbins; ii++ )
    for( jj=0; jj<boxdim->ybins; jj++ )
      for( kk=0; kk<boxdim->zbins; kk++ )
        gridarray[ii][jj][kk] = 0;
  
  for( ii=0; ii<moldata1->natoms; ii++ )
  {
    xcoor = moldata1->coords[0][ii];
    ycoor = moldata1->coords[1][ii];
    zcoor = moldata1->coords[2][ii];
    rad2 = moldata1->radii[ii] * moldata1->radii[ii];
    
    /*
    // no optimization: brute approach
    for( jj=0; jj<boxdim->xbins; jj++ )
    {
      //fprintf( stdout, "DEBUG - jj is now %d\n", jj);
      thisx = xgridorg + jj*step;
      for( kk=0; kk<boxdim->ybins; kk++ )
      {
        //fprintf( stdout, "DEBUG - kk is now %d\n", kk);
        thisy = ygridorg + kk*step;
        for( ll=0; ll<boxdim->zbins; ll++ )
        {
          //fprintf( stdout, "DEBUG - ll is now %d\n", ll);
          if( gridarray[jj][kk][ll] == 1 )
            continue;
          thisz = zgridorg + ll*step;
          //fprintf( stdout, "DEBUG - comparing %8.3f%8.3f%8.3f with %8.3f%8.3f%8.3f -> ", thisx,thisy,thisz,xcoor,ycoor,zcoor);
          //dist = (thisx-moldata1->coords[0][ii] * thisx-moldata1->coords[0][ii]) +
                 //(thisy-moldata1->coords[1][ii] * thisy-moldata1->coords[1][ii]) +
                 //(thisz-moldata1->coords[2][ii] * thisz-moldata1->coords[2][ii]);
          dist = ((thisx-xcoor) * (thisx-xcoor)) +
                 ((thisy-ycoor) * (thisy-ycoor)) +
                 ((thisz-zcoor) * (thisz-zcoor));
          //fprintf( stdout, "DEBUG - dist is %8.3f (rad2: %8.3f) %8.3f %8.3f\n", dist, rad2, (thisx-xcoor * thisx-xcoor), thisx-xcoor);
          if( dist<rad2 )
          {
            gridarray[jj][kk][ll] = (short) 1;
          }
        }
      }
    }
    */
    
    // finding the approx place in the grid to save on cycles
    aaa = xcoor - xgridorg;
    bbb = ycoor - ygridorg;
    ccc = zcoor - zgridorg;
    aaa -= 2 * moldata1->radii[ii];
    bbb -= 2 * moldata1->radii[ii];
    ccc -= 2 * moldata1->radii[ii];
    xgridstart = (int)(aaa/step);
    ygridstart = (int)(bbb/step);
    zgridstart = (int)(ccc/step);
    if( xgridstart < 0 ) xgridstart = 0;
    if( ygridstart < 0 ) ygridstart = 0;
    if( zgridstart < 0 ) zgridstart = 0;
    aaa += 4 * moldata1->radii[ii];
    bbb += 4 * moldata1->radii[ii];
    ccc += 4 * moldata1->radii[ii];
    xgridstop = (int)(aaa/step);
    ygridstop = (int)(bbb/step);
    zgridstop = (int)(ccc/step);
    if( xgridstop > boxdim->xbins ) xgridstop = boxdim->xbins;
    if( ygridstop > boxdim->ybins ) ygridstop = boxdim->ybins;
    if( zgridstop > boxdim->zbins ) zgridstop = boxdim->zbins;

    for( jj=xgridstart; jj<xgridstop; jj++ )
    {
      thisx = xgridorg + jj*step;
      for( kk=ygridstart; kk<ygridstop; kk++ )
      {
        thisy = ygridorg + kk*step;
        for( ll=zgridstart; ll<zgridstop; ll++ )
        {
          if( gridarray[jj][kk][ll] == 1 )
            continue;
          thisz = zgridorg + ll*step;
          dist = ((thisx-moldata1->coords[0][ii]) * (thisx-moldata1->coords[0][ii])) +
                 ((thisy-moldata1->coords[1][ii]) * (thisy-moldata1->coords[1][ii])) +
                 ((thisz-moldata1->coords[2][ii]) * (thisz-moldata1->coords[2][ii]));
          if( dist<rad2 )
          {
            gridarray[jj][kk][ll] = (short) 1;
          }
        }
      }
    }
    
  }
  
  return(0);
}
// -------
float * assignRadii( int natoms, char **types )
{
  int         ii;
  float      *radii;
  radii = (float *)calloc( natoms, sizeof(float));
  
  for( ii=0; ii<natoms; ii++ )
  {
    if( types[ii][0] == 'C' && types[ii][0] != 'L' && types[ii][0] != 'l' )
      radii[ii] = 1.70;
    else if( types[ii][0] == 'C' && (types[ii][0] != 'L'|| types[ii][0] != 'l') )
      radii[ii] = 1.75;
    else if( types[ii][0] == 'N' )
      radii[ii] = 1.55;
    else if( types[ii][0] == 'O' )
      radii[ii] = 1.52;
    else if( types[ii][0] == 'H' )
      radii[ii] = 1.00;
    else if( types[ii][0] == 'S' )
      radii[ii] = 1.80;
    else if( types[ii][0] == 'F' )
      radii[ii] = 1.47;
    else if( types[ii][0] == 'P' )
      radii[ii] = 1.80;
    else if( types[ii][0] == 'I' )
      radii[ii] = 1.98;
    else if( types[ii][0] == 'B' && (types[ii][0] != 'R'|| types[ii][0] != 'r') )
      radii[ii] = 1.85;
    else if( types[ii][1] == 'H' )
      radii[ii] = 1.00;
    //else if( types[ii][0] == '' )
    //  radii[ii] = ;
    else
    {
      fprintf( stderr, "Could not assign radius to atom %s\n", types[ii]);
      exit(10);
    }
    //fprintf( stdout, "DEBUG: %3d %8.3f\n", ii, radii[ii]);
  }
  
  return radii;
}
// -------
moldata * getMolData( char *molfilename )
{
  int         ii;
  Molecule   *mol1;
  moldata    *moldata1;
  
  moldata1 = (moldata *)calloc(1, sizeof(moldata));
  
  mol1 = ReadMolecule( molfilename, wrd_whichmoltype(molfilename));
  moldata1->natoms = mol1->nato;
  
  moldata1->coords = (float **)calloc( 3, sizeof(float *));
  moldata1->coords[0] = (float *)calloc( 3*moldata1->natoms, sizeof(float));
  for( ii=0; ii<3; ii++ )
    moldata1->coords[ii] = moldata1->coords[0] + moldata1->natoms*ii;
  
  for( ii=0; ii<mol1->nato; ii++ )
  {
    moldata1->coords[0][ii] = mol1->coor.xcoor[ii];
    moldata1->coords[1][ii] = mol1->coor.ycoor[ii];
    moldata1->coords[2][ii] = mol1->coor.zcoor[ii];
  }
  
  moldata1->radii = assignRadii( moldata1->natoms, mol1->rawmol.atmtype);
  
  return(moldata1);
}
// -------
//int calcVolumes(struct inp_vol *inpvol, moldata *moldata_l, float ***gridarray_l )
//{
  //int       ii, jj, kk;
  
  
  //for( ll=0; ll<liglist.ntrj; ll+=inpvol.nt )
  //{
    //scanGrid( moldata_l[ll], gridarray_l, boxdim );
    ////fprintf( stdout, "| Molname: %20s        |\n", liglist.trjname[ll]); fflush(stdout);
    //fprintf( stdout, "%-30s ", liglist.trjname[ll]);
    //mol_vol = 0.0;
    //for( ii=0; ii<boxdim->xbins; ii++ )
    //{
      //for( jj=0; jj<boxdim->ybins ; jj++ )
      //{
        //for( kk=0; kk<boxdim->zbins ; kk++ )
        //{
          //if( gridarray_s[ii][jj][kk] && gridarray_l[ii][jj][kk] )
          //{
            //mol_vol += 1.0;
          //}
        //}
      //}
    //}
    //mol_vol = 0.0;
    //for( ii=0; ii<boxdim->totalbins; ii++ )
      //mol_vol += (float) (gridarray_l[0][0][ii]);
    ////fprintf( stdout, "| Molecular Volume:    %8.3f        |\n", mol_vol/8.0);
    //fprintf( stdout, " %8.3f", mol_vol/factor);
    
    //com_vol = 0.0;
    //for( ii=0; ii<boxdim->totalbins; ii++ )
      //com_vol += (float) (gridarray_s[0][0][ii] & gridarray_l[0][0][ii]);
    ////fprintf( stdout, "| Common Volume:       %8.3f        |\n", com_vol/8.0);
    //fprintf( stdout, " %8.3f", com_vol/factor);
    
    //fprintf( stdout, "%9.4f\n", (2*com_vol-mol_vol)/sup_vol);
  //}
  ////fprintf( stdout, "\\--------------------------------------/ \n"); fflush(stdout);
  
  //return(0);
//}
// -------
int CalcVolumes( char **input, int inp_index, FILE *output )
{
  int               ii, jj, kk, ll;
  struct inp_vol    inpvol;
  char              buffer[1024];
  int               gotit=0;
  moldata          *moldata_s;
  moldata         **moldata_l;
  TrjList           liglist;
  BoxDim           *boxdim;
  short          ***gridarray_s;
  short          ***gridarray_l;
  //short          ***gridarray_o;
  float             mol_vol, sup_vol, com_vol, factor;
  float           **results;
  
  inpvol.boxflag = 0;
  inpvol.stepflag = 0;
  inpvol.step = 0.5;
  inpvol.mthread = 0;
  
  while( strncmp (buffer, "END", 3))
  {
    gotit = 0;
    sprintf( buffer, "%s", input[inp_index]);
    if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
      gotit = 1;
    else if ( !strncmp(buffer, "--SUPERMOL", 10))
    {
      sscanf(buffer, "--SUPERMOL %[^\n]%*c ", inpvol.supermol);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--MOLLIST", 9))
    {
      sscanf(buffer, "--MOLLIST %[^\n]%*c ", inpvol.mollist);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--BOXFILE", 9))
    {
      inpvol.boxflag = 1;
      sscanf(buffer, "--BOXFILE %[^\n]%*c ", inpvol.boxfile);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CENTERBOX", 11))
    {
      inpvol.boxflag = 1;
      inpvol.centerboxflag = 1;
      sscanf(buffer, "--CENTERBOX %f %f %f %f %f %f ", &inpvol.box_xcenter, &inpvol.box_ycenter, &inpvol.box_zcenter, &inpvol.box_xsize, &inpvol.box_ysize, &inpvol.box_zsize);
      gotit = 1;
    }
    else if ( !strncmp(buffer, "--CORNERBOX", 11))
    {
      inpvol.boxflag = 1;
      inpvol.cornerboxflag = 1;
      sscanf(buffer, "--CORNERBOX %f %f %f %f %f %f ", &inpvol.box_xcorner, &inpvol.box_ycorner, &inpvol.box_zcorner, &inpvol.box_xsize, &inpvol.box_ysize, &inpvol.box_zsize);
      gotit = 1;
    }
     else if ( !strncmp(buffer, "--STEP", 6))
    {
      inpvol.stepflag = 1;
      sscanf(buffer, "--STEP %f ", &inpvol.step);
      gotit = 1;
    }
     else if ( !strncmp(buffer, "--NT", 4))
    {
      inpvol.mthread = 1;
      sscanf(buffer, "--NT %d ", &inpvol.nt);
      gotit = 1;
    }
   if( gotit==0 )
    {
     fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
     exit(5);
    }
    inp_index++;
  }
  
  factor = 1/(inpvol.step*inpvol.step*inpvol.step);
  moldata_s = getMolData( inpvol.supermol );
  
  ReadTrjList( inpvol.mollist, &liglist );
  moldata_l = (moldata **)calloc( liglist.ntrj, sizeof( moldata * ));
  for( ii=0; ii<liglist.ntrj; ii++ )
    moldata_l[ii] = getMolData( liglist.trjname[ii] );
  
  fprintf( stdout, "#--------------------------------------------------------#\n"); fflush(stdout);
  if( inpvol.boxflag )
    boxdim = deriveBox( &inpvol );
  else
    boxdim = buildBox( moldata_s, moldata_l, liglist.ntrj, inpvol.step );
  
  gridarray_s = mkGrid( boxdim );
  gridarray_l = mkGrid( boxdim );
  //gridarray_o = mkGrid( boxdim );
  
  scanGrid( moldata_s, gridarray_s, boxdim );
  sup_vol = 0.0;
  for( ii=0; ii<boxdim->totalbins; ii++ )
    if( gridarray_s[0][0][ii] )
      sup_vol += 1.0;
  fprintf( stdout, "# Supermol Volume: %8.3f                              #\n", sup_vol/factor);
  fprintf( stdout, "#--------------------------------------------------------#\n" );
  fprintf( stdout, "# Molname                        Volume     Vin     Vdiff \n" );
  //fprintf( stdout, "# ------------------------------------ #\n" );

  results = (float **) calloc( liglist.ntrj, sizeof(float *));
  for( ii=0; ii<liglist.ntrj; ii++)
    results[ii] = (float *) calloc( 2, sizeof(float));
  
  for( ll=0; ll<liglist.ntrj; ll++)
  {
    scanGrid( moldata_l[ll], gridarray_l, boxdim );
    //fprintf( stdout, "| Molname: %20s        |\n", liglist.trjname[ll]); fflush(stdout);
    fprintf( stdout, "%-30s ", liglist.trjname[ll]);
    mol_vol = 0.0;
    for( ii=0; ii<boxdim->xbins; ii++ )
    {
      for( jj=0; jj<boxdim->ybins ; jj++ )
      {
        for( kk=0; kk<boxdim->zbins ; kk++ )
        {
          if( gridarray_s[ii][jj][kk] && gridarray_l[ii][jj][kk] )
          {
            mol_vol += 1.0;
          }
        }
      }
    }
    mol_vol = 0.0;
    for( ii=0; ii<boxdim->totalbins; ii++ )
      mol_vol += (float) (gridarray_l[0][0][ii]);
    //fprintf( stdout, "| Molecular Volume:    %8.3f        |\n", mol_vol/8.0);
    fprintf( stdout, " %8.3f", mol_vol/factor);
    
    com_vol = 0.0;
    for( ii=0; ii<boxdim->totalbins; ii++ )
      com_vol += (float) (gridarray_s[0][0][ii] & gridarray_l[0][0][ii]);
    //fprintf( stdout, "| Common Volume:       %8.3f        |\n", com_vol/8.0);
    fprintf( stdout, " %8.3f", com_vol/factor);
    
    fprintf( stdout, "%9.4f\n", (2*com_vol-mol_vol)/sup_vol);
  }
  //fprintf( stdout, "\\--------------------------------------/ \n"); fflush(stdout);
  
  return(0);

}
