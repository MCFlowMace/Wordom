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
/*! \file kga.h
 \brief Kinetic Grouping Analysis module
 
 Headers and structures for the Kinetic Grouping Analysis (KGA) and
 Free Energy Profile (FEP) module.
*/
#ifndef KGA
#define KGA
#pragma once
// ------------------------------------------------------------------
// ** logbin Structure
struct inp_lb
{
  int            BpD;
  int            target;
  char           clusfile[256];
};
// ------------------------------------------------------------------
// ** ka Structure
struct inp_ka
{
  int            nnodes;
  int            tcomm;
  char           clusfile[256];
};
// ------------------------------------------------------------------
// ** basin Structure
struct inp_basin
{
  int            target;
  int            tcomm;
  char           clusfile[256];
};
// ------------------------------------------------------------------
// ** mfpt Structure
struct inp_mfpt
{
  int            fixed;
//int            nit;
//int            target;
//int            nonsymm;
//float          temp;
//char           clusfile[256];
};
// ------------------------------------------------------------------
// ** mfptnet Structure
struct inp_mfptnet
{
  int            nit;
  int            target;
  int            nonsymm;
  float          temp;
  char           linkfile[256];
};
// ------------------------------------------------------------------
// ** pfoldf Structure
struct inp_pfoldf
{
  int            nit;
  int            target;
  int            target2;
  int            nonsymm;
  float          temp;
  float          lambda;
  char           clusfile[256];
};
// ------------------------------------------------------------------
// ** pfoldfnet Structure
struct inp_pfoldfnet
{
  int            nit;
  int            target;
  int            target2;
  int            nonsymm;
  float          temp;
  float          lambda;
  char           linkfile[256];
};
// ------------------------------------------------------------------
// ** equil Structure
struct inp_equil
{
  int            nit;
  char           linkfile[256];
};
// ------------------------------------------------------------------
int LogBin ( char **input, int inp_index, FILE *output ) ;
int Ka ( char **input, int inp_index, FILE *output ) ;
int Basin ( char **input, int inp_index, FILE *output ) ; // I called the function "basin" because one isolates one basin. Just a suggestion
int MFPT ( char **input, int inp_index, FILE *output ) ; //Input is the clustered timeseries, output the mfpt profile
int MFPTnet ( char **input, int inp_index, FILE *output ) ; //Input is the network, output the mfpt profile
int PFoldf ( char **input, int inp_index, FILE *output ); //Input is the clustered timeseries, output the mfpt profile
int PFoldfnet ( char **input, int inp_index, FILE *output ); //Input is the network, output the mfpt profile
int Equil ( char **input, int inp_index, FILE *output );  //Input is the (out of equilibrium) linkfile, output the equilibrated network file
// ==================================================================
// ** Tarjan reduction Structures
typedef struct _intlist
{
  int     ndata;
  int *   data;
} intlistd;
// ------------------------------------------------------------------
typedef struct _intdlist
{
  int     ndata;
  int *   keys;
  int *   values;
} intdlist;
// ------------------------------------------------------------------
typedef struct _graph
{
  int               nnodes;
  int     *         nodes;
  intlistd *         leenklist;
} graph;
// ------------------------------------------------------------------
// ** Tarjan reduction Functions
intlistd * readSymbTrj( char *filename );
graph * sTrj2Graph( intlistd *strj );
int imin( int num1, int num2 );
int whereinlist( intlistd *list, int number );
int addtolist( intlistd *list, int number );
int popfromlist( intlistd *list );
void printlist( intlistd *ilist );
int addtodlist( intdlist *dlist, int key, int value );
int getfromdlist( intdlist *dlist, int key );
void printdlist( intdlist *dlist );
int whereingraph( graph *graph1, int number );
int addtograph1( graph *lgraph, int key, int value); // this adds to the current value of the 'value'
int addtograph2( graph *lgraph, int key, int value); // this changes the current value of the 'value'
intlistd * getgraphval( graph *lgraph, int key );
void printgraph( graph *grapho );
graph * Tarjan_ite( graph *grapho );
int ELEC( char **input, int inp_index );

#endif
