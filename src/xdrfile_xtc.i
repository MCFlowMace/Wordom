%module xdrfile_xtc

%{
# include "xdrfile.h"
# include "xdrfile_xtc.h"
%}

%include xdrfile.h
%include xdrfile_xtc.h

extern int read_xtc_natoms(char *fn,int *natoms);
extern int read_xtc(XDRFILE *xd,int natoms,int *step,float *time, matrix box,rvec *x,float *prec);
extern int write_xtc(XDRFILE *xd,int natoms,int step,float time, matrix box,rvec *x,float prec);
