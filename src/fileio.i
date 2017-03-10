%module fileio

%{
# include "fileio.h"
# include <assert.h>
static int myErr = 0; // flag to save error state
%}

%include "numpy.i"

%inline %{
// helper functions (object creators, data extractor)
Molecule * NewMolecule()
{
   Molecule *molecule;
   molecule = calloc(1, sizeof( Molecule));
   return molecule;
};
void DesMolecule(Molecule *molecule)
{
   // eventually free all stuff inside
   DelMolecule(molecule);
   free(molecule);
};
// ------------------------------------------------------------------
Selection * NewSele()
{
   Selection     *sele;
   sele = calloc( 1, sizeof(Selection));
   return sele;
};
int geti( int *intvec, int index)
{
   return intvec[index];
};
// ------------------------------------------------------------------
Traj * NewTraj()
{
   Traj     *trj;
   trj = calloc( 1, sizeof(Traj));
   return trj;
};
void DesTraj(Traj *trj)
{
   free(trj);
};
// ------------------------------------------------------------------
struct trjh * NewTrjh()
{
   Trjh    *trjh;
   trjh = calloc( 1, sizeof(Trjh));
   return trjh;
};
// ------------------------------------------------------------------
CoorSet * NewCoor()
{
  CoorSet *coords;
  coords = calloc( 1, sizeof(CoorSet));
  return coords;
};
void DesCoor(CoorSet *coor)
{
   free(coor);
};
float getf( float *floatvec, int index)
{
   return floatvec[index];
};
%}
// ==================================================================
//%module coor{
//struct coor {
//   //! vector of x coordinates
//   float        *xcoor;   // x coordinates
//   //! vector of y coordinates
//   float        *ycoor;   // y coordinates
//   //! vector of z coordinates
//   float        *zcoor;   // z coordinates
//}
//
%extend coorset {
   float x(int ii)
   {
     return self->xcoor[ii];
   };
   float y(int ii)
   {
     return self->ycoor[ii];
   };
   float z(int ii)
   {
     return self->zcoor[ii];
   };
   void setx(int ii, float value)
   {
     self->xcoor[ii] = value;
   }
   void sety(int ii, float value)
   {
     self->ycoor[ii] = value;
   }
   void setz(int ii, float value)
   {
     self->zcoor[ii] = value;
   }
   void setcoor(int ii, float xvalue, float yvalue, float zvalue)
   {
     self->xcoor[ii] = xvalue;
     self->ycoor[ii] = yvalue;
     self->zcoor[ii] = zvalue;
   }
   void copycoor(CoorSet *coorset2)
   {
     int ii;
     if( coorset2->nato != self->nato )
     {
       printf("Can't do: number of atoms doesn't match\n");
       return;
     }
     for( ii=0; ii<3*self->nato; ii++)
       self->cords[ii] = coorset2->cords[ii];
     return;
   }
   CoorSet *getselecoor( Selection *sele1 )
   {
     CoorSet *newcoorset;
     newcoorset = InitCoor(sele1->nselatm);
     GetSeleCoor( self, newcoorset, sele1 );
     return newcoorset;
   }
   
}
%extend traj {
  int nato()
  {
    return self->crdset->nato;
  }
}
%extend selection {
   int get_selatm(int ii)
   {
      return self->selatm[ii];
   }
   void set_selatm(int ii, int atm_index)
   {
      self->selatm[ii] = atm_index;
      return;
   }
   int SelAtm(int ii)
   {
      return self->selatm[ii];
   }
   void SetSelAtm(int ii, int atm_index)
   {
      self->selatm[ii] = atm_index;
      return;
   }
//   char SeleString(int ii)
}
%extend rawmol
{
   int Atomn( int index )
   {
     return self->atomn[index];
   };
   int Resn( int index )
   {
     return self->resn[index];
   };
   int Presn( int index )
   {
     return self->presn[index];
   };
   int Ipresn( int index )
   {
     return self->presn[index];
   };
   char * Restype( int index )
   {
     return self->restype[index];
   };
   char Chain( int index )
   {
     //return 'B';
     return self->chainId[index];
   };
   char * SegId( int index )
   {
     return self->segId[index];
   };
   float BFac( int index )
   {
     return self->bFac[index];
   };
   float Occup( int index )
   {
     return self->occup[index];
   };
   float Weight( int index )
   {
     return self->weight[index];
   };
/*   float ( int index )
   {
     return self->[index];
   };*/
}

// ===== MOLECULE helper functions =====
%exception MoleculeIter::next
{
  assert(!myErr);
  $action
  if (myErr)
  {
    myErr = 0; // clear flag for next time
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}
%inline
%{
  struct MoleculeIter
  {
    struct wrd_chain *ptr;
    size_t len;
  };
%}
%extend MoleculeIter
{
  struct MoleculeIter *__iter__()
  {
    return self;
  }

  struct wrd_chain * next()
  {
    if (self->len--)
    {
      return self->ptr++;
    }
    myErr = 1;
    return NULL;
  }
}

%extend _molecule
{
  struct MoleculeIter __iter__() {
    struct MoleculeIter ret = { self->chains, self->nChains };
    return ret;
  }
  int nAtomPerSegment(int ii)
  {
     return self->nApS[ii];
  };
  int nResPerSegment(int ii)
  {
     return self->nRpS[ii];
  };
  int nAtomPerRes(int ii)
  {
     return self->nApR[ii];
  };
  struct segm * Segment(int ii)
  {
     if( ii < self->nSeg )
       return &self->segment[ii];
     else
     {
       printf("Segnumber range is 0-%d; ", self->nSeg);
       printf("sending you back to segment 0\n");
       return &self->segment[0];
     }
  };
  struct res * get_pRes( int index )
  {
    return self->pRes[index];
  };
  void Info()
  {
     int ii;
     printf("Segments(%d):\n", self->nSeg);
     for( ii=0; ii<self->nSeg; ii++ )
     {
       printf("#0 (%s) - nResidues: %d\n", self->segment[ii].segName, self->segment[ii].nRpS);
     }
  }
  
  void read( char *filename )
  {
     GetMolecule( filename, self, 0 );
     return;
  }
  void write( char *filename )
  {
     WriteMolecule_unstr( self, filename, 0);
     return;
  }
  Selection * select( char *selestring )
  {
    Selection *sele1;
    sele1 = calloc( 1, sizeof(Selection));
    GetSele( selestring, sele1, self);
    return sele1;
  }
  Molecule * getselemol( Selection *sele1 )
  {
    Molecule *newmol;
    
    return MakeSeleMol( sele1, self);
  }
  float distance( int nato1, int nato2 )
  {
    return sqrt((self->coor.xcoor[nato1] - self->coor.xcoor[nato2])*(self->coor.xcoor[nato1] - self->coor.xcoor[nato2]) +
                (self->coor.ycoor[nato1] - self->coor.ycoor[nato2])*(self->coor.ycoor[nato1] - self->coor.ycoor[nato2]) +
                (self->coor.zcoor[nato1] - self->coor.zcoor[nato2])*(self->coor.zcoor[nato1] - self->coor.zcoor[nato2]));
  }
  //void np_copycoor(npdata)
  //{
    //int ii;
    //for( ii=0; ii<self->nato; ii++)
    //{
      //npdata[ii] = (float64) self->xcoor[ii];
      //npdata[ii] = (float64) self->ycoor[ii];
      //npdata[ii] = (float64) self->zcoor[ii];
    //}
  //}
  // additional helper functions in python:
  %pythoncode %{
  def get_coords(self, atmindex):
    xx = self.coor.x(atmindex)
    yy = self.coor.y(atmindex)
    zz = self.coor.z(atmindex)
    return xx, yy, zz
  
  def get_coords_np(self):
    import numpy
    pino = numpy.zeros((3,self.nato), dtype=numpy.float)
    for ii in range(self.nato):
      pino[0][ii] = self.coor.x(ii)
      pino[1][ii] = self.coor.x(ii)
      pino[2][ii] = self.coor.x(ii)
    return pino
  
  def sample( self ):
    return self.nato
  
  def selefromone( self, selection ):
    counter = 1
    reslist = []
    sele1 = self.select(selection)
    for ii in range(self.nRes):
      thisres = self.get_pRes(ii)
      for atmx in thisres:
        for jj in range(sele1.nselatm):
          if atmx.atomn == sele1.get_selatm(jj):
            reslist.append(counter)
            break
      counter += 1
    return reslist
  
  def describe(self, verbose=False):
    print '#============================#'
    print '# Number of Chains:', self.nChains
    print '# Number of Segments:', self.nSeg
    print '# Number of Residues:', self.nRes
    print '# Number of Atoms:', self.nato
    print '# Sequence:', self.seq1
    if verbose == True:
      for chainX in self:
        print '# Chain', chainX.chainId
        print '# n. of segments:', chainX.nSeg
        for segX in chainX:
          print '#  - Segment', segX.segName
          print '#  - n. of residues', segX.nRpS
          for resX in segX:
            print '#   - Residue', resX.resType
            print '#   - n. of atoms', resX.nApR
            for atomX in resX:
              print '#    - Atom', atomX.atmType
              if verbose == True:
                print '#      ', atomX.xCoor, atomX.yCoor, atomX.zCoor
    print '#============================#'
  
  def something(self):
    for chainX in self:
      for segX in chainX:
        for resX in segX:
          for atomX in resX:
            pass
  
  # end of helper functions
  %}  
}

// ===== CHAIN helper functions =====
%exception ChainIter::next
{
  assert(!myErr);
  $action
  if (myErr)
  {
    myErr = 0; // clear flag for next time
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}
%inline
%{
  struct ChainIter
  {
    struct segm *ptr;
    size_t len;
  };
%}
%extend ChainIter
{
  struct ChainIter *__iter__()
  {
    return self;
  }

  struct segm * next()
  {
    if (self->len--)
    {
      return self->ptr++;
    }
    myErr = 1;
    return NULL;
  }
}
%extend wrd_chain
{
  struct ChainIter __iter__() {
    struct ChainIter ret = { self->pSegm, self->nSeg };
    return ret;
  }
}
// ===== SEGMENT helper functions =====
%exception SegmentIter::next
{
  assert(!myErr);
  $action
  if (myErr)
  {
    myErr = 0; // clear flag for next time
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}
%inline
%{
  struct SegmentIter
  {
    struct res *ptr;
    size_t len;
  };
%}
%extend SegmentIter
{
  struct SegmentIter *__iter__()
  {
    return self;
  }

  struct res * next()
  {
    if (self->len--)
    {
      return self->ptr++;
    }
    myErr = 1;
    return NULL;
  }
}

%extend segm
{
  struct SegmentIter __iter__() {
    struct SegmentIter ret = { self->pRes, self->nRpS };
    return ret;
  }
  struct res * Res(int ii)
  {
    if( ii < self->nRpS )
      return &self->pRes[ii];
    else
    {
      printf("Resnumber range is 0-%d; ", self->nRpS);
      printf("sending you back to residue 0\n");
      return &self->pRes[0];
    }
  };
}
// -----------------
%exception ResIter::next
{
  assert(!myErr);
  $action
  if (myErr)
  {
    myErr = 0; // clear flag for next time
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}
%inline
%{
  struct ResIter
  {
    struct atom *ptr;
    size_t len;
  };
%}
%extend ResIter
{
  struct ResIter *__iter__()
  {
    return self;
  }

  struct atom * next()
  {
    if (self->len--)
    {
      return self->ptr++;
    }
    myErr = 1;
    return NULL;
  }
}
%extend res
{
  struct ResIter __iter__()
  {
    struct ResIter ret = { self->pAto, self->nApR };
    return ret;
  }
  struct atom * Atom(int ii)
  {
    if( ii < self->nApR )
      return &self->pAto[ii];
    else
    {
      printf("Atom-in-residue range is 0-%d; ", self->nApR);
      printf("sending you back to atom 0\n");
      return &self->pAto[0];
    }
  }
}
// from GMX vec.h //
/*
static inline void clear_mat(matrix a)
{
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=0.0;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=0.0;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=0.0;
}
*/
/////////////////////


typedef enum Moltype { OOO, PDB, CRD, COR, XYZ, GRO } moltype;       /*!< \brief enum for molecule file type (pdb/crd/txt/gro(?))*/
typedef enum Trjtype { ooo, dcd, xtc, txt, pdb, crd, cor } trjtype;  /*!< \brief enum for trajectory file type (dcd, xtc, txt, pdb, crd)*/
typedef enum Headertype {LITE, FULL } headertype;
/*typedef enum resType {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL} __restype;*/
/* ------------------------------------------------------------------ */
// ** Periodic Boundary Info
//! Periodic Boundary Structure
/*!
  Info about periodic boundary conditions - size, angles, CPT.
*/
typedef struct pbc
{
   //! Constant Pressure and Temperature flag (size may vary)
   int           cpt;
   //! PBC data: lenght along first axis
   double        a_size;
   //! PBC data: lenght along second axis
   double        b_size;
   //! PBC data: lenght along third axis
   double        c_size;
   //! PBC data: first angle
   double        angle1;
   //! PBC data: second angle
   double        angle2;
   //! PBC data: third angle
   double        angle3;
   //! crystal space group - Z value
   int           zval;
} Pbc   ;
/* ------------------------------------------------------------------ */
// ** Coordinates set
//! Coordinate Set Structure
/*!
  The coordinate set from a trj file, stored in three vectors(x, y and z)
*/
typedef struct coorset
{
   //! number of atoms
   int           nato;
   //! single array for all coordinates
   float        *cords;
   //! vector of x coordinates
   float        *xcoor;   // x coordinates (pointer to cords)
   //! vector of y coordinates
   float        *ycoor;   // y coordinates (pointer to cords+nato)
   //! vector of z coordinates
   float        *zcoor;   // z coordinates (pointer to cords+2*nato)
   //! pointers array for all coordinates (coors[0] = xcoor)
   float       **coors;
   //! PBC data: PBC flag
   int           pbc_flag;
   //! PBC data: PBC data structure
   Pbc          *pbc;
} CoorSet;
/* ------------------------------------------------------------------ */
// ** TRJ header - generic trajectory header
//! Trajectory (trj) info structure
/*!
  Structure with info about a trajectory, mostly collected from the
  headers in the case of dcd files.
*/
typedef struct trjh
{		
   char	      type[5];  /*!< \brief trj type can be CORD (coordinate trj) or VELO (velocities trj) */
   int         nframe;  /*!< \brief _real_ number of frames present in trajectories*/
   int         begstp;  /*!< \brief step number at beginning of trj*/
   int         skpstp;  /*!< \brief number of steps between coordinate frames*/
   unsigned int numstp; /*!< \brief total _real_ number of steps*/
   float        tmstp;  /*!< \brief time between steps (time between frames = tmstp*skpstp)*/
   int         varpbc;  /*!< \brief whether (0/1) variable Periodic Boundary Conditions are employed*/


   int           nato;  /*!< \brief number of atoms in the trj*/
   int         fixatm;  /*!< \brief number of (optionallly) fixed atoms in the trj*/
   int   *fixatm_list;  /*!< \brief list of (optionallly) fixed atoms in the trj*/
   CoorSet   *fixcoor;  /*!< \brief coordinates of (optionallly) fixed atoms in the trj*/

   int            sex;  /*!< \brief whether byte order of trj file is the same as the running machine (0=y, 1=n)*/
   int           size;  /*!< \brief size of headers: used to compute offset for random-access reading of frames*/
   int      sixtyfour;  /*!< \brief whether the file was made by a 64-bit machine*/

   
   int         comnum;  /*!< \brief CHARMM only - number of 80-bytes long comment string*/
   char  comm[12][91];  /*!< \brief CHARMM only - comment strings*/

   int       c_nframe;  /*!< \brief CHARMM only - _claimed_ number of frames; incorrect if trj is not complete*/
   int       c_numstp;  /*!< \brief CHARMM only - _claimed_ number of steps; incorrect if trj is not complete*/

   int	    extra1;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
   int	    extra2;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
   int	    extra3;     /*!< \brief CHARMM only - (unused) extra field in the header section*/
} Trjh;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** general TRJ infos
//! Overall Trajectory (trj) structure
/*!
  Structure with overall info about a trajectory, pointers to the file,
  to the header info,  trjtype, etc.
*/
typedef struct traj
{
   //! pointer to the file (if open)
   FILE        *file;
   //! the file name
   char         filename[1024];
   //! whether it's open or not
   short        open;
   //! if it is an xtc (gromacs) file, this points to the file
   XDRFILE *    xtc;
   //! type (dcd/xtc) of the trj file
   trjtype      filetype;
   //! pointer to hdr structure
   struct trjh *hdr;
   //! the detail level at which info from an xtc file has to be gained
   headertype   hdetail;
   //! while reading a trj, this is the last frame read
   int          previousframe;
   //! a container for a coordinate set extracted from the trajectory
   CoorSet     *crdset;
} Traj	;
/* ------------------------------------------------------------------ */
// ** raw (or row) molecule structure
//! Molecule file (pdb/crd) file row-style structure
/*!
  The info from a molecule file (crd/pdb) are gather here in an
  unstructured form, ie in arrays of variables, each relative to a 
  in a row in the file. The info from each atom is thus contained
  across all the arrays at the relevant index. Easy and fast but crude
  access to data. To enable selections and structured access, this data
  is usually processed and organized in a Molecule structure
*/
typedef struct rawmol
{
   int          *presn;  /*!< presumed residue number (for CRD) */
   int        b_cryst;   /*!< bool: if have crystal unit */
   float            a;   /*!< crystal box element - lenght #1 */
   float            b;   /*!< crystal box element - lenght #2 */
   float            c;   /*!< crystal box element - lenght #3 */
   float        alpha;   /*!< crystal box element - angle #1 */
   float         beta;   /*!< crystal box element - angle #2 */
   float        gamma;   /*!< crystal box element - angle #3 */
   char   s_group[11];   /*!< crystal space group - space group */
   int              z;   /*!< crystal space group - Z value */

   int           nato;   /*!< number of atoms */
   char         *name;   /*!< molecule name */
   char        **atom;   /*!< first column label (in a pdb, usually ATOM) */
   int         *atomn;   /*!< atom number (similar to atmId, which is a char) */
   char        *chain;   /*!< chain name (a single character) */
   char      *chainId;   /*!< chain ID (PDB only) */
   char       *altLoc;   /*!< altLoc field (PDB only) */
   char        *iCode;   /*!< iCode field (PDB only) */
   float        *bFac;   /*!< beta factor (temperature) (PDb only) */
   float       *occup;   /*!< occupancy (PDb only) */
   int          *resn;   /*!< residue number */
   float        *weight; /*!< weight (CRD only) */
   
   short      *hetatm;   /*!< whether it is an hetatm or not */
   char       **atmId;   /*!< atom ID */
   char       **segId;   /*!< segment ID */
   char     **element;   /*!< element (PDB only) */
   char      **charge;   /*!< charge (PDB only) */
   char     **atmtype;   /*!< atom type */
   char     **restype;   /*!< residue type (amino acid) */
   char        **tail;   /*!< eventually, what's left after the known fields */
   short      *segend;   /*!< whether it is the last atom in a segment */
   short      *segbeg;   /*!< whether it is the first atom in a segment */

}  RawMol;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Atom
//! Single Atom structure
/*!
  All info regarding a single atom are stored in an instance of this
  structure. ato instances are tipically grouped in res instances.
  ato are allocated and filled by processing data in _pdbs instances
*/
typedef struct atom
{
    int        atomn;      // atom number
    char     atmType[5];   // atom type - actually char[4], but with 5 works

    float       xCoor;      // x coordinates
    float       yCoor;      // y coordinates
    float       zCoor;      // z coordinates

    short      hetatm;     // whether it is an hetatm or not
    char       atmId[6];   // atom ID number (char) new from PDB.org!!!
    char       segId[5];   // segment ID (char)     new from PDB.org!!!
    char     element[2];   // new from PDB.org!!!
    char      charge[2];   // new from PDB.org!!!
    char     chainId;      // new from PDB.org!!!
    char      altLoc;      // new from PDB.org!!!
    char       iCode;      // new from PDB.org!!!
    float       bFac;      // beta Factor, usable as index
    float      occup;      // occupancy

    float     weight;      // added for CRD
    int        presn;      // added for CRD
    int           bb;      // backbone (1) or side chain (0) - not widely used
} Atom;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Residue
//! Single Residue structure
/*!
  info regarding a single residue, like residue type, residue number,
  the number of atoms belonging to it and pointers to the relative 
  _ato instances. res are grouped in segm instances and are allocated
  and filled by processing data in _pdbs instances
*/
typedef struct res
{
    int           nApR;         // number of Atoms Per Residue (nAPR)
    char       resType[5];      // residue type
    int           resn;         // residue number
    int          presn;         // progressive residue number
    Atom         *pAto;         // Atoms 
} Res;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Segment
//! Single Segment structure
/*!
  info regarding a single segment, like segment name, number of atoms
  belonging to it, number of residues belonging to it and pointers to
  the relative _res instances
*/
typedef struct segm
{
    int           nApS;         // number of Atoms (per Segment)
    char       segName[5];      // segment name
    int           nRpS;         // n Residues per Segment (nRpS) 
    Res          *pRes;         // Residues 
} Segm;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Chain
//! Single Chain structure
/*!
  info regarding a single chain, like chain id, number of atoms
  belonging to it, number of segments belonging to it, number of 
  residues belonging to it and pointers to the relative Segm instances
*/
typedef struct wrd_chain
{
    int           nApC;         // number of Atoms (per Chain)
    int           nRpC;         // number of Atoms (per Chain)
    int           nSeg;         // n Segment per Chain (nRpS) 
    char       chainId;         // chain Id
    Segm        *pSegm;         // Segments 
} wrd_Chain;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Model
//! Single Model structure
/*!
  info regarding a single model, actually a whole molecule
  Number of atoms (total, per residue, per chain and per segment),
  number of residues (total, per chain and per segment), number
  of segments, number of chains and pointers to appropriate wrd_Chain
  instances, coordinates and a _coor instance and a _pdbs instance with 
  original data. Also a note about the origin and if it comes from an 
  extended format.
  * NOT USED
*/
typedef struct wrd_model
{
    int           nApM;         // number of Atoms (per Model)
    int        modelId;         // model Id (number)
    int           nCpM;         // n Chains per Model (nCpM) 
    wrd_Chain  *pChain;         // Chains 
} wrd_Model;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Overall Molecular Structure File
//! Single Molecule structure
/*!
  info regarding a molecule. Number of atoms (total, per residue and
  per segment), number of residues (total and per segment), number
  of segments and pointers to appropriate _segm instances, coordinates ad
  a _coor instance and a _pdbs instance with original data. Also a
  note about the origin and if comes from an extended format.
*/
typedef struct _molecule
{
    int             nato;       // number of atoms
    int          nModels;       // number of models
    int          nChains;       // number of chains
    int             nSeg;       // number of segments
    int             nRes;       // number of residues
    int            *nApC;       // number of Atoms per Chain
    int            *nRpC;       // number of Residues per Chain
    int            *nSpC;       // number of Segments per Chain
    int            *nApS;       // number of Atoms per Segment
    int            *nRpS;       // number of Residues per Segment
    int            *nApR;       // number of Atoms per Residue in Molecule - not filled

    wrd_Chain    *chains;       // array of chains
    Segm        *segment;       // array of segments
    Res           **pRes;       // array of residues (pointer to residues inside segments)
    CoorSet        coor;        // coordinates
    CoorSet        coor2;       // additional coordinates (see altLoc)
    CoorSet        coor3;       // additional coordinates (see altLoc)
    CoorSet       *cords;       // used if more than 1 model is present
    RawMol        rawmol;       // raw molecule structural data
    char            *seq1;      // sequence in 1-letter format
    char           **seq3;      // sequence in 3-letter format
    
    int           origin;       // 0 if from pdb, 1 if from crd, 2 if filled for both;
    int           ext_flag;     // whether a crd is EXT(ended)
    
    char         *filled;       // whether the structure has already been used (and thus allocated)
    
    short        altLocB;       // whether a second, alternative set of coordinates is present
    short        altLocC;       // whether a third, alternative set of coordinates is present
   
   /*taken as is from rawmol structure - better have it handy */
   int        b_cryst;   /*!< bool: if have crystal unit */
   float            a;   /*!< crystal box element */
   float            b;   /*!< crystal box element */
   float            c;   /*!< crystal box element */
   float        alpha;   /*!< crystal box element */
   float         beta;   /*!< crystal box element */
   float        gamma;   /*!< crystal box element */
   char   s_group[11];   /*!< crystal space group */
   int              z;   /*!< crystal Z */
} Molecule;
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// ** Selection Structure
//! Atom Selection structure
/*!
  info about an atom selection, such as the selection string which is
  used to scan the molecule file and the arrays where the atom numbers
  are stored once selected. Selestring can also be a file name, with the
  file being a list of atom numbers
*/
typedef struct selection
{
    int   ii;
    //! number of selected atoms
    int   nselatm;
    //! selected atoms' id numbers
    int  *selatm;
    //! selection string
    char  selestring[20480];
} Selection;
/* ------------------------------------------------------------------ */
// ** Frames to extract structure
//! Frame List Structure
/*!
  a list of frames, usually taken from a file
*/
typedef struct intlist
{
   int	   nframes;   // number of frames to process
   int	   *nframe;   // ids of frame to process
} IntList;
/* ------------------------------------------------------------------ */
// ** List of Trajectory files
//! List of Trajectory files
/*!
  a list of trajectory files, usually taken from a file
*/
typedef struct trjlist {
   //! number of trajectories to handle
   int	      ntrj;
   //! names of trajectories to handle
   char	 **trjname;   
} TrjList;
/* ------------------------------------------------------------------ */
// ** List of radii data
//! List of radii data
/*!
  a list of trajectory files, usually taken from a file
*/
typedef struct radlist {
   //! number of radii data
   int	        nradii;
   //! 
   char	      **resname;
   char       **atomname;
   float       *radius;
} RadList;
/* ------------------------------------------------------------------ */
//! prints info about trajectory file
/*!
  prints to stdout (usually terminal) basic info about a trajectory file
  mostly taken from the headers. Complete support for dcd files, more
  lacking for xtc files. Info have to be already been gathered.
*/
void ShowTrjHeader     (const Traj *traj) ;
/* ------------------------------------------------------------------ */
//! copies the content of a dcd header
/*!
 copies the content of the source trajectory header structure to the
 target one.
*/
void HdrCp ( struct trjh *source,   /*!< source header */
             struct trjh *target    /*!< target header */
           ) ;
/* ------------------------------------------------------------------ */
//! reads a molecule
/*!
 general molecule-reading function - calls format-specific functions to
 read in different file formats
*/
Molecule * ReadMolecule ( char          * name, /*!< molecule file name */
                          moltype       filetype/*!< molecule file type */
                        ) ;
/* ------------------------------------------------------------------ */
//! reads a molecule
/*!
 general molecule-reading function - calls format-specific functions to
 read in different file formats
*/
void GetMolecule        ( char          *name,  /*!< molecule file name */
                          Molecule      *molecule ,/*!< molecule structure */
                          moltype        filetype /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! writes a molecule
/*!
 general molecule-writing function - calls format-specific functions to
 write in different file formats

*/
void WriteMolecule      ( Molecule      *molecule,      /*!< molecule data structure */
                          char          *name,          /*!< molecule file name */
                          moltype        filetype       /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! writes a molecule
/*!
 general molecule-writing function - calls format-specific functions to
 write in different file formats. Uses rawmol (order as original file)

*/
void WriteMolecule_unstr( Molecule      *molecule,      /*!< molecule data structure */
                          char          *name,          /*!< molecule file name */
                          moltype        filetype       /*!< molecule file type */
                        );
/* ------------------------------------------------------------------ */
//! updates data inside molecule structure
/*!
 updates data inside molecule structure to complete the informations 
 about file origin
*/
void UpDateMol          ( Molecule      *molecule       /*!< molecule data structure */);  
/* ------------------------------------------------------------------ */
//! Allocate memory for a molecule
/*!
*/
void AllocMol           ( int           nato,           /*!< number of atoms */
                          Molecule      *molecule       /*! molecule data structure */
);
/* ------------------------------------------------------------------ */
//! Create a Molecule from a selection of another molecule
/*!
*/
Molecule * MakeSeleMol    ( Selection *selstr,  /*!< Selection structure to be used (input)*/
                            Molecule *inmol     /*!< Parent molecule (input) */
) ;
//void MakeSeleMol        ( Selection *selstr,  /*!< Selection structure to be used (input)*/
//                          Molecule *inmol,    /*!< Parent molecule (input) */
//                          Molecule **selemol  /*!< New Molecule made from selection (output) */
//) ;
/* ------------------------------------------------------------------ */
//! clean molecule
/*!
 frees all allocated memory inside a molecule structure
*/
void CleanMolecule      ( Molecule      *molecule );
/* ------------------------------------------------------------------ */
//! deletes molecule
/*!
 frees all allocated memory inside a molecule structure and frees molecule pointer
*/
void DelMolecule        ( Molecule      *molecule );
/* ------------------------------------------------------------------ */
//! opens a trajectory file for reading/writing
/*!
 opens the stream, allocates memory inside the traj structure (for trjh)
 annotates filename
*/
int OpenTrj            ( char  *filename,
                         Traj *traj,
                         char  *mode
                        ) ;
/* ------------------------------------------------------------------ */
//! opens a trajectory file for reading/writing
/*!
 opens the stream, allocates memory inside the traj structure (for trjh)
 annotates filename
*/
Traj * InitTrj          ( char  *filename,
                          char  *mode
                        ) ;
/* ------------------------------------------------------------------ */
//! closes a trajectory file
/*!
 free allocated memory (for trjh), erases filename and closes the stream
*/
void CloseTrj           ( Traj  *traj ) ;
/* ------------------------------------------------------------------ */
//! reads info about a trajectory file
/*!
 
*/
void ReadTrjHeader      ( Traj  *traj);
/* ------------------------------------------------------------------ */
//! reads a coordinate set from a trajectory file
/*!
 general header-reading functions; calls format-specific functions to 
 gather info about the trajectory
*/
void ReadTrjCoor        ( Traj          *traj,  /*!< trajectory to read */
                          CoorSet       *crd,   /*!< coordinates (to be filled) */
                          int            frn    /*!< frame to get */
                        );
/* ------------------------------------------------------------------ */
//! copy headers to a (new) trajectory header structure
/*!
 general trajectory headers copying function: calls format-specific
 functions to copy headers
*/
void CopyTrjHeader ( Trjh *source,         // source header
                     Trjh *target  );      // target header
/* ------------------------------------------------------------------ */
//! write trajectory headers to a trajectory file
/*!
 general trajectory headers writing function: calls format-specific
 functions to write headers to a trajectory file
*/
void WriteTrjHeader     ( Traj  *traj,
                          char  *field );
/* ------------------------------------------------------------------ */
//! write a coordinate set to a trajectory file
/*!
 general trajectory coordinates set writing function: calls
 format-specific functions to write a coordinates set to a trajectory
 file
*/
void WriteTrjCoor ( CoorSet     *crd,   /*!< coordinate set*/
                    Traj        *traj   /*!< trajectory file structure*/
                   );
/* ------------------------------------------------------------------ */
//! reads a list of trajectories
/*!
 takes a text file (.txt) and reads a list of trajectories to be handled
 by the program (different procedures can use this)
*/
void ReadTrjList        ( char          *,      /*!< list file */
			  TrjList        *       /*!< list of trj file structure */
		 	) ;
/* ------------------------------------------------------------------ */
//! extracts selected atoms' id numbers
/*!
 high-level selection function: 
*/
int GetSele             ( char          *,      /*!< selection string */
                          Selection     *,      /*!< selection structure */
                          Molecule      *       /*!< molecule file */
                        );
/* ------------------------------------------------------------------ */
//! extracts selected atoms' id numbers
/*!
 selection function: 
*/
int GetSeleAtoms        ( char          *,      /*!< selection string */
                          Selection     *,      /*!< selection structure */
                          Molecule      *       /*!< molecule file */
                        );
/* ------------------------------------------------------------------ */
//! reads a radius list file
/*!
 radius reading function: 
*/
RadList * ReadRadius    (char *filename);
/* ------------------------------------------------------------------ */
//! assign radius based on atomname and residue type
/*!
 radius assigning function: 
*/
int RadAssign           ( RadList *radlist,
                          Molecule *molecule, 
                          double *pdAtomsRadii
                        );

//==============================================================================
//   MEMORY ALLOCATION
//=============================================================================
/* ------------------------------------------------------------------ */
//! creates a coorset and calls CalloCoor
CoorSet * InitCoor          ( int nato);
// -----------------------------------------------------------------
//! frees memory allocated inside a coorset
void CleanCoor           ( CoorSet *);
// =============================================================================
//! copies data from a buffer to another
/*< 
 copies data from buffer1 to buffer2 ; from start position to end position ;
 buffer2 is csize long
*/
int copyvar             ( char *buffer1, int start, int end, char *buffer2, int csize );
/* ------------------------------------------------------------------ */
void GetSeleCoor ( CoorSet *source, CoorSet *target, Selection *sele );
void CopyCoor ( CoorSet *source, CoorSet *target );
/* ------------------------------------------------------------------ */

