// ------------------------------------------------------------------
// ** Sample Computing Module data Structure
struct inp_sample
{
  char           title[64];
  Selection       sele;
  short          someflag;
  short          extout;
  float        **reference;
  float        **moving;
};
// ------------------------------------------------------------------
int Read_iSample(char **input, int inp_index, struct inp_sample * , char * , Molecule * );
// ------------------------------------------------------------------
int Compute_Sample (struct inp_sample *, struct sopt *, CoorSet *, char *  );
//----------------------------------------------------------
int Post_Sample ( struct inp_sample *inp_sample, struct sopt *OPT, int nframe, Molecule *molecule );
