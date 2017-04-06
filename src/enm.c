#include "wordom.h"
#include "tools.h"
#include "fileio.h"
#include "datahandler.h"
#include "kga.h"
#include "enm.h"
#include "volumes.h"

//
# define k 0.001966667
//# define kT 0.59
# define pi 3.14159
# define delta 0.1 
//
/* This code is public domain -- Will Hartung 4/9/09 */
size_t lgetline(char **lineptr, size_t *n, FILE *stream)
{
    char *bufptr = NULL;
    char *p = bufptr;
    size_t size;
    int c;

    if (lineptr == NULL) {
        return -1;
    }
    if (stream == NULL) {
        return -1;
    }
    if (n == NULL) {
        return -1;
    }
    bufptr = *lineptr;
    size = *n;

    c = fgetc(stream);
    if (c == EOF) {
        return -1;
    }
    if (bufptr == NULL) {
        bufptr = malloc(128);
        if (bufptr == NULL) {
                return -1;
        }
        size = 128;
    }
    p = bufptr;
    while(c != EOF) {
        if ((p - bufptr) > (size - 1)) {
                size = size + 128;
                bufptr = realloc(bufptr, size);
                if (bufptr == NULL) {
                        return -1;
                }
        }
        *p++ = c;
        if (c == '\n') {
                break;
        }
        c = fgetc(stream);
    }

    *p++ = '\0';
    *lineptr = bufptr;
    *n = size;

    return p - bufptr - 1;
}
// -------
int getSeleIndex( int element, Selection *sele)
{
  int       ii;
  
  for( ii=0; ii<sele->nselatm; ii++ )
    if( sele->selatm[ii] == element )
    {
      //printf( "DEBUG element %d has index %d\n", element, ii );
      return ii;
    }
  
  return -1;
}
// -------
int Read_iEnm( char **input, int inp_index, struct inp_enm *inp_enm, char *printout, Molecule *molecule )
{
   FILE         *vecmatrix;
   char         buffer[20000];
   char         title[1024];
   char         matrix[1024], *pch;
   int          gotit, ii=0, jj=0;
      
   extern short int     no_frame_par;
   no_frame_par = 1;
   
   inp_enm->sele1_flag = 0;
   inp_enm->mol2_flag = 0;
   inp_enm->sele2_flag = 0;
   inp_enm->vecfile2_flag=0;
   inp_enm->beta_flag = 0;
   inp_enm->perturb_flag = 0;
   inp_enm->corr_flag = 0;
   inp_enm->print=0;
   inp_enm->matrix_chk=0; //ruvido
   inp_enm->network_print=0; //ruvido
   inp_enm->vsa_flag=0;
   inp_enm->proj_flag=0;
   inp_enm->dcd_flag=0;
   inp_enm->rtb_flag=0;
   inp_enm->defen_flag=0;
   inp_enm->mass_flag=0;
   inp_enm->rtblevan_flag=0;
   inp_enm->enm_level_flag=0;
   inp_enm->distmat_flag=0;
   inp_enm->pertmat_flag=0;
     
   memset ( title, '\0', sizeof(title));
   inp_enm->molecule = molecule;
   while( strncmp (buffer, "END", 3))
   {
     gotit = 0;
     memset ( buffer, '\0', sizeof(buffer));
     sprintf( buffer, "%s", input[inp_index]);
     if( !strncmp(buffer, "BEGIN", 5) || !strncmp(buffer, "END", 3) || buffer[0] == '#')
       gotit = 1;
     else if ( !strncmp(buffer, "--TITLE", 7))
     {
       sscanf( buffer, "--TITLE  %s", title);
       sscanf( title, "%s", inp_enm->title);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--INTYPE", 8))
     {
       sscanf(buffer, "--INTYPE %s ", inp_enm->type);
       if( strncmp(inp_enm->type,"linear", 64) && strncmp(inp_enm->type,"kovacs", 64))
       {
         fprintf( stderr, "--INTYPE valid options are \"linear\" or \" kovacs\"\n");
         exit(0);
       }
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CUTOFF", 8))
     {
       sscanf(buffer, "--CUTOFF %f ", &inp_enm->cutoff);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE ", 7))
     {
       sscanf(buffer, "--SELE %[^\n]%*c ", inp_enm->sele1.selestring);
       inp_enm->sele1_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--SELE2 ", 8))
     {
       sscanf(buffer, "--SELE2 %[^\n]%*c ", inp_enm->sele2.selestring);
       inp_enm->sele2_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--VSA", 5))
     {
       sscanf(buffer, "--VSA %[^\n]%*c ", inp_enm->vsasele.selestring);
       inp_enm->vsa_flag = 1;
       gotit = 1;
     }
     //else if ( !strncmp(buffer, "--LEVEL ", 8 ) && inp_enm->rtb_flag == 0)
     //{
     //  sscanf (buffer, "--LEVEL %s", inp_enm->enm_level);
     //  if (strncmp (inp_enm->enm_level, "RES", 64)==0)
     //  {
     //    inp_enm->enm_level_flag=1;
     //  }
     //  gotit = 1;
     //}     
     else if ( !strncmp(buffer, "--RTB ", 6 ))
     {
       sscanf (buffer, "--RTB %s", inp_enm->rtb_method);
       inp_enm->rtb_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--RTBLEVAN", 10) && inp_enm->rtb_flag == 1)
     {
       sscanf(buffer, "--RTBLEVAN %[^\n]%*c ", inp_enm->rtblevan.selestring);
       inp_enm->rtblevan_flag=1;
       gotit = 1; 
     }      
     else if ( !strncmp(buffer, "--MOL2", 6))
     {
       sscanf(buffer, "--MOL2 %[^\n]%*c ", inp_enm->mol2);
       inp_enm->mol2_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--VECFILE2", 10))
     {
       sscanf(buffer, "--VECFILE2 %[^\n]%*c ", inp_enm->vecfile2);
       inp_enm->vecfile2_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NMODES", 8))
     {
       inp_enm->nmodes = calloc( 256, sizeof(char));
       sscanf(buffer, "--NMODES %s", inp_enm->nmodes);
       gotit = 1;
     }          
     else if ( !strncmp(buffer, "--BETA", 6))
     {
       inp_enm->beta = calloc( 64, sizeof(char));
       sscanf (buffer, "--BETA %s", inp_enm->beta);
       inp_enm->beta_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--PERTURB", 9))
     {
       inp_enm->perturb = calloc( 64, sizeof(char));
       sscanf (buffer, "--PERTURB %s", inp_enm->perturb);
       inp_enm->perturb_norm = 0;
       inp_enm->perturb_flag = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NORM_PERTURB", 14) && inp_enm->perturb_flag == 1)
     {
       inp_enm->perturb_norm = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--PERT_MAT", 10) && inp_enm->perturb_flag == 1)
     {
       inp_enm->pertmat_flag = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--DEFEN ", 8))
     {
       inp_enm->defen = calloc( 64, sizeof(char));
       sscanf (buffer, "--DEFEN  %s", inp_enm->defen);
       inp_enm->defen_flag=1;
       gotit = 1;
     }           
     else if ( !strncmp(buffer, "--DIST_MAT", 10) && inp_enm->defen_flag == 1)
     {
       inp_enm->distmat_flag = 1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--CORREL", 8))
     {
       inp_enm->corr = calloc( 64, sizeof(char));
       sscanf (buffer, "--CORREL %s", inp_enm->corr);
       inp_enm->corr_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--PRINT", 7))
     {
       inp_enm->print=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--MATRIX", 7)) //ruvido
     {
       sscanf ( buffer, "--MATRIX  %s", matrix);
       //sprintf( printout, " %8s ", title);     ???? ruvido
       sscanf ( matrix, "%s", inp_enm->matrix);
       inp_enm->matrix_chk=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--NETWORK_PRINT", 7)) //ruvido
     {
       inp_enm->network_print=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--PROJ", 6))
     {
       inp_enm->proj = calloc( 64, sizeof(char));
       sscanf (buffer, "--PROJ %s", inp_enm->proj);
       inp_enm->proj_flag=1;
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--DCDOUT", 8))
     {
       inp_enm->dcd_flag = 1;
       sscanf( buffer, "--DCDOUT %[^\n]%*c ", inp_enm->projdcd);
       if ( inp_enm->proj_flag == 0 )
       { 
         printf("--DCDOUT option needs --PROJ\n");
         exit(0);  
       }
       gotit = 1;
     }
     else if ( (!strncmp(buffer, "--TEMP", 6) && inp_enm->beta_flag == 1) || (!strncmp(buffer, "--TEMP", 6) && inp_enm->proj_flag == 1))
     {
       sscanf (buffer, "--TEMP %s", inp_enm->temp);
       gotit = 1;
     }
     else if ( !strncmp(buffer, "--MASS", 6))
     {
       inp_enm->mass_flag=1;
       gotit = 1;
     }
     if( gotit==0 )
     {
       fprintf( stderr, "Sorry, could not understand option: %s\n", buffer);
       exit(5);
     }
     inp_index++;
   }
 

    GetSele ( inp_enm->sele1.selestring, &inp_enm->sele1, molecule);
    inp_enm->selmol = MakeSeleMol (&inp_enm->sele1, molecule);
   
   /*Allocating memory for system ENM arrays*/     
   /*Obtaining internal numeration for the system*/   
   inp_enm->sys_list = calloc (inp_enm->sele1.nselatm, sizeof(int*));
   for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
   {
     inp_enm->sys_list[ii] = calloc(2, sizeof(int));     
   }   
   
   for (ii=0;ii<inp_enm->sele1.nselatm; ii++)
   {
     inp_enm->sys_list[ii][0] = ii;
     inp_enm->sys_list[ii][1] = inp_enm->sele1.selatm[ii]; 
   }    
   
   //If a selections contains more atoms of a residue, the geometrical centre is considered for the computation

   inp_enm->coord = calloc (inp_enm->sele1.nselatm, sizeof(float*));
   for(ii=0;ii<inp_enm->sele1.nselatm;ii++)
   {
     inp_enm->coord[ii] = calloc(3, sizeof(float));
   }
   
     
   inp_enm->dr = calloc(inp_enm->sele1.nselatm*inp_enm->sele1.nselatm, sizeof( _dr));
   
   inp_enm->He = calloc( 3*inp_enm->sele1.nselatm, sizeof(float *) );
   inp_enm->He[0] = calloc( 9*inp_enm->sele1.nselatm*inp_enm->sele1.nselatm, sizeof(float) );
   for ( ii=0; ii < 3*inp_enm->sele1.nselatm ; ii++)
   {
     inp_enm->He[ii] = inp_enm->He[0] + 3*inp_enm->sele1.nselatm*ii;
   }
   
   inp_enm->eigen = calloc ( 1, sizeof( _eigen));
   
   if(inp_enm->rtb_flag == 1)
   {
     inp_enm->rtb2full_eigen = calloc ( 1, sizeof( _eigen));
   }

 
   if(inp_enm->vsa_flag == 1)
   {
     /*Getting subselections*/

     GetSele ( inp_enm->vsasele.selestring, &inp_enm->vsasele, molecule);
     
     //memset(selestring,'\0',sizeof(selestring));
     //sprintf(selestring,"%s; del %s",inp_enm->sele1.selestring, inp_enm->vsasele.selestring);
     //GetSele ( selestring, &inp_enm->env, molecule);
     SeleDiff( &inp_enm->sele1, &inp_enm->vsasele, &inp_enm->env );
     
     printf("\nHe=%d, Hss=%d, Hee=%d \n",inp_enm->sele1.nselatm, inp_enm->vsasele.nselatm, inp_enm->env.nselatm);
     
     /*Obtaining internal numeration for the subsystem and the environment*/
         
     inp_enm->subsys_list = calloc (inp_enm->vsasele.nselatm, sizeof(int*));
     for(ii=0;ii<inp_enm->vsasele.nselatm;ii++)
     {
       inp_enm->subsys_list[ii] = calloc(2, sizeof(int));    
     }
    
     inp_enm->env_list = calloc (inp_enm->env.nselatm, sizeof(int*));
     for(ii=0;ii<inp_enm->env.nselatm;ii++)
     {
       inp_enm->env_list[ii] = calloc(2, sizeof(int));    
     }
     
     for (ii=0;ii<inp_enm->vsasele.nselatm;ii++)
     {
       inp_enm->subsys_list[ii][1] = inp_enm->vsasele.selatm[ii];     
     }
     for (ii=0;ii<inp_enm->env.nselatm;ii++)
     {
       inp_enm->env_list[ii][1] = inp_enm->env.selatm[ii];     
     }
    	 
     for (ii=0;ii<inp_enm->sele1.nselatm;ii++)
     {
       for (jj=0;jj<inp_enm->vsasele.nselatm;jj++)
       {
         if(inp_enm->sys_list[ii][1] == inp_enm->subsys_list[jj][1])
         {
           inp_enm->subsys_list[jj][0] = inp_enm->sys_list[ii][0];
         }
       }
       for (jj=0;jj<inp_enm->env.nselatm;jj++)
       {
         if(inp_enm->sys_list[ii][1] == inp_enm->env_list[jj][1])
         {
           inp_enm->env_list[jj][0] = inp_enm->sys_list[ii][0];
         }
       }
     }    
     
     if (inp_enm->vsasele.nselatm > inp_enm->sele1.nselatm)
     {
       printf ("Error: --VSA selected atoms must be contained in --SELE selection!\n");
       exit(0);
     }

        
     inp_enm->Hss = calloc( 3*inp_enm->vsasele.nselatm, sizeof(float *) );
     inp_enm->Hss[0] = calloc( 9*inp_enm->vsasele.nselatm*inp_enm->vsasele.nselatm, sizeof(float) );
     for ( ii=0; ii < 3*inp_enm->vsasele.nselatm ; ii++)
     {
       inp_enm->Hss[ii] = inp_enm->Hss[0] + 3*inp_enm->vsasele.nselatm*ii;
     }
     
     if ( inp_enm->env.nselatm > 0)
     {
       inp_enm->Hee = calloc( 3*inp_enm->env.nselatm, sizeof(float *) );
       inp_enm->Hee[0] = calloc( 9*inp_enm->env.nselatm*inp_enm->env.nselatm, sizeof(float) );
       for ( ii=0; ii < 3*inp_enm->env.nselatm ; ii++)
       {
         inp_enm->Hee[ii] = inp_enm->Hee[0] + 3*inp_enm->env.nselatm*ii;
       }
    
       inp_enm->invHee = calloc( 3*inp_enm->env.nselatm, sizeof(float *) );
       inp_enm->invHee[0] = calloc( 9*inp_enm->env.nselatm*inp_enm->env.nselatm, sizeof(float) );
       for ( ii=0; ii < 3*inp_enm->env.nselatm ; ii++)
       {
         inp_enm->invHee[ii] = inp_enm->invHee[0] + 3*inp_enm->env.nselatm*ii;
       }
    
       inp_enm->Hes = calloc( 3*inp_enm->env.nselatm, sizeof(float *) );
       for ( ii=0; ii < 3*inp_enm->env.nselatm ; ii++)
       {
         inp_enm->Hes[ii] = calloc (3*inp_enm->vsasele.nselatm, sizeof (float));
       }
       inp_enm->subeigen = calloc( 1, sizeof( _eigen));
     }
     
     inp_enm->drvsa = calloc(inp_enm->vsasele.nselatm*inp_enm->vsasele.nselatm, sizeof( _dr));
   }
 
   
   if(inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
   {   
     Enm_Getblock(inp_enm, molecule);
     if (inp_enm->rtblevan_flag==0)
     {
       inp_enm->rtblevan = inp_enm->sele1;
     }
     else if (inp_enm->rtblevan_flag==1)
     {
       //GetSele ( inp_enm->rtblevan.selestring, &inp_enm->rtblevan, inp_enm->selmol);
       // rtblevan may be expressed as atom numbers on original molecule, selmol might be smaller
       // following code accounts for all rtblevan while preserving mapping on selmol
      //GetSele ( inp_enm->rtblevan.selestring, &inp_enm->rtblevan, molecule ); //inp_enm->selmol);
      GetSele ( inp_enm->rtblevan.selestring, &inp_enm->rtblevan_org, molecule ); //inp_enm->selmol);
      printf("Level of Analysis %d\n", inp_enm->rtblevan_org.nselatm );
      //fprintf( stdout, "DEBUG rtblevan.nselatm: %d\n", inp_enm->rtblevan.nselatm );
      inp_enm->rtblevan.nselatm = inp_enm->rtblevan_org.nselatm;
      inp_enm->rtblevan.selatm = (int *)calloc( inp_enm->rtblevan.nselatm, sizeof( int) );
      for( ii=0; ii<inp_enm->rtblevan_org.nselatm; ii++ )
      {
        //fprintf( stdout, "DEBUG %d -> %d\n", inp_enm->rtblevan.selatm[ii], getSeleIndex( inp_enm->rtblevan.selatm[ii], &inp_enm->sele1 ) );
        inp_enm->rtblevan.selatm[ii] = getSeleIndex( inp_enm->rtblevan_org.selatm[ii], &inp_enm->sele1 )+1;
      }
      //for(ii=0;ii<inp_enm->sele1.nselatm;ii++)
      //{
      //  printf("%d\n",inp_enm->sele1.selatm[ii]);  
      //}
      //printf("##############");  
      //for(ii=0;ii<inp_enm->rtblevan.nselatm;ii++)
      //{
      //  printf("%d\n",inp_enm->rtblevan.selatm[ii]);  
      //}
     }
   }

   
   if( inp_enm->enm_level_flag==1 && inp_enm->rtb_flag == 0 && inp_enm->vsa_flag == 0 && inp_enm->mol2_flag == 0 )
   {
     GeomCentre( inp_enm, &molecule->coor, molecule);
     printf("   Calculating the geometrical centre from coordinates\n");
   }


   /*Allocating memory for ANALYSIS routines arrays*/
   /*!!Either the SYSTEM or the SUBSYSTEM normal modes will be considered!!*/

   if( inp_enm->vsa_flag == 1 )
   {
     inp_enm->selean = inp_enm->vsasele;
   }
   else if (inp_enm->enm_level_flag ==1 )
   {
      inp_enm->selean = inp_enm->sele1;
      inp_enm->selean.nselatm = inp_enm->nselpart;
   }
   else
   {
     inp_enm->selean = inp_enm->sele1;
   }
   
   if( inp_enm->corr_flag == 1 )
   {
     inp_enm->cov = calloc (3*inp_enm->selean.nselatm, sizeof (float *));
     inp_enm->cov[0] = calloc(9*inp_enm->selean.nselatm*inp_enm->selean.nselatm, sizeof (float));
     for(ii=0;ii<3*inp_enm->selean.nselatm;ii++)
     {
       inp_enm->cov[ii]= inp_enm->cov[0] + 3*inp_enm->selean.nselatm*ii;
     }
     
     inp_enm->cov2 = (float **) calloc(3*inp_enm->selean.nselatm, sizeof(float *));
     for(ii=0;ii<3*inp_enm->selean.nselatm;ii++)
     {
       inp_enm->cov2[ii] = (float *) calloc(3*inp_enm->selean.nselatm, sizeof(float));
       for(jj=0;jj<3*inp_enm->selean.nselatm;jj++)
       {
         inp_enm->cov2[ii][jj] = 0.0;
       }
     }
   }
   
   //if (inp_enm->beta_flag==1 && strncmp(inp_enm->beta,"off",64) != 0 )
   //{
   //  inp_enm->bfact = calloc (inp_enm->selean.nselatm, sizeof (float*));
   //  for (ii=0;ii<inp_enm->selean.nselatm;ii++)
   //  {
   //	    inp_enm->bfact[ii] = calloc (inp_enm->selean.nselatm+1, sizeof (float));	 
   //  }
   //  inp_enm->fluc = calloc (3*inp_enm->selean.nselatm, sizeof (float*));
   //  for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
   //  {
   //     inp_enm->fluc[ii] = calloc (3*inp_enm->selean.nselatm, sizeof (float)); 
   //	 }
  // }

   
   if (inp_enm->mol2_flag==1 && strncmp(inp_enm->mol2, "off",128) !=0 ) 
   {
     inp_enm->reference = calloc (inp_enm->selean.nselatm, sizeof (float * ));
     for ( ii=0; ii<inp_enm->selean.nselatm; ii++ )
     {
       inp_enm->reference[ii] = calloc(3, sizeof(float));
     }     
     inp_enm->moving_sys = calloc (molecule->nato, sizeof (float * ));
     for ( ii=0; ii<molecule->nato; ii++ )
     {
       inp_enm->moving_sys[ii] = calloc(3, sizeof(float));
     }
     inp_enm->moving_selean = calloc (inp_enm->selean.nselatm, sizeof (float * ));
     for ( ii=0; ii<inp_enm->selean.nselatm; ii++ )
     {
       inp_enm->moving_selean[ii] = calloc(3, sizeof(float));
     }    
     inp_enm->deviation = calloc ( 3*inp_enm->selean.nselatm, sizeof (float ) );
     inp_enm->overlap = calloc ( 3*inp_enm->selean.nselatm, sizeof(float) );
     inp_enm->norm = calloc ( 3*inp_enm->selean.nselatm, sizeof(float) );        
     inp_enm->molecule2 = ReadMolecule ( inp_enm->mol2, wrd_whichmoltype(inp_enm->mol2) );
   }
   
   if(inp_enm->vecfile2_flag==1)
   {
     vecmatrix = fopen ( inp_enm->vecfile2, "r"); 
     inp_enm->vecmat2 = calloc ( 3*inp_enm->selean.nselatm, sizeof(float*) );
     for (ii=0;ii< 3*inp_enm->selean.nselatm;ii++)
     {
       inp_enm->vecmat2[ii]= calloc ( 3*inp_enm->selean.nselatm, sizeof (float));
     }     
     inp_enm->matvec_size=0;
     while ( !feof(vecmatrix))
     {
        //memset ( buffer, '\0', sizeof(buffer));
       if(fgets(buffer, 20000, vecmatrix)==NULL){
          fprintf(stderr, "Error! Premature end of the stream reached!\n");
       }
     }
     rewind(vecmatrix);     
     pch = strtok(buffer,".");
     while (pch != NULL)
     {
       if (strncmp (buffer, ".",1)>0 || strncmp (buffer, ".",1)<0)
       {
         inp_enm->matvec_size++;
       }
       pch = strtok(NULL,".");
     }
     rewind(vecmatrix);
     inp_enm->matvec_size--;
     printf("%d\n",inp_enm->matvec_size);   
     for ( ii=0; ii<3*inp_enm->selean.nselatm; ii++)
     {
       for ( jj=0; jj<inp_enm->matvec_size; jj++)
       {
         fscanf(vecmatrix, "%f", &inp_enm->vecmat2[ii][jj]);
       }
       fscanf(vecmatrix,"%*[^\n]%*c");
     }     
     inp_enm->mat_overlap = calloc ( 3*inp_enm->selean.nselatm, sizeof(float*));
     for(ii=0;ii<3*inp_enm->selean.nselatm;ii++)
     {
       inp_enm->mat_overlap[ii] = calloc ( 3*inp_enm->selean.nselatm, sizeof(float));
     } 
   }
 
   if( inp_enm->perturb_flag == 1 )
   {
      //inp_enm->Hs = calloc( 3*inp_enm->selean.nselatm, sizeof(float *) );
      //inp_enm->Hs[0] = calloc( 9*inp_enm->selean.nselatm*inp_enm->selean.nselatm, sizeof(float) );
      //for ( ii=0; ii < 3*inp_enm->selean.nselatm ; ii++)
      //{
      //  inp_enm->Hs[ii] = inp_enm->Hs[0] + 3*inp_enm->selean.nselatm*ii;
      //}
      if( inp_enm->vsa_flag == 1 )
      {     
        inp_enm->Hp = calloc( 3*inp_enm->sele1.nselatm, sizeof(float *) );
        inp_enm->Hp[0] = calloc( 9*inp_enm->sele1.nselatm*inp_enm->sele1.nselatm, sizeof(float) );
        for ( ii=0; ii < 3*inp_enm->sele1.nselatm ; ii++)
        {
          inp_enm->Hp[ii] = inp_enm->Hp[0] + 3*inp_enm->sele1.nselatm*ii;
        }
      }
      else 
      {
        inp_enm->Hp = calloc( 3*inp_enm->selean.nselatm, sizeof(float *) );
        inp_enm->Hp[0] = calloc( 9*inp_enm->selean.nselatm*inp_enm->selean.nselatm, sizeof(float) );
        for ( ii=0; ii < 3*inp_enm->selean.nselatm ; ii++)
        {
          inp_enm->Hp[ii] = inp_enm->Hp[0] + 3*inp_enm->selean.nselatm*ii;
        }
      }
     inp_enm->scalar = calloc ( 3*inp_enm->selean.nselatm, sizeof(float)); 
   }
   
   if( inp_enm->proj_flag == 1 )
   { 
    if( inp_enm->proj_flag == 1 && inp_enm->print == 1 )
    {
      memset (inp_enm->projdcd, '\0', sizeof(inp_enm->projdcd));
      sprintf (inp_enm->projdcd, "%s_proj_mode_comb.dcd",inp_enm->title);
      inp_enm->projdcdfile = O_File( inp_enm->projdcd, "w" );
      FillDcdHeader( molecule, &inp_enm->projdcdhdr );
      WriteDcdHeader( inp_enm->projdcdfile, &inp_enm->projdcdhdr, "all");
      inp_enm->coor = InitCoor( molecule->nato );
    }
    else if ( inp_enm->proj_flag == 1 && inp_enm->print == 0 )
    {
      memset (inp_enm->projdcd, '\0', sizeof (inp_enm->projdcd));
      sprintf (inp_enm->projdcd, "proj_mode_comb.dcd"); 
      inp_enm->projdcdfile = O_File( inp_enm->projdcd, "w" );
      FillDcdHeader( molecule, &inp_enm->projdcdhdr );
      WriteDcdHeader( inp_enm->projdcdfile, &inp_enm->projdcdhdr, "all");
      inp_enm->coor = InitCoor( molecule->nato );
    }
   } 
   return 0;
}
// ------------------------------------------------------------------
int Compute_Enm( struct inp_enm *inp_enm, Molecule *molecule, CoorSet *trj_crd,  char *outprint, int intex )
{
  FILE *eigval;
  FILE *eigvec;

  int ii=0, kk=0, jj=0, size;
  float aa=0, bb=0, cc=0;
  float rotmatrix[3][3], transvec[3], **reference, **moving_sys, **moving_selean;
  char  filename[1280];       

  inp_enm->current_frame = intex;
  Molecule  *mol2;   
  Molecule  *selmol;
  Molecule  *vsaselmol;  
  _dr *dr;
  _dr *drvsa;
  dr = inp_enm->dr;
  drvsa = inp_enm->drvsa;
  //CoorSet *trj_crd;
     
  selmol = MakeSeleMol (&inp_enm->sele1, molecule);   
  WriteMolecule( selmol, "selection_mol.pdb", 0);
  inp_enm->distances = calloc( 3, sizeof(float));  
  size = 3*inp_enm->sele1.nselatm;
  //ssize = size*size;   

  printf("***ENM MODULE***\n"); fflush(stdout);
  if(inp_enm->sele1_flag == 0)
  {
    printf("\n>>  Selection needed!\n    Exiting\n");fflush(stdout);
    exit(0);  
  }
  
  if( inp_enm->enm_level_flag==1 && inp_enm->rtb_flag == 0 && inp_enm->mol2_flag == 0 )
  {
    GeomCentre(inp_enm, trj_crd, molecule);
    printf("   Calculating the geometrical coordinates\n");
    //for (ii=0;ii<inp_enm->nselpart;ii++)
    //{
    //  printf(" %f %f %f \n", inp_enm->virt_coor[ii][0], inp_enm->virt_coor[ii][1], inp_enm->virt_coor[ii][2]); 
    //}
  }
   
  if ( inp_enm->mol2_flag == 1 )
  {
    mol2 = inp_enm->molecule2;
    GetSele ( inp_enm->sele2.selestring, &inp_enm->sele2, mol2 );
    //printf("\n\n%c --> %s \n",molecule->rawmol.name, inp_enm->mol2);

    if ( inp_enm->sele2_flag == 0 && inp_enm->mol2_flag == 1 )
    { 
      printf("\nError: --MOL2 option needs --SELE2!!\n");
      exit(0);
    }
 
    if (inp_enm->selean.nselatm != inp_enm->sele2.nselatm)
    {
      printf(" Different atom selections !! \n");
      printf("ipdb: %4d i2pdb: %4d \n", inp_enm->selean.nselatm, inp_enm->sele2.nselatm ); 
      exit(0);
    }
    //nn=3;
    
    reference = inp_enm->reference;
    moving_sys = inp_enm->moving_sys;
    moving_selean = inp_enm->moving_selean;

    for ( ii=0; ii<inp_enm->selean.nselatm; ii++ )
    {
      reference[ii][0]=0;
      reference[ii][1]=0;
      reference[ii][2]=0;
    } 
    for ( ii=0; ii<inp_enm->selean.nselatm; ii++ )
    {
      moving_selean[ii][0]=0;
      moving_selean[ii][1]=0;
      moving_selean[ii][2]=0;
    }    
    for ( ii=0; ii<molecule->nato; ii++ )
    {
      moving_sys[ii][0]=0;
      moving_sys[ii][1]=0;
      moving_sys[ii][2]=0;
    }    

    for ( ii=0; ii<inp_enm->selean.nselatm; ii++)
    {
      reference[ii][0]  = mol2->coor.xcoor[inp_enm->sele2.selatm[ii]-1];
      reference[ii][1]  = mol2->coor.ycoor[inp_enm->sele2.selatm[ii]-1];
      reference[ii][2]  = mol2->coor.zcoor[inp_enm->sele2.selatm[ii]-1];
      moving_selean[ii][0] = trj_crd->xcoor[inp_enm->selean.selatm[ii]-1];
      moving_selean[ii][1] = trj_crd->ycoor[inp_enm->selean.selatm[ii]-1];
      moving_selean[ii][2] = trj_crd->zcoor[inp_enm->selean.selatm[ii]-1];
    }

    for ( ii=0; ii<molecule->nato; ii++)
    {
      moving_sys[ii][0]    = trj_crd->xcoor[ii];
      moving_sys[ii][1]    = trj_crd->ycoor[ii];
      moving_sys[ii][2]    = trj_crd->zcoor[ii];
      //printf("%f %f %f \n", moving[ii][0], moving[ii][1], moving[ii][2]);
    }
    
    //*Getting roto-translation indexes for the selected regions*/
    superimpose ( inp_enm->selean.nselatm, reference, moving_selean, rotmatrix, transvec );

    //*Apply roto-translation to the whole molecule*/     
    for ( ii=0; ii<molecule->nato; ii++)
    {
      aa = moving_sys[ii][0]*rotmatrix[0][0] + moving_sys[ii][1]*rotmatrix[0][1] + moving_sys[ii][2]*rotmatrix[0][2] ;
      bb = moving_sys[ii][0]*rotmatrix[1][0] + moving_sys[ii][1]*rotmatrix[1][1] + moving_sys[ii][2]*rotmatrix[1][2] ;
      cc = moving_sys[ii][0]*rotmatrix[2][0] + moving_sys[ii][1]*rotmatrix[2][1] + moving_sys[ii][2]*rotmatrix[2][2] ;
      moving_sys[ii][0] = aa + transvec[0]; 
      moving_sys[ii][1] = bb + transvec[1]; 
      moving_sys[ii][2] = cc + transvec[2]; 
    }
    /* Fitting of the Atom coordinates used for Hessian construction */
    for ( ii=0; ii<inp_enm->sele1.nselatm; ii++ )
    {
      inp_enm->coord[ii][0] = moving_sys[inp_enm->sele1.selatm[ii]-1][0];
      inp_enm->coord[ii][1] = moving_sys[inp_enm->sele1.selatm[ii]-1][1];
      inp_enm->coord[ii][2] = moving_sys[inp_enm->sele1.selatm[ii]-1][2];
    }
  } 
  else if (inp_enm->mol2_flag != 1)  
  {
    for (ii=0;ii<inp_enm->sele1.nselatm;ii++)
    {
      inp_enm->coord[ii][0] = trj_crd->xcoor[inp_enm->sele1.selatm[ii]-1];
      inp_enm->coord[ii][1] = trj_crd->ycoor[inp_enm->sele1.selatm[ii]-1];
      inp_enm->coord[ii][2] = trj_crd->zcoor[inp_enm->sele1.selatm[ii]-1];
    }
  }
  

/////*Determination of the INTERACTING PAIRS*/////
/////*Computing 2nd derivative of potential energy function*/////

  Get_Dr (inp_enm, dr,trj_crd, molecule);    
  
	       
/////*HESSIAN GENERATION*/////
  if ( inp_enm->enm_level_flag==0)
  {
    Build_Hessian (inp_enm, inp_enm->He, dr, 3*inp_enm->sele1.nselatm);  
  }
  else if ( inp_enm->enm_level_flag==1)
  {
    Build_Hessian (inp_enm, inp_enm->He, dr, 3*inp_enm->nselpart);    
  }
/////*HESSIAN PROJECTION in the T/R subspace using the RTB method*/
  if( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->enm_level_flag==0)
  {
    Enm_Rtb ( inp_enm, molecule ); 
  }
  else if ( inp_enm->rtb_flag == 0 && inp_enm->enm_level_flag==0)
  {
    inp_enm->eigen->size   = 3*inp_enm->sele1.nselatm; 
    inp_enm->eigen->inpmat = inp_enm->He;
  }
  else if ( inp_enm->enm_level_flag==1 && inp_enm->rtb_flag == 0 && inp_enm->mol2_flag == 0 )
  {
    inp_enm->eigen->size   = 3*inp_enm->nselpart; 
    inp_enm->eigen->inpmat = inp_enm->He;        
  }

/////*VSA METHOD*////

  if  ( inp_enm->vsa_flag==1 &&  inp_enm->vsasele.nselatm < inp_enm->sele1.nselatm && inp_enm->env.nselatm > 0) 
  {
    Enm_GetVsaMat (inp_enm, inp_enm->He, inp_enm->Hss);
    Get_DrVSA (inp_enm, dr, drvsa);
    //Enm_DiagVsaMat( inp_enm, inp_enm->subeigen );
    inp_enm->eigen->size = 3*inp_enm->vsasele.nselatm; 
    inp_enm->eigen->inpmat = inp_enm->Hss; 
  }
  else if (inp_enm->env.nselatm == 0 && inp_enm->vsasele.nselatm == inp_enm->sele1.nselatm)
  {
    printf("Skipping VSA, system equals subsystem!\n");
    inp_enm->vsa_flag = 0;
  }
   
/////*HESSIAN DIAGONALIZATION*////////
  printf("%d\n",inp_enm->eigen->size);fflush(stdout);
  printf("\n>> diagonalization\n");fflush(stdout);
  DiagMatrix ( inp_enm->eigen );
 /////*If RTB method is called, normal modes in the full space are obtained from the R/T space*/
 
  if( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {
    Enm_Rtb2full ( inp_enm, molecule );    
  }

////*Variance calculation*////
  Enm_Variance(inp_enm, inp_enm->eigen);
   
/////*Printing output*////
////*Printing all eig-values*////

  printf("\n>  Printing output\n");
  if(inp_enm->print==0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "eigval.txt");    
  }
  else if (inp_enm->print==1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-eigval_fr%d.txt", inp_enm->title, inp_enm->current_frame);
  }
  eigval = fopen (filename, "w"); 
  for (ii=0; ii < inp_enm->eigen->size; ii++)
  {
    fprintf( eigval, "%15.10f \n", inp_enm->eigen->eigval[inp_enm->eigen->size-ii-1]);
  }
  fclose(eigval);

////*Printing eig-vectors excluding the first 6 zero frequency modes*////      
  if(inp_enm->print==0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "eigvec.txt");    
  }
  else if (inp_enm->print==1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-eigvec_fr%d.txt", inp_enm->title,inp_enm->current_frame);
  }     
  eigvec = fopen (filename, "w");
  if( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {
    for (ii=0; ii<size; ii++)
    {
      for (jj=0;jj<inp_enm->eigen->size/3;jj++)
      {
        fprintf(eigvec,"%15.10f", inp_enm->rtb2full_eigen->eigvec[inp_enm->rtb2full_eigen->size-jj-6-1][ii]);
      }
      fprintf(eigvec,"\n");
    }
  }
  if( inp_enm->rtb_flag != 1 )
  {
    for (ii=0; ii<inp_enm->eigen->size; ii++)
    {
      for (jj=0;jj<inp_enm->eigen->size/3;jj++)
      {
        fprintf(eigvec,"%15.10f", inp_enm->eigen->eigvec[inp_enm->eigen->size-jj-6-1][ii]);
      }
      fprintf(eigvec,"\n");
    }
  }
  fclose(eigvec);
      
////////Printing eigenvectors in a dcd file

  if(inp_enm->vsa_flag==0 && inp_enm->rtb_flag==0)
  {
    if(inp_enm->print==0)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "eigvec.dcd");    
    }
    else if (inp_enm->print==1)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "%s-eigvec_fr%d.dcd", inp_enm->title,inp_enm->current_frame);
    }
    
    inp_enm->evectrj = O_File( filename, "w" );
    FillDcdHeader( selmol, &inp_enm->evechdr );
    WriteDcdHeader( inp_enm->evectrj, &inp_enm->evechdr, "all");
    inp_enm->evecoor = InitCoor( inp_enm->eigen->size/3 );
    
    kk=0;
    while (kk < inp_enm->eigen->size-5)
    {
      for (ii=0; ii<inp_enm->eigen->size/3; ii++)
      {
        if ( kk == 0 )
        {
          inp_enm->evecoor->xcoor[ii] = inp_enm->coord[ii][0];
          inp_enm->evecoor->ycoor[ii] = inp_enm->coord[ii][1];
          inp_enm->evecoor->zcoor[ii] = inp_enm->coord[ii][2];          
        }
        else if ( kk > 0 )
        {
          inp_enm->evecoor->xcoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3];
          inp_enm->evecoor->ycoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3+1];
   	      inp_enm->evecoor->zcoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3+2];
         //fprintf(eigvec,"%15.10f", inp_enm->eigen->eigvec[inp_enm->eigen->size-jj-1][ii]);
        }
      }
    
      WriteDcdCoor( inp_enm->evectrj, inp_enm->evecoor, &inp_enm->evechdr);
      kk++;      
    }

    inp_enm->evechdr.nframe = kk;
    WriteDcdHeader( inp_enm->evectrj, &inp_enm->evechdr, "nframe" );
    fclose(inp_enm->evectrj);
    
  } 
  else if (inp_enm->vsa_flag==1)
  {
    ////Printing eigenvectors in a trajectory file
    vsaselmol = MakeSeleMol (&inp_enm->vsasele, molecule);
    if(inp_enm->print==0)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "eigvec.dcd");
    }
    else if (inp_enm->print==1)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "%s-eigvec_fr%d.dcd", inp_enm->title,inp_enm->current_frame);
    }     
    
    inp_enm->evectrj = O_File( filename, "w" );
    FillDcdHeader( vsaselmol, &inp_enm->evechdr );
    WriteDcdHeader( inp_enm->evectrj, &inp_enm->evechdr, "all");
    inp_enm->evecoor = InitCoor( inp_enm->eigen->size/3 );
    
    kk=0;
    while (kk < inp_enm->eigen->size-5)
    {

      for (ii=0; ii<inp_enm->eigen->size/3; ii++)
      {
        if ( kk == 0 )
        {
          inp_enm->evecoor->xcoor[ii] = trj_crd->xcoor[inp_enm->vsasele.selatm[ii]-1];
          inp_enm->evecoor->ycoor[ii] = trj_crd->ycoor[inp_enm->vsasele.selatm[ii]-1];
          inp_enm->evecoor->zcoor[ii] = trj_crd->zcoor[inp_enm->vsasele.selatm[ii]-1];
        }
        else if ( kk > 0 )
        {
          inp_enm->evecoor->xcoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3];
          inp_enm->evecoor->ycoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3+1];
   	      inp_enm->evecoor->zcoor[ii] = inp_enm->eigen->eigvec[inp_enm->eigen->size-kk-6][ii*3+2];
        }
      }
      //printf("%d\n",kk);fflush(stdout);
      WriteDcdCoor( inp_enm->evectrj, inp_enm->evecoor, &inp_enm->evechdr);
      kk++;      
    }
    
    inp_enm->evechdr.nframe = kk;
    WriteDcdHeader( inp_enm->evectrj, &inp_enm->evechdr, "nframe" );
    fclose(inp_enm->evectrj);
  }
   
   
   /////////////////////////////////////////////////////////////////////////
   //////////////////////////////*ENM ANALYSES MODULES*/////////////////////
   /////////////////////////////////////////////////////////////////////////
   
/////*CORRELATIONS*/////

   if( inp_enm->corr_flag == 1  && strncmp(inp_enm->corr, "off", 64) != 0 && inp_enm->rtb_flag == 0 )
   {
     Enm_Correlation ( inp_enm, inp_enm->eigen,  molecule );    
   }
   else if ( inp_enm->corr_flag == 1 && strncmp(inp_enm->corr, "off", 64) != 0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
   {
      Enm_Correlation ( inp_enm, inp_enm->rtb2full_eigen,  molecule );    
   }   
        
/////*RMSF*/////

   if ( inp_enm->beta_flag == 1 && strncmp(inp_enm->beta,"off",64) != 0 && inp_enm->rtb_flag == 0 )
   {
     Enm_Fluc ( inp_enm, molecule, inp_enm->eigen );
   }
   else if ( inp_enm->beta_flag == 1  && strncmp(inp_enm->beta,"off",64) != 0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
   {
     Enm_Fluc ( inp_enm, molecule, inp_enm->rtb2full_eigen);
   }

/////*INVOLVEMENT COEFFICIENT*/////

   if  (  inp_enm->mol2_flag == 1  && strncmp(inp_enm->mol2, "off",128) !=0 && inp_enm->rtb_flag == 0 ) 
   {
     Enm_Compare(inp_enm, molecule, trj_crd, inp_enm->eigen);
   }
   else if (  inp_enm->mol2_flag == 1  && strncmp(inp_enm->mol2,"off",64) != 0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0)
   {
     Enm_Compare (inp_enm, molecule, trj_crd, inp_enm->rtb2full_eigen);
   }

/////*STRUCTURAL PERTURBATION METHOD*/

   if (  inp_enm->perturb_flag == 1  && strncmp(inp_enm->perturb,"off",64) != 0 && inp_enm->rtb_flag == 0)
   {
     Enm_Perturb(inp_enm, molecule, inp_enm->eigen, dr, trj_crd);
   }
   else if (  inp_enm->perturb_flag == 1  && strncmp(inp_enm->perturb,"off",64) != 0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0)
   {
     Enm_Perturb(inp_enm, molecule, inp_enm->rtb2full_eigen, dr, trj_crd);
   }   
 
/////*DEFORMATION ENERGY*/


   if ( inp_enm->defen_flag == 1 && strncmp(inp_enm->defen,"off",64) != 0 && inp_enm->rtb_flag == 0 )
   {
     Enm_Defen ( inp_enm, molecule, inp_enm->eigen, trj_crd );
   }
   else if ( inp_enm->defen_flag == 1  && strncmp(inp_enm->defen,"off",64) != 0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0)
   {
     Enm_Defen ( inp_enm, molecule, inp_enm->rtb2full_eigen, trj_crd);
   }


/////*PROJECTIONS*/////

   if  (  inp_enm->proj_flag == 1 && strncmp(inp_enm->proj, "off",128) !=0 && inp_enm->rtb_flag == 0 ) 
   {
     Enm_Proj(inp_enm, inp_enm->eigen, molecule, trj_crd);
   }
   else if  ( inp_enm->proj_flag == 1 && strncmp(inp_enm->proj, "off",128) !=0  && inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 ) 
   {
     Enm_Proj(inp_enm, inp_enm->rtb2full_eigen, molecule, trj_crd);
   }

  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////*COMPUTING SUBSYSTEM ENM MODES USING VSA (Vibrational Subsystem Analaysis) METHOD*///////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*Building the Subsystem Hessian*/

int Enm_GetVsaMat ( struct inp_enm *inp_enm, float **He, float **Hss )
{
  int ii=0, jj=0, kk=0;
  float **invHee_Hes, **Hsub;
  
  printf("\n>> VIBRATIONAL SUBSYSTEM ANALYSIS\n");  
    
  invHee_Hes = calloc( 3*inp_enm->env.nselatm, sizeof(float *) );
  for ( ii=0; ii < 3*inp_enm->env.nselatm ; ii++)
  {
    invHee_Hes[ii] = calloc (3*inp_enm->vsasele.nselatm, sizeof (float));
  }    
  Hsub = calloc( 3*inp_enm->vsasele.nselatm, sizeof (float *) );
  for ( ii=0; ii < 3*inp_enm->vsasele.nselatm ; ii++ )
  {
    Hsub[ii] = calloc (3*inp_enm->vsasele.nselatm, sizeof (float));
  }
  
  
  for ( ii=0; ii < 3*inp_enm->env.nselatm ; ii++)
  {
    for ( jj=0; jj < 3*inp_enm->env.nselatm ; jj++)
    {
      inp_enm->Hee[ii][jj] = 0;
    }
    for ( jj=0; jj < 3*inp_enm->vsasele.nselatm; jj++ )
    {
      inp_enm->Hes[ii][jj] = 0;
    }
  }  
 
  for (ii=0;ii<3*inp_enm->vsasele.nselatm;ii++)
  {
    for (jj=0;jj<3*inp_enm->vsasele.nselatm;jj++)
    {
      Hss[ii][jj] = 0;
    }
  }
 
  //if(inp_enm->mass_flag == 1)
  //{
  //}
  
   /*Filling the Hss matrix*/  
  for (ii=0;ii<inp_enm->vsasele.nselatm;ii++)
  {
    for (jj=0;jj<inp_enm->vsasele.nselatm;jj++)
    {
      Hss[3*ii][3*jj]     = He[3*inp_enm->subsys_list[ii][0]][3*inp_enm->subsys_list[jj][0]];
      Hss[3*ii][3*jj+1]   = He[3*inp_enm->subsys_list[ii][0]][3*inp_enm->subsys_list[jj][0]+1];
      Hss[3*ii][3*jj+2]   = He[3*inp_enm->subsys_list[ii][0]][3*inp_enm->subsys_list[jj][0]+2];
      Hss[3*ii+1][3*jj]   = He[3*inp_enm->subsys_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]];
      Hss[3*ii+1][3*jj+1] = He[3*inp_enm->subsys_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]+1];
      Hss[3*ii+1][3*jj+2] = He[3*inp_enm->subsys_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]+2];
      Hss[3*ii+2][3*jj]   = He[3*inp_enm->subsys_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]];
      Hss[3*ii+2][3*jj+1] = He[3*inp_enm->subsys_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]+1];
      Hss[3*ii+2][3*jj+2] = He[3*inp_enm->subsys_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]+2];
    }
  }   
  /*Filling the Hee matrix*/
  for (ii=0;ii<inp_enm->env.nselatm;ii++)
  {
    for (jj=0;jj<inp_enm->env.nselatm;jj++)
    {
      inp_enm->Hee[3*ii][3*jj]     = He[3*inp_enm->env_list[ii][0]][3*inp_enm->env_list[jj][0]];
      inp_enm->Hee[3*ii][3*jj+1]   = He[3*inp_enm->env_list[ii][0]][3*inp_enm->env_list[jj][0]+1];
      inp_enm->Hee[3*ii][3*jj+2]   = He[3*inp_enm->env_list[ii][0]][3*inp_enm->env_list[jj][0]+2];
      inp_enm->Hee[3*ii+1][3*jj]   = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->env_list[jj][0]];
      inp_enm->Hee[3*ii+1][3*jj+1] = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->env_list[jj][0]+1];
      inp_enm->Hee[3*ii+1][3*jj+2] = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->env_list[jj][0]+2];
      inp_enm->Hee[3*ii+2][3*jj]   = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->env_list[jj][0]];
      inp_enm->Hee[3*ii+2][3*jj+1] = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->env_list[jj][0]+1];
      inp_enm->Hee[3*ii+2][3*jj+2] = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->env_list[jj][0]+2];
    }
    /*Filling the Hes matrix*/
    for (jj=0;jj<inp_enm->vsasele.nselatm;jj++)
    {
      inp_enm->Hes[3*ii][3*jj]     = He[3*inp_enm->env_list[ii][0]][3*inp_enm->subsys_list[jj][0]];
      inp_enm->Hes[3*ii][3*jj+1]   = He[3*inp_enm->env_list[ii][0]][3*inp_enm->subsys_list[jj][0]+1];
      inp_enm->Hes[3*ii][3*jj+2]   = He[3*inp_enm->env_list[ii][0]][3*inp_enm->subsys_list[jj][0]+2];
      inp_enm->Hes[3*ii+1][3*jj]   = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]];
      inp_enm->Hes[3*ii+1][3*jj+1] = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]+1];
      inp_enm->Hes[3*ii+1][3*jj+2] = He[3*inp_enm->env_list[ii][0]+1][3*inp_enm->subsys_list[jj][0]+2];
      inp_enm->Hes[3*ii+2][3*jj]   = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]];
      inp_enm->Hes[3*ii+2][3*jj+1] = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]+1];
      inp_enm->Hes[3*ii+2][3*jj+2] = He[3*inp_enm->env_list[ii][0]+2][3*inp_enm->subsys_list[jj][0]+2];
    }
  }

  /*Obtaining the inverse of Hee matrix, Hee^-1*/  
  InvMatrix(3*inp_enm->env.nselatm, inp_enm->Hee, inp_enm->invHee);


  /* Matrix Identity check!*/
  //for (ii=0;ii<3*inp_enm->env.nselatm;ii++)
  //{
  //  for (jj=0;jj<3*inp_enm->env.nselatm;jj++)
  //  {
  //    temp=0;
  //    for (kk=0;kk<3*inp_enm->env.nselatm;kk++)
  //    {
  //      temp+=inp_enm->Hee[ii][kk]*inp_enm->invHee[kk][jj];
  //      //temp+=inp_enm->invHee[ii][kk]*inp_enm->Hee[kk][jj];
  //    }
  //    printf("%f ",temp);
  //  }
  //  printf("\n");
  //}
  
  /*dot(Hee^-1, Hes)*/
  for (ii=0;ii<3*inp_enm->env.nselatm;ii++)
  {
    for (jj=0;jj<3*inp_enm->vsasele.nselatm;jj++)
    {
      //temp=0;
      for (kk=0;kk<3*inp_enm->env.nselatm;kk++)
      {
        invHee_Hes[ii][jj] += inp_enm->invHee[ii][kk]*inp_enm->Hes[kk][jj];
      }
    }
  }
  /*dot (Hes^T, invHee_Hes[ii][jj])*/ 
  for (ii=0;ii<3*inp_enm->vsasele.nselatm;ii++)
  {
    for (jj=0;jj<3*inp_enm->vsasele.nselatm;jj++)
    {
      //temp=0;
      for (kk=0;kk<3*inp_enm->env.nselatm;kk++)
      {
        Hsub[ii][jj] += (inp_enm->Hes[kk][ii]*invHee_Hes[kk][jj]);
      }
    }
  }
  /*Hss-Hsub*/
  for (ii=0;ii<3*inp_enm->vsasele.nselatm;ii++)
  {
    for (jj=0;jj<3*inp_enm->vsasele.nselatm;jj++)
    {
      Hss[ii][jj] -= Hsub[ii][jj];
    }
  }
   
  return (0);
}

/*Diagonalizing the Subsystem Hessian*/
int Enm_DiagVsaMat ( struct inp_enm *inp_enm, _eigen *subeigen)
{
  FILE *subeigval;
  FILE *subeigvec;
  int ii=0, jj=0;
  char  filename[1280];  

  
  inp_enm->subeigen->size = 3*inp_enm->vsasele.nselatm; 
  inp_enm->subeigen->inpmat = inp_enm->Hss; 
  DiagMatrix( inp_enm->subeigen );
       
  if(inp_enm->print==0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "subeigval.txt");    
  }
  else if (inp_enm->print==1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-subeigval_fr%d.txt", inp_enm->title, inp_enm->current_frame);
  }
  subeigval = fopen (filename, "w"); 
  for (ii=0; ii<3*inp_enm->vsasele.nselatm; ii++)
  {
    fprintf( subeigval, "%15.10f \n", inp_enm->subeigen->eigval[subeigen->size-ii-1]);
  }
  fclose(subeigval);

  if(inp_enm->print==0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "subeigvec.txt");    
  }
  else if (inp_enm->print==1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-subeigvec_fr%d.txt", inp_enm->title,inp_enm->current_frame);
  }     

  subeigvec = fopen (filename, "w");
  for (ii=0; ii<subeigen->size; ii++)
  {
    for (jj=0;jj<subeigen->size/3;jj++)
    {
      fprintf(subeigvec,"%15.10f", inp_enm->subeigen->eigvec[subeigen->size-jj-6-1][ii]);
    }
    fprintf(subeigvec,"\n");
  }

  fclose(subeigvec);
 

 
    
  return 0;
}  

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////*ROTATION-TRANSLATION BLOCKS ENM METHOD*/////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int Enm_Getblock( struct inp_enm *inp_enm, Molecule *molecule )
{
  FILE *block_file;
  int ii=0, jj=0, kk=0, res=0, totatom=0;  
  
  inp_enm->nblock=0;
  inp_enm->atm2blk = calloc(inp_enm->sele1.nselatm, sizeof(int));

  printf("\n>> ROTATION TRANSLATION BLOCK\n");    
  if(strncmp(inp_enm->rtb_method,"residue",64)==0)
  {     // RTB residue-level
    inp_enm->nblock = inp_enm->selmol->nRes;  
    inp_enm->blk = calloc (inp_enm->nblock, sizeof (_blk));
    /*Building an array of "block" structures with lists of atoms and masses corresponding to each block*/  
    /*Current implementation uses unitary masses for each atom !*/  
    res=0;
    for ( ii=0; ii<inp_enm->selmol->nSeg; ii++ )
    {
      for ( jj=0; jj<inp_enm->selmol->segment[ii].nRpS; jj++ )
      {
        inp_enm->blk[res].atom = calloc ( inp_enm->selmol->segment[ii].pRes[jj].nApR, sizeof(int) );
        inp_enm->blk[res].sindex = calloc ( inp_enm->selmol->segment[ii].pRes[jj].nApR, sizeof(int) );
        inp_enm->blk[res].atom_mass = calloc ( inp_enm->selmol->segment[ii].pRes[jj].nApR, sizeof(float) );
        inp_enm->blk[res].natom = 0;
	      inp_enm->blk[res].block_mass = 0.0;
        inp_enm->blk[res].resn = inp_enm->selmol->segment[ii].pRes[jj].resn;
        inp_enm->blk[res].blockid = res;
        for ( kk=0; kk<inp_enm->selmol->segment[ii].pRes[jj].nApR; kk++ )
        { 
	        inp_enm->blk[res].natom++;
          inp_enm->blk[res].atom[kk] = inp_enm->selmol->segment[ii].pRes[jj].pAto[kk].atomn;
          inp_enm->blk[res].sindex[kk] = inp_enm->selmol->segment[ii].pRes[jj].pAto[kk].atomn-1;
          //WARNING : should check usage of blk.atom in the future !!!
          //inp_enm->blk[res].atom[kk] = getSeleIndex( inp_enm->selmol->segment[ii].pRes[jj].pAto[kk].atomn, &inp_enm->sele1 );
          inp_enm->atm2blk[inp_enm->blk[res].atom[kk]-1] = res;
          if(inp_enm->mass_flag == 1)
	        {
            /*Masses supplied in PDB's B-factor column are considered in the computation */
            inp_enm->blk[res].atom_mass[kk] = inp_enm->selmol->segment[ii].pRes[jj].pAto[kk].bFac;
          }
          else if (inp_enm->mass_flag == 0)
          {
	          inp_enm->blk[res].atom_mass[kk] = 1.0;
	        }
          inp_enm->blk[res].block_mass += inp_enm->blk[res].atom_mass[kk];
          //printf("atom= %d , atom_mass= %f, block_mass= %f\n", inp_enm->blk[res].atom[kk], inp_enm->blk[res].atom_mass[kk], inp_enm->blk[res].block_mass );
        }
        if( inp_enm->blk[res].natom < 3 )
        {
          printf("\n  ERROR: some residues appear to have less than 3 atoms, RTB approximation could not be applied!\n");
          printf("  Use selections file to incorporate these residues into bigger blocks.\n\n");
          exit(0);  
        }
        res++;
      }
    }
  }
  else if ( strncmp(inp_enm->rtb_method,"residue",64) !=0 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {     // RTB with block file
    printf("\n  Reading Block Selection file\n");
    block_file = fopen (inp_enm->rtb_method, "r");
    
    inp_enm->nblock = linecount( block_file );
    printf("Nr. of blocks %d\n", inp_enm->nblock);
    
    //Each atom is assigned to the specific block
    inp_enm->blk = calloc (inp_enm->nblock, sizeof (_blk));

    totatom = 0;
    ii = 0;
    while ( !feof(block_file) && ii < inp_enm->nblock)
    {
      //memset(inp_enm->blk[ii].sele.selestring,'\0', sizeof(inp_enm->blk[ii].sele.selestring));
      memset(inp_enm->blk[ii].sele.selestring,'\0',1024 );
      //fscanf(block_file,"%s", inp_enm->blk[ii].sele.selestring);
      if(fgets(inp_enm->blk[ii].sele.selestring, 1024, block_file)==NULL && ( !feof(block_file) || ferror(block_file) ))
      {
        fprintf(stderr, "Warning! Premature end of file reached!\n");
      }
      GetSele(inp_enm->blk[ii].sele.selestring, &inp_enm->blk[ii].sele, molecule);
      totatom += inp_enm->blk[ii].sele.nselatm;     
      ii++;
    }
    for(ii=0;ii<inp_enm->nblock;ii++)
    { 
      inp_enm->blk[ii].atom = calloc ( inp_enm->blk[ii].sele.nselatm, sizeof(int) );
      inp_enm->blk[ii].sindex = calloc ( inp_enm->blk[ii].sele.nselatm, sizeof(int) );
      inp_enm->blk[ii].atom_mass = calloc ( inp_enm->blk[ii].sele.nselatm, sizeof(float) );
      inp_enm->blk[ii].natom = inp_enm->blk[ii].sele.nselatm;
      inp_enm->blk[ii].block_mass = 0.0;
      inp_enm->blk[ii].blockid = ii;
      for ( kk=0; kk<inp_enm->blk[ii].sele.nselatm; kk++ )
      { 
        inp_enm->blk[ii].atom[kk] = inp_enm->blk[ii].sele.selatm[kk];
        inp_enm->blk[ii].sindex[kk] = getSeleIndex( inp_enm->blk[ii].sele.selatm[kk], &inp_enm->sele1 );
        // added to handle defen
        inp_enm->atm2blk[inp_enm->blk[ii].sindex[kk]] = ii;
        //printf("%d\n", inp_enm->blk[ii].atom[kk]);fflush(stdout);
        if(inp_enm->mass_flag == 1)
	      {
          /*Masses supplied in PDB's B-factor column are considered in the computation*/
          inp_enm->blk[ii].atom_mass[kk] = molecule->rawmol.bFac[inp_enm->blk[jj].atom[ii]-1];
        }
        else if (inp_enm->mass_flag == 0)
        {
	        inp_enm->blk[ii].atom_mass[kk] = 1.0;
	      }
        inp_enm->blk[ii].block_mass += inp_enm->blk[ii].atom_mass[kk];
      }
      //printf("%d %d %d\n",ii, inp_enm->blk[ii].sele.nselatm, inp_enm->nblock);fflush(stdout); 
    }
    
    if( totatom != inp_enm->sele1.nselatm )
    {
      printf("%d %d\n", totatom, inp_enm->sele1.nselatm);
      printf("  ERROR: not all selected atoms have been assigned to specified blocks!\n");
      printf("  Please check your Block Selection file\n");
      exit(0);
    }
  }  
  else
  {
    printf("\nError: could not understand option\n");
    exit(0);
  }
  
  return (0);
}
  
int Enm_Rtb ( struct inp_enm *inp_enm, Molecule *molecule )
{
  int ii=0, jj=0, kk=0, bc=0, bc2=0, indexbc=0, indexbc2=0,  mmm=0;
  float  **D, **Hrt, **HD, s[6][6], **crd;

  
  /*Projection matrix inp_enm->RT defining the translation/rotation blocks*/
  
  crd = calloc (inp_enm->selean.nselatm, sizeof (float *));
  for(ii=0; ii<inp_enm->selean.nselatm;ii++)
    crd[ii]=calloc (3, sizeof(float));
  
  for( ii=0; ii<inp_enm->selean.nselatm; ii++)
  {
    crd[ii][0] = inp_enm->coord[ii][0];
    crd[ii][1] = inp_enm->coord[ii][1];
    crd[ii][2] = inp_enm->coord[ii][2];
  }
   
  
  inp_enm->RT = calloc (inp_enm->nblock, sizeof (float**));
  for (ii=0;ii<inp_enm->nblock;ii++)
  {
    inp_enm->RT[ii] = calloc(3*inp_enm->blk[ii].natom, sizeof (float*));
    for(jj=0;jj<3*inp_enm->blk[ii].natom; jj++)
    {
      inp_enm->RT[ii][jj] = calloc (6, sizeof(float));
    }
  }

  /*Hessian reduced in the RT subspace*/
  
  Hrt = calloc (6*inp_enm->nblock,sizeof(float*));
  for(ii=0;ii<6*inp_enm->nblock;ii++)
  {
    Hrt[ii] = calloc (6*inp_enm->nblock,sizeof(float));
  }
  
  /*Origin change inside each block*/  
  //printf("%d  %d\n",res,inp_enm->nblock);
  
  for (bc=0; bc<inp_enm->nblock; bc++)
  { 
    inp_enm->blk[bc].x = 0.0;
    inp_enm->blk[bc].y = 0.0;
    inp_enm->blk[bc].z = 0.0;
 
    D = calloc ( 3*inp_enm->blk[bc].natom, sizeof (float*) );
    for (jj=0; jj<3*inp_enm->blk[bc].natom; jj++)
    {
      D[jj] = calloc(6, sizeof (float));
      D[jj][0] = 0.0;
      D[jj][1] = 0.0;
      D[jj][2] = 0.0;
      D[jj][3] = 0.0;
      D[jj][4] = 0.0;
      D[jj][5] = 0.0;
    }
    
    //printf("###############################################\n");
    for ( jj=0; jj<inp_enm->blk[bc].natom; jj++ )
    {
      //inp_enm->blk[bc].x = inp_enm->blk[bc].x + crd[inp_enm->blk[bc].atom[jj]-1][0]*inp_enm->blk[bc].atom_mass[jj];
      //inp_enm->blk[bc].y = inp_enm->blk[bc].y + crd[inp_enm->blk[bc].atom[jj]-1][1]*inp_enm->blk[bc].atom_mass[jj];
      //inp_enm->blk[bc].z = inp_enm->blk[bc].z + crd[inp_enm->blk[bc].atom[jj]-1][2]*inp_enm->blk[bc].atom_mass[jj];
      inp_enm->blk[bc].x = inp_enm->blk[bc].x + crd[inp_enm->blk[bc].sindex[jj]][0]*inp_enm->blk[bc].atom_mass[jj];
      inp_enm->blk[bc].y = inp_enm->blk[bc].y + crd[inp_enm->blk[bc].sindex[jj]][1]*inp_enm->blk[bc].atom_mass[jj];
      inp_enm->blk[bc].z = inp_enm->blk[bc].z + crd[inp_enm->blk[bc].sindex[jj]][2]*inp_enm->blk[bc].atom_mass[jj];
      //printf("%f %f %f\n",crd[inp_enm->blk[bc].atom[jj]-1][0], crd[inp_enm->blk[bc].atom[jj]-1][1], crd[inp_enm->blk[bc].atom[jj]-1][2] );
      //printf("blk>>%f %f %f\n",inp_enm->blk[bc].x, inp_enm->blk[bc].y, inp_enm->blk[bc].z );
    }
    
    inp_enm->blk[bc].x = (inp_enm->blk[bc].x)/(inp_enm->blk[bc].block_mass); 
    inp_enm->blk[bc].y = (inp_enm->blk[bc].y)/(inp_enm->blk[bc].block_mass); 
    inp_enm->blk[bc].z = (inp_enm->blk[bc].z)/(inp_enm->blk[bc].block_mass); 
    
    
        
    for ( jj=0; jj<inp_enm->blk[bc].natom; jj++ )
    {
      //crd[inp_enm->blk[bc].atom[jj]-1][0] = crd[inp_enm->blk[bc].atom[jj]-1][0] - inp_enm->blk[bc].x;
      //crd[inp_enm->blk[bc].atom[jj]-1][1] = crd[inp_enm->blk[bc].atom[jj]-1][1] - inp_enm->blk[bc].y;
      //crd[inp_enm->blk[bc].atom[jj]-1][2] = crd[inp_enm->blk[bc].atom[jj]-1][2] - inp_enm->blk[bc].z;   
      crd[inp_enm->blk[bc].sindex[jj]][0] = crd[inp_enm->blk[bc].sindex[jj]][0] - inp_enm->blk[bc].x;
      crd[inp_enm->blk[bc].sindex[jj]][1] = crd[inp_enm->blk[bc].sindex[jj]][1] - inp_enm->blk[bc].y;
      crd[inp_enm->blk[bc].sindex[jj]][2] = crd[inp_enm->blk[bc].sindex[jj]][2] - inp_enm->blk[bc].z;   
      //printf("%f %f %f\n",crd[inp_enm->blk[bc].atom[jj]-1][0], crd[inp_enm->blk[bc].atom[jj]-1][1], crd[inp_enm->blk[bc].atom[jj]-1][2] );   
    }
    
    for ( jj=0; jj<inp_enm->blk[bc].natom; jj++ )
    {
      D[jj*3][0]   =  sqrt(inp_enm->blk[bc].atom_mass[jj]);
      D[jj*3+1][1] =  sqrt(inp_enm->blk[bc].atom_mass[jj]);
      D[jj*3+2][2] =  sqrt(inp_enm->blk[bc].atom_mass[jj]); 
      //D[jj*3+1][3] = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][2]; 
      //D[jj*3+2][3] =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][1];
      //D[jj*3][4]   =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][2];
      //D[jj*3+2][4] = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][0];
      //D[jj*3][5]   = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][1];
      //D[jj*3+1][5] =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].atom[jj]-1][0];
      D[jj*3+1][3] = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][2]; 
      D[jj*3+2][3] =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][1];
      D[jj*3][4]   =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][2];
      D[jj*3+2][4] = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][0];
      D[jj*3][5]   = -sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][1];
      D[jj*3+1][5] =  sqrt(inp_enm->blk[bc].atom_mass[jj])*crd[inp_enm->blk[bc].sindex[jj]][0];
    }
   
    /*Ortonormalization using Graham-Schmidt procedure*/    

    mmm=6;
    OrtoNorm (mmm, 3*inp_enm->blk[bc].natom, D);
    
    for(kk=0;kk<6;kk++)
    {
      for (jj=0;jj<3*inp_enm->blk[bc].natom;jj++)
      {
        inp_enm->RT[bc][jj][kk] = D[jj][kk];
      }
    }

    for (jj=0; jj<3*inp_enm->blk[bc].natom; jj++)
    {
      free(D[jj]);
    }
    free(D);

  }
  
  indexbc=0; 
  indexbc2=0;
  for(bc=0; bc<inp_enm->nblock; bc++)
  {
    for(bc2=bc; bc2<inp_enm->nblock; bc2++)
    {
      HD = calloc ( 3*inp_enm->blk[bc].natom, sizeof (float*));
      for (ii=0; ii<3*inp_enm->blk[bc].natom; ii++)
      {
        HD[ii]    = calloc(6, sizeof (float));
        HD[ii][0] = 0.0;
        HD[ii][1] = 0.0;
        HD[ii][2] = 0.0;
        HD[ii][3] = 0.0;
        HD[ii][4] = 0.0;
        HD[ii][5] = 0.0;   
      }
      for (ii=0; ii<inp_enm->blk[bc].natom; ii++)
      {
        for(jj=0; jj<6; jj++)
        {
          for(kk=0; kk<inp_enm->blk[bc2].natom; kk++)
          {
	          //HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3][(inp_enm->blk[bc2].atom[kk]-1)*3]*inp_enm->RT[bc2][kk*3][jj];
	          //HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3][(inp_enm->blk[bc2].atom[kk]-1)*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
	          //HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3][(inp_enm->blk[bc2].atom[kk]-1)*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
            //HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+1][(inp_enm->blk[bc2].atom[kk]-1)*3]*inp_enm->RT[bc2][kk*3][jj];
            //HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+1][(inp_enm->blk[bc2].atom[kk]-1)*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
            //HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+1][(inp_enm->blk[bc2].atom[kk]-1)*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
            //HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+2][(inp_enm->blk[bc2].atom[kk]-1)*3]*inp_enm->RT[bc2][kk*3][jj];
            //HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+2][(inp_enm->blk[bc2].atom[kk]-1)*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
            //HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].atom[ii]-1)*3+2][(inp_enm->blk[bc2].atom[kk]-1)*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
	          HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3][(inp_enm->blk[bc2].sindex[kk])*3]*inp_enm->RT[bc2][kk*3][jj];
	          HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3][(inp_enm->blk[bc2].sindex[kk])*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
	          HD[ii*3][jj]   += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3][(inp_enm->blk[bc2].sindex[kk])*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
            HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+1][(inp_enm->blk[bc2].sindex[kk])*3]*inp_enm->RT[bc2][kk*3][jj];
            HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+1][(inp_enm->blk[bc2].sindex[kk])*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
            HD[ii*3+1][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+1][(inp_enm->blk[bc2].sindex[kk])*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
            HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+2][(inp_enm->blk[bc2].sindex[kk])*3]*inp_enm->RT[bc2][kk*3][jj];
            HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+2][(inp_enm->blk[bc2].sindex[kk])*3+1]*inp_enm->RT[bc2][kk*3+1][jj];
            HD[ii*3+2][jj] += inp_enm->He[(inp_enm->blk[bc].sindex[ii])*3+2][(inp_enm->blk[bc2].sindex[kk])*3+2]*inp_enm->RT[bc2][kk*3+2][jj];
	        }
        }
      }
      for (ii=0; ii<6; ii++)
      {
        for(jj=0; jj<6; jj++)
        {
          s[ii][jj] = 0.0;
          for(kk=0; kk<inp_enm->blk[bc].natom; kk++)
          {
	          s[ii][jj] += inp_enm->RT[bc][kk*3][ii]*HD[kk*3][jj];
	          s[ii][jj] += inp_enm->RT[bc][kk*3+1][ii]*HD[kk*3+1][jj];
	          s[ii][jj] += inp_enm->RT[bc][kk*3+2][ii]*HD[kk*3+2][jj];
	        }
	        Hrt[ii+indexbc][jj+indexbc2] = s[ii][jj];
        }
      }
      for (jj=0; jj<3*inp_enm->blk[bc].natom; jj++)
      {
        free(HD[jj]);
      }
      free(HD);
     
      indexbc2+=6;  
    }
    indexbc+=6;
    indexbc2 = indexbc;
  }

 //for(ii=0;ii<6*(inp_enm->nblock);ii++)
 //{
 //  for(jj=0;jj<ii;jj++)
 //  {
 //    Hrt[ii][jj] = Hrt[jj][ii];
 //  }
 //}
 
 /*Matrix traces check!*/ 
 //Hrt_trace=0;
 //for(ii=0;ii<6*inp_enm->nblock;ii++)
 //{
 //  for(jj=0;jj<6*inp_enm->nblock;jj++)
 //  {
 //    if( Hrt[ii][jj] != 0.0 )
 //    {
 //      if( ii == jj )
 //      {
 //        Hrt_trace += Hrt[ii][jj];
 //      }
 //    }
 //  }
 //}
 //
 //He_trace = 0;
 //for(ii=0;ii<3*inp_enm->sele1.nselatm;ii++)
 //{
 //  for(jj=0;jj<3*inp_enm->sele1.nselatm;jj++)
 //  {
 //    if(ii == jj)
 //    {
 //      He_trace += inp_enm->He[ii][jj];
 //    }
 //  }
 //}
 //printf("He  trace = %f\n", He_trace);  
 //printf("RTB trace = %f\n", Hrt_trace);

  inp_enm->eigen->size	 = 6*inp_enm->nblock; 
  inp_enm->eigen->inpmat = Hrt;
  
  
 //for (jj=0; jj<6*inp_enm->nblock; jj++)
 //{
 //  free(Hrt[jj]);
 //}
 //free(Hrt);
  


  return(0) ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
int Enm_Rtb2full (struct inp_enm *inp_enm, Molecule *molecule  )
{
  int ii=0, jj=0, kk=0, ll=0, index=0, ntot=0;

  //norm = calloc (inp_enm->eigen->size, sizeof(float));
 
  //???Using selean.nselatm instead of sele1.nselatm ???
    
  inp_enm->rtb2full_eigen->eigvec = calloc (inp_enm->eigen->size, sizeof(float *));
  for (ii=0;ii<inp_enm->eigen->size;ii++)
  {
    inp_enm->rtb2full_eigen->eigvec[ii] = calloc (3*inp_enm->sele1.nselatm, sizeof(float));
  }

 /*Reduced eigenvectors norm check*/   
 //for (ii=0;ii<inp_enm->eigen->size;ii++)
 //{
 //  norm[ii] = 0.0;
 //  for (jj=0;jj<inp_enm->eigen->size;jj++)
 //  {
 //    dvec = inp_enm->eigen->eigvec[ii][jj]*inp_enm->eigen->eigvec[ii][jj];
 //    norm[ii] += dvec;
 //  }
 //  norm[ii] = sqrt(norm[ii]);
 //  //printf("Norm of eigenvectors in projected coordinates: %f\n",norm[ii]);
 //}

  for( ii=0; ii<inp_enm->eigen->size; ii++ )
  { 
    ntot=0;
    index=0;
    for( jj=0; jj<inp_enm->nblock; jj++ )
    {
      for( kk=0; kk<3*inp_enm->blk[jj].natom; kk++ )
      { 
        inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][kk+ntot] = 0.0;
        for (ll=0; ll<6; ll++)
	      {
	        inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][kk+ntot] += inp_enm->RT[jj][kk][ll]*inp_enm->eigen->eigvec[inp_enm->eigen->size-ii-1][ll+index];
	      }
      }
      index += 6;
      ntot += 3*inp_enm->blk[jj].natom;
    }
  }

 /*Full eigenvectors norm check*/ 
 //dvec=0;
 //for(ii=0;ii<inp_enm->eigen->size;ii++)
 //{
 //  norm[ii] = 0.0;
 //  for(jj=0;jj<inp_enm->sele1.nselatm;jj++)
 //  {
 //    dvec = inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3]*inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3] + inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3+1]*inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3+1] + inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3+2]*inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-ii-1][jj*3+2];
 //    norm[ii] += dvec;
 //  }
 //  norm[ii] = sqrt(norm[ii]);
 // // printf("Norm of eigenvectors in Cartesian coordinates: %f\n",norm[ii]);
 //}
 
  //for (ii=0; ii<3*inp_enm->sele1.nselatm; ii++)
  //{
  //  for (jj=0;jj<inp_enm->eigen->size/2;jj++)
  //  {
  //    printf("%15.10f", inp_enm->rtb2full_eigen->eigvec[inp_enm->eigen->size-jj-1][ii]);
  //  }
  //  printf("\n");
  //}
 
 
  inp_enm->rtb2full_eigen->eigval = inp_enm->eigen->eigval;
  inp_enm->rtb2full_eigen->size   = inp_enm->eigen->size;

  return(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*CORRELATIONS of ATOMIC FLUCTUATIONS*/

int Enm_Correlation (struct inp_enm *inp_enm, _eigen *eig, Molecule *molecule)
{
 
  FILE *correlation, *fhCovMatFile;
  FILE *corr4psn;
  int   len_mode=0, size=0, ssize=0, ii=0, jj=0, kk=0, *mode;
  int   mm, nn, n1, n2;
  float eignorm1=0, eignorm2=0, euclidean=0, **cov;
  float **ppfTemp;
  char  filename1[1280], filename2[1280], filename3[1280], *tmp;
  char  cLabelA[20], cLabelB[20], cResCode1[3];
  
  if(strncmp(inp_enm->corr,"all",64)==0)
  {
    printf("\n>> Computation of Correlation from Covariance Matrix using all modes\n");
    len_mode = eig->size-6;
    mode = calloc (len_mode, sizeof(int));
    for(ii=0; ii<len_mode; ii++)
    {
      mode[ii] = ii+1;
    }
  }
  else if(strncmp(inp_enm->corr," ",64) !=0)
  {
    printf("\n>> Computation of Correlation from Covariance Matrix using modes: %s\n", inp_enm->corr );
    tmp = expandrange ( inp_enm->corr, &len_mode );
    mode = readintstring( tmp, len_mode );
    free (tmp);
  }
  
  size = inp_enm->selean.nselatm;
  ssize = size*size;
  
  cov = inp_enm->cov;
  for( ii=0; ii<ssize; ii++)
  {
    cov[0][ii] = 0.0;
  }
   
  for(ii=0;ii<size;ii++)
  {
    for(jj=0;jj<size;jj++)
    {
      eignorm1 = 0;
      eignorm2 = 0;
      for(kk=0; kk<len_mode; kk++)
      {
        eignorm1    += (eig->eigvec[eig->size-mode[kk]-6][3*ii]*eig->eigvec[eig->size-mode[kk]-6][3*ii])/eig->eigval[eig->size-mode[kk]-6];
        eignorm2    += (eig->eigvec[eig->size-mode[kk]-6][3*jj]*eig->eigvec[eig->size-mode[kk]-6][3*jj])/eig->eigval[eig->size-mode[kk]-6];
        cov[ii][jj] += (eig->eigvec[eig->size-mode[kk]-6][3*ii]*eig->eigvec[eig->size-mode[kk]-6][3*jj])/eig->eigval[eig->size-mode[kk]-6];
        eignorm1    += (eig->eigvec[eig->size-mode[kk]-6][3*ii+1]*eig->eigvec[eig->size-mode[kk]-6][3*ii+1])/eig->eigval[eig->size-mode[kk]-6];
        eignorm2    += (eig->eigvec[eig->size-mode[kk]-6][3*jj+1]*eig->eigvec[eig->size-mode[kk]-6][3*jj+1])/eig->eigval[eig->size-mode[kk]-6];
        cov[ii][jj] += (eig->eigvec[eig->size-mode[kk]-6][3*ii+1]*eig->eigvec[eig->size-mode[kk]-6][3*jj+1])/eig->eigval[eig->size-mode[kk]-6];
        eignorm1    += (eig->eigvec[eig->size-mode[kk]-6][3*ii+2]*eig->eigvec[eig->size-mode[kk]-6][3*ii+2])/eig->eigval[eig->size-mode[kk]-6];
        eignorm2    += (eig->eigvec[eig->size-mode[kk]-6][3*jj+2]*eig->eigvec[eig->size-mode[kk]-6][3*jj+2])/eig->eigval[eig->size-mode[kk]-6];
        cov[ii][jj] += (eig->eigvec[eig->size-mode[kk]-6][3*ii+2]*eig->eigvec[eig->size-mode[kk]-6][3*jj+2])/eig->eigval[eig->size-mode[kk]-6];
      }
      /*printf("%9.5f %9.5f %9.5f \n",cov[ii][jj],eignorm1, eignorm2);*/
      cov[ii][jj] = cov[ii][jj]/(sqrt(eignorm1)*sqrt(eignorm2));
    }
  }

  if (inp_enm->print==0)
  {
    memset ( filename1, '\0', sizeof(filename1));
    sprintf( filename1, "corrmat.txt");
    memset ( filename2, '\0', sizeof(filename2));
    sprintf( filename2, "corrpairs.txt");
  }
  else if (inp_enm->print==1)
  {
    memset ( filename1, '\0', sizeof(filename1));
    sprintf( filename1, "%s-corrmat_fr%d.txt", inp_enm->title,inp_enm->current_frame);
    memset ( filename2, '\0', sizeof(filename2));
    sprintf( filename2, "%s-corrpairs_fr%d.txt", inp_enm->title,inp_enm->current_frame);
  }

  // === *-* ANF *-* ================================================ //
  if (inp_enm->print==0)
  {
    memset ( filename3, '\0', sizeof(filename3));
    sprintf( filename3, "covarmatrix.txt");
  }
  else if (inp_enm->print==1)
  {
    memset ( filename3, '\0', sizeof(filename3));
    sprintf( filename3, "%s-covarmatrix_fr%d.txt", inp_enm->title,inp_enm->current_frame);
  }
  
  ppfTemp = (float **) calloc(len_mode, sizeof(float *));
  for(mm=0; mm<len_mode; ++mm)
    ppfTemp[mm] = (float *) calloc(inp_enm->selean.nselatm*3, sizeof(float));
  
  for(mm=0; mm<len_mode; ++mm)
    for(nn=0; nn<inp_enm->selean.nselatm*3; ++nn)
      ppfTemp[mm][nn] = (eig->eigval[eig->size-mode[mm]-6] * eig->eigvec[eig->size-mode[mm]-6][nn]);
  
  for(n1=0; n1<inp_enm->selean.nselatm*3; ++n1)
    for(n2=0; n2<inp_enm->selean.nselatm*3; ++n2)
      for(mm=0; mm<len_mode; ++mm)
        inp_enm->cov2[n1][n2] += (eig->eigvec[eig->size-mode[mm]-6][n1] * ppfTemp[mm][n2]);
  
  fhCovMatFile = fopen(filename3, "w");
  for(n1=0; n1<inp_enm->selean.nselatm*3; ++n1)
  {
    for(n2=0; n2<inp_enm->selean.nselatm*3; ++n2)
    {
      fprintf(fhCovMatFile,"%9.5f ",(inp_enm->cov2[n1][n2]));
    }
    fprintf(fhCovMatFile,"\n");
  }
  fclose(fhCovMatFile);
  // === *-* ANF *-* ================================================ //

  correlation = fopen(filename1,"w");
  if (inp_enm->rtb_flag==0 || inp_enm->rtblevan_flag == 0)
  {
    for(ii=0;ii<size;ii++)
    { 
     for(jj=0;jj<size;jj++)
     {
       //fprintf(correlation,"%9.5f ",cov[ii][jj]);/*Printing R values from -1 to 1*/
       fprintf(correlation,"%9.5f ",(cov[ii][jj]/2)+0.5);/*Printing R values from 0 to 1*/
       euclidean += cov[ii][jj]*cov[ii][jj];
     }
     fprintf(correlation,"\n");
    }
  }
  else if (inp_enm->rtblevan_flag == 1)
  {
    for(ii=0; ii<inp_enm->rtblevan.nselatm; ii++)
    { 
      for(jj=0; jj<inp_enm->rtblevan.nselatm; jj++)
      {
        //fprintf(correlation,"%9.5f ",cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
        fprintf(correlation,"%9.5f ",(cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]/2)+0.5);/*Printing R values from 0 to 1*/
        //euclidean += cov[inp_enm->selean.selatm[ii]-1][inp_enm->selean.selatm[jj]-1]*cov[inp_enm->selean.selatm[ii]-1][inp_enm->selean.selatm[jj]-1];
      }
      fprintf(correlation,"\n");
    }
  }

  fclose(correlation);

  corr4psn = fopen(filename2,"w");
  if (inp_enm->rtb_flag==0 || inp_enm->rtblevan_flag == 0)
  {
    for(ii=0; ii<size; ii++)
    { 
      for(jj=0; jj<size; jj++)
      {
        {
          Res3ToRes1(molecule->rawmol.restype[inp_enm->selean.selatm[ii]-1], cResCode1);
          sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_enm->selean.selatm[ii]-1], cResCode1,
                                         molecule->rawmol.resn[inp_enm->selean.selatm[ii]-1]);        
          Res3ToRes1(molecule->rawmol.restype[inp_enm->selean.selatm[jj]-1], cResCode1);
          sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_enm->selean.selatm[jj]-1], cResCode1,
                                         molecule->rawmol.resn[inp_enm->selean.selatm[jj]-1]);
        }
        // WARNING: modded to comply with raw2struct in fileio.c
        fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,cov[ii][jj]);/*Printing R values from -1 to 1*/
        //fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", ii+1, jj+1, cLabelA, cLabelB,cov[ii][jj]);/*Printing R values from -1 to 1*/
        //fprintf(corr4psn,"%9.5f ",(cov[ii][jj]/2)+0.5);/*Printing R values from 0 to 1*/
      }
    }
  }
  else if (inp_enm->rtblevan_flag == 1)
  {
    for(ii=0; ii<inp_enm->rtblevan.nselatm; ii++)
    { 
      for(jj=0; jj<inp_enm->rtblevan.nselatm; jj++)
      {
        {
          Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->rtblevan.selatm[ii]-1], cResCode1);
          sprintf(cLabelA, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->rtblevan.selatm[ii]-1], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->rtblevan.selatm[ii]-1]);        
          Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->rtblevan.selatm[jj]-1], cResCode1);
          sprintf(cLabelB, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->rtblevan.selatm[jj]-1], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->rtblevan.selatm[jj]-1]);
        }
        // WARNING: modded to comply with raw2struct in fileio.c
        //fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", ii+1, jj+1, cLabelA, cLabelB,cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
        //fprintf(corr4psn,"%9.5f ",(cov[ii][jj]/2)+0.5);/*Printing R values from 0 to 1*/
        //fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
        
        // WARNING: further modded to account for rtblevan correction (ca line 450)
        //fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
        fprintf(corr4psn," %8d   %8d   %15s   %15s   %5.2f\n", molecule->rawmol.presn[inp_enm->rtblevan_org.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan_org.selatm[jj]-1], cLabelA, cLabelB,cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
      }
    }
      
  }
  fclose(corr4psn);

  return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /*B-FACTORS*/

int Enm_Fluc (struct inp_enm *inp_enm, Molecule *molecule, _eigen *eig)
{

  FILE *bfactors;
  int len_mode=0, ii=0, jj=0, size=0, *mode;
  float **fluc, **bfact, bsumcalc=0, bsum=0, meanb=0, meanc=0,corr=0, corrf2=0, corrb2=0, temp=0, kT=0;
  char  filename[1280], *tmp; 
  # define eigencut 0.00001  

  temp=atof(inp_enm->temp); 
  if ( temp == 0)
  {
    temp=300;
  } 
  
  kT=k*temp;

  if(strncmp(inp_enm->beta,"all",64)==0)
  {
    printf("\n>> Isotropic Mean Square Displacement at %5.2fK using all modes\n", temp);
    len_mode = eig->size-6;
    mode = calloc (len_mode, sizeof(int));
    for(ii=0; ii<len_mode; ii++)
    {
      mode[ii] = ii+1;
    }
  }
  else if(strncmp(inp_enm->beta," ",64) !=0)
  {
    printf("\n>> Isotropic Mean Square Displacement at %5.2fK using %s modes\n", temp, inp_enm->beta);
    tmp = expandrange ( inp_enm->beta, &len_mode );
    mode = readintstring( tmp, len_mode );
    free (tmp);
  }
   
  size  = inp_enm->selean.nselatm;
  //ssize = size*size;
  //bfact = inp_enm->bfact;
  //fluc  = inp_enm->fluc; 

  bfact = calloc (inp_enm->selean.nselatm, sizeof (float*));
  for (ii=0;ii<inp_enm->selean.nselatm;ii++)
  {
	bfact[ii] = calloc (len_mode+1, sizeof (float));	 
  }
  fluc = calloc (3*inp_enm->selean.nselatm, sizeof (float*));
  for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
  {
    fluc[ii] = calloc (len_mode+1, sizeof (float)); 
  }


  for (ii=0;ii<inp_enm->selean.nselatm;ii++)
  {
	for (jj=0;jj<len_mode+1;jj++)
    {
      bfact[ii][jj]=0;
      fluc[(ii)*3][jj]=0;
      fluc[(ii)*3+1][jj]=0;
      fluc[(ii)*3+2][jj]=0;
    }
  }
  
  //printf("%d\n", eig->size);
  for (ii=0;ii<3*size;ii++)
  {
    for (jj=0; jj<len_mode; jj++)
    {
      if(eig->eigval[eig->size-mode[jj]-6] > eigencut)
      {
        fluc[ii][jj] = eig->eigvec[eig->size-mode[jj]-6][ii]*eig->eigvec[eig->size-mode[jj]-6][ii]/eig->eigval[eig->size-mode[jj]-6];
        fluc[ii][len_mode] += eig->eigvec[eig->size-mode[jj]-6][ii]*eig->eigvec[eig->size-mode[jj]-6][ii]/eig->eigval[eig->size-mode[jj]-6];
      }
    }
  }
  
  for (ii=0;ii<size;ii++)
  {
    for (jj=0; jj<=len_mode; jj++)
    {
	  bfact[ii][jj] = bfact[ii][jj] + (fluc[(ii)*3][jj]+fluc[(ii)*3+1][jj]+fluc[(ii)*3+2][jj]);
      bfact[ii][jj] = bfact[ii][jj]*8*pi*pi*(kT)/3;
    }
  }

  //*Computation of correlation between experimental and calculated beta factors*/
  
  if (inp_enm->rtb_flag==0 || inp_enm->rtblevan_flag == 0)
  {
    for(ii=0;ii<size;ii++)
    {
      bsumcalc += bfact[ii][len_mode];
      bsum += molecule->rawmol.bFac[inp_enm->selean.selatm[ii]-1];
      corr     += (bfact[ii][len_mode])*(molecule->rawmol.bFac[inp_enm->selean.selatm[ii]-1]);
      corrf2   += (bfact[ii][len_mode])*(bfact[ii][len_mode]);
      corrb2   += (molecule->rawmol.bFac[inp_enm->selean.selatm[ii]-1])*(molecule->rawmol.bFac[inp_enm->selean.selatm[ii]-1]);
    } 
    meanb = bsum;
    meanc = bsumcalc/size;
    meanb = meanb/size;
    corrf2=sqrt((corrf2/inp_enm->selean.nselatm)-meanc*meanc);
    corrb2=sqrt((corrb2/inp_enm->selean.nselatm)-meanb*meanb);
 
    if (corrb2 > eigencut)
    {
      corr=corr/size-meanc*meanb;
    }
    if ( corrf2 > eigencut )
    {
      corr=corr/(corrf2*corrb2);
    }  
  }
  else if (inp_enm->rtblevan_flag == 1)
  {
    for(ii=0;ii<inp_enm->rtblevan.nselatm;ii++)
    {
      bsumcalc += bfact[inp_enm->rtblevan.selatm[ii]-1][len_mode];
      bsum     += molecule->rawmol.bFac[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1];
      corr     += (bfact[inp_enm->rtblevan.selatm[ii]-1][len_mode])*(molecule->rawmol.bFac[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1]);
      corrf2   += (bfact[inp_enm->rtblevan.selatm[ii]-1][len_mode])*(bfact[inp_enm->rtblevan.selatm[ii]-1][len_mode]);
      corrb2   += (molecule->rawmol.bFac[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1])*(molecule->rawmol.bFac[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1]);
    } 
    meanb = bsum;
    meanc = bsumcalc/inp_enm->rtblevan.nselatm;
    meanb = meanb/inp_enm->rtblevan.nselatm;
    corrf2=sqrt((corrf2/inp_enm->rtblevan.nselatm)-meanc*meanc);
    corrb2=sqrt((corrb2/inp_enm->rtblevan.nselatm)-meanb*meanb);
 
    if (corrb2 > eigencut)
    {
      corr=corr/inp_enm->rtblevan.nselatm-meanc*meanb;
    }
    if ( corrf2 > eigencut )
    {
      corr=corr/(corrf2*corrb2);
    }        
  }
  
  
  if (bsum == 0)
  {
    corr=0;
  }


  if(inp_enm->print == 0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "bfactors.txt");
  }
  else if (inp_enm->print == 1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-bfactors_fr%d.txt", inp_enm->title, inp_enm->current_frame);
  } 
  bfactors = fopen (filename, "w");

  fprintf(bfactors,"#  Cross Correlation between calculated and experimental B factors CC = %9.5f\n",corr);

  if (inp_enm->rtb_flag==0 || inp_enm->rtblevan_flag == 0)
  { 
    for (ii=0;ii<size;ii++)
    {
      fprintf(bfactors,"%3d %9.5f",molecule->rawmol.resn[inp_enm->selean.selatm[ii]-1],molecule->rawmol.bFac[inp_enm->selean.selatm[ii]-1]);
      for (jj=0;jj<=len_mode;jj++)
      {
	    fprintf(bfactors," %9.5f ", bfact[ii][jj]); 
	  }
	  fprintf(bfactors,"\n");
    }    
  }
  else if (inp_enm->rtblevan_flag == 1)
  {
    for (ii=0;ii<inp_enm->rtblevan.nselatm;ii++)
    {
      fprintf(bfactors,"%3d %9.5f",molecule->rawmol.resn[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1],molecule->rawmol.bFac[inp_enm->sys_list[inp_enm->rtblevan.selatm[ii]-1][1]-1]);
      for (jj=0;jj<=len_mode;jj++)
      {
	    fprintf(bfactors," %9.5f ", bfact[inp_enm->rtblevan.selatm[ii]-1][jj]); 
	  }
	  fprintf(bfactors,"\n");
    }
    
  }
  fclose(bfactors);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /*ELASTIC NETWORK STRUCTURAL PERTURBATION */

int Enm_Perturb(struct inp_enm *inp_enm, Molecule *molecule, _eigen *eig, _dr *dr, CoorSet *trj_crd)
{ 
  FILE *perturbation;
  FILE *pertmat_file;
  int resid=0, len_mode=0, *mode, ii=0, jj=0, kk=0,  nn=0, size=0;
  int nres=0, *pert_list, pert_list_size=0, cont=0;
  float  *scalar, **Hp, **Hssp, *max, *max_pertmat;
  char  filename_plot[1280], **filename_pdb, *tmp, **pertmat_name;
  _dr *drpert; 
  Molecule *frame;
  scalar = inp_enm->scalar;
 // Hs = inp_enm->Hs;
  
  printf("\n>> Structural Perturbation Method\n");

  if (strncmp(inp_enm->perturb,"all",64) == 0)
  {
    printf (">  Calculating perturbations over all modes\n"); 
    printf (">  Calculation may take a lot... \n");
    len_mode = eig->size-6;
    mode = calloc (len_mode, sizeof(int));
    for(ii=0; ii<len_mode; ii++)
    {
      mode[ii] = ii+1;
    }
  }
  else if (strncmp(inp_enm->perturb," ",64) != 0)
  {
    printf (">  Calculating perturbations using modes: %s\n", inp_enm->perturb);
    tmp = expandrange ( inp_enm->perturb, &len_mode );
    mode = readintstring( tmp, len_mode );
    free (tmp);
  }

  if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0  )
  {
    nres  = inp_enm->nblock;
    size  = 3*inp_enm->nblock;
    //ssize = 9*inp_enm->selean.nselatm*inp_enm->selean.nselatm;
  }
  else
  {
    nres  = inp_enm->selean.nselatm;
    size  = 3*inp_enm->selean.nselatm;
    //ssize = size*size;
  }

  drpert = calloc (inp_enm->nint, sizeof( _dr));    
  inp_enm->delta_omega = calloc (size, sizeof(float *));
  for ( ii=0; ii < size; ii++)
  {
    inp_enm->delta_omega[ii] = calloc ( len_mode, sizeof(float));
  } 
  inp_enm->delta_omega_plot = calloc (nres, sizeof(float *));
  for ( ii=0; ii < nres; ii++)
  {
    inp_enm->delta_omega_plot[ii] = calloc (len_mode+1, sizeof(float));
  } 
  max = calloc ( len_mode+1, sizeof(float));  
  
  filename_pdb = calloc (len_mode+1, sizeof(char*));
  for (ii=0; ii<=len_mode; ii++)
  {
    filename_pdb[ii] = calloc (1280, sizeof(char));
  }
  
  pert_list = calloc (inp_enm->nint, sizeof(int)); 
  
  if (inp_enm->pertmat_flag == 1)
  {
    inp_enm->pertmat = calloc(nres, sizeof(float **) );
    for(ii=0;ii<nres;ii++)
    {
      inp_enm->pertmat[ii] = calloc(nres, sizeof(float*));
      for (jj=0;jj<nres;jj++)
      {
        inp_enm->pertmat[ii][jj] = calloc(len_mode+1,sizeof(float)); 
      } 
    }  
    pertmat_name = calloc (len_mode+1, sizeof(char*));
    for (ii=0; ii<=len_mode; ii++)
    {
      pertmat_name[ii] = calloc (1280, sizeof(char));
    }
    max_pertmat = calloc ( len_mode+1, sizeof(float) );
  }


  
  for( kk=0; kk<inp_enm->nint; kk++)
  {
    drpert[kk].drdx  = dr[kk].drdx;  
    drpert[kk].drdy  = dr[kk].drdy;  
    drpert[kk].drdz  = dr[kk].drdz;  
    drpert[kk].part1 = dr[kk].part1;
    drpert[kk].part2 = dr[kk].part2;
    drpert[kk].cost  = dr[kk].cost;
  }


  resid=1;
  printf("Perturbing residue nr.");fflush(stdout); 
  while ( resid <= nres )
  { 
    printf(" %d,", resid);fflush(stdout);
    for(kk=0; kk<inp_enm->nint; kk++)
    {
      pert_list[kk] = 0;       
      drpert[kk].cost = dr[kk].cost;
    }   
    pert_list_size = 0;

    if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
    { 
      for(ii=0; ii<inp_enm->nint; ii++)
      {
        for ( jj=0; jj<inp_enm->blk[resid-1].natom; jj++ )
        {
          if(drpert[ii].part1 == inp_enm->blk[resid-1].atom[jj] || drpert[ii].part2 == inp_enm->blk[resid-1].atom[jj] )    
          { 
            drpert[ii].cost = delta*dr[ii].cost;
            pert_list[pert_list_size] = ii; 
            pert_list_size++;
          }
        }
      }
    }
    else if ( inp_enm->vsa_flag == 1 )
    {
      for(ii=0; ii<inp_enm->nint; ii++)
      {        
        if( drpert[ii].part1-1 == inp_enm->subsys_list[resid-1][0] || drpert[ii].part2-1 == inp_enm->subsys_list[resid-1][0] )    
        {   
          drpert[ii].cost = delta*dr[ii].cost;
          pert_list[pert_list_size] = ii; 
          pert_list_size++;
        }
        else if(inp_enm->subsys_list[resid-1][0] != drpert[ii].part1-1 && inp_enm->subsys_list[resid-1][0] != drpert[ii].part2-1)
        {
          drpert[ii].cost = dr[ii].cost;
        }
      }        
    }
    else if (inp_enm->rtb_flag == 0 && inp_enm->vsa_flag == 0)
    {
      for(ii=0;ii<inp_enm->nint;ii++)
      {        
        if(drpert[ii].part1 == resid || drpert[ii].part2 == resid)    
        {   
          drpert[ii].cost = delta*dr[ii].cost;
          pert_list[pert_list_size] = ii; 
          pert_list_size++;
        }
        else if(resid != drpert[ii].part1 && resid != drpert[ii].part2)
        {
          drpert[ii].cost = dr[ii].cost;
        }
      }
    }
   
    if  ( inp_enm->vsa_flag == 1 )
    {
      Hssp = inp_enm->Hss;
      Build_Perthessian (inp_enm, inp_enm->Hp, drpert, 3*inp_enm->sele1.nselatm, pert_list, pert_list_size);      
      Enm_GetVsaMat (inp_enm, inp_enm->Hp, Hssp);
      Hp = 0;
      Hp = Hssp;
    }
    else
    {
      Build_Perthessian (inp_enm, inp_enm->Hp, drpert, 3*inp_enm->selean.nselatm, pert_list, pert_list_size);
      Hp = inp_enm->Hp;
    }
  
    for( nn=0; nn<len_mode; nn++ )
    {
      for( ii=0; ii<3*inp_enm->selean.nselatm; ii++)
      {
        inp_enm->scalar[ii] = 0.0;   
      }
          
      for( jj=0; jj<3*inp_enm->selean.nselatm; jj++)
      {
        cont=0;
        for(ii=0;ii<jj;ii++)
        {
          scalar[cont] += eig->eigvec[eig->size-(mode[nn]+6)][jj] * Hp[ii][jj];
          cont++;  
        }
        for( ii=jj; ii<3*inp_enm->selean.nselatm; ii++)
        { 
          scalar[cont] += eig->eigvec[eig->size-(mode[nn]+6)][jj] * Hp[jj][ii];
          cont++;
        }
      }
            
      if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->pertmat_flag==1)
      {
        for (ii=0;ii<inp_enm->selean.nselatm;ii++)
        {
          //printf("%d\n",inp_enm->atm2blk[ii]);
          inp_enm->pertmat[resid-1][inp_enm->atm2blk[ii]][nn] = (inp_enm->scalar[3*ii]+scalar[3*ii+1]+scalar[3*ii+2]) * (nres/(2*eig->eigval[eig->size-(mode[nn]+6)]));
          inp_enm->pertmat[resid-1][inp_enm->atm2blk[ii]][len_mode] += (inp_enm->scalar[3*ii]+scalar[3*ii+1]+scalar[3*ii+2]) * (nres/(2*eig->eigval[eig->size-(mode[nn]+6)]));
        }
      }
      else if( inp_enm->pertmat_flag==1 )
      {
        for (ii=0;ii<inp_enm->selean.nselatm;ii++)
        {        
          inp_enm->pertmat[resid-1][ii][nn] = (inp_enm->scalar[3*ii]+scalar[3*ii+1]+scalar[3*ii+2]) * (nres/(2*eig->eigval[eig->size-(mode[nn]+6)]));
          inp_enm->pertmat[resid-1][ii][len_mode] += (inp_enm->scalar[3*ii]+scalar[3*ii+1]+scalar[3*ii+2]) * (nres/(2*eig->eigval[eig->size-(mode[nn]+6)]));
        }
      }
      
      for(ii=0;ii<inp_enm->selean.nselatm;ii++)
      {
        inp_enm->delta_omega[3*(resid-1)][nn]   += inp_enm->scalar[3*ii]   * eig->eigvec[eig->size-(mode[nn]+6)][3*ii];
        inp_enm->delta_omega[3*(resid-1)+1][nn] += inp_enm->scalar[3*ii+1] * eig->eigvec[eig->size-(mode[nn]+6)][3*ii+1];
        inp_enm->delta_omega[3*(resid-1)+2][nn] += inp_enm->scalar[3*ii+2] * eig->eigvec[eig->size-(mode[nn]+6)][3*ii+2];
      }
     
      if (inp_enm->perturb_norm == 1)
      {
        inp_enm->delta_omega[3*(resid-1)][nn]   = (nres/(2*eig->eigval[eig->size-(mode[nn]+6)])) * inp_enm->delta_omega[3*(resid-1)][nn];
        inp_enm->delta_omega[3*(resid-1)+1][nn] = (nres/(2*eig->eigval[eig->size-(mode[nn]+6)])) * inp_enm->delta_omega[3*(resid-1)+1][nn];
        inp_enm->delta_omega[3*(resid-1)+2][nn] = (nres/(2*eig->eigval[eig->size-(mode[nn]+6)])) * inp_enm->delta_omega[3*(resid-1)+2][nn];
      }
      else
      {
        inp_enm->delta_omega[3*(resid-1)][nn]   = inp_enm->delta_omega[3*(resid-1)][nn];
        inp_enm->delta_omega[3*(resid-1)+1][nn] = inp_enm->delta_omega[3*(resid-1)+1][nn];
        inp_enm->delta_omega[3*(resid-1)+2][nn] = inp_enm->delta_omega[3*(resid-1)+2][nn];
      }
    }
    resid++;
  }
  
  if (inp_enm->pertmat_flag==1)
  {
    for ( ii=0; ii<=len_mode; ii++)
    {
      max_pertmat[ii] = 0;
      for (jj=0;jj<nres;jj++)
      {
        for (kk=0;kk<nres;kk++)
        {
          if (inp_enm->pertmat[jj][kk][ii] > max_pertmat[ii])
          {
            max_pertmat[ii] = inp_enm->pertmat[jj][kk][ii];
          }
        }
      }
    }
  }

  
  if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {  
    for(ii=0;ii<nres;ii++)
    {
      for(nn=0;nn<len_mode;nn++)
      {
        inp_enm->delta_omega_plot[ii][nn] = inp_enm->delta_omega[3*ii][nn]+inp_enm->delta_omega[3*ii+1][nn]+inp_enm->delta_omega[3*ii+2][nn];
        inp_enm->delta_omega_plot[ii][len_mode] += inp_enm->delta_omega[3*ii][nn]+inp_enm->delta_omega[3*ii+1][nn]+inp_enm->delta_omega[3*ii+2][nn];
      }
    }   
    for ( nn=0; nn<=len_mode; nn++ )
    {
      max[nn]=0.0;
      for ( ii=0;ii<nres;ii++ )
      {
        if( inp_enm->delta_omega_plot[ii][nn] > max[nn] )
        {
          max[nn] = inp_enm->delta_omega_plot[ii][nn];
        }
      }
    }
  }
  else
  {
    for(ii=0;ii<inp_enm->selean.nselatm;ii++)
    {
//      inp_enm->delta_omega_plot[ii][0] = molecule->rawmol.resn[inp_enm->selean.selatm[ii]-1];
      for(nn=0;nn<len_mode;nn++)
      {
        inp_enm->delta_omega_plot[ii][nn] = inp_enm->delta_omega[3*ii][nn]+inp_enm->delta_omega[3*ii+1][nn]+inp_enm->delta_omega[3*ii+2][nn];
        inp_enm->delta_omega_plot[ii][len_mode] += inp_enm->delta_omega[3*ii][nn]+inp_enm->delta_omega[3*ii+1][nn]+inp_enm->delta_omega[3*ii+2][nn];
      }
    }   
    for ( nn=0; nn<=len_mode; nn++ )
    {
      max[nn]=0.0;
      for ( ii=0;ii<inp_enm->selean.nselatm;ii++ )
      {
        if( inp_enm->delta_omega_plot[ii][nn] > max[nn] )
        {
          max[nn] = inp_enm->delta_omega_plot[ii][nn];
        }
      }
    }    
  }
 
   /*Printing output*/
  
  if(inp_enm->print == 0)
  {
    memset ( filename_plot, '\0', sizeof(filename_plot));
    sprintf( filename_plot, "perturbation.txt");
    for(ii=0; ii<len_mode; ii++)
    {
      memset ( filename_pdb[ii], '\0', sizeof(filename_pdb));
      sprintf( filename_pdb[ii], "perturbation_eig%d.pdb", mode[ii]);
      if (inp_enm->pertmat_flag==1)
      {
        memset ( pertmat_name[ii], '\0', sizeof(pertmat_name));
        sprintf( pertmat_name[ii], "pertmat_eig%d.txt", mode[ii]);
      }
    }
    memset ( filename_pdb[len_mode], '\0', sizeof(filename_pdb));
    sprintf( filename_pdb[len_mode], "perturbation_eig%d_comb.txt", mode[len_mode]);
    if(inp_enm->pertmat_flag==1)
    {
      memset ( pertmat_name[len_mode], '\0', sizeof(pertmat_name));
      sprintf( pertmat_name[len_mode], "pertmat_eig%d_comb.txt", mode[len_mode]);
    }
  }
  else if (inp_enm->print == 1)
  {
    memset ( filename_plot, '\0', sizeof(filename_plot));
    sprintf( filename_plot, "%s-perturbation_fr%d.txt", inp_enm->title, inp_enm->current_frame);
    for(ii=0; ii<len_mode; ii++)
    {
      memset ( filename_pdb[ii], '\0', sizeof(filename_pdb));
      sprintf( filename_pdb[ii], "%s-perturbation_eig%d_fr%d.pdb", inp_enm->title, mode[ii], inp_enm->current_frame);    
      if(inp_enm->pertmat_flag==1)
      {
        memset ( pertmat_name[ii], '\0', sizeof(pertmat_name));
        sprintf( pertmat_name[ii], "%s-pertmat_eig%d_fr%d.txt", inp_enm->title, mode[ii], inp_enm->current_frame);
      }
    }
    memset ( filename_pdb[len_mode], '\0', sizeof(filename_pdb));
    sprintf( filename_pdb[len_mode], "%s-perturbation_eig%d_comb_fr%d.pdb", inp_enm->title, mode[ii], inp_enm->current_frame);    
    if(inp_enm->pertmat_flag==1)
    {
      memset ( pertmat_name[len_mode], '\0', sizeof(pertmat_name));
      sprintf( pertmat_name[len_mode], "%s-pertmat_eig%d_comb_fr%d.txt", inp_enm->title, len_mode, inp_enm->current_frame);
    }
  } 
  
  /*Writing file with delta omega values for plotting*/
  perturbation = fopen(filename_plot,"w");    
  for(ii=0; ii<nres; ii++)
  {
    for(nn=0;nn<=len_mode;nn++)
    {
      fprintf(perturbation,"%15.11f",inp_enm->delta_omega_plot[ii][nn]);
    }
    fprintf(perturbation,"\n");
  } 
  fclose(perturbation);
  printf("\n>  Done!\n");  

/*Writing perturbation matrices*/
  if (inp_enm->pertmat_flag==1)
  {
    for (nn=0; nn<len_mode+1; nn++)
    {
	    pertmat_file = fopen(pertmat_name[nn], "w");
	    for (ii=0; ii<nres; ii++)
	    {
	      for (jj=0; jj<nres; jj++)
	      {
	        //fprintf(distance," %11.8f ", pertmat[ii][jj][nn]);  
	        fprintf(pertmat_file," %11.8f ", (inp_enm->pertmat[ii][jj][nn]/max_pertmat[nn])/2+0.5);  
	      }
        fprintf(pertmat_file, "\n"); 	
      }
      fclose(pertmat_file);  
    }
  }

  /*Writing pdb files with delta omega values in B-factor field*/
  frame = inp_enm->molecule;
  if ( inp_enm->mol2_flag == 1 )
  {
    for (ii=0; ii<molecule->nato; ii++)
    {
      frame->coor.xcoor[ii] = inp_enm->moving_sys[ii][0];
      frame->coor.ycoor[ii] = inp_enm->moving_sys[ii][1];
      frame->coor.zcoor[ii] = inp_enm->moving_sys[ii][2];
    }
  }
  else
  {
    for (ii=0; ii<molecule->nato; ii++)
    {
      frame->coor.xcoor[ii] = trj_crd->xcoor[ii];
      frame->coor.ycoor[ii] = trj_crd->ycoor[ii];
      frame->coor.zcoor[ii] = trj_crd->zcoor[ii];
    }
  }
    
  if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {
    for( jj=0; jj< molecule->nato; jj++)
    {
      molecule->rawmol.bFac[jj] = 0.0;
    } 
    for (nn=0;nn<=len_mode;nn++)
    {
      for( jj=0; jj<nres; jj++)
      {
        for (ii=0;ii<inp_enm->blk[jj].natom;ii++)
        {
          molecule->rawmol.bFac[inp_enm->blk[jj].atom[ii]-1] = 100 * inp_enm->delta_omega_plot[jj][nn]/max[nn] ;//fabs(100*(A[ii+jj*(inp_pca->sele.nselatm*3)]));
        }
      }  
      WritePdb_unstr( filename_pdb[nn], frame );
    }
  }
  else 
  {
    for( jj=0; jj< molecule->nato; jj++)
    {
      molecule->rawmol.bFac[jj] = 0.0;
    }
 
    for (nn=0; nn<=len_mode; nn++)
    {
      for( jj=0; jj<inp_enm->selean.nselatm; jj++)
      {
        molecule->rawmol.bFac[inp_enm->selean.selatm[jj]-1] = 100 * inp_enm->delta_omega_plot[jj][nn]/max[nn] ;//fabs(100*(A[ii+jj*(inp_pca->sele.nselatm*3)]));
      }  
      WritePdb_unstr( filename_pdb[nn], frame );      
    }
  }

 
  for ( ii=0; ii < size; ii++)
  {
    free(inp_enm->delta_omega[ii]);
  }
  free(inp_enm->delta_omega);

  for ( ii=0; ii < nres; ii++)
  {
    free(inp_enm->delta_omega_plot[ii]);
  }
  free(inp_enm->delta_omega_plot);
   
  free(drpert);
  
  for(ii=0;ii<len_mode;ii++)
  {
    free (filename_pdb[ii]);
  }
  free (filename_pdb);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /*COMPARISON BETWEEN DEFORMATION MODES AND OBSERVED CONFOMATIONAL CHANGE*/
   /*POSITION DIFFERENCE VECTOR BETWEEN THE 1st AND 2nd PDB*/
 
   
int Enm_Compare (struct inp_enm *inp_enm, Molecule *molecule, CoorSet *trj_crd, _eigen *eig )
{

  FILE *difference_vector;
  FILE *overlap_list;
  FILE *CSO;
  FILE *displacement;
  FILE *vec_compare;
  FILE *deforVectFile;

  int ii=0, jj=0, nn=0, meaningfulmodes=0, *mode, len_mode=0, kk=0;
  float *deviation, *norm, rmsd=0, *overlap, **difference_matrix,  *distances;
  float **mat_overlap, cumulative=0, rmsip=0;
  char  filename[1280],filename1[1280], filename2[1280], deforVectFileName[1280], *tmp;
 
  Molecule *mol2;
  mol2 = inp_enm->molecule2 ;
 
  //if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0  )
  //{
    //nres  = inp_enm->nblock;
    //size  = 3*inp_enm->nblock;
    //ssize = 9*inp_enm->selean.nselatm*inp_enm->selean.nselatm;
  //}
  //else
  //{ 
    //size  = 3*inp_enm->selean.nselatm;
    //ssize = size*size;
  //}
  
  distances = calloc(3, sizeof(float));
 
  if(inp_enm->vecfile2_flag != 1)
  {  
    overlap = inp_enm->overlap;
    norm = inp_enm->norm;
    deviation = inp_enm->deviation;
    for ( ii=0; ii<inp_enm->selean.nselatm; ii++ )
    {
      deviation[3*ii]   = 0.0;
      deviation[3*ii+1] = 0.0;
      deviation[3*ii+2] = 0.0;
    }
  
    if( trj_crd->pbc == NULL )
    {
      for (ii=0;ii<inp_enm->selean.nselatm;ii++)
      {
        deviation[3*ii]   = (inp_enm->moving_sys[inp_enm->selean.selatm[ii]-1][0] - mol2->coor.xcoor[inp_enm->sele2.selatm[ii]-1]);
        deviation[3*ii+1] = (inp_enm->moving_sys[inp_enm->selean.selatm[ii]-1][1] - mol2->coor.ycoor[inp_enm->sele2.selatm[ii]-1]);
        deviation[3*ii+2] = (inp_enm->moving_sys[inp_enm->selean.selatm[ii]-1][2] - mol2->coor.zcoor[inp_enm->sele2.selatm[ii]-1]); 
      }
    }
    else
    {
      for(ii=0;ii<inp_enm->selean.nselatm;ii++)
      {
        DistanceAxes( distances, inp_enm->coord[ii][0], inp_enm->coord[ii][1], inp_enm->coord[ii][2], mol2->coor.xcoor[inp_enm->sele2.selatm[ii]-1], mol2->coor.ycoor[inp_enm->sele2.selatm[ii]-1], mol2->coor.zcoor[inp_enm->sele2.selatm[ii]-1], trj_crd->pbc);
        deviation[3*ii]   = distances[0];
        deviation[3*ii+1] = distances[1];
        deviation[3*ii+2] = distances[2];
      }
    }
  
    displacement = fopen("displacement.txt", "w");
    for ( ii=0; ii<3*inp_enm->selean.nselatm; ii++ )
    {
      rmsd += deviation[ii]*deviation[ii];
      fprintf (displacement,"%9.5f\n",deviation[ii]);
    }
    printf("RMSD = %9.5f\n", sqrt(rmsd/inp_enm->selean.nselatm) );
  
    /*========================================================================*/
    /*-Overlap between the position difference vector and the ith eigvector-*/
    /*========================================================================*/  
  
    if (inp_enm->print == 0)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "overlap.txt");
      memset ( filename1, '\0', sizeof(filename1));
      sprintf( filename1, "cso.txt");
      memset ( filename2, '\0', sizeof(filename2));
      sprintf( filename2, "overlap_list.txt");
    }
    else if (inp_enm->print == 1)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "%s-overlap_fr%d.txt", inp_enm->title,inp_enm->current_frame);
      memset ( filename1, '\0', sizeof(filename1));
      sprintf( filename1, "%s-cso_fr%d.txt", inp_enm->title,inp_enm->current_frame);     
      memset ( filename2, '\0', sizeof(filename2));
      sprintf( filename2, "%s-overlap_list_fr%d.txt", inp_enm->title,inp_enm->current_frame);
    }
  
    for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
    {
      overlap[ii]=0.0;
      norm[ii]=0.0;
    }

    for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
    {
      overlap[ii]=0.0;
      norm[ii]=0.0;
    }
    difference_vector = fopen(filename,"w");
    overlap_list = fopen(filename2,"w");
    printf("\n>  Mode	 Overlap \n\n");  
    //fprintf(difference_vector,"\n#>  Mode	 Overlap \n\n");  
    meaningfulmodes=0;
    for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
    {
      for (jj=0;jj<3*inp_enm->selean.nselatm;jj++)
      {
        norm[ii] += (eig->eigvec[ii][jj]*eig->eigvec[ii][jj]);   
        overlap[ii] += (eig->eigvec[ii][jj]*deviation[jj]);
      }
      norm[ii]=sqrt(norm[ii])*sqrt(rmsd); 
      fprintf(overlap_list,"Mode%d:%3.2f  \n",eig->size-ii-1, overlap[ii]/norm[ii]);   
      if ((fabs(overlap[ii])/norm[ii]) > 0.2)
      {
        meaningfulmodes++;
        printf("   %3d      %9.5f \n",(eig->size-ii-1),overlap[ii]/norm[ii]); 
        //fprintf(difference_vector,"#   %3d      %9.5f \n",(eig->size-ii-1),overlap[ii]/norm[ii]); 
        fprintf(difference_vector,"Mode%d:%3.2f  ",eig->size-ii-1, overlap[ii]/norm[ii]);   
      }
    }
    fprintf(difference_vector," Avg_10       Disp.        Res.num.\n\n");

    // === *-* ANF *-* =================================================
    memset ( deforVectFileName, '\0', sizeof(deforVectFileName));
    sprintf( deforVectFileName, "%s-deforvect.txt", inp_enm->title);
    deforVectFile = fopen(deforVectFileName,"w");
    for (jj=0;jj<3*inp_enm->selean.nselatm;jj++)
    {
      fprintf(deforVectFile, "%.2f\n", deviation[jj]); fflush(stdout);
    }
    fclose(deforVectFile);
    printf("RMSD -> %.f\n", rmsd);
    // === *-* ANF *-* =================================================
    
    /*********************************/	  
    /****Cumulative Square Overlap****/
    /*********************************/
  
    CSO = fopen (filename1, "w");
    cumulative=0;
    for(ii=1;ii<=eig->size;ii++)
    {
      cumulative += (overlap[eig->size-ii-6]*overlap[eig->size-ii-6])/(norm[eig->size-ii-6]*norm[eig->size-ii-6]);	 
      fprintf(CSO,"%9.5f\n", sqrt(cumulative));
    }
  
    /***************************************************************************/     
    /****Printing the contribution of each mode to the conformational change****/ 
    /***************************************************************************/     

    difference_matrix = calloc(inp_enm->selean.nselatm, sizeof(float*));
    for(ii=0;ii<inp_enm->selean.nselatm;ii++)
    {
      difference_matrix[ii] = calloc((meaningfulmodes+3), sizeof(float));
    }

    for(ii=0;ii<inp_enm->selean.nselatm;ii++)
    {
      for(jj=0;jj<meaningfulmodes+2;jj++)
      {
        difference_matrix[ii][jj]=0;
      }
    }
  
    for(jj=0;jj<inp_enm->selean.nselatm;jj++)
    {
      difference_matrix[jj][meaningfulmodes+2] = molecule->rawmol.resn[inp_enm->selean.selatm[jj]-1];
      difference_matrix[jj][meaningfulmodes+1] = sqrt(deviation[3*jj]*deviation[3*jj]+deviation[3*jj+1]*deviation[3*jj+1]+deviation[3*jj+2]*deviation[3*jj+2]);
    }
  
    for (ii=0;ii<inp_enm->selean.nselatm;ii++)
    {
      for (jj=0;jj<10;jj++)
      {
        difference_matrix[ii][meaningfulmodes] = eig->eigvec[eig->size-jj-7][3*ii]*eig->eigvec[eig->size-jj-7][3*ii]+eig->eigvec[eig->size-jj-7][3*ii+1]*eig->eigvec[eig->size-jj-7][3*ii+1]+eig->eigvec[eig->size-jj-7][3*ii+2]*eig->eigvec[eig->size-jj-7][3*ii+2];
      }
      difference_matrix[ii][meaningfulmodes] = sqrt(difference_matrix[ii][meaningfulmodes]*rmsd);
    }
    
    nn=0;
    for (ii=0;ii<3*inp_enm->selean.nselatm;ii++)
    {
      if ((fabs(overlap[ii])/norm[ii]) > 0.2)
      {
        for (jj=0;jj<inp_enm->selean.nselatm;jj++)
        {
          difference_matrix[jj][nn] = sqrt(eig->eigvec[ii][3*jj]*eig->eigvec[ii][3*jj]*rmsd+eig->eigvec[ii][3*jj+1]*eig->eigvec[ii][3*jj+1]*rmsd+eig->eigvec[ii][3*jj+2]*eig->eigvec[ii][3*jj+2]*rmsd);
        }
        nn++;
      }
    }
    for(ii=0;ii<inp_enm->selean.nselatm;ii++)
    {
      for(jj=0;jj<meaningfulmodes+3;jj++)
      {
        fprintf(difference_vector," %-9.5f   ", difference_matrix[ii][jj]);
      }
      fprintf(difference_vector,"\n");
    }

    for( ii=0; ii<inp_enm->selean.nselatm; ii++ )
    {
      free( difference_matrix[ii]);
    }
    free( difference_matrix );
    fclose(difference_vector);
    fclose(CSO);
    fclose(displacement);
  }

///*Comparing two eigenvector sets*/
 
  if(inp_enm->vecfile2_flag==1)
  {
    if (inp_enm->print == 0)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "vec_compare.txt");     
    }
    else if (inp_enm->print == 1)
    {
      memset ( filename, '\0', sizeof(filename));
      sprintf( filename, "%s-vec_compare_fr%d.txt", inp_enm->title,inp_enm->current_frame);     
    }  
    
    mat_overlap = inp_enm->mat_overlap;
    for(ii=0;ii<3*inp_enm->selean.nselatm;ii++)
    {
      for(jj=0;jj<3*inp_enm->selean.nselatm;jj++)
      {
        mat_overlap[ii][jj] = 0.0;
      }
    }
    
    if (strncmp(inp_enm->nmodes,"essential",64)==0)
    {
      len_mode = inp_enm->essential;
      mode = calloc (len_mode, sizeof(int));
      for(ii=0; ii<len_mode; ii++)
      {
        mode[ii] = ii+1;
      }
    }
    else
    {
      tmp  = expandrange ( inp_enm->nmodes, &len_mode );
      mode = readintstring( tmp, len_mode );
      free (tmp);
    }
    if( inp_enm->matvec_size < len_mode )
    {
      printf("Error: too many vectors requested!\n");
      exit(0);
    }
    
    vec_compare = fopen(filename, "w");
    printf ("NMODES=%s\n",inp_enm->nmodes);
    printf ("Overlap between ENM and VECFILE subspace\n");
    fprintf(vec_compare,"%s : %s\n", inp_enm->title, inp_enm->vecfile2);
    fprintf(vec_compare,"\n>  Mode	 Overlap \n\n");  
    //printf("\n>  Mode	 Overlap \n\n");  
    meaningfulmodes=0;
    

    for (ii=0;ii<len_mode;ii++)
    {
      for (jj=0;jj<len_mode;jj++)
      {
        for (kk=0;kk<3*inp_enm->selean.nselatm;kk++)
        {
          mat_overlap[ii][jj] += (eig->eigvec[eig->size-mode[ii]-6][kk]*inp_enm->vecmat2[kk][mode[jj]-1]);
          //printf("%f\n", inp_enm->vecmat2[kk][mode[jj]]);
        }
      }
    }
  
    for(ii=0;ii<len_mode;ii++)
    {
      for(jj=0;jj<len_mode;jj++)
      {
        if ((fabs(mat_overlap[ii][jj])) > 0.2)
        {
          meaningfulmodes++;
          fprintf(vec_compare,"%3d:%3d     %9.5f \n",mode[ii]+5,mode[jj],mat_overlap[ii][jj]); 
          //printf("%3d:%3d     %9.5f \n",ii+1,jj+1,mat_overlap[ii][jj]); 
        }
      }
    }
 
    rmsip = 0;  
    for(ii=0;ii<len_mode;ii++)
    {
      cumulative = 0;
      for(jj=0;jj<len_mode;jj++)
      {
        cumulative += mat_overlap[jj][ii]*mat_overlap[jj][ii];/*CSO between selected normal modes and jjth file vector*/
        rmsip += mat_overlap[jj][ii]*mat_overlap[jj][ii];/*RMSIP between normal modes and file vector sets*/
      }
      fprintf(vec_compare,"CSO = %3s:%3d     %9.5f \n",inp_enm->nmodes,mode[ii],sqrt(cumulative)); 
    }
    fprintf(vec_compare,"RMSIP = %3s:%3s     %9.5f \n",inp_enm->nmodes,inp_enm->nmodes,sqrt(rmsip/len_mode));
    fclose (vec_compare);
  }
  
  //free(distances);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /*DEFORMATION ENERGY*/

int Enm_Defen (struct inp_enm *inp_enm, Molecule *molecule, _eigen *eig, CoorSet *trj_crd)
{
  FILE *defen;
  FILE *distance;
  int len_mode=0, ii=0, jj=0, kk=0, size=0, *mode, nres=0, nn=0 ;
  float **deform, ***distmat, *max, *max_distmat, d1=0, normdist=0;
  char  filename[1280], **distmat_name, **distpairs_name, *tmp, **filename_pdb; 
  char     cLabelA[20], cLabelB[20], cResCode1[3];
  
  Molecule *frame;
 
  if(strncmp( inp_enm->defen,"all",64)==0)
  {
    printf("\n>> Computation of Deformation Energy using all modes\n");
    len_mode = eig->size-6;
    mode = calloc (len_mode, sizeof(int));
    for(ii=0; ii<len_mode; ii++)
    {
      mode[ii] = ii+1;
    }
  }
  else if(strncmp(inp_enm->defen," ",64) !=0)
  {
    printf("\n>> Computation of Deformation Energy using %s modes\n", inp_enm->defen);
    tmp  = expandrange ( inp_enm->defen, &len_mode );
    mode = readintstring( tmp, len_mode );
    free (tmp);
  }
  
  if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {
    nres = inp_enm->nblock;
    size = inp_enm->nblock;
  }
  //else if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->rtblevan_flag == 1 )
  //{
  //  nres = inp_enm->nblock;
  //  size = inp_enm->rtblevan.nselatm;
  //}
  else
  {  
    nres = inp_enm->selean.nselatm;
    size = inp_enm->selean.nselatm;
  }
  
  deform = calloc(size, sizeof(float*));
  for (ii=0;ii<nres;ii++)
  {
    deform[ii] = calloc (len_mode+1, sizeof(float));
  }

  filename_pdb = calloc (len_mode+1, sizeof(char*));
  for (ii=0; ii<=len_mode; ii++)
  {
    filename_pdb[ii] = calloc (1280, sizeof(char));
  }
    
  max = calloc ( len_mode+1, sizeof(float) );

  for (ii=0; ii<nres; ii++)
  {
    for(jj=0; jj<len_mode+1; jj++)
    {
      deform[ii][jj]=0;
    }
  }
   
  if( inp_enm->distmat_flag == 1 )
  {
    distmat = calloc (nres, sizeof (float**));
    for (ii=0;ii<nres;ii++)
    {
      distmat[ii] = calloc (nres, sizeof(float*));
      for(jj=0; jj<nres; jj++)
      {
        distmat[ii][jj] = calloc (len_mode+1, sizeof(float));  
      }
    }  
    distmat_name = calloc (len_mode+1, sizeof(char*));
    distpairs_name = calloc (len_mode+1, sizeof(char*));
    for (ii=0; ii<=len_mode; ii++)
    {
      distpairs_name[ii] = calloc (1280, sizeof(char));
      distmat_name[ii] = calloc (1280, sizeof(char));
    }
    max_distmat = calloc ( len_mode+1, sizeof(float) );
  }
  
  for(nn=0; nn<len_mode;nn++)
  {
    //printf("mode %d\n", mode[nn]);
    if (inp_enm->vsa_flag == 0 && inp_enm->rtb_flag == 0)
    {
      for (kk=0; kk<inp_enm->nint; kk++) 
      {
        d1=0;      
        d1 =  pow(((inp_enm->dr[kk].drdx*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3]   - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3]), 2)  +
              pow(((inp_enm->dr[kk].drdy*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3+1] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3+1]), 2)+
              pow(((inp_enm->dr[kk].drdz*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3+2] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3+2]), 2);
              
        d1 = sqrt(d1);
      
        if(inp_enm->distmat_flag==1)
        {
          distmat[(inp_enm->dr[kk].part1-1)][(inp_enm->dr[kk].part2-1)][nn] = (d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist);
          distmat[(inp_enm->dr[kk].part2-1)][(inp_enm->dr[kk].part1-1)][nn] = (d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist);
        }
        deform[(inp_enm->dr[kk].part1-1)][nn] += (0.5*inp_enm->dr[kk].cost*((d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
        deform[(inp_enm->dr[kk].part2-1)][nn] += (0.5*inp_enm->dr[kk].cost*((d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
      }      
    }
    else if (inp_enm->vsa_flag == 1 && inp_enm->rtb_flag == 0)
    {
      for (kk=0; kk<inp_enm->nintvsa; kk++) 
      {
        d1=0;      
        d1 =  pow(((inp_enm->drvsa[kk].drdx*inp_enm->drvsa[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part2-1)*3]   - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part1-1)*3]), 2)  +
              pow(((inp_enm->drvsa[kk].drdy*inp_enm->drvsa[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part2-1)*3+1] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part1-1)*3+1]), 2)+
              pow(((inp_enm->drvsa[kk].drdz*inp_enm->drvsa[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part2-1)*3+2] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->drvsa[kk].part1-1)*3+2]), 2);

        d1 = sqrt(d1);

        if(inp_enm->distmat_flag==1)
        {
          distmat[(inp_enm->drvsa[kk].part1-1)][(inp_enm->drvsa[kk].part2-1)][nn] = (d1-inp_enm->drvsa[kk].dist)*(d1-inp_enm->drvsa[kk].dist);
          distmat[(inp_enm->drvsa[kk].part2-1)][(inp_enm->drvsa[kk].part1-1)][nn] = (d1-inp_enm->drvsa[kk].dist)*(d1-inp_enm->drvsa[kk].dist);
        }
        deform[(inp_enm->drvsa[kk].part1-1)][nn] += (0.5*inp_enm->drvsa[kk].cost*((d1-inp_enm->drvsa[kk].dist)*(d1-inp_enm->drvsa[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
        deform[(inp_enm->drvsa[kk].part2-1)][nn] += (0.5*inp_enm->drvsa[kk].cost*((d1-inp_enm->drvsa[kk].dist)*(d1-inp_enm->drvsa[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
      }
    }
    else if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
    {
      ////printf("Assigning each interacting residue to its block...\n");
      for (kk=0; kk<inp_enm->nint; kk++) 
      {
        inp_enm->dr[kk].block1 = inp_enm->atm2blk[inp_enm->dr[kk].part1-1];
        inp_enm->dr[kk].block2 = inp_enm->atm2blk[inp_enm->dr[kk].part2-1]; 
      }     
      for (kk=0; kk<inp_enm->nint; kk++) 
      {
        d1=0;      
        d1 = pow(((inp_enm->dr[kk].drdx*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3]   - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3]), 2)  +
             pow(((inp_enm->dr[kk].drdy*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3+1] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3+1]), 2)+
             pow(((inp_enm->dr[kk].drdz*inp_enm->dr[kk].dist) + eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part2-1)*3+2] - eig->eigvec[eig->size-mode[nn]-6][(inp_enm->dr[kk].part1-1)*3+2]), 2);
              
        d1 = sqrt(d1);
        
        if(inp_enm->distmat_flag==1)
        {
          distmat[(inp_enm->dr[kk].block1)][(inp_enm->dr[kk].block2)][nn] += (d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist);
          distmat[(inp_enm->dr[kk].block2)][(inp_enm->dr[kk].block1)][nn] += (d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist);
        }
        deform[(inp_enm->dr[kk].block1)][nn] += (0.5*inp_enm->dr[kk].cost*((d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
        deform[(inp_enm->dr[kk].block2)][nn] += (0.5*inp_enm->dr[kk].cost*((d1-inp_enm->dr[kk].dist)*(d1-inp_enm->dr[kk].dist)))/(nres*eig->eigval[eig->size-mode[nn]-6]*eig->eigval[eig->size-mode[nn]-6]);
      }
    }
  } 
  
  for (ii=0; ii<nres; ii++)
  {
    for (jj=0;jj<len_mode;jj++)
    {
      //deform[ii][len_mode] += deform[ii][jj]*norm[jj];
      deform[ii][len_mode] += deform[ii][jj];
    }
    /*Averaging the deformation energies of each residue over the considered modes*/
    //deform[ii][len_mode] /= len_mode; 
  }
  
  if(inp_enm->distmat_flag==1)
  {
    if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
    {
      /*Getting an average distance fluctuation for each block pair*/
      for (ii=0; ii<nres; ii++)
      {
        for (jj=0; jj<nres; jj++)
        {
          normdist = inp_enm->blk[ii].natom*inp_enm->blk[jj].natom;        
          for (kk=0; kk<len_mode; kk++)
          {
            distmat[ii][jj][kk] = distmat[ii][jj][kk]/normdist;
          }
        }
      }      
    }
    for (ii=0; ii<nres; ii++)
    {
      for (jj=0; jj<nres; jj++)
      {
        for (kk=0; kk<len_mode; kk++)
        {
          //deform[ii][len_mode] += deform[ii][jj]*norm[jj];
          distmat[ii][jj][len_mode] += distmat[ii][jj][kk]/eig->eigval[eig->size-mode[kk]-6];
        }
      }
    }
  }
  
  for ( ii=0; ii<=len_mode; ii++)
  {
    max[ii] = 0;
    for (jj=0;jj<nres;jj++)
    {
      if (deform[jj][ii] > max[ii])
      {
        max[ii] = deform[jj][ii];
      }
    }
  }
  
  if(inp_enm->distmat_flag==1)
  {
    for ( ii=0; ii<=len_mode; ii++)
    {
      max_distmat[ii] = 0;
      for (jj=0;jj<nres;jj++)
      {
        for (kk=0;kk<nres;kk++)
        {
          if (distmat[jj][kk][ii] > max_distmat[ii])
          {
            max_distmat[ii] = distmat[jj][kk][ii];
          }
        }
      }
    }
  }

  if(inp_enm->print == 0)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "deformation_energy.txt");
    for(ii=0; ii<len_mode; ii++)
    {
      memset ( filename_pdb[ii], '\0', sizeof(filename_pdb));
      sprintf( filename_pdb[ii], "deformation_eig%d.pdb", mode[ii]);
      if(inp_enm->distmat_flag==1)
      {
        memset ( distmat_name[ii], '\0', sizeof(distmat_name));
        sprintf( distmat_name[ii], "distmat_eig%d.txt", mode[ii]);
        memset ( distpairs_name[ii], '\0', sizeof(distpairs_name));
        sprintf( distpairs_name[ii], "distpairs_eig%d.txt", mode[ii]);
      } 
    }
    memset ( filename_pdb[len_mode], '\0', sizeof(filename_pdb));
    sprintf( filename_pdb[len_mode], "deformation_eig%d_eig%d_comb.pdb",mode[0], mode[len_mode-1]);
    if(inp_enm->distmat_flag==1)
    {
      memset ( distmat_name[len_mode], '\0', sizeof(distmat_name));
      sprintf( distmat_name[len_mode], "distmat_eig%d_eig%d_comb.txt",mode[0], mode[len_mode-1]);
      memset ( distpairs_name[len_mode], '\0', sizeof(distpairs_name));
      sprintf( distpairs_name[len_mode], "distpairs_eig%d_eig%d_comb.txt",mode[0], mode[len_mode-1]);
    }
  }
  else if (inp_enm->print == 1)
  {
    memset ( filename, '\0', sizeof(filename));
    sprintf( filename, "%s-deformation_energy_fr%d.txt", inp_enm->title, inp_enm->current_frame);
    for(ii=0; ii<len_mode; ii++)
    {
      memset ( filename_pdb[ii], '\0', sizeof(filename_pdb));
      sprintf( filename_pdb[ii], "%s-deformation_eig%d_fr%d.pdb", inp_enm->title, mode[ii], inp_enm->current_frame);    
      if(inp_enm->distmat_flag==1)
      {
        memset ( distmat_name[ii], '\0', sizeof(distmat_name));
        sprintf( distmat_name[ii], "%s-distmat_eig%d_fr%d.txt", inp_enm->title, mode[ii], inp_enm->current_frame);
        memset ( distpairs_name[ii], '\0', sizeof(distpairs_name));
        sprintf( distpairs_name[ii], "%s-distpairs_eig%d_fr%d.txt", inp_enm->title, mode[ii], inp_enm->current_frame);
      }
    }
    memset ( filename_pdb[len_mode], '\0', sizeof(filename_pdb));
    sprintf( filename_pdb[len_mode], "%s-deformation_eig%d_eig%d_comb_fr%d.pdb", inp_enm->title, mode[0], mode[len_mode-1], inp_enm->current_frame);
    if(inp_enm->distmat_flag==1)
    {
      memset ( distmat_name[len_mode], '\0', sizeof(distmat_name));
      sprintf( distmat_name[len_mode], "%s-distmat_eig%d_eig%d_comb_fr%d.txt", inp_enm->title,mode[0], mode[len_mode-1], inp_enm->current_frame);
      memset ( distpairs_name[len_mode], '\0', sizeof(distmat_name));
      sprintf( distpairs_name[len_mode], "%s-distpairs_eig%d_eig%d_comb_fr%d.txt", inp_enm->title,mode[0], mode[len_mode-1], inp_enm->current_frame);
    }
  } 
  
  /*Writing deformation energies*/
  defen = fopen (filename, "w");  
  for (jj=0; jj<nres; jj++)
  {
    for (ii=0; ii<len_mode+1; ii++)
    { 
      fprintf(defen," %10.6f ", deform[jj][ii]);
    }
    fprintf(defen,"\n");
  }
  
///////#########################################/////
///////*Writing distance fluctuations matrices*//////
///////#########################################/////
  if(inp_enm->distmat_flag==1)
  {
    for (nn=0; nn<len_mode+1; nn++)
    {
	    distance = fopen(distmat_name[nn], "w");
	    for (ii=0; ii<nres; ii++)
	    {
	      for (jj=0; jj<nres; jj++)
	      {
	        fprintf(distance," %11.8f ", distmat[ii][jj][nn]);  
	        //fprintf(distance," %11.8f ", distmat[ii][jj][nn]/max_distmat[nn]);  
	      }
        fprintf(distance, "\n"); 	
      }
      fclose(distance);
      if (inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method,"residue",64)==0) //( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
      {     // RTB residue-level
        distance = fopen(distpairs_name[nn], "w");
        for(ii=0; ii<nres; ii++)
        { 
          for(jj=0; jj<nres; jj++)
          {
            //printf("%d %d\n", inp_enm->blk[ii].atom[0], inp_enm->blk[jj].atom[0]);
            Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->blk[ii].atom[0]], cResCode1);
            sprintf(cLabelA, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->blk[ii].atom[0]], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->blk[ii].atom[0]]);        
            Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->blk[jj].atom[0]], cResCode1);
            sprintf(cLabelB, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->blk[jj].atom[0]], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->blk[jj].atom[0]]);
            // WARNING: modded to comply with raw2struct in fileio.c
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", ii+1, jj+1, cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", inp_enm->selmol->rawmol.presn[inp_enm->selean.selatm[ii]-1], inp_enm->selmol->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
            fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
          }
        }
        fclose(distance); 
      }
      else if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method,"residue",64) !=0 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->rtblevan_flag == 1) //( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
      {     // RTB with block file & rtblevan selection
        distance = fopen(distpairs_name[nn], "w");
        //fprintf( distance, "DEBUG %d %d\n", inp_enm->rtblevan.nselatm, inp_enm->rtblevan.selatm[inp_enm->rtblevan.nselatm]-1 );
        for(ii=0; ii<inp_enm->rtblevan.nselatm; ii++)
        { 
          for(jj=0; jj<inp_enm->rtblevan.nselatm; jj++)
          {
            //printf("%d %d\n", inp_enm->blk[ii].atom[0], inp_enm->blk[jj].atom[0]);
            Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->rtblevan.selatm[ii]-1], cResCode1);
            sprintf(cLabelA, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->rtblevan.selatm[ii]-1], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->rtblevan.selatm[ii]-1]);        
            Res3ToRes1(inp_enm->selmol->rawmol.restype[inp_enm->rtblevan.selatm[jj]-1], cResCode1);
            sprintf(cLabelB, "%s:%s%d", inp_enm->selmol->rawmol.segId[inp_enm->rtblevan.selatm[jj]-1], cResCode1,
                                         inp_enm->selmol->rawmol.resn[inp_enm->rtblevan.selatm[jj]-1]);
            // WARNING: modded to comply with raw2struct in fileio.c
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", ii+1, jj+1, cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
            
            // further mods to comply with fixed rtblevan sele (line 451) and to run along rtblevan->nselatm
            // first commented line is same as last of above lines ... 
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[inp_enm->atm2blk[ inp_enm->rtblevan.selatm[ii]-1]][inp_enm->rtblevan.selatm[jj]-1][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], inp_enm->selmol->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[inp_enm->atm2blk[ inp_enm->rtblevan.selatm[ii]-1]][inp_enm->atm2blk[ inp_enm->rtblevan.selatm[jj]-1]][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);
            fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan_org.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan_org.selatm[jj]-1], cLabelA, cLabelB,distmat[inp_enm->atm2blk[ inp_enm->rtblevan.selatm[ii]-1]][inp_enm->atm2blk[ inp_enm->rtblevan.selatm[jj]-1]][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);
          }
        }
        fclose(distance); 
      }
      /* backup copy - before modding to run along rtblevan->nselatm
      else if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method,"residue",64) !=0 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->rtblevan_flag == 1) //( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
      {     // RTB with block file & rtblevan selection
        distance = fopen(distpairs_name[nn], "w");
        fprintf( distance, "%d %d\n", inp_enm->rtblevan.nselatm, inp_enm->rtblevan.selatm[inp_enm->rtblevan.nselatm]-1 );
        for(ii=0; ii<nres; ii++)
        { 
          for(jj=0; jj<nres; jj++)
          {
            //printf("%d %d\n", inp_enm->blk[ii].atom[0], inp_enm->blk[jj].atom[0]);
            Res3ToRes1(molecule->rawmol.restype[inp_enm->blk[ii].sindex[0]], cResCode1);
            sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_enm->blk[ii].sindex[0]], cResCode1,
                                         molecule->rawmol.resn[inp_enm->blk[ii].sindex[0]]);        
            Res3ToRes1(molecule->rawmol.restype[inp_enm->blk[jj].sindex[0]], cResCode1);
            sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_enm->blk[jj].sindex[0]], cResCode1,
                                         molecule->rawmol.resn[inp_enm->blk[jj].sindex[0]]);
            // WARNING: modded to comply with raw2struct in fileio.c
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", ii+1, jj+1, cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
            // further mods to comply with fixed rtblevan sele (line 451) and to get run along rtblevan->nselatm
            fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);
            //fprintf(distance," DEBUG %8d   %8d   %15s   %15s   %11.8f\n", inp_enm->rtblevan.selatm[ii], inp_enm->rtblevan.selatm[jj], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);
          }
        }
        fclose(distance); 
      }
      */
      else if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method,"residue",64) !=0 && strncmp(inp_enm->rtb_method, "off", 64) != 0 && inp_enm->rtblevan_flag == 0 ) //( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
      {     // RTB with block file - no rtblevan WARNING : still to fix!!!
        distance = fopen(distpairs_name[nn], "w");
        for(ii=0; ii<nres; ii++)
        { 
          for(jj=0; jj<nres; jj++)
          {
            //printf("%d %d\n", inp_enm->blk[ii].atom[0], inp_enm->blk[jj].atom[0]);
            Res3ToRes1(molecule->rawmol.restype[inp_enm->blk[ii].sindex[0]], cResCode1);
            sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_enm->blk[ii].sindex[0]], cResCode1,
                                         molecule->rawmol.resn[inp_enm->blk[ii].sindex[0]]);        
            Res3ToRes1(molecule->rawmol.restype[inp_enm->blk[jj].sindex[0]], cResCode1);
            sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_enm->blk[jj].sindex[0]], cResCode1,
                                         molecule->rawmol.resn[inp_enm->blk[jj].sindex[0]]);
            // WARNING: modded to comply with raw2struct in fileio.c
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", ii+1, jj+1, cLabelA, cLabelB,distmat[ii][jj][nn]);
            //fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
            
            fprintf(distance," %8d   %8d   %15s   %15s   %11.8f\n", molecule->rawmol.presn[inp_enm->rtblevan.selatm[ii]-1], molecule->rawmol.presn[inp_enm->rtblevan.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
            //fprintf(distance," DEBUG %8d   %8d   %15s   %15s   %11.8f\n", inp_enm->rtblevan.selatm[ii], inp_enm->rtblevan.selatm[jj], cLabelA, cLabelB,distmat[ii][jj][nn]); //cov[inp_enm->rtblevan.selatm[ii]-1][inp_enm->rtblevan.selatm[jj]-1]);/*Printing R values from -1 to 1*/
          }
        }
        fclose(distance); 
      }
      else
      {
        distance = fopen(distpairs_name[nn], "w");
        for(ii=0; ii<nres; ii++)
        { 
          for(jj=0; jj<nres; jj++)
          {
            
            Res3ToRes1(molecule->rawmol.restype[inp_enm->selean.selatm[ii]-1], cResCode1);
            sprintf(cLabelA, "%s:%s%d", molecule->rawmol.segId[inp_enm->selean.selatm[ii]-1], cResCode1,
                                         molecule->rawmol.resn[inp_enm->selean.selatm[ii]-1]);        
            Res3ToRes1(molecule->rawmol.restype[inp_enm->selean.selatm[jj]-1], cResCode1);
            sprintf(cLabelB, "%s:%s%d", molecule->rawmol.segId[inp_enm->selean.selatm[jj]-1], cResCode1,
                                         molecule->rawmol.resn[inp_enm->selean.selatm[jj]-1]);
            // WARNING: modded to comply with raw2struct in fileio.c
            //fprintf(distance," %8d   %8d   %15s   %15s   %7.4f\n", ii+1, jj+1, cLabelA, cLabelB,distmat[ii][jj][nn]);
            fprintf(distance," %8d   %8d   %15s   %15s   %7.4f\n", molecule->rawmol.presn[inp_enm->selean.selatm[ii]-1], molecule->rawmol.presn[inp_enm->selean.selatm[jj]-1], cLabelA, cLabelB,distmat[ii][jj][nn]);
          }
        }
        fclose(distance); 
      }
    }
    
  }
    
  /*Writing pdb files with deformation energy values in B-factor field*/
  
  frame = inp_enm->molecule;
  if ( inp_enm->mol2_flag == 1 )
  {
    for (ii=0; ii<molecule->nato; ii++)
    {
      frame->coor.xcoor[ii] = inp_enm->moving_sys[ii][0];
      frame->coor.ycoor[ii] = inp_enm->moving_sys[ii][1];
      frame->coor.zcoor[ii] = inp_enm->moving_sys[ii][2];
    }
  }
  else
  {
    for (ii=0; ii<molecule->nato; ii++)
    {
      frame->coor.xcoor[ii] = trj_crd->xcoor[ii];
      frame->coor.ycoor[ii] = trj_crd->ycoor[ii];
      frame->coor.zcoor[ii] = trj_crd->zcoor[ii];
    }
  }
    
  if ( inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0 )
  {
    for( jj=0; jj< molecule->nato; jj++)
    {
      molecule->rawmol.bFac[jj] = 0.0;
    } 
    for (nn=0;nn<=len_mode;nn++)
    {
      for( jj=0; jj<nres; jj++)
      {
        for (ii=0;ii<inp_enm->blk[jj].natom;ii++)
        {
          molecule->rawmol.bFac[inp_enm->blk[jj].atom[ii]-1] = 100 * deform[jj][nn]/max[nn] ;//fabs(100*(A[ii+jj*(inp_pca->sele.nselatm*3)]));
        }
      }  
      WritePdb_unstr( filename_pdb[nn], frame );
    }
  }
  else 
  {
    for(jj=0; jj< molecule->nato; jj++)
    {
      molecule->rawmol.bFac[jj] = 0.0;
    }
 
    for(nn=0; nn<=len_mode; nn++)
    {
      for( jj=0; jj<inp_enm->selean.nselatm; jj++)
      {
        molecule->rawmol.bFac[inp_enm->selean.selatm[jj]-1] = 100 * deform[jj][nn]/max[nn] ;//fabs(100*(A[ii+jj*(inp_pca->sele.nselatm*3)]));
      }  
      WritePdb_unstr( filename_pdb[nn], frame );      
    }
  }
  
  //for(ii=0;ii<len_mode+1;ii++)
  for(ii=0;ii<nres;ii++)
  {
    free(deform[ii]);
  }
  free(deform);
 
  fclose(defen);
  
  return (0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*NORMAL MODES PROJECTIONS*/

int Enm_Proj ( struct inp_enm *inp_enm, _eigen *eig, Molecule *molecule, CoorSet *trj_crd )
{
  FILE *projsinglefile;
  int ii=0, kk=0, jj=0, aa=0;
  int size=0, esize=0, len_mode=0, *mode, scale=0;
  float amp=0.0, *amp_max_list, amp_max=0, spacing=0.0, projection=0.0, temp=0, kT=0, move=0;
  double nframe=0;
  char *tmp, **filename;
  
  Selection selec;

  if( inp_enm->vsa_flag == 1 )
  {
    printf("\n>  Projecting the sub-system along modes: %s\n", inp_enm->proj);
  }
  else
  {
    printf("\n>  Projecting the system along modes: %s \n", inp_enm->proj);
  }

  if (strncmp( inp_enm->proj, "all", 64) == 0)
  {
    len_mode = eig->size-6;
    mode = calloc (len_mode, sizeof(int));
    for(ii=0; ii<len_mode; ii++)
    {
      mode[ii] = ii+1;
    }
  }
  else if (strncmp( inp_enm->proj, " ", 64) != 0)
  {
    tmp = expandrange ( inp_enm->proj, &len_mode );
    mode = readintstring( tmp, len_mode );
    free (tmp);
  }  
   
  filename = calloc (len_mode, sizeof (char*));
  for (ii=0; ii<len_mode; ii++)
  {
    filename[ii] = calloc (1280, sizeof(char));
    memset  ( filename[ii], '\0', sizeof(filename));
    if(inp_enm->print == 0)
    {
      sprintf ( filename[ii], "proj_mode_%d.dcd", mode[ii]);
    }
    else if(inp_enm->print == 1)
    {
      sprintf ( filename[ii], "%s_proj_mode_%d.dcd",inp_enm->title, mode[ii]);
    }
  }
  amp_max_list = calloc (len_mode, sizeof(float));

  if (inp_enm->rtb_flag == 1 && strncmp(inp_enm->rtb_method, "off", 64) != 0  )
  {
    size  = inp_enm->sele1.nselatm; 
    selec = inp_enm->sele1;
    esize = inp_enm->rtb2full_eigen->size;
    scale = 10;
  }
  else if ( inp_enm->vsa_flag == 1)
  {
    size  = inp_enm->vsasele.nselatm; 
    selec = inp_enm->vsasele;
    esize = eig->size; 
    scale =2;   
  }
  else 
  {
    size  = inp_enm->sele1.nselatm;
    selec = inp_enm->sele1;
    esize = eig->size; 
    scale = 2;   
  }

  temp=atof(inp_enm->temp); 
  if ( temp == 0)
  {
    temp=300;
  } 
  kT=k*temp;
  
   
 /*Determining the largest fluctuations observed for the chosen modes*/
 /*and the spacing between each frame*/

  for (ii=0; ii<len_mode; ii++)
  {
    for(jj=0; jj<size; jj++)
    {
      amp = eig->eigvec[esize-mode[ii]-6][jj*3+0]*eig->eigvec[esize-mode[ii]-6][jj*3+0]+eig->eigvec[esize-mode[ii]-6][jj*3+1]*eig->eigvec[esize-mode[ii]-6][jj*3+1]+eig->eigvec[esize-mode[ii]-6][jj*3+2]*eig->eigvec[esize-mode[ii]-6][jj*3+2];
      amp = amp*1/(eig->eigval[esize-mode[ii]-6]);
      if( amp > amp_max_list[ii])
      {
        amp_max_list[ii] = amp; 
      }
      if(amp > amp_max)
      {
        amp_max = amp;
      }   
    }
  }
  
  /*Computing projection along single normal modes*/
  for (aa=0;aa<len_mode;aa++)
  {
    projsinglefile = O_File( filename[aa], "w" );
    WriteDcdHeader( projsinglefile, &inp_enm->projdcdhdr, "all");
    amp_max_list[aa] = amp_max_list[aa]*8*pi*pi*kT/3;
    projection = 0;
    spacing = 1;
    nframe = 50;
    //spacing = 1/sqrt(eig->eigval[esize-mode[len_mode-1]-6]);
    //fracpart = modf((2*amp_max_list[aa]/spacing), &nframe);
    //spacing = nframe/50;
    //printf("%lf  %f\n", nframe, spacing);
    jj=0;
    for( projection = -nframe; projection<nframe; projection+=spacing)
    {

      if( inp_enm->mol2_flag == 1 )
      {
        for(ii=0; ii<molecule->nato; ii++)
        {   
          inp_enm->coor->xcoor[ii] = inp_enm->moving_sys[ii][0];
          inp_enm->coor->ycoor[ii] = inp_enm->moving_sys[ii][1];
          inp_enm->coor->zcoor[ii] = inp_enm->moving_sys[ii][2];   
        }
      }
      else
      {
        for(ii=0; ii<molecule->nato; ii++)
        {
          inp_enm->coor->xcoor[ii] = trj_crd->xcoor[ii];
          inp_enm->coor->ycoor[ii] = trj_crd->ycoor[ii];
          inp_enm->coor->zcoor[ii] = trj_crd->zcoor[ii];		  
        }
      }	 
      
      move = projection*pi/nframe;
      move = cos(move)/sqrt(eig->eigval[esize-mode[aa]-6]);
      for( ii=0; ii<size; ii++ )
      {
        inp_enm->coor->xcoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[aa]-6][ii*3+0];
        inp_enm->coor->ycoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[aa]-6][ii*3+1];
        inp_enm->coor->zcoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[aa]-6][ii*3+2];
      }
 
      WriteDcdCoor( projsinglefile, inp_enm->coor, &inp_enm->projdcdhdr);

      for( ii=0; ii<size; ii++ )
      { 
        inp_enm->coor->xcoor[selec.selatm[ii]-1] = 0;
        inp_enm->coor->ycoor[selec.selatm[ii]-1] = 0;
        inp_enm->coor->zcoor[selec.selatm[ii]-1] = 0;
      }  
      jj++; 
    }
    inp_enm->projdcdhdr.nframe = jj;
    //printf("%d\n",jj);
    WriteDcdHeader( projsinglefile, &inp_enm->projdcdhdr, "nframe" );
    fclose(projsinglefile);    
  }

  /*Computing projection along combination of normal modes*/    
  //spacing = 1/(scale*sqrt(eig->eigval[esize-mode[len_mode-1]-6]));
  //spacing = 1/(eig->eigval[esize-mode[len_mode-1]-6]);
  //printf("nframes=%lf\n", nframe);
 // for( projection = -1*amp_max; projection<=amp_max; projection+=spacing)
 //fracpart = modf((2*amp_max/spacing), &nframe);
 //spacing = nframe/50;
    
   amp_max = amp_max*8*pi*pi*kT/3;
   projection = 0;
   spacing = 1;
   nframe = 50;
   jj=0;
   for( projection = -nframe; projection<=nframe; projection+=spacing)
   {
    if( inp_enm->mol2_flag == 1 )
    {
      for(ii=0; ii<molecule->nato; ii++)
      {   
        inp_enm->coor->xcoor[ii] = inp_enm->moving_sys[ii][0];
        inp_enm->coor->ycoor[ii] = inp_enm->moving_sys[ii][1];
        inp_enm->coor->zcoor[ii] = inp_enm->moving_sys[ii][2];   
       //printf("%f %f %f\n",inp_enm->coor->xcoor[ii],inp_enm->coor->ycoor[ii],inp_enm->coor->zcoor[ii]);
      }
    }
    else
    {
      for(ii=0; ii<molecule->nato; ii++)
      {
       inp_enm->coor->xcoor[ii] = trj_crd->xcoor[ii];
       inp_enm->coor->ycoor[ii] = trj_crd->ycoor[ii];
       inp_enm->coor->zcoor[ii] = trj_crd->zcoor[ii];		  
      }
    }	 
    for ( kk=0; kk<len_mode; kk++ ) 
    {
      move = projection*pi/nframe;
      move = cos(move)/sqrt(eig->eigval[esize-mode[kk]-6]);
      for( ii=0; ii<size; ii++ )
      {
        //inp_enm->coor->xcoor[selec.selatm[ii]-1] += projection*eig->eigvec[esize-mode[kk]-6][ii*3+0]*1/sqrt(eig->eigval[esize-mode[kk]-6]);
        //inp_enm->coor->ycoor[selec.selatm[ii]-1] += projection*eig->eigvec[esize-mode[kk]-6][ii*3+1]*1/sqrt(eig->eigval[esize-mode[kk]-6]);
        //inp_enm->coor->zcoor[selec.selatm[ii]-1] += projection*eig->eigvec[esize-mode[kk]-6][ii*3+2]*1/sqrt(eig->eigval[esize-mode[kk]-6]);
        inp_enm->coor->xcoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[kk]-6][ii*3+0];
        inp_enm->coor->ycoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[kk]-6][ii*3+1];
        inp_enm->coor->zcoor[selec.selatm[ii]-1] += scale*move*eig->eigvec[esize-mode[kk]-6][ii*3+2];
      }
    }
 
    WriteDcdCoor( inp_enm->projdcdfile, inp_enm->coor, &inp_enm->projdcdhdr);
    jj++;   
  
    for( ii=0; ii<size; ii++ )
    { 
      inp_enm->coor->xcoor[selec.selatm[ii]-1] = 0;
      inp_enm->coor->ycoor[selec.selatm[ii]-1] = 0;
      inp_enm->coor->zcoor[selec.selatm[ii]-1] = 0;
      //printf("%f %f %f\n",inp_enm->coor->xcoor[ii],inp_enm->coor->ycoor[ii],inp_enm->coor->zcoor[ii]);
    }  

  }
  inp_enm->projdcdhdr.nframe = jj;
  WriteDcdHeader( inp_enm->projdcdfile, &inp_enm->projdcdhdr, "nframe" );
  fclose(inp_enm->projdcdfile);
  
  free(amp_max_list);
  for (ii=0; ii<len_mode; ii++)
  {
    free (filename[ii]);
  }
  free (filename);
  
  return (0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////ENM related FUNCTIONS///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Subroutine for calculating the 2nd derivative of the Potential energy function*/

int Get_Dr (struct inp_enm *inp_enm, _dr *dr,  CoorSet *trj_crd, Molecule *molecule)
{  

  FILE *m_file;
  int ii=0, jj=0, kk=0, m_ok, m_ato1, m_ato2;
  int kforce=40, nselatm=0; /*Constant force for Kovacs algorithm (Kcal/molxA2)*/
  float cutoff=0, distance=0 , m_cost, **coord;
  char *m_buff, *m_stmp; 
  size_t  m_size=1024;  
  
  /*****Using a linear cutoff*****/
  
  kk = 0;
  if( strncmp (inp_enm->enm_level, "RES", 64)==0 && inp_enm->rtb_flag == 0 && inp_enm->mol2_flag == 0)
  {
    coord = inp_enm->virt_coor;
    nselatm = inp_enm->nselpart;
  }
  else
  {
    coord = inp_enm->coord;
    nselatm = inp_enm->sele1.nselatm; 
  }

  if(strncmp(inp_enm->type,"linear",64)==0)
  {
    printf ("\n>> LINEAR CUTOFF ");
    cutoff = inp_enm->cutoff;
    if( trj_crd->pbc == NULL )
    {
      for (ii=0;ii<nselatm-1;ii++)  
  	  for (jj=ii+1;jj<nselatm;jj++)
  	  {
  	    distance = sqrt((coord[ii][0] - coord[jj][0])*
  			     (coord[ii][0] - coord[jj][0])+
  			     (coord[ii][1] - coord[jj][1])*
  			     (coord[ii][1] - coord[jj][1])+
  			     (coord[ii][2] - coord[jj][2])*
  			     (coord[ii][2] - coord[jj][2]));
  	  
  	    if ((distance < cutoff) && (molecule->rawmol.resn[inp_enm->sele1.selatm[ii]-1] == molecule->rawmol.resn[inp_enm->sele1.selatm[jj-1]-1]) 
  				  && (molecule->rawmol.chainId[inp_enm->sele1.selatm[jj]-1] == molecule->rawmol.chainId[inp_enm->sele1.selatm[ii]-1]) && inp_enm->rtb_flag == 0)     
  	    { 
  	      dr[kk].drdx  = (coord[ii][0] - coord[jj][0])/distance;
  	      dr[kk].drdy  = (coord[ii][1] - coord[jj][1])/distance;
  	      dr[kk].drdz  = (coord[ii][2] - coord[jj][2])/distance;
          dr[kk].dist  = distance;
  	      dr[kk].part1 = ii+1;
  	      dr[kk].part2 = jj+1;
  	      //dr[kk].cost  = 1;
          dr[kk].cost  = 10;
          dr[kk].block1 = 0;
          dr[kk].block2 = 0;
          dr[kk].Id = kk;
  	      kk++;
  	      /*printf("%c\n", molecule.rawmol.chainId[inp_enm->sele1.selatm[ii]-1]);*/
  	    }	 
  	    else if ((distance < cutoff) /*&& (molecule.rawmol.resn[inp_enm->sele1.selatm[jj]-1] > molecule.rawmol.resn[inp_enm->sele1.selatm[ii+1]-1])*/)
  	    {
  	      dr[kk].drdx  = (coord[ii][0] - coord[jj][0])/distance;
  	      dr[kk].drdy  = (coord[ii][1] - coord[jj][1])/distance;
  	      dr[kk].drdz  = (coord[ii][2] - coord[jj][2])/distance;
          dr[kk].dist  = distance;
  	      dr[kk].part1 = ii+1;  
  	      dr[kk].part2 = jj+1;  
  	      dr[kk].cost  = 1;	  
          dr[kk].block1 = 0;
          dr[kk].block2 = 0;
          dr[kk].Id = kk;
           kk++;		  
  	    } 
  	    else
  	    {
  	      dr[kk].drdx  = (coord[ii][0] - coord[jj][0])/distance;
  	      dr[kk].drdy  = (coord[ii][1] - coord[jj][1])/distance;
  	      dr[kk].drdz  = (coord[ii][2] - coord[jj][2])/distance;
          dr[kk].dist  = distance;
  	      dr[kk].part1 = ii+1;  
  	      dr[kk].part2 = jj+1;  
  	      dr[kk].cost  = 0;	  
          dr[kk].block1 = 0;
          dr[kk].block2 = 0;
          dr[kk].Id = 0;
          kk++;
  	    } 
  	  }
    }
    else
    {
      for (ii=0;ii<nselatm-1;ii++)  
  	    for (jj=ii+1;jj<nselatm;jj++)
  	    {
  	      distance = DistanceCoor(coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	  
          if ((distance < cutoff) && (molecule->rawmol.resn[inp_enm->sele1.selatm[ii]-1] == molecule->rawmol.resn[inp_enm->sele1.selatm[jj-1]-1]) 
  				  && (molecule->rawmol.chainId[inp_enm->sele1.selatm[jj]-1] == molecule->rawmol.chainId[inp_enm->sele1.selatm[ii]-1]))     
  	      { 
  	        DistanceAxes( inp_enm->distances, coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	        dr[kk].drdx = inp_enm->distances[0]/distance;
  	        dr[kk].drdy = inp_enm->distances[1]/distance;
  	        dr[kk].drdz = inp_enm->distances[2]/distance;
            dr[kk].dist  = distance;
  	        dr[kk].part1 = ii+1;
  	        dr[kk].part2 = jj+1;
  	        dr[kk].cost = 10;
  	        //dr[kk].cost = 1;
            dr[kk].block1 = 0;
            dr[kk].block2 = 0;
            dr[kk].Id = kk;
  	        kk++;
  	        /*printf("%c\n", molecule.rawmol.chainId[inp_enm->sele1.selatm[ii]-1]);*/
  	      }	 
  	      else if (distance < cutoff)
  	      {
  	        DistanceAxes( inp_enm->distances, coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	        dr[kk].drdx = inp_enm->distances[0]/distance;
  	        dr[kk].drdy = inp_enm->distances[1]/distance;
  	        dr[kk].drdz = inp_enm->distances[2]/distance;
            dr[kk].dist  = distance;
  	        dr[kk].part1 = ii+1;  
  	        dr[kk].part2 = jj+1;  
  	        dr[kk].cost = 1;	  
            dr[kk].block1 = 0;
            dr[kk].block2 = 0;
            dr[kk].Id = kk;
  	        kk++;		  
  	      }
          else
          {
  	        DistanceAxes( inp_enm->distances, coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	        dr[kk].drdx  = inp_enm->distances[0]/distance;
  	        dr[kk].drdy  = inp_enm->distances[1]/distance;
  	        dr[kk].drdz  = inp_enm->distances[2]/distance;
            dr[kk].dist  = distance;
  	        dr[kk].part1 = ii+1;  
  	        dr[kk].part2 = jj+1;  
  	        dr[kk].cost  = 0;	  
            dr[kk].block1 = 0;
            dr[kk].block2 = 0;
	          dr[kk].Id = 0;
            kk++;
          }
  	    }
    }
  }
   /****Using KOVACS' definition of spring constant for interacting pairs*****/
  else if (strncmp(inp_enm->type,"kovacs",64)==0)
  {
    printf ("\n>> KOVACS METHOD");
    for (ii=0;ii<nselatm-1;ii++)  
      for (jj=ii+1;jj<nselatm;jj++)
      {
  	    if( trj_crd->pbc == NULL )
  	    {
  	      distance = sqrt((coord[ii][0] - coord[jj][0])*
  			   (coord[ii][0] - coord[jj][0])+
  			   (coord[ii][1] - coord[jj][1])*
  			   (coord[ii][1] - coord[jj][1])+
  			   (coord[ii][2] - coord[jj][2])*
  			   (coord[ii][2] - coord[jj][2]));
  	      dr[kk].drdx = (coord[ii][0] - coord[jj][0])/distance;
  	      dr[kk].drdy = (coord[ii][1] - coord[jj][1])/distance;
  	      dr[kk].drdz = (coord[ii][2] - coord[jj][2])/distance;
          dr[kk].dist  = distance;
        }
  	    else
  	    {
  	      distance = DistanceCoor(coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	      DistanceAxes( inp_enm->distances, coord[ii][0], coord[ii][1], coord[ii][2], coord[jj][0], coord[jj][1], coord[jj][2], trj_crd->pbc);
  	      dr[kk].drdx = inp_enm->distances[0]/distance;
  	      dr[kk].drdy = inp_enm->distances[1]/distance;
  	      dr[kk].drdz = inp_enm->distances[2]/distance;
          dr[kk].dist  = distance;
  	    }
  	    dr[kk].part1 = ii+1;
  	    dr[kk].part2 = jj+1;
  	    dr[kk].cost = ( pow( (3.8/distance), 6) ) * kforce ;
        dr[kk].block1 = 0;
        dr[kk].block2 = 0;
  	    kk++;
      } 
  }
   
  inp_enm->nint = kk;

    //--MATRIX : added by Rao (ca 2010/02/12)
     //  Module to arbitrarely change pairwise interactions previously created.
     //  Use the keyword MATRIX in the input file to activate this part
     if(inp_enm->matrix_chk==1)
     {
       m_buff = (char *) malloc (m_size + 1);
       m_stmp = (char *) malloc (m_size + 1);

       //--read-in the modified interaction pairs
       m_file =  O_File ( inp_enm->matrix, "r" );
       while ( !feof( m_file ))
       {
         if ( lgetline (&m_buff, &m_size, m_file)!=-1  )
         {
           m_stmp=strtok(m_buff," ");
           m_ato1=atoi(m_stmp);

           m_stmp=strtok(NULL," ");
           m_ato2=atoi(m_stmp);

           m_stmp=strtok(NULL," ");
           m_cost=atof(m_stmp);

           printf ( "MODIFIED INTERACTION-- %5d %5d %10.2f",m_ato1,m_ato2,m_cost);
           //--searching for the correct pairs (prev. created above)
           m_ok=0;
           for (ii=0;ii<kk;ii++)
           {
             if ( inp_enm->sele1.selatm[ dr[ii].part1-1 ] == m_ato1 &&
                     inp_enm->sele1.selatm[ dr[ii].part2-1 ] == m_ato2     ) 
             {
               m_ok=1;
               dr[ii].cost=m_cost;
               printf (" ... done!\n");
             }
           }

           //--warning
           if ( m_ok==0 )
             printf (" ... WARNING!\n    these two atoms are not supposed to interact. Try a larger cutoff!\n");
         }
       }
     }

     // PRINTOUT the contact MATRIX
     if( inp_enm->network_print==1 )
     {
       printf     ( "#%10s %10s %10s %11s %14s\n", "Atom_ID", "Atom_ID" , "K-force", "NAME", "NAME" );
       for (ii=0;ii<kk;ii++)
       {
         printf ( "%10d %10d %10.2f %10s_%-3d %10s_%-3d\n", molecule->rawmol.atomn  [inp_enm->sele1.selatm[dr[ii].part1-1]-1],
                                                  molecule->rawmol.atomn  [inp_enm->sele1.selatm[dr[ii].part2-1]-1], dr[ii].cost,
                                                  molecule->rawmol.atmtype[inp_enm->sele1.selatm[dr[ii].part1-1]-1],
                                                  molecule->rawmol.resn   [inp_enm->sele1.selatm[dr[ii].part1-1]-1],
                                                  molecule->rawmol.atmtype[inp_enm->sele1.selatm[dr[ii].part2-1]-1],
                                                  molecule->rawmol.resn   [inp_enm->sele1.selatm[dr[ii].part2-1]-1] );
       }
     }
     // end --MATRIX Rao mod
     
  printf("\n>  cutoff		 = %-3.1f", cutoff);
  printf("\n>  nodes		 = %-5d", nselatm);
  printf("\n>  pair interactions = %-7d\n", inp_enm->nint);	    

  return 0;
}

int Get_DrVSA (struct inp_enm *inp_enm, _dr *dr,  _dr *drvsa)
{  

  int ii=0, jj=0, kk=0, cont=0;
  
  inp_enm->nintvsa=0;
  cont=0;
  for (ii=0;ii<inp_enm->nint; ii++)
  {
    for (jj=0;jj<inp_enm->vsasele.nselatm;jj++)
    {
      if (dr[ii].part1-1 == inp_enm->subsys_list[jj][0])
      {
        for (kk=0;kk<inp_enm->vsasele.nselatm;kk++)
        {
          if(dr[ii].part2-1 == inp_enm->subsys_list[kk][0])
          {
            drvsa[cont].dist = dr[ii].dist;
            drvsa[cont].cost = dr[ii].cost;
            drvsa[cont].drdx = dr[ii].drdx;
            drvsa[cont].drdy = dr[ii].drdy;
            drvsa[cont].drdz = dr[ii].drdz;
            drvsa[cont].part1 = jj+1;
            drvsa[cont].part2 = kk+1;
            cont++; 
          }
        }
      }
    }
  }
  inp_enm->nintvsa = cont;  

  return 0;
}




/*Subroutine for building the HESSIAN matrix*/

int Build_Hessian (struct inp_enm *inp_enm, float **hessian, _dr *dr, int size)
{
 int  ssize=0, ii=0, jj=0;
 int  indx_11=0, indx_12=0, indx_13=0, indx_21=0, indx_22=0, indx_23=0; 
 float dr1=0, dr2=0, dr3=0, dr4=0, dr5=0, dr6=0, dr7=0, dr8=0, dr9=0;
 ssize = size*size;

 for( ii=0; ii<ssize; ii++)
   hessian[0][ii] = 0.0;
  
 for (ii=0;ii<inp_enm->nint;ii++)
 {
   indx_11 = 3*dr[ii].part1-1;
   indx_12 = 3*dr[ii].part1-2;
   indx_13 = 3*dr[ii].part1-3;
   indx_21 = 3*dr[ii].part2-1;
   indx_22 = 3*dr[ii].part2-2;
   indx_23 = 3*dr[ii].part2-3;
   dr1 = dr[ii].drdx*dr[ii].drdx*dr[ii].cost;
   dr2 = dr[ii].drdx*dr[ii].drdy*dr[ii].cost;
   dr3 = dr[ii].drdx*dr[ii].drdz*dr[ii].cost;
   dr4 = dr[ii].drdy*dr[ii].drdy*dr[ii].cost;
   dr5 = dr[ii].drdy*dr[ii].drdz*dr[ii].cost;
   dr6 = dr[ii].drdz*dr[ii].drdz*dr[ii].cost;
   dr7 = dr[ii].drdy*dr[ii].drdx*dr[ii].cost;
   dr8 = dr[ii].drdz*dr[ii].drdx*dr[ii].cost;
   dr9 = dr[ii].drdz*dr[ii].drdy*dr[ii].cost;
   
   hessian[indx_13][indx_13] +=  dr1;
   hessian[indx_13][indx_12] +=  dr2;
   hessian[indx_13][indx_11] +=  dr3;
   hessian[indx_12][indx_12] +=  dr4;
   hessian[indx_12][indx_11] +=  dr5;
   hessian[indx_11][indx_11] +=  dr6;
   
   hessian[indx_23][indx_23] +=  dr1;
   hessian[indx_23][indx_22] +=  dr2;
   hessian[indx_23][indx_21] +=  dr3;
   hessian[indx_22][indx_22] +=  dr4;
   hessian[indx_22][indx_21] +=  dr5;
   hessian[indx_21][indx_21] +=  dr6;
   
   hessian[indx_13][indx_23] -=  dr1;
   hessian[indx_13][indx_22] -=  dr2;
   hessian[indx_13][indx_21] -=  dr3;
   hessian[indx_12][indx_23] -=  dr7;
   hessian[indx_12][indx_22] -=  dr4;
   hessian[indx_12][indx_21] -=  dr5;
   hessian[indx_11][indx_23] -=  dr8;
   hessian[indx_11][indx_22] -=  dr9;
   hessian[indx_11][indx_21] -=  dr6;
 }
    
     
 for ( ii=1; ii<size; ii++)
 {
   for ( jj=0; jj<ii; jj++)
   { 
     hessian[ii][jj] = hessian[jj][ii];
     //printf("%f ", hessian[jj][ii]);
   }
   //printf("\n");
 }
    
 return 0;
}

/*Subroutine for building a perturbed HESSIAN matrix*/

int Build_Perthessian (struct inp_enm *inp_enm, float **hessian, _dr *dr, int size, int *pert_list, int pert_list_size)
{
  int  ssize=0, ii=0;
  int  indx_11=0, indx_12=0, indx_13=0, indx_21=0, indx_22=0, indx_23=0; 
  float dr1=0, dr2=0, dr3=0, dr4=0, dr5=0, dr6=0, dr7=0, dr8=0, dr9=0;
  ssize = size*size;

  for( ii=0; ii<ssize; ii++)
    hessian[0][ii] = 0.0;
  
  for (ii=0;ii<pert_list_size;ii++)
  {
    indx_11 = 3*dr[pert_list[ii]].part1-1;
    indx_12 = 3*dr[pert_list[ii]].part1-2;
    indx_13 = 3*dr[pert_list[ii]].part1-3;
    indx_21 = 3*dr[pert_list[ii]].part2-1;
    indx_22 = 3*dr[pert_list[ii]].part2-2;
    indx_23 = 3*dr[pert_list[ii]].part2-3;
    dr1 = dr[pert_list[ii]].drdx*dr[pert_list[ii]].drdx*dr[pert_list[ii]].cost;
    dr2 = dr[pert_list[ii]].drdx*dr[pert_list[ii]].drdy*dr[pert_list[ii]].cost;
    dr3 = dr[pert_list[ii]].drdx*dr[pert_list[ii]].drdz*dr[pert_list[ii]].cost;
    dr4 = dr[pert_list[ii]].drdy*dr[pert_list[ii]].drdy*dr[pert_list[ii]].cost;
    dr5 = dr[pert_list[ii]].drdy*dr[pert_list[ii]].drdz*dr[pert_list[ii]].cost;
    dr6 = dr[pert_list[ii]].drdz*dr[pert_list[ii]].drdz*dr[pert_list[ii]].cost;
    dr7 = dr[pert_list[ii]].drdy*dr[pert_list[ii]].drdx*dr[pert_list[ii]].cost;
    dr8 = dr[pert_list[ii]].drdz*dr[pert_list[ii]].drdx*dr[pert_list[ii]].cost;
    dr9 = dr[pert_list[ii]].drdz*dr[pert_list[ii]].drdy*dr[pert_list[ii]].cost;
   
    hessian[indx_13][indx_13] += dr1;
    hessian[indx_13][indx_12] += dr2;
    hessian[indx_13][indx_11] += dr3;
    hessian[indx_12][indx_12] += dr4;
    hessian[indx_12][indx_11] += dr5;
    hessian[indx_11][indx_11] += dr6;
   
    hessian[indx_23][indx_23] += dr1;
    hessian[indx_23][indx_22] += dr2;
    hessian[indx_23][indx_21] += dr3;
    hessian[indx_22][indx_22] += dr4;
    hessian[indx_22][indx_21] += dr5;
    hessian[indx_21][indx_21] += dr6;
   
    hessian[indx_13][indx_23] -= dr1;
    hessian[indx_13][indx_22] -= dr2;
    hessian[indx_13][indx_21] -= dr3;
    hessian[indx_12][indx_23] -= dr7;
    hessian[indx_12][indx_22] -= dr4;
    hessian[indx_12][indx_21] -= dr5;
    hessian[indx_11][indx_23] -= dr8;
    hessian[indx_11][indx_22] -= dr9;
    hessian[indx_11][indx_21] -= dr6;
  }
    
    /*TESTING: Comment the part below to gain in performance!!!!*/
    //for ( ii=1; ii<size; ii++)
    //{
    //  for ( jj=0; jj<ii; jj++)
    // { 
      // hessian[ii][jj] = hessian[jj][ii];
      //printf("%f ", hessian[jj][ii]);
      //}
    //printf("\n");
    //}
    
  return 0;
}

/*Computing Variance from ENM eigenvalues*/
int Enm_Variance (struct inp_enm *inp_enm, _eigen *eig)
{
  int ii=0, ess_perc=90;
  float total_variance=0, ess_sbsp, sum=0;
  
  for (ii=1; ii<eig->size-7;ii++)
  {
    total_variance += 1/pow(eig->eigval[eig->size-6-ii],2);
  }
  
  ess_sbsp = (total_variance*ess_perc)/100;
  
  inp_enm->essential=0;
  sum=0;
  ii=1;
  while (sum < ess_sbsp)
  {
    sum += 1/pow(eig->eigval[eig->size-6-ii],2);
    ii++;  
  }
  
  inp_enm->essential = ii;
  //printf("   Total variance from ENM calculation = %f\n", total_variance);
  //printf("   Essential subspace described by the first %d eigenvectors", inp_enm->essential);

  return 0;  
}


/*Calling LAPACK sgetrf and sgetri routines to compute the inverse of a matrix*/

void InvMatrix ( int size, float **inpmat, float **outmat)
{
   int ii, jj;
   //SGETRF and SGETRI variables

   static int M;   
   static int N;
   static int LDA;
   static int INFO;
   static int LWORK;
   static int *IPIV;
   static float *A;
   static float *WORK;
      
  #ifdef LAPACK
    M     = size;
    N     = size;
    LDA   = size;
    LWORK = LDA*N;
    IPIV  = calloc(N, sizeof (int));
    WORK  = calloc(LWORK,sizeof(float));
    A     = (float *)walloc ( size*size , sizeof(float));
    INFO  = 0;
   
   for( ii=0; ii<size; ii++)
   {
    for ( jj=0; jj<size; jj++)
    {
     A[ii*size+jj] = inpmat[jj][ii];
    }
   }
   
   sgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
   if(INFO < 0)
   {
    fprintf(stderr,"The i-th argument has an illegal value\n");
   }
   else if (INFO > 0)
   {   
    fprintf(stderr,"The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.\n");
   }

   sgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
   if(INFO < 0)
   {
    fprintf(stderr,"The i-th argument has an illegal value\n");
   }
   else if (INFO > 0)
   {   
    fprintf(stderr,"The matrix is singular and its inverse could not be computed\n");
   }

   for( ii=0; ii<size; ii++)
   {
    for ( jj=0; jj<size; jj++)
    {
     outmat[jj][ii] = A[ii*size+jj];
    }
   }
 free(A);

  #else
   fprintf( stderr, "Sorry, all modules requiring matrix inverse with lapack libraries\n have been disabled in this executable\n");
   exit(0);
  #endif
   
   return;
}


void OrtoNorm ( int deg, int num, float **mat)
{
  int ii=0, jj=0, kk=0;
  float aaa=0, anorm=0, rec[6][6];

  anorm = 0.0;
  for (ii=0;ii<num;ii++)
  {
    anorm += mat[ii][0]*mat[ii][0];
  }
  
  anorm = 1/(sqrt(anorm));
  
  for (ii=0;ii<num;ii++)
  {
    mat[ii][0] = anorm*mat[ii][0];
  }  
  
  for (ii=1; ii<deg; ii++)
  {
    for (jj=0; jj<ii; jj++)
    {
      rec[jj][ii] = 0;
      for ( kk=0; kk<num; kk++)
      {
        rec[jj][ii] += mat[kk][ii]*mat[kk][jj];
      }
    }
    for (kk=0; kk<num; kk++)
    {
      aaa = 0.0;
      for (jj=0; jj<ii; jj++)
      {
        aaa += mat[kk][jj]*rec[jj][ii];
      }
      mat[kk][ii] -= aaa;
    }
    
    anorm = 0.0;
    for(kk=0;kk<num;kk++)
    {
      anorm += mat[kk][ii]*mat[kk][ii];
    }
    anorm = 1/sqrt(anorm);
    for(kk=0;kk<num;kk++)
    {
      mat[kk][ii] = anorm*mat[kk][ii];
    }
  }
    
  return;
}

int GeomCentre (struct inp_enm *inp_enm, CoorSet *crd, Molecule *molecule)
{
  int ii=0, NumOfSelRes=0, MultiAtomFlag=0, iResCont=0, iResNum=0;//, iZeroMassWarning=0 ;
  int *piSelResLen, *piProgNum, iVirtAtom=0;
  float fVirtAtomXCoord=0, fVirtAtomYCoord=0, fVirtAtomZCoord=0, fTotalWeight=0, *pfMasses;
  char cTmpString[10], cLastRes[50]; 
  
  for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
  {
    sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_enm->sele1.selatm[ii]-1], 
                                molecule->rawmol.resn[inp_enm->sele1.selatm[ii]-1]);
    if(strcmp(cTmpString, cLastRes)!=0)
    {
      NumOfSelRes++;
      strcpy(cLastRes, cTmpString);
    }
    else
      MultiAtomFlag=1;                                     // More than 1 atom per residue
  }
  
  if(inp_enm->mass_flag==1)
  {
    pfMasses = (float *) calloc(inp_enm->sele1.nselatm, sizeof(float));
    for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
      pfMasses[ii] = molecule->rawmol.bFac[inp_enm->sele1.selatm[ii]-1];
      
    //iZeroMassWarning=0;
    for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
    {
      pfMasses[ii] = molecule->rawmol.bFac[inp_enm->sele1.selatm[ii]-1];
      //if(pfMasses[ii]==0.0)
        //iZeroMassWarning=1;
    }
  }
    
    
  if( MultiAtomFlag==1 )                                    // Calculates the number of atoms per selected residue
  {
    piSelResLen = (int   *) calloc(NumOfSelRes, sizeof(int   ));
    piProgNum   = (int   *) calloc(NumOfSelRes, sizeof(int   ));

    cTmpString[0]='\0';
    cLastRes[0]='\0';
    sprintf(cLastRes, "%s%d", molecule->rawmol.segId[inp_enm->sele1.selatm[0]-1], 
                              molecule->rawmol.resn[inp_enm->sele1.selatm[0]-1]);    
    iResNum=-1;
    iResCont=-1;
                                                      
   piProgNum[0]=molecule->rawmol.presn[inp_enm->sele1.selatm[0]-1];
   
   inp_enm->virt_coor = (float **) calloc(NumOfSelRes, sizeof(float *));
   for(ii=0; ii<NumOfSelRes; ii++)                            // Stores Geometrical/Mass centre coordinates
   {
     inp_enm->virt_coor[ii] = (float *) calloc(3, sizeof(float));
   }
   inp_enm->nselpart = NumOfSelRes;
   for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
   {
     iResCont++;
     sprintf(cTmpString, "%s%d", molecule->rawmol.segId[inp_enm->sele1.selatm[ii]-1], 
                                 molecule->rawmol.resn[inp_enm->sele1.selatm[ii]-1]);
     
     if(strcmp(cTmpString, cLastRes)!=0)
     {
       iResNum++;
       piSelResLen[iResNum]=iResCont;
       strcpy(cLastRes, cTmpString);
       piProgNum[iResNum+1]=molecule->rawmol.presn[inp_enm->sele1.selatm[ii]-1];
       iResCont=0;
     }
   }
    iResNum++;
    piSelResLen[iResNum]=iResCont+1;
    
    // === Calculates geo/mass centres of reference molecule res =======
    iResNum=0;
    iResCont=0;
    iVirtAtom=-1;
    
    fVirtAtomXCoord = 0.0;
    fVirtAtomYCoord = 0.0;
    fVirtAtomZCoord = 0.0;
    
    for(ii=0; ii<inp_enm->sele1.nselatm; ii++)
    {
      if(iResCont==piSelResLen[iResNum])
      {
        if(inp_enm->mass_flag==0)
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/(float)piSelResLen[iResNum]);
          fVirtAtomYCoord = (fVirtAtomYCoord/(float)piSelResLen[iResNum]);
          fVirtAtomZCoord = (fVirtAtomZCoord/(float)piSelResLen[iResNum]);
        }
        else
        {
          fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
          fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
          fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
          
          fTotalWeight = 0.0;
        }
        
        iResNum++;
        iResCont=0;
        iVirtAtom++;
        
        inp_enm->virt_coor[iVirtAtom][0]=fVirtAtomXCoord;
        inp_enm->virt_coor[iVirtAtom][1]=fVirtAtomYCoord;
        inp_enm->virt_coor[iVirtAtom][2]=fVirtAtomZCoord;
        

        fVirtAtomXCoord = 0.0;
        fVirtAtomYCoord = 0.0;
        fVirtAtomZCoord = 0.0;
        printf("%d %d %f %f %f \n",iResCont, piSelResLen[iResNum], inp_enm->virt_coor[iVirtAtom][0], inp_enm->virt_coor[iVirtAtom][1], inp_enm->virt_coor[iVirtAtom][2]);
     }

      if(inp_enm->mass_flag==0)
      {
        // Geometrical centre
        fVirtAtomXCoord += crd->xcoor[inp_enm->sele1.selatm[ii]-1];
        fVirtAtomYCoord += crd->ycoor[inp_enm->sele1.selatm[ii]-1];
        fVirtAtomZCoord += crd->zcoor[inp_enm->sele1.selatm[ii]-1];
      }
      else
      {
        // Mass centre
        fVirtAtomXCoord += (crd->xcoor[inp_enm->sele1.selatm[ii]-1] * pfMasses[ii]);
        fVirtAtomYCoord += (crd->ycoor[inp_enm->sele1.selatm[ii]-1] * pfMasses[ii]);
        fVirtAtomZCoord += (crd->zcoor[inp_enm->sele1.selatm[ii]-1] * pfMasses[ii]);
        
        fTotalWeight += pfMasses[ii];
      }
      iResCont++;
 
    }

    // Process the last residue
    if(inp_enm->mass_flag==0)
    {
      printf("%d %d %f %f %f \n",iResCont, piSelResLen[iResNum], fVirtAtomXCoord, fVirtAtomYCoord, fVirtAtomZCoord);
      fVirtAtomXCoord = (fVirtAtomXCoord/(float)piSelResLen[iResNum]);
      fVirtAtomYCoord = (fVirtAtomYCoord/(float)piSelResLen[iResNum]);
      fVirtAtomZCoord = (fVirtAtomZCoord/(float)piSelResLen[iResNum]);
      

    }
    else
    {
      fVirtAtomXCoord = (fVirtAtomXCoord/fTotalWeight);
      fVirtAtomYCoord = (fVirtAtomYCoord/fTotalWeight);
      fVirtAtomZCoord = (fVirtAtomZCoord/fTotalWeight);
      
      fTotalWeight = 0.0;
    }
    
    iVirtAtom++;
    inp_enm->virt_coor[iVirtAtom][0]=fVirtAtomXCoord;
    inp_enm->virt_coor[iVirtAtom][1]=fVirtAtomYCoord;
    inp_enm->virt_coor[iVirtAtom][2]=fVirtAtomZCoord;
    
    // =================================================================
    
  }
  
  return(0);  
}
