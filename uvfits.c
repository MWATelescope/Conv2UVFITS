/************************************
 ! MODULE: uvfits
 !         utilities for writing/reading UVFITS files in AIPS/MIRIAD format
 !         Author: Randall Wayth. Sep, 2006.
 !         Feb, 2007. Added uvfits reader.
 !         2013: Updates for interator mode operation (single timestep to conserve memory)
 !         2014: Re-sync with RTS version, including multiple IF support originally by Stew Gleadow
 !              Note that multi-IF support is inherently limited to the same number of polarisations
 !              and same pol type etc in the data as discussed in the memo.
 !  The reference for the AIPS uvfits format is AIPS Memo 117. This code was meant to
 !  implement a subset of that for single source, single pointing data
 !  and has been tested with Miriad, AIPS and CASA.
 !  Note: this is not FITS-IDI.
 !  For the FITS-IDI format, see: http://www.aips.nrao.edu/FITS-IDI.html
 ************************************/

/* a note about sign conventions etc:
  The de-facto standard for UVFITS files follows the AIPS convention:
  - a baseline is defined as location(ant1)-location(ant2) where the antenna indices are such that ant2 > ant1
  - using this definition of a baseline, a visibility (for a point source) is V(u,v) = I*exp(-2*pi*(ul+vm))
  - this means that if a baseline is defined as (ant2-ant1),
    then the exponent in the exponential must be +2pi, not -2pi

  This writer does NOT do any reformatting of data to/from the de-facto standard. It is up to the user of this
  code to encode/decode visibilities as appropriate before the reader/writer is called.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <float.h>
#include <assert.h>
#include <star/pal.h>
#include "uvfits.h"

#define MAX_CSIZE 8


/* private function prototypes */
static int writeAntennaData(fitsfile *fptr, uvdata *data);
static int writeFQData(fitsfile *fptr, uvdata *data);
void printHeader(fitsfile *fptr);
int readAntennaTable(fitsfile *fptr,array_data *array);
int writeSourceData(fitsfile *fptr, uvdata *data);
int readFQTable(fitsfile *fptr,uvdata *obj);
int readUVFITSOpenInit(char *filename, uvdata **data, fitsfile **outfptr);
void addGroupRow(uvdata *obj, uvReadContext *iter);

/* private global vars */
static int debug=0;

/******************************
 ! NAME:      uvfitsSetDebugLevel(
 ! PURPOSE:   set the module global debug level.
 ! ARGUMENTS: in_debug - integer global debug setting. 0= no debugging messages.
 ! RETURNS:   void
******************************/
void uvfitsSetDebugLevel(int in_debug) {
    debug = in_debug;
}


/******************************
 ! NAME:      readUVFITSOpenInit
 ! PURPOSE:   open a UVFITS file and create a uvdata data structure
 !            This also reads the antenna table so you know how many antennas ther are.
 !            This sets up for iterative reading of data via readUVFITSInitIterator
 ! ARGUMENTS: filename: name of UVFITS file to read
 !            data: pointer to pointer data struct. The struct will be allocated by this function.
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes if a CFITSIO error
 !            occurred. returns -1 for malloc failure.
******************************/
int readUVFITSOpenInit(char *filename, uvdata **data, fitsfile **outfptr) {
    fitsfile *fptr;
    int i,status=0,num_hdu=0,hdu_type=0;
    int bitpix,naxis;
    long naxes[MAX_CSIZE];
    char temp[84];
    uvdata *obj=NULL;

    /* open the file */
    fits_open_file(&fptr,filename, READONLY, &status);
    if(status!=0) {
        fprintf(stderr,"readUVFITSOpenInit: failed to open file <%s>\n",filename);
        return status;
    }
    *outfptr = fptr;

    /* get data size, number of extensions etc */
    fits_get_num_hdus(fptr,&num_hdu,&status);
    if(debug) fprintf(stdout,"readUVFITSOpenInit: there are %d HDUs\n",num_hdu);
    if (num_hdu > 2) {
        fprintf(stderr,"readUVFITSOpenInit WARNING: %s appears to be a multi-source or multi-IF file. This is poorly supported.\n",filename);
        //    return 1;
    }
    fits_get_hdu_type(fptr,&hdu_type,&status);
    if (hdu_type != IMAGE_HDU) {
        fprintf(stderr,"readUVFITSOpenInit: first HDU in file is not an image. type: %d\n",hdu_type);
    }

    /* get dimensions of vis data */
    fits_get_img_param(fptr,MAX_CSIZE,&bitpix,&naxis,naxes,&status);
    if(status!=0) {
        fits_report_error(stderr,status);
        fprintf(stderr,"%s: failed to get img params.\n",__func__);
        return status;
    }

    if (debug) {
        fprintf(stdout,"readUVFITSOpenInit: NAXIS is %d\n",naxis);
        for(i=0;i<naxis;i++) fprintf(stdout,"\tNAXIS%d is %d\n",i+1,(int)naxes[i]);
    }
    if (naxis > 7) {
        fprintf(stderr,"ERROR: NAXIS is %d. Extra image dimensions aren't supported yet\n",(int)naxis);
        return -1;
    }

    /* basic stuff seems to be in order, allocate space for data structure (not visibilites yet) */
    *data = calloc(1,sizeof(uvdata));
    if(*data==NULL) {
        fprintf(stderr,"readUVFITSOpenInit: no malloc for main struct\n");
        return -1;
    }
    obj = *data;
    obj->visdata=NULL;
    obj->weightdata=NULL;
    obj->baseline=NULL;
    obj->u=NULL;
    obj->v=NULL;
    obj->w=NULL;
    obj->n_baselines=NULL;
    obj->source = calloc(1,sizeof(source_table)); /* no source table yet... */
    obj->array = calloc(1,sizeof(array_data)); /* data filled in in readAntennaTable */
    obj->date = calloc(1,sizeof(double));
    if (obj->source==NULL || obj->date==NULL) {
        fprintf(stderr,"readUVFITSOpenInit: no malloc for source table\n");
        return -1;
    }
    obj->n_pol = naxes[2];    // copy across number of pols
    obj->n_freq = naxes[3];   // copy across number of freqs
    obj->n_IF = 1;            // set the number of IFs
    if (naxis==7) {
      // extra image axis means there are multiple IFs to support
      obj->n_IF = naxes[4];
      if (obj->n_IF < 1 || obj->n_IF > 128) {    // max is arbitrary here. Just a sanity check
          fprintf(stderr,"ERROR: bad number of IFs: %d\n",(int)obj->n_IF);
          return -1;
      }
    }

    /* examine headers and read the antenna table */
    for(i=1;i<=num_hdu; i++) {

        fits_movabs_hdu(fptr, i, NULL, &status);
        if (status) {
            fits_report_error(stderr,status);
            return status;
        }

        if (debug) {
            fprintf(stdout," *** Printing header index %d\n",i);
            printHeader(fptr);
        }

        temp[0] = '\0';
        fits_read_key(fptr,TSTRING,"EXTNAME",temp,NULL,&status);
        if (status == 0 && strncmp(temp,"AIPS AN",7)==0) {
            if(debug) printf("AN table is in HDU %d\n",i);
            status = readAntennaTable(fptr,obj->array);
            if (status) return status;
        }
        if (status == 0 && strncmp(temp,"AIPS FQ",7)==0) {
            if(debug) printf("FQ table is in HDU %d\n",i);
            status = readFQTable(fptr,obj);
            if (status) return status;
        }
        /* clear errors from trying to read keywords that possibly aren't there */
        fits_clear_errmsg();
        status=0;
    }
    fits_movabs_hdu(fptr, 1, NULL, &status);  // go back to start
    return 0;
}


/******************************
 ! NAME:      readUVFITSInitIterator
 ! PURPOSE:   open a UVFITS file and create a uvdata data structure for successive data reads
 !            each data read fetches a single timestep of data for all baselines,freqs,pols
 !            This also reads the antenna table
 ! ARGUMENTS: filename: name of UVFITS file to read
 !            data: pointer to pointer data struct. The struct will be allocated by this function.
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes if a CFITSIO error
 !            occurred. returns -1 for malloc failure.
******************************/
int readUVFITSInitIterator(char *filename, uvdata **data, uvReadContext **iterator) {
    fitsfile *fptr=NULL;
    int i,status=0,max_bl,grp_row_size,found_obsra=0;
    char temp[84]="";
    uvReadContext *iter=NULL;

    /* init */

    /* open the file and get basic info*/
    status = readUVFITSOpenInit(filename, data, &fptr);
    if(status!=0) {
        return status;
    }
    iter = calloc(1,sizeof(uvReadContext));
    assert(iter != NULL);
    *iterator = iter;
    iter->fptr = fptr;
    iter->base_date = 0;   // set this sentinel for check below.

    /* find out how many groups there are. There is one group for each source,time,baseline combination.
       we don't know how many baselines there are yet- have to work this out while reading the data.
       The rows must be time ordered or everything will break.  */
    fits_read_key(fptr,TSTRING,"GROUPS",temp,NULL,&status);
    if (temp[0] != 'T') {
        fprintf(stderr,"%s: GROUPS fits keyword not set. This is not a uvfits file.\n",__func__);
        return -1;
    }
    fits_read_key(fptr,TINT,"PCOUNT",&(iter->pcount),NULL,&status);
    fits_read_key(fptr,TINT,"GCOUNT",&(iter->gcount),NULL,&status);
    if (iter->pcount==0 || iter->gcount==0) {
        fprintf(stderr,"ERROR: pcount: %d, gcount %d. Cannot be zero\n",iter->pcount,iter->gcount);
        return -1;
    }
    fits_read_key(fptr,TSTRING,"OBJECT",(*data)->source->name,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSDEC",&((*data)->source->dec),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSRA",&((*data)->source->ra),NULL,&status);
    if (status !=0) {
        fprintf(stderr,"%s WARNING: status %d while getting object keywords.\n",__func__,status);
        fits_clear_errmsg();
        status=0;
    }
    else {
        found_obsra=1;
    }

    /* get the details of the observations PTYPE/CTYPE keyword from header*/
    /* set defaults fisrt */
    for (i=0; i<3; i++) {
        iter->pscal[i] = 1.0;
        iter->pzero[i] = 0.0;
    }
    for(i=0; i<iter->pcount; i++) {
        char typedesc[84];
        double temp_pscal=1.0,temp_pzero=0.0;

        sprintf(temp,"PTYPE%d",i+1);
        fits_read_key(fptr,TSTRING,temp,typedesc,NULL,&status);
        if (status) {
            fprintf(stderr,"ERROR: Didn't find %s keyword in header\n",typedesc);
            return -1;
        }
        sprintf(temp,"PSCAL%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,&temp_pscal,NULL,&status);
        sprintf(temp,"PZERO%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,&temp_pzero,NULL,&status);
        if(debug) fprintf(stdout,"PTYPE%d: '%s', PSCAL: %g, PZERO: %g\n",i+1,typedesc,temp_pscal,temp_pzero);
        fits_clear_errmsg(); status=0;// clear fits errors in case PSCAL and PZERO are missing since they have defaults

        /* remember the stuff we care about */
        if (strncmp("DATE",typedesc,4)==0) {
            iter->base_date += temp_pzero;
            // handle 2nd occurance of DATE field which can be in non-standard CASA format
            // in this case, the first DATE is what should be in the CRVALs in the header, which is the
            // base date that the offsets in the visibilities are added to to get the JD
            // so copy over the index of the previous occurrance of DATE for date0
            if (iter->ptype_map.date != 0) {
                iter->ptype_map.date0 = i;
            }
            else {
                iter->ptype_map.date = i;
            }
        }
        if (strncmp("UU",typedesc,2)==0) {
            iter->pscal[0] = temp_pscal;
            iter->pzero[0] = temp_pzero;
            iter->ptype_map.u = i;
        }
        if (strncmp("VV",typedesc,2)==0) {
            iter->pscal[1] = temp_pscal;
            iter->pzero[1] = temp_pzero;
            iter->ptype_map.v = i;
        }
        if (strncmp("WW",typedesc,2)==0) {
            iter->pscal[2] = temp_pscal;
            iter->pzero[2] = temp_pzero;
            iter->ptype_map.w = i;
        }
        if (strncmp("BASELINE",typedesc,8)==0) {
            iter->ptype_map.bl = i;
        }
        if (strncmp("FREQSEL",typedesc,7)==0) {
            iter->ptype_map.fq = i;
        }
    }
    if(iter->base_date <= 0.0) {
        fprintf(stderr,"ERROR: no JD base date read from PTYPE keywords in header.\n");
    }

    /* get the CTYPEs, which start at CRVAL2 */
    for(i=0; i<=7; i++) {
        char typedesc[84]="";
        double crval=0,crpix=1,cdelt=0;

        sprintf(temp,"CTYPE%d",i+2);
        fits_read_key(fptr,TSTRING,temp,typedesc,NULL,&status);
        if (status != 0) {
            if(debug) fprintf(stdout,"no more C keys. status: %d\n",status);
            status =0;
            fits_clear_errmsg();
            break;
        }
        sprintf(temp,"CRVAL%d",i+2);
        fits_read_key(fptr,TDOUBLE,temp,&crval,NULL,&status);
        sprintf(temp,"CRPIX%d",i+2);
        fits_read_key(fptr,TDOUBLE,temp,&crpix,NULL,&status);
        sprintf(temp,"CDELT%d",i+2);
        fits_read_key(fptr,TDOUBLE,temp,&cdelt,NULL,&status);
        if(debug) fprintf(stdout,"CRVAL%d: %g,\tCRPIX%d: %g,\tCDELT%d: %g,\tCTYPE%d: %s\n",
                                 i+2,crval,i+2,crpix,i+2,cdelt,i+2,typedesc);

        /* check for the stuff we're interested in */
        if (strncmp("FREQ",typedesc,4)==0) {
            /* reference freq is CRVAL at pixel CRPIX, BW is CDELT */
            (*data)->cent_freq = crval + ( ((*data)->n_freq)/2 +1-crpix)*cdelt;
            (*data)->freq_delta = cdelt;
        }
        if (strncmp("STOKES",typedesc,6)==0) {
            (*data)->pol_type = crval;
        }
        /* RA and DEC might be in the OBSRA and OBSDEC fields already */
        if (strncmp("RA",typedesc,2)==0) {
            if (found_obsra) {
                fprintf(stderr,"%s: WARNING: overwriting OBSRA (%g) with CTYPE (%g)\n",__func__,(*data)->source->ra,crval);
            }
            (*data)->source->ra = crval;
        }
        if (strncmp("DEC",typedesc,3)==0) {
            if (found_obsra) {
                fprintf(stderr,"%s: WARNING: overwriting OBSDEC (%g) with CTYPE (%g)\n",__func__,(*data)->source->dec,crval);
            }
            (*data)->source->dec = crval;
        }
        if (strncmp("IF",typedesc,2)==0) {
        }
    }
    /* allocate space for visibility data. In iterator mode, we only keep one timestep (n_vis=1), so the
       arrays of pointers in the uvdata structure are overkill but it doesn't hurt */
    (*data)->visdata = calloc(1,sizeof(float *));
    (*data)->weightdata = calloc(1,sizeof(float *));
    (*data)->baseline = calloc(1,sizeof(float *));
    (*data)->n_baselines = calloc(1,sizeof(int));
    (*data)->u = calloc(1,sizeof(double *));
    (*data)->v = calloc(1,sizeof(double *));
    (*data)->w = calloc(1,sizeof(double *));
    (*data)->n_vis = 1;

    if ((*data)->weightdata==NULL || (*data)->visdata==NULL || (*data)->baseline==NULL || 
        (*data)->u==NULL || (*data)->v==NULL || (*data)->w==NULL || (*data)->date==NULL ) {
        fprintf(stderr,"%s: no malloc for vis data\n",__func__);
        return -1;
    }

    /* now alloc space for the visibilities and u,v,weight values. We already know the max number
       of baselines from how many antennas are in the antenna table */
    max_bl = (*data)->array->n_ant * ((*data)->array->n_ant+1)/2;
    (*data)->visdata[0] = malloc(sizeof(float)*(*data)->n_freq*(*data)->n_pol*max_bl*2);
    (*data)->weightdata[0] = malloc(sizeof(float)*(*data)->n_freq*(*data)->n_pol*max_bl);
    (*data)->baseline[0]=malloc(max_bl*sizeof(float));
    (*data)->u[0] = calloc(max_bl,sizeof(double));
    (*data)->v[0] = calloc(max_bl,sizeof(double));
    (*data)->w[0] = calloc(max_bl,sizeof(double));
    if ((*data)->visdata[0]==NULL || (*data)->weightdata[0]==NULL || (*data)->baseline[0]==NULL) {
        fprintf(stderr,"%s: No malloc for vis data\n",__func__);
      return -1;
    }

    grp_row_size = 3*((*data)->n_pol)*((*data)->n_freq); // 3*npols*nfreqs
    iter->grp_row = malloc(sizeof(float)*grp_row_size);
    iter->grp_par = malloc(sizeof(float)*iter->pcount);
    if (iter->grp_row==NULL || iter->grp_par==NULL) {
        fprintf(stderr,"%s: no malloc for group row or params\n",__func__);
        return -1;
    }

    return 0;
}


/******************************
 ! NAME:      addGroupRow
 ! PURPOSE:   add a row of parameters (one baseline/time) to the output data structure.
 ! ARGUMENTS: obj: data structure returned by InitIterator for resulting data
 !            iter: a uvReadContext iterator context structure
 ! RETURNS:   void
******************************/
void addGroupRow(uvdata *obj, uvReadContext *iter) {
    int k,l,m,index=0,j,time_index=0;
    long out_ind;

    j = obj->n_baselines[time_index];
    obj->baseline[time_index][j] = iter->grp_par[iter->ptype_map.bl];
    obj->date[time_index] = iter->grp_par[iter->ptype_map.date];
    if (iter->ptype_map.date0 > 0) obj->date[time_index] += iter->grp_par[iter->ptype_map.date0];   // add optional second DATE
    obj->date[time_index] += iter->base_date;
    // handle special non-standard CASA dates
    if (iter->ptype_map.date0 != 0) {
        obj->date[time_index] = (double) iter->grp_par[iter->ptype_map.date] + (double) iter->grp_par[iter->ptype_map.date0];
        //printf("Date: %f, date0: %f\n",iter->grp_par[iter->ptype_map.date],iter->grp_par[iter->ptype_map.date0]);
    }
    obj->u[time_index][j] = iter->grp_par[iter->ptype_map.u]*iter->pscal[0] + iter->pzero[0];
    obj->v[time_index][j] = iter->grp_par[iter->ptype_map.v]*iter->pscal[1] + iter->pzero[1];
    obj->w[time_index][j] = iter->grp_par[iter->ptype_map.w]*iter->pscal[2] + iter->pzero[2];

    if (debug>1) fprintf(stdout,"Adding vis. (u,v,w): %g,%g,%g. Bl: %f. Time: %f. N_baselines: %d\n",
                obj->u[time_index][j], obj->v[time_index][j], obj->w[time_index][j], obj->baseline[time_index][j],
                iter->grp_par[4],j);

    for(m=0; m<obj->n_IF; m++) {
        for(k=0; k<obj->n_freq; k++) {
            for(l=0; l< obj->n_pol; l++) {
//                out_ind = (j*obj->n_freq*obj->n_pol + k*obj->n_pol + l);
                out_ind = l + obj->n_pol*(k + obj->n_freq*(j + obj->n_IF*m));       // include multiple IFs
                obj->visdata[time_index][out_ind*2+0] = iter->grp_row[index++];   // real
                obj->visdata[time_index][out_ind*2+1] = iter->grp_row[index++];   // imag
                obj->weightdata[time_index][out_ind]  = iter->grp_row[index++];   // weight
            }
        }
    }
    obj->n_baselines[time_index]++;
}


/******************************
 ! NAME:      readUVFITSnextIter
 ! PURPOSE:   get the next chunk of time from an open uvfits file and fill uvdata structure.
 !            requires prior call to readUVFITSInitIterator
 ! ARGUMENTS: obj: data structure returned by InitIterator for resulting data
 !            iter: a uvReadContext iterator context structure
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes if a CFITSIO error
 !            occurred. returns -1 for malloc failure.
 !            return code 1 means normal completion end of file
******************************/
int readUVFITSnextIter(uvdata *obj, uvReadContext *iter) {
    int grp_row_size,status=0,done=0;
    int max_bl,n_bl=0;
    double currtime=FLT_MAX,this_time;
    fitsfile *fptr;
 
    fptr = (fitsfile *)iter->fptr;
    max_bl = obj->array->n_ant*(obj->array->n_ant + 1)/2;
    grp_row_size = 3*(obj->n_pol)*(obj->n_freq)*(obj->n_IF); // 3*npols*nfreqs*n_ifs

    obj->n_baselines[0]=0;  // this is a new time, so reset the number of baselines for this time before adding new data
    if (debug) fprintf(stdout,"grp_index is %d. gcount is %d\n",iter->grp_index, iter->gcount);
    if (iter->grp_index >= iter->gcount) {
        // we've reached the end
        done=1;
        status=1;
        if (debug) fprintf(stdout,"grp_index is %d. All done\n",iter->grp_index);
    }

    while (!done) {
        
        /* read a row from the parameters. contains U,V,W,baseline,time */
        fits_read_grppar_flt(fptr,iter->grp_index+1,1,iter->pcount, iter->grp_par, &status);
        if (status != 0) {
            fprintf(stderr,"readUVFITSnextIter: error reading group %d. status: %d\n",iter->grp_index+1,status);
            fits_report_error(stderr,status);
            return status;
        }

        if (debug > 1) fprintf(stdout,"Read group %d of %d. Time: %f\n",iter->grp_index+1,iter->gcount,iter->grp_par[4]);

        // calculate the time, which might be the sum of two "DATE" fields
        this_time = iter->grp_par[iter->ptype_map.date];
        if (iter->ptype_map.date0 > 0) {
            this_time += iter->grp_par[iter->ptype_map.date0];  // add optional second date field
        }

        /* if the time changes, we're done for this chunk */
        if (currtime != this_time && currtime != FLT_MAX) {

            if (debug) printf("*** new time *** %f\n",iter->grp_par[4]);
            /* check that time is increasing only */
            if (this_time < currtime) {
                fprintf(stderr,"readUVFITSnextIter ERROR: data is not time ordered. current: %f, new: %f\n",currtime,iter->grp_par[4]);
                return -1;
            }
           done=1;
        }
        else {
            // read the vis/weight values for all freqs for this time/baseline
            fits_read_img_flt(fptr, iter->grp_index+1, 1, grp_row_size, 0.0, iter->grp_row, NULL, &status);
            if (status) {
                fprintf(stderr,"readUVFITSnextIter: Error reading image row for time index %d\n",iter->grp_index+1);
                fits_report_error(stderr,status);
                return status;
            }

            addGroupRow(obj,iter);
            n_bl++; // count how many baselines we've read and do some sanity checks
            if (n_bl > max_bl) {
                fprintf(stderr,"ERROR: counted %d baselines for time %g. Expected max: %d\n",n_bl,iter->grp_par[4],max_bl);
            }
            iter->grp_index++;
            currtime = this_time;
        }
        if (iter->grp_index >= iter->gcount) {
            // we've reached the end
            done=1;
            if (debug) fprintf(stdout,"grp_index is %d. All done\n",iter->grp_index);
        }
    }
    return status;
}


/******************************
 ! NAME:      readUVFITSCloseIter
 ! PURPOSE:   close the FITS file and free the iterator structure
 ! ARGUMENTS: iter: pointer to pointer data struct
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes for a CFITSIO error
******************************/
int readUVFITSCloseIter(uvReadContext *iter) {
    int status=0;

    /* tidy up */
    fits_close_file(iter->fptr, &status);            /* close the file */
    fits_report_error(stderr, status);  /* print out any error messages */

    if (iter->grp_par != NULL) free(iter->grp_par);
    if (iter->grp_row != NULL) free(iter->grp_row);
    free(iter);
    return status;
}


/******************************
 ! NAME:      readUVFITS
 ! PURPOSE:   read an entire UVFITS file and create/populate a uvdata data structure
 !            NOTE this method is deprecated since it uses a LOT of memory. Use the iterator instead.
 ! ARGUMENTS: filename: name of UVFITS file to read
 !            data: pointer to pointer data struct. The struct will be allocated by this function.
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes if a CFITSIO error
 !            occurred. returns -1 for malloc failure.
******************************/
int readUVFITS(char *filename, uvdata **data) {
    fitsfile *fptr=NULL;
    int i,j,k,l,status=0,retcode=0,grp_row_size=0,index=0;
    int pcount,gcount,num_bl=0,time_index=0,bl_index=0;
    int date_pindex=4,uscale_pindex=0,vscale_pindex=1,wscale_pindex=2,freq_cindex=4;
    int if_cindex=3;
    char temp[84];
    char *ptype[MAX_CSIZE],*ctype[MAX_CSIZE];
    float crval[MAX_CSIZE],crpix[MAX_CSIZE],cdelt[MAX_CSIZE],*grp_row=NULL,*grp_par=NULL;
    float currtime=FLT_MAX,mintime=FLT_MAX;
    double pscal[MAX_CSIZE],pzero[MAX_CSIZE];
    uvdata *obj;

    /* init */
    memset(ptype,0,MAX_CSIZE*sizeof(char *));
    memset(ctype,0,MAX_CSIZE*sizeof(char *));

    /* open the file and get basic info*/
    status = readUVFITSOpenInit(filename, data, &fptr);
    if(status!=0) {
        fprintf(stderr,"readUVFITS: failed to read file <%s>\n",filename);
        goto EXIT;
    }
    obj = *data;

    /* find out how many groups there are. There is one group for each source,time,baseline combination.
       we don't know how many baselines there are yet- have to work this out while reading the data.
       The rows must be time ordered or everything will break.  */
    fits_read_key(fptr,TSTRING,"GROUPS",temp,NULL,&status);
    if (temp[0] != 'T') {
      fprintf(stderr,"readUVFITS: GROUPS fits keyword not set. This is not a uvfits file.\n");
      status = -1;
      goto EXIT;
    }
    fits_read_key(fptr,TINT,"PCOUNT",&pcount,NULL,&status);
    fits_read_key(fptr,TINT,"GCOUNT",&gcount,NULL,&status);
    fits_read_key(fptr,TSTRING,"OBJECT",(*data)->source->name,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSRA",&((*data)->source->ra),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSDEC",&((*data)->source->dec),NULL,&status);
    if (status !=0) {
        fprintf(stderr,"readUVFITS WARNING: status %d while getting basic keywords.\n",status);
        fits_clear_errmsg();
        status=0;
    }

    /* get the details of the observations PTYPE/CTYPE keyword from header*/
    for(i=0; i<pcount; i++) {
        sprintf(temp,"PTYPE%d",i+1);
        ptype[i] = malloc(sizeof(char)*84);
        if (ptype[i]==NULL) {
            fprintf(stderr,"readUVFITS: no malloc for ptype char array\n");
            status = -1; goto EXIT;
        }
        ptype[i][0]='\0';
        fits_read_key(fptr,TSTRING,temp,ptype[i],NULL,&status);
        sprintf(temp,"PSCAL%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,pscal+i,NULL,&status);
        sprintf(temp,"PZERO%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,pzero+i,NULL,&status);

        /* remember the stuff we care about */
        if (strncmp("DATE",ptype[i],4)==0) {
            date_pindex = i;
        }
        if (strncmp("UU",ptype[i],2)==0) {
            uscale_pindex = i;
        }
        if (strncmp("VV",ptype[i],2)==0) {
            vscale_pindex = i;
        }
        if (strncmp("WW",ptype[i],2)==0) {
            wscale_pindex = i;
        }
        if (strncmp("FREQSEL",ptype[i],7)==0) {
        /* FIXME: this is unfinished for the multi IF case */
//            fscale_pindex = i;
        }

        if(debug) fprintf(stdout,"PZERO%d: %f, PSCAL%d: %g, PTYPE%d: %s\n",i+1,pzero[i],i+1,pscal[i],i+1,ptype[i]);
    }

    /* get the CRs, which start at CRVAL2. Should go up to 7 max */
    for(i=0; i<MAX_CSIZE; i++) {
        sprintf(temp,"CTYPE%d",i+2);
        ctype[i] = malloc(sizeof(char)*84); /* strings can be 80 chars not including terminator */
        if (ctype[i]==NULL) {
            fprintf(stderr,"readUVFITS: no malloc for ctype char array\n");
            status = -1; goto EXIT;
        }
        ctype[i][0] = '\0';
        fits_read_key(fptr,TSTRING,temp,ctype[i],NULL,&status);
        if (status != 0) {
            if(debug) fprintf(stdout,"no more C keys. status: %d\n",status);
            status =0;
            fits_clear_errmsg();
            break;
        }
        sprintf(temp,"CRVAL%d",i+2);
        fits_read_key(fptr,TFLOAT,temp,crval+i,NULL,&status);
        sprintf(temp,"CRPIX%d",i+2);
        fits_read_key(fptr,TFLOAT,temp,crpix+i,NULL,&status);
        sprintf(temp,"CDELT%d",i+2);
        fits_read_key(fptr,TFLOAT,temp,cdelt+i,NULL,&status);
        if(debug) fprintf(stdout,"CRVAL%d: %g,\tCRPIX%d: %g,\tCDELT%d: %g,\tCTYPE%d: %s\n",
		      i+2,crval[i],i+2,crpix[i],i+2,cdelt[i],i+2,ctype[i]);

        /* check for the stuff we're interested in */
        if (strncmp("FREQ",ctype[i],4)==0) {
            /* centre freq is CRVAL, BW is CDELT */
            freq_cindex = i+1;
            (*data)->cent_freq = crval[i];
            (*data)->freq_delta = cdelt[i];
        }
        if (strncmp("IF", ctype[i],2)==0) {
            if_cindex = i+1;
        }
    }

    // Added by DAM. If the IF key isn't present, if_cindex ends up set to freq_cindex.
    if ( if_cindex == freq_cindex ) if_cindex++;

/*
    if(debug) fprintf(stdout,"NAXES: Complex(%d)=%ld, Pol(%d)=%ld, IF(%d)=%ld, Freq(%d)=%ld\n",complex_cindex,
                    naxes[complex_cindex], pol_cindex, naxes[pol_cindex], if_cindex, naxes[if_cindex],
                    freq_cindex, naxes[freq_cindex]);
*/

    grp_row_size = 3*((*data)->n_pol)*((*data)->n_freq)*(*data)->n_IF; // 3*npol*nfreq*n_IF
    /* each row in the group table has this size */
    grp_row = malloc(sizeof(float)*grp_row_size*gcount);
    grp_par = malloc(sizeof(float)*pcount);
    if (grp_row==NULL || grp_par==NULL) {
        fprintf(stderr,"readUVFITS: no malloc for big arrays\n");
        status = -1; goto EXIT;
    }

    /* copy data into data structure. We know how many groups there are, but we don't know
     how many different baselines or time units there are. The data must be ordered by increasing
     time, so we keep track of the time until it changes. This then tells us how many baselines
     there are. If there are missing baselines in the first time step, it will be bad.... */
    for (i=0; i< gcount; i++) {

        fits_read_img_flt(fptr, i+1, 1, grp_row_size, 0.0, grp_row+i*grp_row_size, NULL, &status);
        if (status) {
            fprintf(stderr,"Error reading image row index %d\n",i);
            fits_report_error(stderr,status);
            goto EXIT;
        }

        /* read a row from the parameters. contains U,V,W,baseline,time */
        fits_read_grppar_flt(fptr,i+1,1,pcount, grp_par, &status);
        if (status != 0) {
            fprintf(stderr,"readUVFITS: error reading group %d. status: %d\n",i+1,status);
            status = -1; goto EXIT;
        }
        /* if the time changes, make some space for more data */
        if (currtime != grp_par[4]) {

            /* new set of baselines */
            currtime = grp_par[4];
            if (debug) printf("new time: %f\n",currtime);
            /* check that time is increasing only */
            if (currtime < mintime) {
                if (mintime != FLT_MAX) {
                    fprintf(stderr,"readUVFITS ERROR: data is not time ordered. current: %f, min: %f\n",currtime,mintime);
                }
                else {
                    /* first time around, set mintime */
                    mintime = currtime;
                }
            }

            /* allocate space for visibility data. Increase
               the size of these arrays to have one for each time sample */
            obj->visdata = realloc(obj->visdata,sizeof(float *)*(time_index+1));
            obj->weightdata = realloc(obj->weightdata,sizeof(float *)*(time_index+1));
            obj->baseline = realloc(obj->baseline,sizeof(float *)*(time_index+1));
            obj->n_baselines = realloc(obj->n_baselines,sizeof(int)*(time_index+1));
            obj->u = realloc(obj->u,sizeof(double *)*(time_index+1));
            obj->v = realloc(obj->v,sizeof(double *)*(time_index+1));
            obj->w = realloc(obj->w,sizeof(double *)*(time_index+1));
            obj->date = realloc(obj->date,sizeof(double)*(time_index+1));

            if (obj->weightdata==NULL || obj->visdata==NULL || obj->baseline==NULL || 
                obj->u==NULL || obj->v==NULL || obj->w==NULL || obj->date==NULL ) {
                fprintf(stderr,"readUVFITS: no malloc in param loop\n");
                status = -1; goto EXIT;
            }

            obj->visdata[time_index]=NULL;
            obj->weightdata[time_index]=NULL;
            obj->baseline[time_index]=NULL;
            obj->u[time_index]=NULL;
            obj->v[time_index]=NULL;
            obj->w[time_index]=NULL;
            obj->n_baselines[time_index]=0;

            /* there may be less than the full number of baselines in the first and last scan from VLA data.
               so keep track of this so we don't run over the array end below */
            bl_index=0;  /* reset for next scan */
            time_index++;
        }

        /* put the UVW data into the data structure now */
        /* first, resize the arrays for the next baseline */
        obj->baseline[time_index-1] = realloc(obj->baseline[time_index-1],sizeof(float)*(bl_index+1));
        obj->u[time_index-1] = realloc(obj->u[time_index-1],sizeof(double)*(bl_index+1));
        obj->v[time_index-1] = realloc(obj->v[time_index-1],sizeof(double)*(bl_index+1));
        obj->w[time_index-1] = realloc(obj->w[time_index-1],sizeof(double)*(bl_index+1));
        if (obj->u[time_index-1]==NULL || obj->v[time_index-1]==NULL || obj->w[time_index-1]==NULL) {
            fprintf(stderr,"readUVFITS: realloc failed for UVW arrays\n");
            status = -1; goto EXIT;
        }
        /* copy the UVW, baseline, time */
        obj->u[time_index-1][bl_index] = grp_par[0]*pscal[uscale_pindex] + pzero[uscale_pindex];
        obj->v[time_index-1][bl_index] = grp_par[1]*pscal[vscale_pindex] + pzero[vscale_pindex];
        obj->w[time_index-1][bl_index] = grp_par[2]*pscal[wscale_pindex] + pzero[wscale_pindex];
        obj->baseline[time_index-1][bl_index] = grp_par[3];
        obj->date[time_index-1] = grp_par[4]+pzero[date_pindex];
        obj->n_baselines[time_index-1] = bl_index+1;

        bl_index++;
        /* keep track of the maximum number of baselines that are listed */
        if (num_bl < bl_index) num_bl = bl_index;
        if(debug) printf("PARAMS %d: %g,\t%g,\t%g,\t%g,\t%g\n",
		     i+1,grp_par[0],grp_par[1],grp_par[2],grp_par[3],grp_par[4]);
    }

    /* now we know how many time steps there are and how many groups there are.
       The number of baselines has been tracked in the loop. */
    obj->n_vis = time_index;
    if (debug) printf("There are %d time steps. Counted maximum %d baselines\n",time_index,num_bl);

    /* now alloc space for the visibilities */
    for (i=0; i< time_index; i++) {
        obj->visdata[i] = malloc(sizeof(float)*obj->n_freq*obj->n_pol*num_bl*obj->n_IF*2);
        obj->weightdata[i] = malloc(sizeof(float)*obj->n_freq*obj->n_pol*num_bl*obj->n_IF);
        if (obj->visdata[i]==NULL || obj->weightdata[i]==NULL) {
            fprintf(stderr,"No malloc in time loop\n");
            status = -1; goto EXIT;
        }
    }

    /* pack the visibility data structure */
    for (i=0; i< time_index; i++) {
        for (j=0; j<obj->n_baselines[i]; j++) {
            int m,out_ind;
            /* vis data and weights. U,V,W already done above. */
            for(m=0; m< obj->n_IF; m++) {
            for(k=0; k<obj->n_freq; k++) {
                for(l=0; l< obj->n_pol; l++) {
                    out_ind = l + obj->n_pol*(k + obj->n_freq*(j + obj->n_IF*m));
                    obj->visdata[i][out_ind*2  ] = grp_row[index++];   // real
                    obj->visdata[i][out_ind*2+1] = grp_row[index++];   // imag
                    obj->weightdata[i][out_ind]  = grp_row[index++];   // weight
                }
            }
            }
        }
    }

 EXIT:
    retcode = status;
    status=0;
    /* all done, close file */
    fits_close_file(fptr, &status);
    if (status!=0) fprintf(stderr,"readUVFITS: problem closing file. status: %d\n",status);
    /* done with raw data, so free everything we have allocated */
    if (grp_row !=NULL) free(grp_row);
    if (grp_par!=NULL) free(grp_par);
    for (i=0; i< MAX_CSIZE; i++) {
        if (ptype[i]!=NULL) free(ptype[i]);
        if (ctype[i]!=NULL) free(ctype[i]);
    }
    return retcode;
}


/********************************
 ! NAME:      writeUVFITS
 ! PURPOSE:   write a UVFITS file (old method where all data is known up front - deprecated)
 !            This function is now just a wrapper for the new iterator method.
 ! ARGUMENTS: filename: name of output file (will be created)
 !            data: pointer to uvdata struct that contains the data
 ! RETURNS:   integer: 0 success, nonzero failure. returns CFITSIO errors in the case of a CFITSIO error.
*********************************/
int writeUVFITS(char *filename, uvdata *data) {
    int status=0,i;
    uvWriteContext *iter;

    /* set up the iterator */
    status= writeUVFITSiterator(filename, data, &iter);
    if (status) return status;

    /* write the acutal UV data to the main array.
       this is done in "random groups" with one group for each
       instant of visibilities. Each instant contains all freqs,
       pols and baselines */
    for (i=0; i<data->n_vis; i++) {
        /*    fprintf(stderr,"doing vis set %d of %d. There are %d baselines\n",i+1,data->n_vis,data->n_baselines);*/
        status = writeUVinstant(iter, data, data->date[i]-(iter->jd_day_trunc),i);
        if (status) {
            fprintf(stderr,"ERROR: code %d when writing time instant %d\n",status,i);
            exit(EXIT_FAILURE);
        }
    }

    /* write other stuff at the end like antenna table etc and close the file */  
    status = writeUVFITSfinalise(iter, data);
    return status;
}

/**************************
 ! NAME:      writeUVFITSiterator
 ! PURPOSE:   create an output uvfits file and set it up for iterative writing of individual time steps.
 ! ARGUMENTS: filename: an open fits file.
 !            data: uvdata struct with existing data.
 !            vfptr: resulting pointer to fitsfile in void format
 !            requires most observation info, expect the actual data, to be existing already in the uvdata including:
 !            n_pol, n_freq, source name ra dec, dates, pol types and freq info.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeUVFITSiterator(char *filename, uvdata *data, uvWriteContext **iter) {

    fitsfile *fptr;
    FILE *fp;
    int status=0;
    long naxes[7],i;
    long naxis=6,n_grp_params=5,n_ctypes=4;   // set defaults for single IF data
    float temp;
    char tempstr[40],*ctypes[5];
    char *ctypes_multiIF[5] = {"STOKES","FREQ","IF","RA","DEC"};
    char *params[6] = {"UU","VV","WW","BASELINE","DATE","FREQSEL"};
    int mon=0,day=0,year=0;
    double jd_day,dtemp=0;

    /* create the iterator context */
    *iter = calloc(1,sizeof(uvWriteContext));
    assert(*iter !=NULL);

    /* set up for cfitsio interface */
    if (data->n_IF ==0) data->n_IF =1;    // sanity check
    // If this is multi-IF data, adjust arrays sizes
    if (data->n_IF > 1) {
        naxis=7;
        n_grp_params=6;
        n_ctypes=5;
    }
    (*iter)->n_grp_params = n_grp_params;

    /* here we set the dimensions of the random group data. We exclude the FQ table index if there
       is only 1 IF to prevent writing redundant indices of 0 in the data and follow miriad's format */
    i=0;
    naxes[i++] = 0;     // indicates this is random group data
    naxes[i++] = 3;     // 
    naxes[i++] = data->n_pol;
    naxes[i++] = data->n_freq;
    if (data->n_IF > 1) naxes[i++] = data->n_IF;
    naxes[i++] = 1;
    naxes[i++] = 1;

    /* set the CTYPE names depending on whether this is multi-IF or not 
       don't use the "IF" keyword for single IF */
    for (i=0; i< 5; i++) ctypes[i] = ctypes_multiIF[i];   // default, copy all over in place
    if (data->n_IF < 2) {
        for (i=2; i< 4; i++) ctypes[i] = ctypes_multiIF[i+1];   // move RA and DEC up one place to remove "IF"
    }

    if (data->n_baselines[0] <1) {
        fprintf(stderr,"There are no baselines.\n");
        return 1;
    }
    if (filename==NULL || filename[0]=='\0') {
        fprintf(stderr,"writeUVFITSiterator: empty file name\n");
        return 1;
    }

    /* open the file after checking to see if it already exists and clobbering */
    if ((fp=fopen(filename,"r"))!=NULL) {
        fclose(fp);
        remove(filename);
    }
    fits_create_file(&fptr,filename,&status);
    if (status !=0) {
        fits_report_error(stderr, status);
        return 1;
    }
    (*iter)->fptr=fptr;
    /* set up for writing fits "random groups". There is one group for each baseline
       of visibilities per time sample including all pols and channels. The group consists of:
       1- the (N_GRP_PARAMS) preamble of: U,V,W (nanoseconds),baseline (see below), time offset (days)
       2- the vis data as real,imag,weight for pol1, real,imag,weight for pol2 etc
       for each pol, then repeated for each freq channel.
       optionally another block of vis data for channels in the Nth IF for multi-IF data */
    /* using this function with a single IF creates the keywords PCOUNT=5,GROUPS=T and GCOUNT=data->n_baselines. */
    fits_write_grphdr(fptr,TRUE,FLOAT_IMG,naxis,naxes,n_grp_params,data->n_baselines[0],TRUE,&status);

    /* Write a keyword; must pass pointers to values */
    fits_update_key(fptr,TSTRING, "OBJECT", data->source->name, NULL, &status);
    dtemp = data->source->ra*15.0;    // FITS uses degrees for all angles
    fits_update_key(fptr,TDOUBLE, "OBSRA", &dtemp, "deprecated. Use CTYPE RA key instead", &status);
    fits_update_key(fptr,TDOUBLE, "OBSDEC", &(data->source->dec), "deprecated. Use CTYPE DEC key instead", &status);
    fits_update_key(fptr,TSTRING, "TELESCOP", data->array->name, NULL, &status);
    fits_update_key(fptr,TSTRING, "INSTRUME", data->array->instrument , NULL, &status);
    temp=2000.0;
    fits_update_key(fptr,TFLOAT, "EPOCH", &temp, NULL, &status);
    /* the following is a fix for an AIPS bug. Technically, shouldn't need it */
    temp=1.0;
    fits_update_key(fptr,TFLOAT, "BSCALE", &temp, NULL, &status);
    temp=0.0;
//  fits_update_key(fptr,TFLOAT, "BZERO", &temp, NULL, &status);

    /* get the JD of the first obs. Then truncate the time within the day.
       The times associated with the visibilities are then increments to
       this truncated start-of-jd time. JDs start at 12:00h on a day. */
    JD_to_Cal(data->date[0],&year,&mon,&day); /* get year, month, day for this JD */
    Cal_to_JD(year,mon,day,&jd_day);          /* gets the JD for the start of this day */
    (*iter)->jd_day_trunc = ((int)jd_day)+0.5;         /* truncate */
    sprintf(tempstr,"%d-%02d-%02dT00:00:00.0",year,mon,day);
    fits_update_key(fptr,TSTRING, "DATE-OBS", tempstr , NULL, &status);

    /* UV parameters */
    /* this sets the UU,VV,WW,BASELINE, DATE and optionally FREQSEL header params */
    for (i=0; i<n_grp_params; i++) {
        sprintf(tempstr,"PTYPE%d",(int)i+1);
        fits_update_key(fptr,TSTRING, tempstr, params[i] , NULL, &status);
        sprintf(tempstr,"PSCAL%d",(int)i+1);
        temp=1.0;
        fits_update_key(fptr,TFLOAT, tempstr,&temp , NULL, &status);
        sprintf(tempstr,"PZERO%d",(int)i+1);
        temp=0.0;
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
    }
    /* except the PZERO for the date, which must match the OBS-DATE in Julian Days
       the following relies on the tempstr from previous loop. */
    fits_update_key(fptr,TDOUBLE,tempstr,&((*iter)->jd_day_trunc), NULL, &status);

    /* write the coord system keywords CRVAL, CRPIX, CDELT, CTYPE.*/
    // CTYPE2 is always fixed as below for uvfits data 
    temp=1.0;
    fits_update_key(fptr,TFLOAT, "CRVAL2", &temp , NULL, &status);
    fits_update_key(fptr,TFLOAT, "CRPIX2", &temp , NULL, &status);
    fits_update_key(fptr,TFLOAT, "CDELT2", &temp , NULL, &status);
    fits_update_key(fptr,TSTRING,"CTYPE2", "COMPLEX" , NULL, &status);

    for (i=0; i<n_ctypes; i++){
        /* this fills in the CRVAL etc arrays with some generic values, then specific cases are updated below. */
        sprintf(tempstr,"CTYPE%d",(int)i+3);
        fits_update_key(fptr,TSTRING, tempstr, ctypes[i] , NULL, &status);
        sprintf(tempstr,"CRVAL%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
        sprintf(tempstr,"CRPIX%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
        sprintf(tempstr,"CDELT%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
    }

    /* set the input pol type. for circular, CRVAL3 should be -1, CDELT3 -1
       for linear CRVAL3 should be -5, CDELT3 -1
       for Stokes CRVAL3 should be  1, CDELT3  1 */
    temp = data->pol_type;
    fits_update_key(fptr,TFLOAT, "CRVAL3", &temp , NULL, &status);
    temp = (data->pol_type < 0 ? -1.0 : 1.0);
    fits_update_key(fptr,TFLOAT, "CDELT3", &temp , NULL, &status);

    /* specific updates for freq channels. FIXME: check this is correct!!
       need to check that delta might want to be negative. */
    fits_update_key(fptr,TFLOAT,"CRVAL4",&(data->cent_freq),NULL,&status);
    fits_update_key(fptr,TFLOAT,"CDELT4",&(data->freq_delta) , NULL, &status);
    /* set the reference pixel for the central freq to be the center of the channels */
    /* FIXME: is there supposed to be a 1/2 channel offset here for an even number of channels? */
    temp = data->n_freq/2 + 1; /* pixel indexes start at 1... */
    fits_update_key(fptr,TFLOAT,"CRPIX4",&temp , NULL, &status);

    /* specific updates for RA/DEC */
    temp = data->source->ra*15;
    fits_update_key(fptr,TFLOAT,"CRVAL5",&temp,NULL,&status);
    temp = data->source->dec;
    fits_update_key(fptr,TFLOAT,"CRVAL6",&temp,NULL,&status);

    /* write the HISTORY AIPS WTSCAL item which seems to be required */
    fits_write_history(fptr,"AIPS WTSCAL =  1.0",&status);
    fits_write_comment(fptr,"written by the UV FITS writer of RBW.",&status);

    return EXIT_SUCCESS;
}


/**************************
 ! NAME:      writeUVFITSfinalise
 ! PURPOSE:   write antenna tables etc at end of uvfits file and close up
 !            after writing a uvfits file with the iterator process
 ! ARGUMENTS: vfptr: an open fits file (as a void pointer)
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeUVFITSfinalise(uvWriteContext *iter, uvdata *data) {
    int status=0;

    if(debug) {
        fprintf(stdout,"Finalising and writing antenna table\n");
      fflush(stdout);
    }

    if (iter==NULL || iter->fptr==NULL) {
        fprintf(stderr,"writeUVFITSfinalise ERROR: invalid uvWriteContext structure\n");
        return EXIT_FAILURE;
    }

    /* create and write the antenna table */
    status = writeAntennaData(iter->fptr, data);
    if (status !=0) {
      fprintf(stderr,"ERROR: writeAntennaData failed with status %d\n",status);
    }

    /* if this is multi-IF data, create and write the FQ table */
    if (data->n_IF > 1) {
        if(debug) {
            fprintf(stdout,"Writing FQ table\n");
            fflush(stdout);
        }
        status=writeFQData(iter->fptr, data);
        if (status !=0) {
            fprintf(stderr,"ERROR: writeAntennaData failed with status %d\n",status);
        }
    }

    /* tidy up */
    fits_close_file(iter->fptr, &status);            /* close the file */
    fits_report_error(stderr, status);  /* print out any error messages */

    free(iter);

    return status;
}


/**************************
 ! NAME:      writeUVinstant
 ! PURPOSE:   write one time instant of uv data into the UVFITS file. Requires prior call to writeUVFITSiterator
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 !            n_elements: number of items in a group for a time instant
 !            i: time index of data in uvdata array. This should normally be zero.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeUVinstant(uvWriteContext *iter, uvdata *data, double jd_frac, int i) {
    double *u_ptr, *v_ptr, *w_ptr;
    float  *vis_ptr, *wt_ptr;
    float *array=NULL;
    int nelements,status=0;
    int l,j,k,m;
    char msg[80];

    /* sanity checks */
    if (data->n_pol <1 || data->n_freq < 1 || data->n_baselines[i] < 1) {
        fprintf(stderr,"%s ERROR: npols %d, nfreqs %d and n_baselines %d all must be > 0 at time index %d\n",
                        __func__,data->n_pol, data->n_freq, data->n_baselines[i],i);
        return EXIT_FAILURE;
    }
    if (iter->fptr==NULL || iter->n_grp_params==0) {
        fprintf(stderr,"%s ERROR: uvWriteContext is invalid\n",__func__);
        return EXIT_FAILURE;
    }

    if (debug) fprintf(stdout,"%s: jd: %f, scan: %d, n_grp_params: %d. So far ntimes %ld, ngroups %ld\n",__func__,
                        jd_frac,i,(int) iter->n_grp_params, iter->ntimes, iter->ngroups);

    /* magic number 3 here matches naxes[1] from above. */
    nelements = 3*data->n_pol*data->n_freq*data->n_IF;
    assert(nelements > 0);
    array = malloc((nelements+iter->n_grp_params)*sizeof(float));
    if(array==NULL) {
        fprintf(stderr,"%s: no malloc for %ld floats\n",__func__,nelements+iter->n_grp_params);
        exit(EXIT_FAILURE);
    }

    u_ptr = data->u[i];
    v_ptr = data->v[i];
    w_ptr = data->w[i];
    vis_ptr = data->visdata[i];
    wt_ptr = data->weightdata[i];
    for(l=0; l<data->n_baselines[i]; l++) {
      long out_ind,inp_ind;

      array[0] = u_ptr[l];
      array[1] = v_ptr[l];
      array[2] = w_ptr[l];
      /*      EncodeBaseline(data->array->baseline_ant1[l],data->array->baseline_ant2[l],array+3); */
      array[3] = data->baseline[i][l];
      /* last, the time offset */
      array[4] = jd_frac;

      for (m=0; m<data->n_IF; m++){
        for (j=0; j<data->n_freq; j++){
          for (k=0; k<data->n_pol; k++) {
//              inp_ind = l*data->n_freq*data->n_pol + j*data->n_pol + k;
//              out_ind = (j*data->n_pol*3)+(k*3) + iter->n_grp_params;
              inp_ind = k + data->n_pol*(j + data->n_freq*(m + data->n_IF*l));  // include IF index
                // we write 1 baseline at a time, so no l index for output
              out_ind = 3*(k+ data->n_pol*(j + data->n_freq*(m))) + iter->n_grp_params; 

              array[out_ind  ] = vis_ptr[(inp_ind)*2  ];    /* real part */
              array[out_ind+1] = vis_ptr[(inp_ind)*2+1];    /* imag part */
              array[out_ind+2] = wt_ptr[inp_ind];           /* weight */
          }
        }
      }
      /* Write each vis batch as a "random group" including n_grp_params extra params at the front */
      /* the starting index must be 1, not 0.  fortran style. also for the vis index */
      fits_write_grppar_flt(iter->fptr,(iter->ngroups)+1,1,nelements + iter->n_grp_params, array, &status);
      /*      fprintf(stderr,"write group number: %d\n",i*data->n_baselines+l+1);*/
      fits_report_error(stderr, status);  /* print out any error messages */
      if (status) return status;
      (iter->ngroups)++;
    }
    sprintf(msg,"%ld timesteps, %d baselines",++(iter->ntimes),data->n_baselines[i]);
    fits_update_key(iter->fptr, TLONG, "GCOUNT", &(iter->ngroups), NULL, &status);
    fits_report_error(stderr, status);  /* print out any error messages */
    free(array);
    return 0;
}

/**************************
 ! NAME:      writeSourceData
 ! PURPOSE:   write the source ("SU") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeSourceData(fitsfile *fptr, uvdata *data) {
    return 0;
}


/**************************
 ! NAME:      writeAntennaData
 ! PURPOSE:   write the antenna ("AIPS AN") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeAntennaData(fitsfile *fptr, uvdata *data) {

    array_data *array;
    int status=0,i,itemp=0,year,mon,day;
    double temp=0,mjd;
    char *ptemp;
    char *col_names[] = {"ANNAME", "STABXYZ", "NOSTA", "MNTSTA", "STAXOF",
                         "POLTYA", "POLAA", "POLCALA", "POLTYB", "POLAB", "POLCALB"};
    char *col_format[] = {"8A","3D","1J","1J","1E",
                          "1A","1E","3E","1A","1E","3E"};
    char *col_units[] = {"","METERS","","","METERS","","DEGREES","","","DEGREES",""};
    char tempstr[80];

    if (debug) fprintf(stdout,"Writing antenna table\n");

    array = data->array;

    /* convert JD date to Calendar date. I don't think sub-day accuracy is needed here. */
    JD_to_Cal(data->date[0],&year,&mon,&day);

    /* create a new binary table */
    fits_create_tbl(fptr,BINARY_TBL,0, 11 ,col_names, col_format,col_units,"AIPS AN", &status);
    if (status) {
        fits_report_error(stderr,status);
        fprintf(stderr,"%s: failed to create antenna table\n",__func__);
        return status;
    }

    /* write table header info */
    fits_update_key(fptr,TDOUBLE,"ARRAYX", &(data->array->xyz_pos[0]), NULL, &status);
    fits_update_key(fptr,TDOUBLE,"ARRAYY", &(data->array->xyz_pos[1]), NULL, &status);
    fits_update_key(fptr,TDOUBLE,"ARRAYZ", &(data->array->xyz_pos[2]), NULL, &status);
    fits_update_key(fptr,TFLOAT,"FREQ", &(data->cent_freq) , "[Hz]", &status);
    if (status) {
        fits_report_error(stderr,status);
        fprintf(stderr,"%s: WARNING: failed to add array/freq keywords to antenna table\n",__func__);
        fprintf(stderr,"%s: vals were: ARRAY X,Y,Z: %g,%g,%g freq: %g\n",__func__, data->array->xyz_pos[0],
                    data->array->xyz_pos[1], data->array->xyz_pos[2], data->cent_freq);
        fits_clear_errmsg(); status=0;
    }

    /* GSTIAO is the GST at zero hours in the time system of TIMSYS (i.e. UTC) */
    mjd = trunc(data->date[0] - 2400000.5);
    temp = palGmst(mjd)*180.0/M_PI;  // technically, palGmst takes UT1, but it won't matter here.
    if (debug) {
      fprintf(stdout,"GMST is: %g degs\n",temp);
    }
    fits_update_key(fptr,TDOUBLE,"GSTIA0",&temp , NULL, &status);
    temp = 3.60985e2;
    fits_update_key(fptr,TDOUBLE,"DEGPDY",&temp, NULL, &status);

    sprintf(tempstr,"%d-%02d-%02dT00:00:00.0",year,mon,day);
    fits_update_key(fptr,TSTRING,"RDATE",tempstr , NULL, &status);

    temp=0.0;
    fits_update_key(fptr,TDOUBLE,"POLARX", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"POLARY", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"UT1UTC", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"DATUTC", &temp, NULL, &status);
    fits_update_key(fptr,TSTRING,"TIMSYS","UTC" , NULL, &status);
    fits_update_key(fptr,TSTRING,"ARRNAM", data->array->name, NULL, &status);
    itemp = 0;
    fits_update_key(fptr,TINT,"NUMORB",&itemp , NULL, &status);
    itemp = 3;
    fits_update_key(fptr,TINT,"NOPCAL",&itemp , NULL, &status);
    itemp = -1;
    fits_update_key(fptr,TINT,"FREQID",&itemp, NULL, &status);
    temp= palDat(mjd);
    fits_update_key(fptr,TDOUBLE,"IATUTC",&temp , NULL, &status);
    if (status) {
        fits_report_error(stderr,status);
        fprintf(stderr,"%s: WARNING: failed to add keywords to antenna table\n",__func__);
        fits_clear_errmsg(); status=0;
    }

    /* write data row by row. CFITSIO automatically adjusts the size of the table */
    for(i=0; i < array->n_ant; i++) {
        ptemp = array->antennas[i].name; /* need this to get a **char in the next line */
        fits_write_col(fptr,TSTRING,1,i+1,1,1,&ptemp, &status);  /* ANNAME */
        fits_write_col_dbl(fptr,2,i+1,1,3,array->antennas[i].xyz_pos, &status); /* STABXYZ */
        itemp = i+1;
        fits_write_col_int(fptr,3,i+1,1,1,&itemp,&status);  /* NOSTA */
        itemp = array->antennas[i].mount_type;
        fits_write_col_int(fptr,4,i+1,1,1,&itemp,&status);  /* MNTSTA */
        ptemp = array->antennas[i].pol_typeA;
        fits_write_col_str(fptr,6,i+1,1,1,&ptemp,&status);  /* POLTYA */
        fits_write_col_flt(fptr,7,i+1,1,1,&(array->antennas[i].pol_angleA),&status); /* POLAA */
        fits_write_col_flt(fptr,8,i+1,1,1,&(array->antennas[i].pol_calA),&status); /* POL calA */
        ptemp = array->antennas[i].pol_typeB;
        fits_write_col_str(fptr,9,i+1,1,1,&ptemp,&status);  /* POLTYB */
        fits_write_col_flt(fptr,10,i+1,1,1,&(array->antennas[i].pol_angleB),&status); /*POLAB */
        fits_write_col_flt(fptr,11,i+1,1,1,&(array->antennas[i].pol_calB),&status); /*POL calB */
    }
    if (status) {
        fprintf(stderr,"Problems writing the AN table\n");
        fits_report_error(stderr,status);
    }
    return status;
}

/**************************
 ! NAME:      JD_to_Cal
 ! PURPOSE:   convert a Julian date to a Gregorian calendar date
 ! ARGUMENTS: jd: input julian date in double format
 !            year, month, day: output calendar date. day/months start at 1.
 ! RETURNS:   void
 **************************/
/* follows from: http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html */
void JD_to_Cal(double jd, int *year, int *month, int *day) {
    int z,w,x,a,b,c,d,e,f;

    z = jd+0.5;
    w = (z-1867216.25)/36524.25;
    x = w/4;
    a = z+1+w-x;
    b = a+1524;
    c = (b-122.1)/365.25;
    d = 365.25*c;
    e = (b-d)/30.6001;
    f = 30.6001*e;
    *day = b-d-f;
    while (e-1 > 12) e-=12;
    *month = e-1;
    *year = c-4715-((e-1)>2?1:0);
    /*  printf("z:%d, w:%d, x: %d, a: %d, b: %d, c: %d, e: %d, f: %d\n",z,w,x,a,b,c,e,f);*/
}


/**************************
 ! NAME:      Cal_to_JD
 ! PURPOSE:   convert a Gregorian calendar date to a Julian date
 ! ARGUMENTS: jd: output julian date in double format
 !            year, month, day: input calendar date. day/months start at 1.
 ! RETURNS:   void
 **************************/
/* follows from: http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html */
void Cal_to_JD(int year, int month, int day, double *jd) {
    int L;

    if (month < 1 || day < 1) {
        fprintf(stderr,"Cal_to_JD WARNING: month (%d) or day (%d) in conversion is less than 1\n",month, day);
    }

    /* this one ripped off from IDL */
    L = (month-14)/12;        //   In leap years, -1 for Jan, Feb, else 0
    *jd = day - 32075 + 1461*(year+4800+L)/4 +  367*(month - 2-L*12)/12 - 3*((year+4900+L)/100)/4 - 0.5;
}


/**************************
 ! NAME:      printUVData
 ! PURPOSE:   debugging utility to print the uvdata structure.
 ! ARGUMENTS: data: input data struct
 !            fp: open file pointer. can also be "stderr" or "stdout".
 ! RETURNS:   void
 **************************/
void printUVData(uvdata *data,FILE *fp) {
  int i,j,k,l;

  fprintf(fp,"Source:  \t<%s>\n",data->source->name);
  fprintf(fp,"OBSRA:   \t%f\n",data->source->ra);
  fprintf(fp,"OBSDEC:  \t%f\n",data->source->dec);
  fprintf(fp,"Num Pols:\t%d\n",data->n_pol);
  fprintf(fp,"Num Freq:\t%d\n",data->n_freq);
  fprintf(fp,"Cen Freq:\t%g\n",data->cent_freq);
  fprintf(fp,"Del Freq:\t%g\n",data->freq_delta);
  fprintf(fp,"Num times:\t%d\n",data->n_vis);
  for (i=0; i< data->n_vis; i++) {
    fprintf(fp,"\nScan %d. Julian date:\t%g, Num BL:\t\t%d ",i,data->date[i],data->n_baselines[i]);
    for (j=0; j<data->n_baselines[i]; j++) {
      fprintf(fp,"\nUVW: %g,%g,%g. Baseline: %g\n",data->u[i][j],data->v[i][j],data->w[i][j],data->baseline[i][j]);
      for(k=0; k<data->n_freq; k++) { 
        fprintf(fp,"    Freq: %d\n",k);
        for(l=0;l<data->n_pol; l++) {
          fprintf(fp,"\tPol %d Val,wt: (%g,%g),%g\n",l,data->visdata[i][(j*data->n_freq*data->n_pol+k*data->n_pol+l)*2  ],
                                                     data->visdata[i][(j*data->n_freq*data->n_pol+k*data->n_pol+l)*2+1],
                                                  data->weightdata[i][j*data->n_freq*data->n_pol+k*data->n_pol+l]);
        }
      }
    }
  }
}


/**************************
 ! NAME:      printAntennaData
 ! PURPOSE:   debugging utility to print the uvdata antenna structure.
 ! ARGUMENTS: data: input data struct
 !            fp: open file pointer. can be "stderr" or "stdout".
 ! RETURNS:   void
 **************************/
void printAntennaData(array_data *array,FILE *fp) {
    int i;

    fprintf(fp,"Array X,Y,Z: %.10g,%.10g,%.10g\n",array->xyz_pos[0],array->xyz_pos[1],array->xyz_pos[2]);
    for (i=0; i< array->n_ant; i++) {
        fprintf(fp,"Index: %d, Name: %s, id: %d, XYZ: (%g,%g,%g), Mount: %d\n",
	        i,array->antennas[i].name, array->antennas[i].station_num, array->antennas[i].xyz_pos[0],
	        array->antennas[i].xyz_pos[1],array->antennas[i].xyz_pos[2],
	        array->antennas[i].mount_type);
    }
}


/**************************
 ! NAME:      printHeader
 ! PURPOSE:   debugging utility to print the header of the current HDU from an open FITS file.
 ! ARGUMENTS: fptr: pointer to open FITS file.
 ! RETURNS:   void
 **************************/
void printHeader(fitsfile *fptr) {
    int i,status=0,nkeys=0;
    char *str=NULL,line[81];

    fits_hdr2str(fptr,TRUE,NULL,0,&str,&nkeys,&status);
    printf("Current header. There are %d keys\n",nkeys);
    for(i=0; i<nkeys; i++) {
        memset(line,'\0',81);
        strncpy(line,str+i*80,80);
        printf("%s\n",line);
    }
    if (str!=NULL) free(str);
}


/**************************
 ! NAME:      readAntennaTable
 ! PURPOSE:   extract data we care about from the antenna table in the FITS file
 ! ARGUMENTS: fptr: pointer to open FITS file.
 !            obj: pointer to existing uvdata struct. Assumes no existing antenna/array data in the struct.
 ! RETURNS:   integer: 0 for success. returns CFITSIO error codes or -1 for memory allocation problems.
 **************************/
int readAntennaTable(fitsfile *fptr,array_data *array) {
    int status=0,ncols,i,curr_col;
    long nrows;
    char **str_names=NULL;

    /* get the properties of the table */
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading number of rows/cols\n",status);
        goto EXIT;
    }
    if (debug) printf("AN table has %d rows, %d columns\n",(int)nrows,ncols);

    /* make space for 'array' struct and an array of antenna structs */
    if (array==NULL) {
        fprintf(stderr,"array is NULL. You forgot to malloc it\n");
        return 1;
    }
    if(array->antennas != NULL) {
        fprintf(stderr,"readAntennaTable WARNING: existing antennas array will be overwritten and lost. This is a memory leak.\n");
    }
    array->antennas = calloc(nrows,sizeof(ant_table));
    str_names = calloc(nrows+1,sizeof(char *));
    if (array->antennas==NULL || str_names==NULL) {
        fprintf(stderr,"readAntennaTable: no malloc\n");
        status = -1; goto EXIT;
    }
    array->n_ant = nrows;

    /* load the array X,Y,Z location */
    fits_read_key(fptr,TDOUBLE,"ARRAYX",&(array->xyz_pos[0]),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"ARRAYY",&(array->xyz_pos[1]),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"ARRAYZ",&(array->xyz_pos[2]),NULL,&status);
    if (status != 0) {
        fprintf(stderr,"readAntennaTable: status %d reading array XYZ\n",status);
        goto EXIT;
    }
  
    /* get ANNAME */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = array->antennas[i].name;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "ANNAME", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no ANNAME column in antenna table\n");
        goto EXIT;
    }

    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d name: %s\n",i+1,array->antennas[i].name);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column ANNAME\n",status);
    }

    /* get STABXYZ */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "STABXYZ", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no STABXYZ column in antenna table\n");
        goto EXIT;
    }

    for (i=0; i<nrows; i++) {
        fits_read_col_dbl(fptr, curr_col ,i+1, 1, 3, 0.0, array->antennas[i].xyz_pos, NULL, &status);
        if (debug) printf("Antenna %d pos: %f,%f,%f\n",i+1,array->antennas[i].xyz_pos[0],
                         array->antennas[i].xyz_pos[1],array->antennas[i].xyz_pos[2]);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for STABXYZ\n",status);
        }
    }

    /* get NOSTA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "NOSTA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no NOSTA column in antenna table\n");
        goto EXIT;
    }

    for (i=0; i<nrows; i++) {
        fits_read_col_int(fptr, curr_col ,i+1, 1, 1, 0, &(array->antennas[i].station_num), NULL, &status);
        if (debug) printf("Antenna %d num: %d\n",i+1,array->antennas[i].station_num);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for NOSTA\n",status);
        }
    }

    /* get MNTSTA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "MNTSTA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no MNTSTA column in antenna table\n");
        goto EXIT;
    }

    for (i=0; i<nrows; i++) {
        fits_read_col_int(fptr, curr_col ,i+1, 1, 1, 0, &(array->antennas[i].mount_type), NULL, &status);
        if (debug) printf("Antenna %d mount type: %d\n",i+1,array->antennas[i].mount_type);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for MNTSTA\n",status);
        }
    }

    /* get POLTYA */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = array->antennas[i].pol_typeA;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLTYA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLTYA column in antenna table\n");
        goto EXIT;
    }

    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d POLTYA: %s\n",i+1,array->antennas[i].pol_typeA);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column POLTYA\n",status);
    }

    /* get POLTYB */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = array->antennas[i].pol_typeB;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLTYB", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLTYB column in antenna table\n");
        goto EXIT;
    }

    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d POLTYB: %s\n",i+1,array->antennas[i].pol_typeB);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column POLTYB\n",status);
    }

    /* get POLAA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLAA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLAA column in antenna table\n");
        goto EXIT;
    }

    for (i=0; i<nrows; i++) {
        fits_read_col_flt(fptr, curr_col ,i+1, 1, 1, 0, &(array->antennas[i].pol_angleA), NULL, &status);
        if (debug) printf("Antenna %d pol angle A: %f\n",i+1,array->antennas[i].pol_angleA);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for POLAA\n",status);
        }
    }

    /* get POLAB */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLAB", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLAB column in antenna table\n");
        goto EXIT;
    }

    for (i=0; i<nrows; i++) {
        fits_read_col_flt(fptr, curr_col ,i+1, 1, 1, 0, &(array->antennas[i].pol_angleB), NULL, &status);
        if (debug) printf("Antenna %d pol angle B: %f\n",i+1,array->antennas[i].pol_angleB);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for POLAB\n",status);
        }
    }

 EXIT:
    if(str_names!=NULL) free(str_names);
    return status;
}


/**************************
 ! NAME:      freeUVFITSdata
 ! PURPOSE:   free space associated with a uvfits data struct
 ! ARGUMENTS: data: pointer to data struct
 ! RETURNS:   void
 **************************/
void freeUVFITSdata(uvdata *data) {
  int i;

  if (data->fq != NULL) free(data->fq);
  if (data->date !=NULL) free(data->date);
  if (data->array->antennas!=NULL) free(data->array->antennas);
  if (data->array !=NULL) free(data->array);
  if (data->source != NULL) free(data->source);
  if (data->n_baselines !=NULL) free(data->n_baselines);
  for(i=0; i<data->n_vis; i++) {
    if (data->visdata!=NULL && data->visdata[i] !=NULL) free(data->visdata[i]);
    if (data->weightdata!=NULL && data->weightdata[i] !=NULL) free(data->weightdata[i]);
    if (data->u!=NULL && data->u[i] !=NULL) free(data->u[i]);
    if (data->v!=NULL && data->v[i] !=NULL) free(data->v[i]);
    if (data->w!=NULL && data->w[i] !=NULL) free(data->w[i]);
    if (data->baseline!=NULL && data->baseline[i]!=NULL) free(data->baseline[i]);
  }
  if (data->visdata!=NULL) free(data->visdata);
  if (data->weightdata!=NULL) free(data->weightdata);
  if (data->u != NULL) free(data->u);
  if (data->v != NULL) free(data->v);
  if (data->w != NULL) free(data->w);
  if (data->baseline !=NULL) free(data->baseline);
  free(data);
}


/* encode the baseline. Use the miriad convention to handle more than 255 antennas (up to 2048).
   this is backwards compatible with the standard UVFITS convention. Antenna indices start at 1 */
void EncodeBaseline(int b1, int b2, float *result) {
  if (b2 > 255) {
    *result = b1*2048 + b2 + 65536;
  } else {
    *result = b1*256 + b2;
  }
}


/* decode a baseline using same (extended) miriad format as above */
void DecodeBaseline(float blcode, int *b1, int *b2) {
  if (blcode > 65536) {
    blcode -= 65536;
    *b2 = (int) blcode % 2048;
    *b1 = (blcode - *b2)/2048;
  }
  else {
    *b2 = (int) blcode % 256;
    *b1 = (blcode - *b2)/256;
  }
}


/**************************
 ! NAME:      writeFQData
 ! PURPOSE:   write the frequency ("FQ") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeFQData(fitsfile *fptr, uvdata *data) {
    double f;
    int status=0,i, id;
    int NUM_COLS = 5;
    char *col_names[] = {"FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"};
    char *col_format[] = {"1I","1D","1E","1E","1I"};
    char *col_units[] = {"","Hz","Hz","Hz",""};

    /* create a new binary table */
    fits_create_tbl(fptr,BINARY_TBL,0, NUM_COLS ,col_names, col_format,col_units,"AIPS FQ", &status);

    /* write table header info: just the number of IFs */
    fits_update_key(fptr,TINT,"NO_IF", &(data->n_IF), NULL, &status);

    // write data row by row. CFITSIO automatically adjusts the size of the table
    // Array index is used to create the identifier (index+1 so starts at 1 not 0)
    // Arg order: fptr, col#, first row, first elem, num elems, value, status
    // Value is passed as an address, as uses arrays for multiple row writing
    for(i=0; i<data->n_IF; i++) {
        // Changed from centre of entire IF to just centre of first channel
        f = data->fq[i].freq + data->fq[i].chbw/2 - data->cent_freq;
        //printf("if #%d, f=%f\n", (i+1), f);
        id = 1;//i+1;

        fits_write_col_int(fptr,1,i+1,1,1,&(id),&status); /* FRQSEL */
        fits_write_col_dbl(fptr,2,i+1,1,1,&(f), &status); /* IF FREQ */
        fits_write_col_flt(fptr,3,i+1,1,1,&(data->fq[i].chbw), &status); /* CH WIDTH */
        fits_write_col_flt(fptr,4,i+1,1,1,&(data->fq[i].bandwidth), &status); /* TOTAL BANDWIDTH */
        fits_write_col_int(fptr,5,i+1,1,1,&(data->fq[i].sideband), &status); /* SIDEBAND */
    }

    return status;
}

/**************************
 ! NAME:      readFQTable
 ! PURPOSE:   read the IF FQ data from the table "AIPS FQ"
 ! ARGUMENTS: fptr: pointer to open FITS file.
 !            obj: pointer to existing uvdata struct. Assumes no existing fq array.
 ! RETURNS:   integer: 0 for success. returns CFITSIO error codes or -1 for memory allocation problems.
 **************************/
int readFQTable(fitsfile *fptr,uvdata *obj) {
    int status=0,ncols,i,j,curr_col, num_IF=0;
    long nrows;
    int NUM_COLS = 5;
    char *col_names[] = {"FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"};

    /* get the properties of the table */
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    if (status !=0) {
        fprintf(stderr,"readFQTable: status %d reading number of rows/cols\n",status);
        goto EXIT;
    }
    if (debug) printf("FQ table has %d rows, %d columns\n",(int)nrows,ncols);
    if (nrows != obj->n_IF) {
        fprintf(stderr,"ERROR: mismatch between number of IFs in primary HDU axes (%d) and FQ table (%d)\n",obj->n_IF,(int)nrows);
        return -1;
    }

    /* make space for fqtable struct and an array of antenna structs */
    if(obj->fq != NULL)
        fprintf(stderr,"readFQTable WARNING: existing fq array will be overwritten and lost. This is a memory leak.\n");

    obj->fq = calloc(nrows,sizeof(fq_table));
    if (obj->fq==NULL) {
        fprintf(stderr,"readFQTable: no malloc for fqtable\n");
        status = -1; goto EXIT;
    }

    /* load the NO_IF parameter */
    fits_read_key(fptr,TINT,"NO_IF",&num_IF,NULL,&status);
    if (status != 0) {
        fprintf(stderr,"readFQTable: status %d reading NO_IF\n",status);
        goto EXIT;
    }
    // Check the number of rows matches the IFs already read, and the IF parameter matches also
    if(nrows != obj->n_IF || num_IF != obj->n_IF) {
        fprintf(stderr,
                "readFQTable: ERROR: n_rows in IF table or NO_IF value does not match number of IFs in main data\n");
        status = -1; goto EXIT;
    }
   
    // Check for the existence of the necessary columns
    for(i=0; i<NUM_COLS; i++) {
        curr_col = 0;
        fits_get_colnum(fptr,CASEINSEN, col_names[i], &curr_col, &status);
        if (curr_col==0 || status !=0) {
            fprintf(stderr,"readFQTable: no %s column in FQ table\n", col_names[i]);
            goto EXIT;
        }
    
        // Read the data out of each row
        for (j=0; j<nrows; j++) {
            if(i==1) {
                fits_read_col_dbl(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].freq), NULL, &status);
                if (debug) printf("FQ%d, %s: %lf\n", j+1, col_names[i], obj->fq[j].freq);
            }
            else if(i==2) {
                fits_read_col_flt(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].chbw), NULL, &status);
                if (debug) printf("FQ%d, %s: %f\n", j+1, col_names[i], obj->fq[j].chbw);
            }
            else if(i==3) {
                fits_read_col_flt(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].bandwidth), NULL, &status);
                if (debug) printf("FQ%d, %s: %f\n", j+1, col_names[i], obj->fq[j].bandwidth);
            }
            else if(i==4) {
                fits_read_col_int(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].sideband), NULL, &status);
                if (debug) printf("FQ%d, %s: %d\n", j+1, col_names[i], obj->fq[j].sideband);
            }

            if (status !=0) {
                fprintf(stderr,"readFQTable: status %d reading column for %s\n",status, col_names[i]);
            }
        }
    }

    // Convert back to absolute frequencies (in file as relative)
    for(i=0; i<obj->n_IF; i++) {
        obj->fq[i].freq = obj->fq[i].freq - obj->fq[i].bandwidth/2 + obj->cent_freq;
    }

EXIT:
    return status;
}


/***************************
  convert Geodetic lat/lon/height to XYZ coords
  uses constants from WGS84 system, defined in header
***************************/
void Geodetic2XYZ(double lat_rad, double lon_rad, double height_meters, double *X, double *Y, double *Z) {
    double s_lat,s_lon,c_lat,c_lon;
    double chi;

    s_lat = sin(lat_rad); c_lat = cos(lat_rad);
    s_lon = sin(lon_rad); c_lon = cos(lon_rad);
    chi = sqrt(1.0 - E_SQUARED*s_lat*s_lat);

    *X = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*c_lon;
    *Y = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*s_lon;
    *Z = (EARTH_RAD_WGS84*(1.0-E_SQUARED)/chi + height_meters)*s_lat;
}


/*****************************
 convert local topocentric East, North, Up units to Earth-centred earth-fixed cartesian
 coords relative to a specified center. Latitude is geodetic, not geocentric
 lat,lon in radian. E,N,U in same units as X,Y,Z
********************************/
void ENH2XYZ_absolute(double E,double N, double H, double lat_rad, double lon_rad, double *X, double *Y, double *Z) {
    double c_lat,c_lon,s_lat,s_lon;
 
    s_lat = sin(lat_rad); c_lat = cos(lat_rad);
    s_lon = sin(lon_rad); c_lon = cos(lon_rad);

    *X = -s_lon*E - s_lat*c_lon*N + c_lat*c_lon*H;
    *Y =  c_lon*E - s_lat*s_lon*N + c_lat*s_lon*H;
    *Z =            c_lat*N       + s_lat*H;
}

/*********************************
  convert coords in local topocentric East, North, Height units to
  'local' XYZ units. Local means Z point north, X points through the equator from the geocenter
  along the local meridian and Y is East.
  This is like the absolute system except that zero lon is now
  the local meridian rather than prime meridian.
  Latitude is geodetic, in radian.
  This is what you want for constructing the local antenna positions in a UVFITS antenna table.
**********************************/
void ENH2XYZ_local(double E,double N, double H, double lat, double *X, double *Y, double *Z) {
  double sl,cl;

  sl = sin(lat);
  cl = cos(lat);
  *X = -N*sl + H*cl;
  *Y = E;
  *Z = N*cl + H*sl;
}
