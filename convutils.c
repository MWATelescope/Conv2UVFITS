#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <ctype.h>
#include <star/pal.h>
#include "convutils.h"
#include "fitsio.h"

// 2022-03 replaced old SLA lib with starlink PAL equivalents
//#include "slalib.h"

int qsort_compar_float(const void *p1, const void *p2);

/* use common debug file handle from main program */
extern FILE *fpd;

/* private globals */
static int debug_flag=0,debug=0;

void setConvDebugLevel(int level) {
    debug=level;
    if (debug>0) {
        fprintf(fpd,"%s: New debug level %d\n",__FILE__,level);
    }
}


/* comparison function for qsort in the autoflagger */
int qsort_compar_float(const void *p1, const void *p2) {
    float val1,val2;

    val1 = *((float *)p1);
    val2 = *((float *)p2);
    if (val1 < val2) return -1;
    if (val1 > val2) return 1;
    return 0;
}


/* simple median-based flagging based on autocorrelations */
/* to avoid making an additional copy of all the data, which needs to be sorted,
   we pass through the data once for each ant/pol combination
   and extract the data vs time for each single channel. This is then sorted by magnitude
   per channel and the variance estimate made. */
int autoFlag(uvdata *uvdata, float sigma, int n_neighbours,int flag_type) {
    int status=0,ant,pol,chan,t,bl,ant1,ant2,visindex;
    float *chan_sigma=NULL, *chan_data=NULL, *chan_median=NULL, *chan_sorted=NULL, sig, med;
    char  *flags=NULL;
    float local_medians[AUTOFLAG_N_NEIGHBOURS*2+1],local_stdevs[AUTOFLAG_N_NEIGHBOURS*2+1];
    
    /* init */
    chan_sigma    = calloc(uvdata->n_freq, sizeof(float));
    chan_median   = calloc(uvdata->n_freq, sizeof(float));
    chan_data     = malloc(uvdata->n_vis*uvdata->n_freq*sizeof(float));
    chan_sorted   = malloc(uvdata->n_vis*uvdata->n_freq*sizeof(float));
    flags         = malloc(uvdata->n_vis*uvdata->n_freq*sizeof(char));
    if(chan_sigma ==NULL || chan_data==NULL || chan_median==NULL || flags ==NULL || chan_sorted==NULL) {
        fprintf(stderr,"ERROR: autoFlag: no malloc\n");
        exit(1);
    }


    /* each pol on each antenna is treated separately, then all baselines for affected visibilities for
       the ant/pol are flagged */
    for(ant=1; ant <= uvdata->array->n_ant; ant++) {
    for(pol=0; pol < (int)ceil(uvdata->n_pol/2.0); pol++) {
        // make pass through data to determine variance estimates
        for(bl=0; bl < uvdata->n_baselines[0]; bl++) {
           // is this the autocorrelation we're looking for?
           DecodeBaseline(uvdata->baseline[0][bl], &ant1, &ant2);
           for(t=0; t<uvdata->n_vis; t++) {
                for(chan=0; chan < uvdata->n_freq; chan++) {
                   if( ant1 == ant2 && ant1 == ant) {
                        // extract data for this time,pol,channel
                        visindex = bl*uvdata->n_pol*uvdata->n_freq + chan*uvdata->n_pol + pol;
                        chan_data[chan*uvdata->n_vis + t] = uvdata->visdata[t][visindex*2];
                    }
                }
            }
        }
        // make a copy of the data that will be sorted
        memcpy(chan_sorted,chan_data,uvdata->n_vis*uvdata->n_freq*sizeof(float));
        for(chan=0; chan < uvdata->n_freq; chan++) {
            /* sort the data by magnitude */
            qsort(chan_sorted+chan*uvdata->n_vis,uvdata->n_vis,sizeof(float),qsort_compar_float);
            /* estimate channel variance based on interquartile range . ((d[0.75*n_vis] - d[0.25*n_vis])/1.35)^2*/
            sig = ( (chan_sorted[(int)(chan*uvdata->n_vis+0.75*uvdata->n_vis)]-chan_sorted[(int)(chan*uvdata->n_vis+0.25*uvdata->n_vis)])/1.35);
            med = chan_sorted[(int)(chan*uvdata->n_vis + 0.5*uvdata->n_vis)];
            chan_sigma[chan]  = fabs(sig);
            chan_median[chan] = med;
            if(debug_flag) fprintf(fpd,"ant: %d, pol: %d, chan: %d. Median: %g, stddev: %g\n",ant,pol,chan,med,sig);
        }

        /* next pass through data to apply flags */
        /* clear existing flags */
        memset(flags,'\0',uvdata->n_vis*uvdata->n_freq);
        for(chan=0; chan < uvdata->n_freq; chan++) {
            int cpc=0;  // channels per coarse channel
            int lower,upper,n_points,n_rej=0;
            float local_median,local_stdev,val;

            // set bounds for lower and upper chans for neighbor comparisions
            lower = chan-n_neighbours;
            if (lower < 0) lower=0;
            upper = chan+n_neighbours;
            if (upper >= uvdata->n_freq) upper = uvdata->n_freq-1;
            switch(flag_type) {

                case 1:
                    // generic case - don't do anything special
                    break;
                case 2: 
                    // MWA 40kHz channels. do not use the n_rej edge channels to compare to neighbours
                    // there are 32 fine 40kHz channels per coarse 1.28MHz channel. So we ignore the edge
                    // channels within a coarse chan
                    n_rej = 6;
                    cpc=32;
                    if (chan%cpc < n_rej || chan%cpc >= cpc-n_rej) lower = upper = chan;
                    if (chan%cpc > n_rej && lower%cpc < n_rej) lower += n_rej-lower%cpc;
                    if (chan%cpc < cpc-n_rej && upper%cpc >= cpc-n_rej) upper -= (upper%cpc)-(cpc-1-n_rej);
                    break;
                case 3:
                    // MWA 10kHz channels. do not use the n_rej edge channels to compare to neighbours
                    n_rej = 20;
                    cpc=128;
                    if (chan%cpc < n_rej || chan%cpc >= cpc-n_rej) lower = upper = chan;
                    if (chan%cpc > n_rej && lower%cpc < n_rej) lower += n_rej-lower%cpc;
                    if (chan%cpc < cpc-n_rej && upper%cpc >= cpc-n_rej) upper -= (upper%cpc)-(cpc-1-n_rej);
                    break;
                default:
                    fprintf(stderr,"Unsupported/unknown autoflag mode %d. Not flagging\n",flag_type);
                    goto EXIT;
            }
            n_points = upper-lower+1;
            // extract median and stdev channel values for channel plus neighbours. By using the neighbours, we avoid
            // being contaminated by long-lived, narrowband RFI.
            memcpy(local_medians,chan_median+lower,sizeof(float)*n_points);
            memcpy(local_stdevs,chan_sigma+lower,sizeof(float)*n_points);
            // sort the local medians and stdevs to get the median of medians and median of stdevs
            qsort(local_medians,n_points,sizeof(float),qsort_compar_float);
            qsort(local_stdevs,n_points,sizeof(float),qsort_compar_float);
            // use the local median,stdev as the median,stdev for this channel
            local_median = local_medians[n_points/2];
            local_stdev  = local_stdevs[n_points/2];
            if(debug_flag) fprintf(fpd,"ant: %d, pol: %d, chan: %d. lower: %d, upper: %d, Local median: %g, local stdev: %g\n",ant,pol,chan,lower,upper,local_median,local_stdev);
            // now scan data to find outliers and flag them
            for(t=0; t<uvdata->n_vis; t++) {
                val = chan_data[chan*uvdata->n_vis + t];
                if ( fabs(val - local_median) > sigma*local_stdev) {
                    if(debug_flag) fprintf(fpd,"flagging ant: %d, pol: %d, chan: %d, time: %d. Val: %g\n",ant,pol,chan,t,val);
                    flag_antenna(uvdata,ant,pol,chan,t);
                    flags[chan + t*uvdata->n_freq] = 1;
                }
            }
        }
        // make an image of the flags for this input (for debugging purposes)

        if (debug_flag) {
            char filename[100];
            long fpixel[2];
            fitsfile *fp=NULL;
            //dump flags
            sprintf(filename,"flags_inp%02d_pol%d.fits",ant,pol);
            remove(filename);
            fits_create_file(&fp,filename,&status);
            if (status !=0) {
                fprintf(stderr,"autoFlag: cannot create file %s\n",filename);
                return 1;
            }
            fprintf(fpd,"created flag dump file %s\n",filename);
            fpixel[0] = uvdata->n_freq;
            fpixel[1] = uvdata->n_vis;
            fits_create_img(fp, BYTE_IMG, 2, fpixel, &status);
            fpixel[0] = 1;
            fpixel[1] = 1;
            fits_write_pix(fp, TBYTE, fpixel, uvdata->n_freq*uvdata->n_vis,flags, &status);
            fits_close_file(fp,&status);
            // dump data
            sprintf(filename,"data_inp%02d_pol%d.fits",ant,pol);
            remove(filename);
            fits_create_file(&fp,filename,&status);
            if (status !=0) {
                fprintf(stderr,"autoFlag: cannot create file %s\n",filename);
                return 1;
            }
            fprintf(fpd,"created flag dump file %s\n",filename);
            fpixel[1] = uvdata->n_freq;
            fpixel[0] = uvdata->n_vis;
            fits_create_img(fp, FLOAT_IMG, 2, fpixel, &status);
            fpixel[0] = 1;
            fpixel[1] = 1;
            fits_write_pix(fp, TFLOAT, fpixel, uvdata->n_freq*uvdata->n_vis,chan_data, &status);
            fits_close_file(fp,&status);
        }

    }   
    }

EXIT:
    /* free working arrays */
    if (chan_sigma!= NULL) free(chan_sigma);
    if (chan_median!=NULL) free(chan_median);
    if (chan_data != NULL) free(chan_data);
    if (chan_sorted!=NULL) free(chan_sorted);
    if (flags != NULL)     free(flags);
    return status;
}


/* flag all antennas/pols for channel chan at time index t */
/* setting a visibility weight to negative means it is flagged */
int flag_all_antennas(uvdata *uvdata, const int chan, const int t) {
    int bl,pol,visindex;
    float weight;

    for (bl=0; bl<uvdata->n_baselines[t]; bl++) {
         for(pol=0; pol < uvdata->n_pol; pol++) {
            visindex = bl*uvdata->n_pol*uvdata->n_freq + chan*uvdata->n_pol + pol;
            weight = uvdata->weightdata[t][visindex];
            if (weight > 0) uvdata->weightdata[t][visindex] = -weight;
         }
    }
    return 0;
}


/* flag all visibilities that are formed by the antenna "ant" on polarisation "pol"
   in channel "chan" at time index "t".
*/
int flag_antenna(uvdata *uvdata, const int ant, const int pol, const int chan, const int t) {
    int bl,ant1,ant2,visindex;
    float weight;

    for (bl=0; bl<uvdata->n_baselines[t]; bl++) {
        // antenna we're looking for?
        DecodeBaseline(uvdata->baseline[t][bl], &ant1, &ant2);
        if (ant1 == ant || ant2==ant) {
            // negative weights mean flagged
            visindex = bl*uvdata->n_pol*uvdata->n_freq + chan*uvdata->n_pol + pol;
            weight = uvdata->weightdata[t][visindex];
            if (weight > 0) {
                uvdata->weightdata[t][visindex] = -weight;
            }
            // and the other product with this pol, if applicable
            if (uvdata->n_pol > 1) {
                visindex += 2;
                weight = uvdata->weightdata[t][visindex];
                if (weight > 0) {
                    uvdata->weightdata[t][visindex] = -weight;
                }
            }
        }  
    }
    return 0;
}


/* apply a simple global flags file to the data */
int applyFlagsFile(char *filename,uvdata *data) {
    FILE *fp=NULL;
    char line[MAX_LINE],key[80];
    int val;

    if ((fp=fopen(filename,"r")) ==NULL) {
        fprintf(stderr,"Cannot open flags file %s\n",filename);
        return 1;
    }
    if(debug) fprintf(fpd,"Applying global flags in flags file: %s\n",filename);

    /* process each line in the file */
    while( fgets(line,MAX_LINE,fp) != NULL) {
        /* skip comment and blank lines */
        if (line[0] == '\0' || line[0] == '\n' || line[0]=='#') continue;

        /* process the line */
        sscanf(line,"%s %d",key,&val);
        if (strncmp(CHAN_ALL_ANT_ALL_TIME,line,strlen(CHAN_ALL_ANT_ALL_TIME))==0) {
            int t;

            /* flag a channel for all antennas and times */
            /* check that channel is valid */
            if (val >= data->n_freq) {
                fprintf(stderr,"ERROR: asked to flag channel %d, but only have %d\n",val,data->n_freq);
                fprintf(stderr,"Offending line: %s\n",line);
            }
            if(debug) fprintf(fpd,"Flagging chan %d on all antennas\n",val);
            for (t=0; t<data->n_vis; t++) flag_all_antennas(data, val, t);
        }
    }

    if (fp !=NULL) fclose(fp);
    return 0;
}


/******************************
 read the station locations from a text
 file and populate the antenna positions in the
 data structure.
*******************************/
int readArray(char *filename, const double lat_radian, const double arr_lon_rad, array_data *array) {
  FILE *fp=NULL;
  char line[MAX_LINE];
  int index=0,nscan;
  double east=0,north=0,height=0;
  ant_table *antennas;

  if( (fp=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: readArray: failed to open array file <%s>\n",filename);
    return 1;
  }

  antennas = array->antennas;
  array->arr_lat_rad = lat_radian;
  array->arr_lon_rad = arr_lon_rad;

  /* scan through lines. convert east,north,height units to XYZ units */
  while((fgets(line,MAX_LINE-1,fp)) !=NULL) {
    if(line[0]=='\n' || line[0]=='#' || line[0]=='\0') continue; // skip blank/comment lines
    nscan = sscanf(line,"%8s %lf %lf %lf",antennas[index].name,&east,&north,&height);
    if(nscan != 4) {
        fprintf(stderr,"Failed scanning antenna file with line: <%s>\n",line);
        return 1;
    }
    ENH2XYZ_local(east,north,height,lat_radian,antennas[index].xyz_pos,antennas[index].xyz_pos+1,antennas[index].xyz_pos+2);
    if (debug) {
      fprintf(fpd,"ant %s. Pos (ENH) (%g,%g,%g).\tPos (XYZ): (%g,%g,%g)\n",antennas[index].name, east,north,height,
          antennas[index].xyz_pos[0],antennas[index].xyz_pos[1],antennas[index].xyz_pos[2]);
    }
    antennas[index].station_num = index;
    index++;
    if (index > MAX_ANT) {
        fprintf(stderr,"ERROR: there are too many antennas. Increase MAX_ANT and recompile.\n");
        exit(1);
    }
  }
  array->n_ant = index;
  fclose(fp);
  return 0;
}


/******************************
 read the mapping between antennas and correlator inputs.
*******************************/
int readInputConfig(char *filename, InpConfig *inp) {
  FILE *fp=NULL;
  char line[MAX_LINE],pol_char,cable_len[128];
  int index=0,dummy,nscan,i,inp_flag;

  if( (fp=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: readInputConfig: failed to open array file <%s>\n",filename);
    return 1;
  }
  
  for (i=0; i<MAX_ANT*2; i++) inp->inpFlag[i] = 0;

  /* scan through lines.  */
  while((fgets(line,MAX_LINE-1,fp)) !=NULL) {
    if(line[0]=='\n' || line[0]=='#' || line[0]=='\0') continue; // skip blank/comment lines

    inp_flag=0;

    nscan = sscanf(line,"%d %d %c %s %d",&dummy,inp->ant_index+index,&pol_char,cable_len,&inp_flag);
    if(nscan < 4) {
        fprintf(stderr,"Failed scanning instr config file with line: <%s>\n",line);
        return 1;
    }

    // decode the string with the cable length. for a prefix of "EL_" this means the value is an electrical length
    // not a physical one, so no velocity factor should be applied.
    if (strncmp(cable_len,"EL_",3)==0) {
        inp->cable_len_delta[index] = atof(cable_len+3);
    }
    else {
        inp->cable_len_delta[index] = atof(cable_len)*VEL_FACTOR;
    }

    inp->inpFlag[index] = inp_flag;
    
    inp->pol_index[index] = decodePolChar(pol_char);
    if (debug) {
        fprintf(fpd,"input: %d is antenna %d with pol %d. Length delta: %g\n",index,inp->ant_index[index],
                    inp->pol_index[index],inp->cable_len_delta[index]);
    }
    index++;
  }
  inp->n_inputs = index;
  fclose(fp);
  return 0;
}


/***************************
 ***************************/
int readHeader(char *header_filename, Header *header) {
    FILE *fp=NULL;
    char line[MAX_LINE],key[MAX_LINE],value[MAX_LINE],junk[MAX_LINE];
    
    if((fp=fopen(header_filename,"r"))==NULL) {
        fprintf(stderr,"ERROR: failed to open obs metadata file <%s>\n",header_filename);
        exit(1);
    }
    memset(header,'\0', sizeof(Header));
    header->corr_type = 'N';
    header->ha_hrs_start = -99.0;
    header->ra_hrs = -99.0;
    header->dec_degs = -99.0;
    header->ref_el = M_PI/2.0;
    header->ref_az = 0.0;
    strcpy(header->telescope,"MWA");
    strcpy(header->instrument,"128T");
    header->invert_freq = 0;    // correlators have now been fixed
    header->conjugate = 0;      // just in case.
    header->geom_correct = 1;   // default, stop the fringes

    while((fgets(line,MAX_LINE-1,fp)) !=NULL) {
        if(line[0]=='\n' || line[0]=='#' || line[0]=='\0') continue; // skip blank/comment lines

        if (sscanf(line,"%s %s %s",key,value,junk) < 2) {
            fprintf(stderr,"WARNING: failed to make 2 conversions on line: %s\n",line);
        }
        if (strncmp(key,"FIELDNAME",MAX_LINE)==0) strncpy(header->field_name,value,SIZE_SOURCE_NAME);
        if (strncmp(key,"TELESCOPE",MAX_LINE)==0) strncpy(header->telescope,value,SIZ_TELNAME);
        if (strncmp(key,"INSTRUMENT",MAX_LINE)==0) strncpy(header->instrument,value,SIZ_TELNAME);
        if (strncmp(key,"POL_PRODS",MAX_LINE)==0) strncpy(header->pol_products,value,SIZ_PRODNAME);
        if (strncmp(key,"N_SCANS",MAX_LINE)==0) header->n_scans = atoi(value);
        if (strncmp(key,"N_INPUTS",MAX_LINE)==0) header->n_inputs = atoi(value);
        if (strncmp(key,"N_CHANS",MAX_LINE)==0) header->n_chans = atoi(value);
        if (strncmp(key,"CORRTYPE",MAX_LINE)==0) header->corr_type = toupper(value[0]);
        if (strncmp(key,"INT_TIME",MAX_LINE)==0) header->integration_time = atof(value);
        if (strncmp(key,"FREQCENT",MAX_LINE)==0) header->cent_freq = atof(value);
        if (strncmp(key,"BANDWIDTH",MAX_LINE)==0) header->bandwidth = atof(value);
        if (strncmp(key,"INVERT_FREQ",MAX_LINE)==0) header->invert_freq = atoi(value);
        if (strncmp(key,"CONJUGATE",MAX_LINE)==0) header->conjugate = atoi(value);
        if (strncmp(key,"GEOM_CORRECT",MAX_LINE)==0) header->geom_correct = atoi(value);
        if (strncmp(key,"REF_AZ",MAX_LINE)==0) header->ref_az = atof(value)*(M_PI/180.0);
        if (strncmp(key,"REF_EL",MAX_LINE)==0) header->ref_el = atof(value)*(M_PI/180.0);
        if (strncmp(key,"HA_HRS",MAX_LINE)==0) header->ha_hrs_start = atof(value);
        if (strncmp(key,"RA_HRS",MAX_LINE)==0) header->ra_hrs = atof(value);
        if (strncmp(key,"DEC_DEGS",MAX_LINE)==0) header->dec_degs = atof(value);
        if (strncmp(key,"DATE",MAX_LINE)==0) {
            header->day = atoi(value+6); value[6]='\0';
            header->month = atoi(value+4); value[4]='\0';
            header->year = atoi(value);
        }
        if (strncmp(key,"TIME",MAX_LINE)==0) {
            header->ref_second = atof(value+4); value[4]='\0';
            header->ref_minute = atoi(value+2); value[2]='\0';
            header->ref_hour = atoi(value);
            if (debug) fprintf(fpd,"Time: %d:%d:%f\n",header->ref_hour,header->ref_minute,header->ref_second);
        }
    
    }

    /* sanity checks and defaults */
    if (header->pol_products[0]=='\0') strcpy(header->pol_products,"XXXYYXYY");
    if (header->n_scans==0) {
        header->n_scans=1;
        fprintf(stderr,"WARNING: N_SCANS unspecified. Assuming: %d\n",header->n_scans);
    }
    if(header->field_name[0] == '\0') strcpy(header->field_name,"unknown");
    if (header->n_inputs==0) {
        fprintf(stderr,"ERROR: N_INPUTS unspecified.\n");
        return EXIT_FAILURE;
    }
    if (header->n_chans==0) {
        fprintf(stderr,"ERROR: N_CHANS unspecified. \n");
        return EXIT_FAILURE;
    }
    if (header->corr_type=='N') {
        header->corr_type='B';
        fprintf(stderr,"WARNING: CORRTYPE unspecified. Assuming: %c\n",header->corr_type);
    }
    if (header->integration_time==0) {
        fprintf(stderr,"ERROR: INT_TIME unspecified.\n");
        return EXIT_FAILURE;
    }
    if (header->bandwidth==0) {
        fprintf(stderr,"ERROR: BANDWIDTH unspecified.\n");
        return EXIT_FAILURE;
    }
    if (header->cent_freq==0) {
        fprintf(stderr,"ERROR: FREQCENT unspecified. There is no default.\n");
        return 1;
    }
    if (header->ra_hrs==-99 && header->ha_hrs_start == -99.0) {
        fprintf(stderr,"ERROR: RA_HRS and HA unspecified. There is no default.\n");
        return 1;
    }
    if (header->dec_degs==-99) {
        fprintf(stderr,"ERROR: DEC_DEGS unspecified. There is no default.\n");
        return 1;
    }
    if (header->year==0 || header->month==0 || header->day==0) {
        fprintf(stderr,"ERROR: DATE unspecified. There is no default.\n");
        return 1;
    }

    if(fp!=NULL) fclose(fp);
    return 0;
}


/**************************
**************************/
void calcUVW(double ha,double dec,double x,double y,double z,double *u,double *v,double *w) {
    double sh,ch,sd,cd;

    sh = sin(ha); sd = sin(dec);
    ch = cos(ha); cd = cos(dec);
    *u  = sh*x + ch*y;
    *v  = -sd*ch*x + sd*sh*y + cd*z;
    *w  = cd*ch*x  - cd*sh*y + sd*z;
}


/***************************************
 examine the ordering of polarisation products from maps and create an
 index to order them the way miriad likes: i.e. XX, YY, XY, YX
 typically they will be XX,XY,YX,YY from MAPS
**************************************/
int createPolIndex(char *polprods, int *index) {
  int p1='\0',p2='\0',i;

  /* find the unique letters representing the pol products. i.e. X,Y or R,L or II */
  for (i=0; i<strlen(polprods); i++) {
    if (p1=='\0' && polprods[i] != '\0') {
      p1 = polprods[i];
      continue;
    }
    if (p2=='\0' && polprods[i] != '\0' && p1 != polprods[i]) {
      p2 = polprods[i];
      continue;
    }
  }
  if (debug) fprintf(fpd,"Found pol keys '%c' and '%c'\n",p1,p2);

  /* find the index of products */
  for (i=0; i<4; i++) {
    if (polprods[i*2]==p1 && polprods[i*2+1]==p1) index[0] = i;
    if (polprods[i*2]==p2 && polprods[i*2+1]==p2) index[1] = i;
    if (polprods[i*2]==p1 && polprods[i*2+1]==p2) index[2] = i;
    if (polprods[i*2]==p2 && polprods[i*2+1]==p1) index[3] = i;
  }
  if (debug) {
    for (i=0; i<4; i++) fprintf(fpd,"polindex: %d ",index[i]);
    fprintf(fpd,"\n");
  }
  return 0;
}


/**************************
***************************/
int decodePolChar(int pol_char) {
    int temp;
    temp = toupper(pol_char);
    if (temp == 'X' || temp == 'R'|| temp=='I') return 0;
    if (temp == 'Y' || temp == 'L') return 1;
    fprintf(stderr,"WARNING: Unknown pol char: <%c>\n",pol_char);
    return 0;
}


/**************************
***************************/
int decodePolIndex(int pol1, int pol2) {

    if (pol1==0 && pol2==0) return 0;
    if (pol1==1 && pol2==1) return 1;
    if (pol1==0 && pol2==1) return 2;
    if (pol1==1 && pol2==0) return 3;
    return 0;
}



/**************************
***************************/
void azel2xyz(double az, double el, double *x, double *y, double *z) {
    double sa,ca,se,ce;
    
    sa = sin(az); se = sin(el);
    ca = cos(az); ce = cos(el);

    *x = sa*ce;
    *y = ca*ce;
    *z = se;
}


/*************************
 count the number of antennas actually present in this data.
 might be less than the number of antennas in the array configuration 
**************************/
int countPresentAntennas(InpConfig *inputs) {
    int i,ant_present[MAX_ANT],total_ants=0;

    memset(ant_present,'\0',sizeof(int)*MAX_ANT);
    
    for(i=0; i<inputs->n_inputs; i++) ant_present[inputs->ant_index[i]] =1;
    for(i=0; i<MAX_ANT; i++) total_ants += ant_present[i];

    if (!ant_present[0]) {
        fprintf(stderr,"WARNING: Antenna column in instr_config file did not appear to have an antenna with index 0\n");
        fprintf(stderr,"This is probably a mistake. Please ensure antenna IDs are zero-indexed in instr_config file.\n");
        fprintf(stderr,"Continuing and hoping for the best...\n");
    }

    return total_ants;

}

/***********************
 check if an antenna is being used in the current config\
 returns 1 (true) or 0 (false)
 ***********************/
int checkAntennaPresent(InpConfig *inputs, int ant_index) {
    int i;

    for (i=0; i<inputs->n_inputs; i++) {
        if (inputs->ant_index[i] == ant_index) return 1;
    }

    return 0;
}

/**************************
***************************/
/* lmst, lmst2000 are the local mean sidereal times in radians
 * for the obs. and J2000 epochs.
 */
void precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
         double *xp, double *yp, double *zp, double lmst2000)
{
  double sep, cep, s2000, c2000;
  double xpr, ypr, zpr, xpr2, ypr2, zpr2;

  sep = sin(lmst);
  cep = cos(lmst);
  s2000 = sin(lmst2000);
  c2000 = cos(lmst2000);

  /* rotate to frame with x axis at zero RA */
  xpr = cep*x - sep*y;
  ypr = sep*x + cep*y;
  zpr = z;

  xpr2 = (rmat[0][0])*xpr + (rmat[0][1])*ypr + (rmat[0][2])*zpr;
  ypr2 = (rmat[1][0])*xpr + (rmat[1][1])*ypr + (rmat[1][2])*zpr;
  zpr2 = (rmat[2][0])*xpr + (rmat[2][1])*ypr + (rmat[2][2])*zpr;

  /* rotate back to frame with xp pointing out at lmst2000 */
  *xp = c2000*xpr2 + s2000*ypr2;
  *yp = -s2000*xpr2 + c2000*ypr2;
  *zp = zpr2;
}


/**************************
***************************/
/* rmat = 3x3 rotation matrix for going from one to another epoch
 * ra1, dec1, ra2, dec2 are in radians
 */

void rotate_radec(double rmat[3][3], double ra1, double dec1,
          double *ra2, double *dec2)
{
   double v1[3], v2[3];

  palDcs2c(ra1,dec1,v1);
  palDmxv(rmat,v1,v2);
  palDcc2s(v2,ra2,dec2);
  *ra2 = palDranrm(*ra2);
}

/**************************
***************************/
/* ra, dec are in radians in this function call
 */

void aber_radec_rad(double eq, double mjd, double ra1, double dec1, double *ra2, double *dec2)
{
  double v1[3], v2[3];

  palDcs2c(ra1,dec1,v1);
  stelaber(eq,mjd,v1,v2);
  palDcc2s(v2,ra2,dec2);
  *ra2 = palDranrm(*ra2);
}


/**************************
***************************/
/* eq = epoch of equinox to be used (e.g., 2000.0 for J2000)
 * mjd = Modified Julian Date (TDB) of correction
 *  will ignore MJD(UTC) vs. MJD(TDB) difference here
 * v1[3] = vector in barycenter frame
 * v2[3] = corresponding vector in Earth-centered frame
 *       = apparent direction from Earth
 */

void stelaber(double eq, double mjd, double v1[3], double v2[3])
{
   double amprms[21], v1n[3], v2un[3], w, ab1, abv[3], p1dv;
   int i;

   palMappa(eq,mjd,amprms);

/* code from mapqk.c (w/ a few names changed): */

/* Unpack scalar and vector parameters */
   ab1 = amprms[11];
   for ( i = 0; i < 3; i++ )
   {
      abv[i] = amprms[i+8];
   }

   palDvn ( v1, v1n, &w );

/* Aberration (normalization omitted) */
   p1dv = palDvdv ( v1n, abv );
   w = 1.0 + p1dv / ( ab1 + 1.0 );
   for ( i = 0; i < 3; i++ ) {
      v2un[i] = ab1 * v1n[i] + w * abv[i];
   }

/* normalize  (not in mapqk.c */
   palDvn ( v2un, v2, &w );

}

/**************************
***************************/
void mat_transpose(double rmat1[3][3], double rmat2[3][3])
{
  int i, j;

  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      rmat2[j][i] = rmat1[i][j];
    }
  }

}

/**************************
***************************/
/* find components of epoch unit vectors in j2000 frame */

void unitvecs_j2000(double rmat[3][3], double xhat[3], double yhat[3], double zhat[3])

{
  int i;

  for(i=0;i<3;++i) {
    xhat[i] = rmat[i][0];
    yhat[i] = rmat[i][1];
    zhat[i] = rmat[i][2];
  }
}



/**************************
***************************/
/* ra, dec, lmst units = radians */

void ha_dec_j2000(double rmat[3][3], double lmst, double lat_rad, double ra2000,
                  double dec2000, double *newha, double *newlat, double *newlmst)
{
  double nwlmst, nwlat;

  rotate_radec(rmat,lmst,lat_rad,&nwlmst,&nwlat);
  *newlmst = nwlmst;
  *newha = palDranrm(nwlmst - ra2000);
  *newlat = nwlat;
}


/**************************
create a lookup table mapping baselines to correlation product index 
***************************/
int makeBaselineLookup(InpConfig *inps, Header *header, array_data *array, int bl_ind_lookup[MAX_ANT][MAX_ANT]) {
    int ant1,ant2,bl_index=0;

    assert(array != NULL);
    assert(header != NULL);
    assert(bl_ind_lookup != NULL);

    /* make a lookup table for which baseline corresponds to a correlation product */
    if(header->corr_type=='A') {
        /* autocorrelations only */
        for(ant1=0; ant1 < array->n_ant; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
          if(debug) fprintf(fpd,"AUTO: bl %d is for ant %d\n",bl_index,ant1);
          bl_ind_lookup[ant1][ant1] = bl_index++;
      }
    }
    else if(header->corr_type=='C') {
        /* this for cross correlations only */
        for (ant1=0; ant1 < array->n_ant-1; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
            for(ant2=ant1+1; ant2 < array->n_ant; ant2++) {
              if (checkAntennaPresent(inps,ant2) == 0) continue;
              if(debug) fprintf(fpd,"CROSS: bl %d is for ants %d-%d\n",bl_index,ant1,ant2);
              bl_ind_lookup[ant1][ant2] = bl_index++;
            }
        }
    }
    else {
        /* this for auto and cross correlations */
        for (ant1=0; ant1 < array->n_ant; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
            for(ant2=ant1; ant2 < array->n_ant; ant2++) {
              if (checkAntennaPresent(inps,ant2) == 0) continue;
              if(debug> 1) fprintf(fpd,"BOTH: bl %d is for ants %d-%d\n",bl_index,ant1,ant2);
              bl_ind_lookup[ant1][ant2] = bl_index++;
            }
        }
    }
    return 0;
}


/*****************
Calculate the phase of each antenna for a given time relative to the array centre.
Applies precession and nutation corrections to calculate the correct u,v,w for
mean (J2000) coords. The phase on the baseline is then just ant1-conj(ant2)
*****************/
int calcAntPhases(double mjd, Header *header, array_data *array,double ant_u[], double ant_v[], double ant_w[]) {
    double ha,ra_app,dec_app,lmst,ra_aber,dec_aber,ha2000,lmst2000,newarrlat;
    double rmatpr[3][3], rmattr[3][3];
    double x,y,z,xprec, yprec, zprec,ant_u_ep, ant_v_ep, ant_w_ep;
    int i;

    lmst = palRanorm(palGmst(mjd) + array->arr_lon_rad);  // local mean sidereal time, given array location

    /* convert mean RA/DEC of phase center to apparent for current observing time. This applies precession,
     nutation, annual abberation. */
    palMap(header->ra_hrs*(M_PI/12.0), header->dec_degs*(M_PI/180.0), 0.0, 0.0, 0.0, 0.0, 2000.0, mjd, &ra_app, &dec_app);
    if (debug) {
        fprintf(fpd,"Mean coords (radian):               RA: %g, DEC: %g\n",header->ra_hrs*(M_PI/12.0),
                      header->dec_degs*(M_PI/180.0));
        fprintf(fpd,"Precessed apparent coords (radian): RA: %g, DEC: %g\n",ra_app,dec_app);
    }
    /* calc apparent HA of phase center, normalise to be between 0 and 2*pi */
    ha = palRanorm(lmst - ra_app);

    /* I think this is correct - it does the calculations in the frame with current epoch
    * and equinox, nutation, and aberrated star positions, i.e., the apparent geocentric
    * frame of epoch. (AML)
    */

    if(debug) fprintf(fpd,"lmst: %g (radian). HA (calculated): %g (radian)\n",lmst,ha);

    /* calc el,az of desired phase centre for debugging */
    if (debug) {
        double az,el;
        palDe2h(ha,dec_app,array->arr_lat_rad,&az,&el);
        fprintf(fpd,"Phase cent ha/dec: %g,%g. az,el: %g,%g\n",ha,dec_app,az,el);
    }

    /* Compute the apparent direction of the phase center in the J2000 coordinate system */
    aber_radec_rad(2000.0,mjd,header->ra_hrs*(M_PI/12.0),header->dec_degs*(M_PI/180.0), &ra_aber,&dec_aber);
    if (debug) {
        fprintf(fpd,"aber_radec_rad apparent coords (radian): RA: %g, DEC: %g\n",ra_aber,dec_aber);
    }

    /* Below, the routines "palPrecl" and "palPreces" do only a precession correction,
     * i.e, they do NOT do corrections for aberration or nutation.
     *
     * We want to go from apparent coordinates at the observation epoch
     * to J2000 coordinates which do not have the nutation or aberration corrections
     * (and since the frame is J2000 no precession correction is needed).
     */

    // slaPrecl(slaEpj(mjd),2000.0,rmatpr);  /* 2000.0 = epoch of J2000 */
    palPrenut(2000.0,mjd,rmattr);
    mat_transpose(rmattr,rmatpr);
    /* rmatpr undoes precession and nutation corrections */
    ha_dec_j2000(rmatpr,lmst,array->arr_lat_rad,ra_aber,dec_aber,&ha2000,&newarrlat,&lmst2000);

    if (debug) {
        fprintf(fpd,"Dec, dec_app, newarrlat (radians): %f %f %f\n", header->dec_degs*(M_PI/180.0),dec_app,newarrlat);
        fprintf(fpd,"lmst, lmst2000 (radians): %f %f\n",lmst,lmst2000);
        fprintf(fpd,"ha, ha2000 (radians): %f %f\n",ha,ha2000);
    }
    /* calc u,v,w at phase center and reference for all antennas relative to center of array */
    for(i=0; i<array->n_ant; i++) {
        // double x,y,z;   /* moved to front of this function (Aug. 12, 2011) */
        x = array->antennas[i].xyz_pos[0];
        y = array->antennas[i].xyz_pos[1];
        z = array->antennas[i].xyz_pos[2];
      /* value of lmst at current epoch - will be changed to effective value in J2000 system
       *
       * To do this, need to precess "ra, dec" (in quotes on purpose) of array center
       * from value at current epoch
       */
        precXYZ(rmatpr,x,y,z,lmst,&xprec,&yprec,&zprec,lmst2000);
        calcUVW(ha,dec_app,x,y,z,&ant_u_ep,&ant_v_ep,&ant_w_ep);
        calcUVW(ha2000,dec_aber,xprec,yprec,zprec,ant_u+i,ant_v+i,ant_w+i);
        if (debug) {
        /* The w value should be the same in either reference frame. */
            fprintf(fpd,"Ant: %d, u,v,w: %g,%g,%g.\n",i,ant_u[i],ant_v[i],ant_w[i]);
            fprintf(fpd,"Ant at epoch: %d, u,v,w: %g,%g,%g.\n",i,ant_u_ep,ant_v_ep,ant_w_ep);
        }
    }
    return 0;
}


/***********************
Apply geometric and/or cable length corrections to the data. Handles various MWA and other
oddities like sign errors and baseline antenna reversals.
************************/
int correctPhases(double mjd, Header *header, InpConfig *inps, array_data *array, int bl_ind_lookup[MAX_ANT][MAX_ANT], float *ac_data, float complex *cc_data,double ant_u[MAX_ANT],double ant_v[MAX_ANT],double ant_w[MAX_ANT]) {
    int chan_ind, baseline_reverse=0, inp1, inp2, ac_chunk_index=0, cc_chunk_index=0,res=0;
    double *k=NULL;
    float vis_weight=1.0;

    k  = malloc(sizeof(double)*header->n_chans);
    assert(k != NULL);

    if(header->integration_time > 0.0) vis_weight = header->integration_time;

    res = calcAntPhases(mjd, header, array,ant_u, ant_v, ant_w);
    if (res) return res;

    /* pre-calculate wavenumber for each channel */
    for(chan_ind=0; chan_ind<header->n_chans; chan_ind++) {
        double freq;

        /* calc wavenumber for this cannel. header freqs are in MHz*/
        freq = (header->cent_freq + (header->invert_freq? -1.0:1.0)*(chan_ind - header->n_chans/2.0)/header->n_chans*header->bandwidth);
        k[chan_ind] = freq/(VLIGHT/1e6);
    }

    for(inp1=0; inp1 < header->n_inputs ; inp1++) {
        for(inp2=inp1; inp2 < header->n_inputs ; inp2++) {
            double w,cable_delay=0.0;
            int pol_ind=0,ant1,ant2,bl_index,pol1,pol2;
            float *visdata=NULL;
            float complex *cvis;

            /* decode the inputs into antennas and pols */
            baseline_reverse=0;
            ant1 = inps->ant_index[inp1];
            ant2 = inps->ant_index[inp2];
            pol1 = inps->pol_index[inp1];
            pol2 = inps->pol_index[inp2];
            /* UVFITS by convention expects the index of ant2 to be greater than ant1, so swap if necessary */
            if (ant1>ant2) {
                int temp;
                temp=ant1;
                ant1=ant2;
                ant2=temp;
                temp=pol1;
                pol1=pol2;
                pol2=temp;
                baseline_reverse=1;
            }
            if (debug) pol_ind = decodePolIndex(pol1, pol2);
            bl_index = bl_ind_lookup[ant1][ant2];

            /* cable delay: the goal is to *correct* for differential cable lengths. The inputs include a delta (offset)
           of cable length relative to some ideal length. (positive = longer than ideal)
           Call the dot product of the baseline (ant2-ant1) and look direction 'phi'.
           Then if ant1 has more delay than ant2, then this is like having phi be positive where
           the visibility is V = Iexp(-j*2*pi*phi)
           Hence we want to add the difference ant2-ant1 (in wavelengths) to phi to correct for the length difference.
            */
            cable_delay = (inps->cable_len_delta[inp2] - inps->cable_len_delta[inp1]);
            if (baseline_reverse) cable_delay = -cable_delay;

            /* only process the appropriate correlations */
            if (header->corr_type=='A' && inp1!=inp2) continue;
            if (header->corr_type=='C' && inp1==inp2) continue;

            /* keep track of which chunk of channels we're up to in autos and crosses.
            Use visdata as a sort of void pointer to point to the start of the block
            of channels. */
            if (inp1 != inp2) {
                /* process a block of cross-correlations */
                visdata = (float *)(cc_data + header->n_chans*cc_chunk_index);
                cc_chunk_index += 1;
            }
            else {
                /* process a block of auto-correlations */
                visdata = (ac_data + header->n_chans*ac_chunk_index);
                ac_chunk_index += 1;
            }
 
            /* calc u,v,w for this baseline in meters */
            w=0.0;
            /* if not correcting for geometry, don't apply w */
            if(header->geom_correct) {
                //u = ant_u[ant1] - ant_u[ant2];
                //v = ant_v[ant1] - ant_v[ant2];
                w = ant_w[ant1] - ant_w[ant2];
            }
            // add the cable delay term and 2pi here so that we don't have to add it many times in the loop below.
            w = (w+cable_delay)*(-2.0*M_PI);

            if (debug>1) {
                fprintf(fpd,"doing inps %d,%d. ants: %d,%d pols: %d,%d, polind: %d, bl_ind: %d, w (m): %g, delay (m): %g, blrev: %d\n",
                    inp1,inp2,ant1,ant2,pol1,pol2,pol_ind,bl_index,w,cable_delay,baseline_reverse);
            }

            /* populate the visibility arrays */
            cvis = (float complex *) visdata;    /* cast this so we can use pointer arithmetic */
            for(chan_ind=0; chan_ind < header->n_chans; chan_ind++) {
                float complex vis,phase;

                if (inp1 != inp2) {
                    phase = cexp(I*w*k[chan_ind]);
                    //vis = visdata[chan_ind*2] + I*(header->conjugate ? -visdata[chan_ind*2+1]: visdata[chan_ind*2+1]);
                    vis = header->conjugate ? conjf(cvis[chan_ind]) : cvis[chan_ind];
/*
                    if(debug && chan_ind==header->n_chans/2) {
                        fprintf(fpd,"Chan %d, w: %g (wavelen), vis: %g,%g. ph: %g,%g. rot vis: %g,%g\n",
                                    chan_ind,w*k[chan_ind],creal(vis),cimag(vis),creal(phase),cimag(phase),
                                    creal(vis*phase),cimag(vis*phase));
                    }
*/
                    // apply input-based flags if necessary

                    if ( (inps->inpFlag[inp1] || inps->inpFlag[inp2]) && vis_weight > 0) {
                        vis=0.0;
                    }

                    if (baseline_reverse) vis = conjf(vis);
                    cvis[chan_ind] = vis*phase; /* update the input data with the phase correction */
                }
            }
        }
    }
    if (k != NULL) free(k);
    return 0;
}

