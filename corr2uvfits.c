/* Program to convert the binary output format of the 32T on-site correlators
 * UV FITS format. Written by Randall Wayth. Feb, 2008.
 *
 * August 12, 2011 - precess u, v, w, data to J2000 frame (Alan Levine)
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <fitsio.h>
#include <assert.h>
#include <star/pal.h>
#include "uvfits.h"
#include "convutils.h"


/* private function prototypes */
void printusage(const char *progname);
int readScan(FILE *fp_ac, FILE *fp_cc,int scan,Header *header, InpConfig *inps,uvdata *data);
void parse_cmdline(const int argc, char * const argv[]);
void initData(uvdata *data);
int applyHeader(Header *header, uvdata *data);
void checkInputs(Header *header,uvdata *data,InpConfig *inputs);

/* allow convutils to use same file handle */
FILE *fpd=NULL;

/* private global vars */
static int bl_ind_lookup[MAX_ANT][MAX_ANT];
static int debug=0,do_flag=0,lock_pointing=0;
static int pol_index[4];
char *stationfilename="antenna_locations.txt",*outfilename=NULL;
char *configfilename="instr_config.txt";
char *header_filename="header.txt";
char *crosscor_filename=NULL;
char *autocorr_filename=NULL;
char *bothcorr_filename=NULL;
char *flagfilename=NULL;
static double arr_lat_rad=MWA_LAT*(M_PI/180.0),arr_lon_rad=MWA_LON*(M_PI/180.0),height=MWA_HGT;

/************************
************************/
int main(const int argc, char * const argv[]) {
  FILE *fpin_ac=NULL,*fpin_cc=NULL;
  int i,scan=0,res=0;
  Header header;
  InpConfig inputs;
  uvdata *data;
  array_data *arraydat;
  source_table *source;
  ant_table *antennas;
  uvWriteContext *iter=NULL;

  fpd=stderr;

  if(argc < 2) printusage(argv[0]);
  parse_cmdline(argc,argv);
  setConvDebugLevel(debug);

  /* initialise some values for the UV data array and antennas*/
  data = calloc(1,sizeof(uvdata));
  source = calloc(1,sizeof(source_table));
  arraydat = calloc(1,sizeof(array_data));
  antennas = calloc(MAX_ANT,sizeof(ant_table));
  memset(&inputs,'\0', sizeof(InpConfig));
  assert(antennas!=NULL && source!=NULL && arraydat!=NULL && data !=NULL);
  data->array    = arraydat;
  data->array->antennas = antennas;
  data->source   = source;
  initData(data);

  /* get the mapping of inputs to anntena numbers and polarisations */
  if ((res = readInputConfig(configfilename, &inputs)) != 0) {
      fprintf(stderr,"readInputConfig failed with code %d. exiting\n",res);
  }

  /* get the number of antennas and their locations relative to the centre of the array */
  if ((res = readArray(stationfilename, arr_lat_rad, arr_lon_rad, arraydat)) != 0) {
      fprintf(stderr,"readArray failed with code %d. exiting\n",res);
  }

  /* read the header/metadata  */
  res = readHeader(header_filename,&header);
  if (res != 0) {
    fprintf(stderr,"Error reading main header. exiting.\n");
    exit(1);
  }

  checkInputs(&header,data,&inputs);

  /* open input files */
  if (bothcorr_filename != NULL) {
    /* this code to support the newer raw data file format that includes both autos and cross
     * correlations in one file. Since we read sequentially and read in the order that is implicit in
     * the new format, we can re-use the same file pointer */
    fpin_cc=fopen(bothcorr_filename,"r");
    if (fpin_cc==NULL) {
        fprintf(stderr,"cannot open combined correlation input file <%s>\n",bothcorr_filename);
        exit(1);
    }
    fpin_ac = fpin_cc; // re-use file pointer, need to ensure subsequent reads are in the implied order
  }
  else {
    if (header.corr_type!='A' && (fpin_cc=fopen(crosscor_filename,"r"))==NULL) {
        fprintf(stderr,"cannot open cross correlation input file <%s>\n",crosscor_filename);
        exit(1);
    }
    if (header.corr_type!='C' && (fpin_ac=fopen(autocorr_filename,"r"))==NULL) {
        fprintf(stderr,"cannot open auto correlation input file <%s>\n",autocorr_filename);
        exit(1);
    }
  }

  /* assign vals to output data structure from inputs */
  res = applyHeader(&header, data);
  if (bothcorr_filename != NULL) {
    fprintf(fpd,"Setting correlation type to S for combined (Single) file input\n");
    header.corr_type='S';
  }

  /* populate antenna info */
  if (debug) fprintf(fpd,"there are %d antennas\n",arraydat->n_ant);
  for (i=0; i<arraydat->n_ant; i++){
    //sprintf(antennas[i].name,"ANT%03d",i+1);
    // FIXME: update this for correct pol type
    sprintf(antennas[i].pol_typeA,"X");
    sprintf(antennas[i].pol_typeB,"Y");
    antennas[i].pol_angleA = 0.0;
    antennas[i].pol_angleB = 90.0;
    antennas[i].pol_calA = 0.0;
    antennas[i].pol_calB = 0.0;
    antennas[i].mount_type = 0;
  }

  /* assign XYZ positions of the array for the site. */
  Geodetic2XYZ(arr_lat_rad,arr_lon_rad,height,&(arraydat->xyz_pos[0]),&(arraydat->xyz_pos[1]),&(arraydat->xyz_pos[2]));
  if (debug) fprintf(fpd,"converted array location to XYZ\n");

  /* create correlator->baseline mapping lookup table */
  assert(makeBaselineLookup(&inputs, &header, data->array, bl_ind_lookup)==0);

  /* open the iterator for writing. Need to have set the date in  */
  res = writeUVFITSiterator(outfilename, data, &iter);
  if (res!=0) {
    fprintf(stderr,"writeUVFITSiterator returned %d.\n",res);
    return res;
  }
  if (debug) fprintf(fpd,"created UVFITS writer iterator\n");

  /* read each scan, populating the data structure. */
  scan=0;
  while ((res = readScan(fpin_ac,fpin_cc,scan, &header, &inputs, data))==0) {
    int status;

    if (flagfilename != NULL) {
      if (debug) fprintf(fpd,"Applying flags...\n");
      res = applyFlagsFile(flagfilename,data);
      if(res!=0) {
        fprintf(stderr,"Problems in applyFlagsFile. exiting\n");
        exit(1);
      }
    }

    /* write this chunk of data to the output file */
    status = writeUVinstant(iter, data, data->date[0]-iter->jd_day_trunc,0);
    if (status) {
        fprintf(stderr,"ERROR: code %d when writing time instant %d\n",status,i);
        exit(EXIT_FAILURE);
    }
    scan++;
    if (scan ==header.n_scans) break;   // don't read more than is specified
  }
  if(res < 0) {
      fprintf(stderr,"Problems in readScan(). exiting\n");
      exit(1);
  }

  if (debug) fprintf(fpd,"Read %d time steps\n",scan);
  if (abs(header.n_scans-scan) > 4) {
    fprintf(stderr,"WARNING: expected to read %d scans, but actually read %d. Are you sure you have the correct number of freqs, inputs and timesteps? Carrying on and hoping for the best...\n",header.n_scans,scan);
  }
  else if (res > 0) {
      fprintf(stderr,"WARNING: Wanted to read %d time steps, but actually read %d. Carrying on...\n",header.n_scans, scan);
  }

  if (do_flag) {
    if (debug) fprintf(fpd,"Auto flagging...\n");
    res = autoFlag(data,5.0,AUTOFLAG_N_NEIGHBOURS,do_flag);
    if(res!=0) {
      fprintf(stderr,"Problems in autoflag. exiting\n");
      exit(1);
    }
  }

  /* finish up and write the antenna table */
  writeUVFITSfinalise(iter, data);

  if(debug) fprintf(fpd,"finished writing UVFITS file\n");
  if(fpin_ac !=NULL) fclose(fpin_ac);
  if(fpin_cc !=NULL) fclose(fpin_cc);

  freeUVFITSdata(data);

  return 0;
}


/***************************
 ***************************/
int readScan(FILE *fp_ac, FILE *fp_cc,int scan_count, Header *header, InpConfig *inps, uvdata *uvdata) {

  static int init=0;
  static float vis_weight=1.0;
  static double date_zero=0.0;  // time of zeroth scan

  array_data *array;    /*convenience pointer */
  double mjd;
  double ant_u[MAX_ANT],ant_v[MAX_ANT],ant_w[MAX_ANT]; //u,v,w for each antenna, in meters
  int res=0,chan_ind,n_read,inp1,inp2;
  int scan=0,ac_chunk_index=0,cc_chunk_index=0;
  size_t size_ac, size_cc, size_comb;
  float *ac_data=NULL,*visdata=NULL;
  float complex *cc_data=NULL,*cc_temp=NULL;

  array = uvdata->array;

  /* allocate space to read 1 time step of binary correlation data.
    There are n_inputs*nchan floats of autocorrelations per time step 
    and n_inp*(n_inp-1)/2*nchan float complex cross correlations*/
  size_ac = uvdata->n_freq*header->n_inputs*sizeof(float);
  size_cc = uvdata->n_freq*header->n_inputs*(header->n_inputs-1)/2*sizeof(float complex);
  /* if using a single combined auto/cross file, then autos have redundant imaginary components */
  size_comb = uvdata->n_freq*header->n_inputs*(header->n_inputs+1)/2*sizeof(float complex);
  
  if (header->corr_type=='S') {
    // handle special case of combined single input auto/cross file, which has redundant imag terms for autos
    cc_temp = calloc(1,size_comb);
  }
  ac_data = calloc(1,size_ac);
  cc_data = calloc(1,size_cc);
  assert(ac_data != NULL);
  assert(cc_data != NULL);

  if (!init) {
    /* count the total number of antennas actually present in the data */
    if (debug) fprintf(fpd,"Init %s.\n",__func__);

    /* increase size of arrays for the new scan */
    // date and n_baselines should already be allocated
    assert(uvdata->date != NULL);
    assert(uvdata->n_baselines[0] > 0);
    uvdata->n_vis=1;
    uvdata->visdata=calloc(1,sizeof(double *));
    uvdata->weightdata=calloc(1,sizeof(double *));
    uvdata->u=calloc(1,sizeof(double *));
    uvdata->v=calloc(1,sizeof(double *));
    uvdata->w=calloc(1,sizeof(double *));
    uvdata->baseline = calloc(1,sizeof(float *));

    /* make space for the actual visibilities and weights */
    if (debug) fprintf(fpd,"%s: callocing array of %d floats\n",__func__,uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol*2);
    uvdata->visdata[scan]    = calloc(uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol*2,sizeof(float));
    uvdata->weightdata[scan] = calloc(uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol  ,sizeof(float));
    uvdata->u[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->v[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->w[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->baseline[scan] = calloc(uvdata->n_baselines[scan],sizeof(float));
    if(uvdata->visdata[scan]==NULL || uvdata->weightdata[scan]==NULL || uvdata->visdata[scan]==NULL
       || uvdata->visdata[scan]==NULL || uvdata->visdata[scan]==NULL || uvdata->baseline[scan]==NULL) {
      fprintf(stderr,"%s: no malloc for BIG arrays\n",__func__);
      exit(1);
    }
    /* set a weight for the visibilities based on integration time */
    if(header->integration_time > 0.0) vis_weight = header->integration_time;

    date_zero = uvdata->date[0];    // this is already initialised in applyHeader

    init=1;
  }

  /* read all the data for this timestep */
  /* need to read autos first for the case of a combined single raw data input file,
   * in which case the two file pointers will point to the same file */
  if (header->corr_type=='S') {
    n_read = fread(cc_temp,size_comb,1,fp_ac);
  }
  else {
    n_read = fread(ac_data,size_ac,1,fp_ac);
  }
  if (n_read != 1) {
    fprintf(stderr,"EOF: reading auto correlations. Wanted to read %d bytes.\n",(int) size_ac);
    return 1;
  }
  if (header->corr_type!='S') {
    n_read = fread(cc_data,size_cc,1,fp_cc);
    if (n_read != 1) {
        fprintf(stderr,"EOF: reading cross correlations. Wanted to read %d bytes.\n",(int) size_cc);
        return 1;
    }
  }
  else {
    float complex *tmpvis,*tmpcc;
    int prod_ind=0;
    /* copy out the combined data into the original separate data structures for exisiting code below */
    for(inp1=0; inp1 < header->n_inputs ; inp1++) {
        for(inp2=inp1; inp2 < header->n_inputs ; inp2++) {
            tmpvis = (cc_temp + header->n_chans*prod_ind);
            if (inp1 != inp2) {
                /* process a block of cross-correlations */
                tmpcc = cc_data + header->n_chans*cc_chunk_index;
                memcpy(tmpcc,tmpvis,header->n_chans*sizeof(float complex));
                cc_chunk_index += 1;
            }
            else {
                /* process a block of auto-correlations */
                visdata = (ac_data + header->n_chans*ac_chunk_index);
                for (int i=0; i<header->n_chans; i++) {
                    visdata[i] = crealf(tmpvis[i]);
                }
                ac_chunk_index += 1;
            }
            prod_ind++;
        }
    }
  }

  /* set time of scan. Note that 1/2 scan time offset already accounted for in date[0]. */
  if (scan_count > 0) uvdata->date[scan] = date_zero + scan_count*header->integration_time/86400.0;
  mjd = uvdata->date[scan] - 2400000.5;  // get Modified Julian date of scan.

  /* set default ha/dec from header, if HA was specified. Otherwise, it will be calculated below */
  if (lock_pointing!=0) {   // special case for RTS output which wants a phase centre fixed at an az/el, not ra/dec
    mjd = date_zero - 2400000.5;  // get Modified Julian date of scan.
  }

  /* apply geometric and/or cable length corrections to the visibilities */
  res = correctPhases(mjd, header, inps, array, bl_ind_lookup, ac_data, cc_data, ant_u, ant_v, ant_w);
  if (res) return res;

  /* copy out the phase rotated data into the uvfits data structure */
  ac_chunk_index=cc_chunk_index=0;
  for(inp1=0; inp1 < header->n_inputs ; inp1++) {
    for(inp2=inp1; inp2 < header->n_inputs ; inp2++) {
        float complex *cvis;
        int ant1,ant2,pol1,pol2,pol_ind,bl_index;
        double u,v,w;

        /* decode the inputs into antennas and pols */
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
        }
        pol_ind = decodePolIndex(pol1, pol2);
        bl_index = bl_ind_lookup[ant1][ant2];


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

        if (header->corr_type=='C' && ant1==ant2) continue;

        /* calc u,v,w for this baseline in meters */
        u=v=w=0.0;
        if(ant1 != ant2) {
            u = ant_u[ant1] - ant_u[ant2];
            v = ant_v[ant1] - ant_v[ant2];
            w = ant_w[ant1] - ant_w[ant2];
        }

        /* populate the baseline info. Antenna numbers start at 1 in UVFITS.  */
        uvdata->baseline[scan][bl_index] = 0; // default: no baseline. useful to catch bugs.
        EncodeBaseline(ant1+1, ant2+1, uvdata->baseline[scan]+bl_index);
        if (debug) fprintf(fpd,"a1: %d, a2: %d. bl_index: %d, polind: %d, uvbl: %f\n",ant1,ant2,bl_index,pol_ind,*(uvdata->baseline[scan]+bl_index));
        /* arcane units of UVFITS require u,v,w in light-seconds */
        uvdata->u[scan][bl_index] = u/VLIGHT;
        uvdata->v[scan][bl_index] = v/VLIGHT;
        uvdata->w[scan][bl_index] = w/VLIGHT;

        cvis = (float complex *)visdata;    /* cast this so we can use pointer arithmetic */
        for(chan_ind=0; chan_ind<uvdata->n_freq; chan_ind++) {
            int visindex;

            visindex = bl_index*uvdata->n_pol*uvdata->n_freq + chan_ind*uvdata->n_pol + pol_ind;

            if (inp1 != inp2) {

                /* cross correlation, use imaginary and real */
                uvdata->visdata[scan][visindex*2   ] = crealf(cvis[chan_ind]);
                uvdata->visdata[scan][visindex*2 +1] = cimagf(cvis[chan_ind]);
            }
            else {
                /* auto correlation, set imag to zero */
                uvdata->visdata[scan][visindex*2   ] = visdata[chan_ind];
                uvdata->visdata[scan][visindex*2 +1] = 0.0;
            }
            uvdata->weightdata[scan][visindex] = vis_weight;
            // apply input-based flags if necessary
            if ( (inps->inpFlag[inp1] || inps->inpFlag[inp2]) && vis_weight > 0) {
                uvdata->weightdata[scan][visindex] = -vis_weight;
            }
        }
    }
  }


  /* sanity check */
  if (debug) {
    fprintf(fpd,"At end of %s. ac_chunk_index: %d, cc_chunk_index: %d\n", __func__,ac_chunk_index, cc_chunk_index);
  }

  if (ac_data != NULL) free(ac_data);
  if (cc_data != NULL) free(cc_data);
  if (cc_temp != NULL) free(cc_temp);
  return 0;
}


/****************************
*****************************/
void parse_cmdline(const int argc,char * const argv[]) {
    const char optstring[] = "vldS:a:c:o:I:H:A:F:f:b:";
    int result=0;
    char arrayloc[80],*lon,*lat;

    arrayloc[0]='\0';

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'S': stationfilename = optarg;
            break;
          case 'o': outfilename = optarg;
            break;
          case 'a': autocorr_filename = optarg;
            break;
          case 'c': crosscor_filename = optarg;
            break;
          case 'b': bothcorr_filename = optarg;
            break;
          case 'd': debug += 1;
            fprintf(fpd,"Debugging level %d\n",debug);
            break;
          case 'I': configfilename = optarg;
            break;
          case 'H': header_filename = optarg;
            break;
          case 'f': do_flag=atoi(optarg);
            break;
          case 'F': flagfilename=optarg;
            break;
          case 'l': lock_pointing=1;
            fprintf(fpd,"Locking phase center to initial HA/DEC\n");
            break;
          case 'A':
            strncpy(arrayloc,optarg,80);
            break;
          case 'v':
            fprintf(stdout,"corr2uvfits revision $Rev: 4135 $\n");
            exit(1);
            break;
          default:
              fprintf(stderr,"unknown option: %c\n",result);
              printusage(argv[0]);
        }
    }

    /* convert array lon/lat */
    if(arrayloc[0]!='\0') {
        lon = arrayloc;
        lat = strpbrk(arrayloc,",");
        if (lat ==NULL) {
            fprintf(stderr,"Cannot find comma separator in lon/lat. Typo?\n");
            printusage(argv[0]);
        }
        /* terminate string for lon, then offset for lat */
        *lat = '\0';
        lat++;
        arr_lat_rad = atof(lat);
        arr_lon_rad = atof(lon);
        fprintf(fpd,"User specified array lon,lat: %g, %g (degs)\n",arr_lon_rad,arr_lat_rad);
        arr_lat_rad *= (M_PI/180.0); // convert to radian
        arr_lon_rad *= (M_PI/180.0);
    }

    /* do some sanity checks */
    if(outfilename==NULL) {
        fprintf(stderr,"ERROR: no output file name specified\n");
        exit(1);
    }
    /* auto flagging requires autocorrelations */
    if (autocorr_filename==NULL && do_flag) {
        fprintf(stderr,"ERROR: auto flagging requires the autocorrelations to be used\n");
        exit(1);
    }
}


/***************************
 ***************************/
void printusage(const char *progname) {
  fprintf(stderr,"Usage: %s [options]\n\n",progname);
  fprintf(stderr,"options are:\n");
  fprintf(stderr,"-a filename\tThe name of autocorrelation data file. no default.\n");
  fprintf(stderr,"-c filename\tThe name of cross-correlation data file. no default.\n");
  fprintf(stderr,"-b filename\tThe name of combined auto/cross data file. no default.\n");
  fprintf(stderr,"-o filename\tThe name of the output file. No default.\n");
  fprintf(stderr,"-S filename\tThe name of the file containing antenna name and local x,y,z. Default: %s\n",stationfilename);
  fprintf(stderr,"-I filename\tThe name of the file containing instrument config. Default: %s\n",configfilename);
  fprintf(stderr,"-H filename\tThe name of the file containing observing metadata. Default: %s\n",header_filename);
  fprintf(stderr,"-A lon,lat \tSpecify East Lon and Lat of array center (degrees). Comma separated, no spaces. Default: MWA\n");
  fprintf(stderr,"-l         \tLock the phase center to the initial HA/DEC\n");
  fprintf(stderr,"-f mode    \tturn on automatic flagging. Requires autocorrelations\n");
  fprintf(stderr,"\t\t0:\tno flagging\n");
  fprintf(stderr,"\t\t1:\tgeneric flagging based on autocorrelation median\n");
  fprintf(stderr,"-F filename\tOptionally apply global flags as specified in filename.\n");
  fprintf(stderr,"-d         \tturn debugging on.\n");
  fprintf(stderr,"-v         \treturn revision number and exit.\n");
  exit(1);
}


/*******************************
*********************************/
void initData(uvdata *data) {
  data->date = calloc(1,sizeof(double));
  data->n_pol=0;
  data->n_baselines = calloc(1,sizeof(int));
  data->n_freq=0;
  data->n_vis=0;
  data->cent_freq=0.0;
  data->freq_delta = 0.0;
  data->u=NULL;
  data->v=NULL;
  data->w=NULL;
  data->baseline=NULL;
  data->weightdata=NULL;
  data->visdata=NULL;
  data->pol_type=1;    /* default is Stokes pol products */
  data->array->n_ant=0;
  data->array->arr_lat_rad = arr_lat_rad;
  data->array->arr_lon_rad = arr_lon_rad;
  strcpy(data->array->name,"UNKNOWN");
}


/***********************
***********************/
int applyHeader(Header *header, uvdata *data) {

  double jdtime_base=0,mjd,lmst;
  int n_polprod=0,i,res;

  data->n_freq = header->n_chans;
  data->cent_freq = header->cent_freq*1e6;
  data->freq_delta = header->bandwidth/(header->n_chans)*1e6*(header->invert_freq ? -1.0: 1.0);
  strncpy(data->array->name,header->telescope,SIZ_TELNAME);
  strncpy(data->array->instrument,header->instrument,SIZ_TELNAME);

  /* discover how many pol products there are */
  createPolIndex(header->pol_products, pol_index);
  for(i=0; i<4; i++) if (header->pol_products[2*i] !='\0') n_polprod++;
  if(n_polprod<1 || n_polprod > 4) {
    fprintf(stderr,"bad number of stokes: %d\n",n_polprod);
    exit(1);
  }
  data->n_pol = n_polprod;

  /* set the polarisation product type. linear, circular or stokes */
  if (toupper(header->pol_products[0]) == 'X' || toupper(header->pol_products[0]) == 'Y') data->pol_type=-5;
  if (toupper(header->pol_products[0]) == 'R' || toupper(header->pol_products[0]) == 'L') data->pol_type=-1;
  if (debug) fprintf(fpd,"Found %d pol products. pol_type is: %d\n",n_polprod,data->pol_type);

  /* calculate the JD of the beginning of the data */
  palCaldj(header->year, header->month, header->day, &jdtime_base, &res); // get MJD for calendar day of obs
  jdtime_base += 2400000.5;  // convert MJD to JD
  jdtime_base += header->ref_hour/24.0+header->ref_minute/1440.0+header->ref_second/86400.0; // add intra-day offset
  data->date[0] = jdtime_base+0.5*(header->integration_time/86400.0);
  if (debug) fprintf(fpd,"JD time is %.2lf\n",jdtime_base);

  memset(data->source->name,0,SIZE_SOURCE_NAME+1);
  strncpy(data->source->name,header->field_name,SIZE_SOURCE_NAME);

  mjd = data->date[0] - 2400000.5;  // get Modified Julian date of scan.
  lmst = palRanorm(palGmst(mjd) + arr_lon_rad);  // local mean sidereal time, given array location

  /* if no RA was specified in the header,then calculate the RA based on lmst and array location
     and update the ra */
  if (header->ra_hrs < -98.0 ) {
    // set the RA to be for the middle of the scan
//    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start + header->n_scans*header->integration_time*1.00274/(3600.0*2);   // include 1/2 scan offset
    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start;  // match existing code. RA defined at start of scan
    if (debug) fprintf(fpd,"Calculated RA_hrs: %g of field centre based on HA_hrs: %g and lmst_hrs: %g\n",
                        header->ra_hrs,header->ha_hrs_start,lmst*(12.0/M_PI));
  }

  /* extract RA, DEC from header. Beware negative dec and negative zero bugs. */
  data->source->ra  = header->ra_hrs;
  data->source->dec = header->dec_degs;

  /* calcualte the number of baselines, required to be constant for all data */
  data->n_baselines[0] = (data->array->n_ant)*(data->array->n_ant+1)/2; //default: both auto and cross
  if (header->corr_type=='A') data->n_baselines[0] = data->array->n_ant;
  if (header->corr_type=='C') data->n_baselines[0] = data->array->n_ant*(data->array->n_ant-1)/2;
  if (debug) fprintf(fpd,"Corr type %c, so there are %d baselines\n",header->corr_type,data->n_baselines[0]);

  return 0;
}


/**************************
***************************/
void checkInputs(Header *header,uvdata *data,InpConfig *inputs) {
    int total_ants;

    if(inputs->n_inputs != header->n_inputs) {
        fprintf(stderr,"ERROR: mismatch between the number of inputs in %s (%d) and header (%d)\n",
                configfilename,inputs->n_inputs,header->n_inputs);
        exit(1);
    }

    total_ants = countPresentAntennas(inputs);

    if (total_ants > data->array->n_ant) {
        fprintf(stderr,"ERROR: mismatch between the number of antennas in %s (%d) and %s (%d)\n",
                stationfilename,data->array->n_ant,configfilename,total_ants);
    }
    
    if (do_flag && header->corr_type == 'C') {
        fprintf(stderr,"ERROR: CORRTYPE must be auto or both for autoflagging\n");
        exit(1);
    }
    if (header->ha_hrs_start == -99.0 && lock_pointing) {
        fprintf(stderr,"ERROR: HA must be specified in header if -l flag is used.\n");
        exit(1);
    }

}

