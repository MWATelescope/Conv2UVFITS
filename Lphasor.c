/* Program to apply geometric and cable length corrections to
 * "L-file" data. Based heavily on corr2uvfits which has had
 * the core code pulled out into convutils.
 *  Randall Wayth. May 2014.
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <star/pal.h>
#include "uvfits.h"
#include "convutils.h"

typedef struct {
    char *outccfilename;
    char *stationfilename;
    char *configfilename;
    char *header_filename;
    char *crosscor_filename;
    char *flagfilename;
    double arr_lat_rad;
    double arr_lon_rad;
    double height;
    int do_flag;
} cl_options;

/* private function prototypes */
void printusage(const char *progname);
int doScan(FILE *fp_cc,FILE *fpout_cc,int scan,Header *header, InpConfig *inps, uvdata *data);
void parse_cmdline(const int argc, char * const argv[], const char *optstring);
int applyHeader(Header *header, uvdata *data);
void checkInputs(Header *header,uvdata *data,InpConfig *inputs);
void initData(uvdata *data);

/* allow convutils to use same file handle */
FILE *fpd=NULL;

/* private global vars */
static int bl_ind_lookup[MAX_ANT][MAX_ANT];
static int debug=0,lock_pointing=0;
static int pol_index[4];
static cl_options options;

/************************
************************/
int main(const int argc, char * const argv[]) {
  const char optstring[] = "vldS:c:o:I:H:A:F:";
  FILE *fpin_cc=stdin,*fpout_cc=stdout;
  int scan=0,res=0;
  Header header;
  InpConfig inputs;
  uvdata *data;
  array_data *arraydat;
  ant_table *antennas;

  fpd=stderr;
  memset(&options,'\0',sizeof(cl_options));
  options.stationfilename="antenna_locations.txt";
  options.configfilename="instr_config.txt";
  options.header_filename="header.txt";
  options.arr_lat_rad=MWA_LAT*(M_PI/180.0);
  options.arr_lon_rad=MWA_LON*(M_PI/180.0);
  options.height=MWA_HGT;
  options.do_flag=0;

  if(argc < 2) printusage(argv[0]);
  parse_cmdline(argc,argv,optstring);
  setConvDebugLevel(debug);

  /* initialise some values for the UV data array and antennas*/
  data = calloc(1,sizeof(uvdata));
  arraydat = calloc(1,sizeof(array_data));
  antennas = calloc(MAX_ANT,sizeof(ant_table));
  assert(antennas!=NULL && arraydat!=NULL && data != NULL);
  memset(&inputs,'\0', sizeof(InpConfig));
  arraydat->antennas = antennas;
  arraydat->arr_lat_rad = options.arr_lat_rad;
  arraydat->arr_lon_rad = options.arr_lon_rad;
  data->array = arraydat;
  initData(data);

  /* get the mapping of inputs to anntena numbers and polarisations */
  if ((res = readInputConfig(options.configfilename, &inputs)) != 0) {
      fprintf(stderr,"readInputConfig failed with code %d. exiting\n",res);
  }

  /* get the number of antennas and their locations relative to the centre of the array */
  if ((res = readArray(options.stationfilename, options.arr_lat_rad, options.arr_lon_rad, arraydat)) != 0) {
      fprintf(stderr,"readArray failed with code %d. exiting\n",res);
  }

  /* read the header/metadata  */
  res = readHeader(options.header_filename,&header);
  if (res != 0) {
    fprintf(stderr,"Error reading main header. exiting.\n");
    exit(1);
  }

  checkInputs(&header,data,&inputs);

  /* open input */
  if (options.crosscor_filename != NULL && (fpin_cc=fopen(options.crosscor_filename,"r"))==NULL) {
    fprintf(stderr,"cannot open cross correlation input file <%s>\n",options.crosscor_filename);
    exit(1);
  }

  /* open output */
  if (options.outccfilename != NULL && (fpout_cc=fopen(options.outccfilename,"wb"))==NULL) {
    fprintf(stderr,"cannot open cross correlation output file <%s>\n",options.outccfilename);
    exit(1);
  }

  /* assign vals to output data structure from inputs */
  res = applyHeader(&header, data);

  /* populate antenna info */
  if (debug) fprintf(fpd,"there are %d antennas\n",arraydat->n_ant);

  /* assign XYZ positions of the array for the site. */
  Geodetic2XYZ(arraydat->arr_lat_rad,arraydat->arr_lon_rad,options.height,
                &(arraydat->xyz_pos[0]),&(arraydat->xyz_pos[1]),&(arraydat->xyz_pos[2]));
  if (debug) fprintf(fpd,"converted array location to XYZ\n");

  /* create correlator->baseline mapping lookup table */
  assert(makeBaselineLookup(&inputs, &header, data->array, bl_ind_lookup)==0);

  /* read each scan, populating the data structure. */
  scan=0;
  while ((res = doScan(fpin_cc,fpout_cc,scan, &header, &inputs,data))==0) {

    if (options.flagfilename != NULL) {
      if (debug) fprintf(fpd,"Applying flags...\n");
      res = applyFlagsFile(options.flagfilename,data);
      if(res!=0) {
        fprintf(stderr,"Problems in applyFlagsFile. exiting\n");
        exit(1);
      }
    }
    scan++;
  }
  if(res < 0) {
      fprintf(stderr,"Problems in readScan(). exiting\n");
      exit(1);
  }

  if (debug) fprintf(fpd,"Read %d time steps\n",scan);

  /* finish up  */
  if(fpin_cc !=NULL && fpin_cc != stdin) fclose(fpin_cc);
  if(fpout_cc !=NULL && fpout_cc != stdout) fclose(fpout_cc);

  freeUVFITSdata(data);

  return 0;
}


/***************************
 ***************************/
int doScan(FILE *fp_cc, FILE *fpout_cc, int scan_count, Header *header, InpConfig *inps, uvdata *uvdata) {

  static int init=0;
  //static float vis_weight=1.0;
  static double date_zero=0.0;  // time of zeroth scan

  double mjd;
  double ant_u[MAX_ANT],ant_v[MAX_ANT],ant_w[MAX_ANT]; //u,v,w for each antenna, in meters
  int res=0,n_read,scan=0;
  size_t size_cc;
  float complex *cc_data=NULL;
  array_data *array;

  array = uvdata->array;
  /* allocate space to read 1 time step of binary correlation data.
    There are n_inputs*nchan floats of autocorrelations per time step 
    and n_inp*(n_inp-1)/2*nchan float complex cross correlations*/
  size_cc = header->n_chans*header->n_inputs*(header->n_inputs-1)/2*sizeof(float complex);
  cc_data = malloc(size_cc);
  assert(cc_data != NULL);

  if (!init) {
    /* count the total number of antennas actually present in the data */
    if (debug) fprintf(fpd,"Init %s.\n",__func__);
    /* set a weight for the visibilities based on integration time */
//    if(header->integration_time > 0.0) vis_weight = header->integration_time;

    date_zero = uvdata->date[0];    // this is already initialised in applyHeader

    init=1;
  }

  /* read all the data for this timestep */
  n_read = fread(cc_data,size_cc,1,fp_cc);
  if (n_read != 1) {
    //fprintf(stderr,"EOF: reading cross correlations. Wanted to read %d bytes.\n",(int) size_cc);
    res=1;
    goto EXIT;
  }

  /* set time of scan. Note that 1/2 scan time offset already accounted for in date[0]. */
  if (scan_count > 0) uvdata->date[scan] = date_zero + scan_count*header->integration_time/86400.0;
  mjd = uvdata->date[scan] - 2400000.5;  // get Modified Julian date of scan.

  /* set default ha/dec from header, if HA was specified. Otherwise, it will be calculated below */
  if (lock_pointing!=0) {   // special case for RTS output which wants a phase centre fixed at an az/el, not ra/dec
    mjd = date_zero - 2400000.5;  // get Modified Julian date of scan.
  }

  /* apply geometric and/or cable length corrections to the visibilities */
  res = correctPhases(mjd, header, inps, array, bl_ind_lookup, NULL, cc_data, ant_u, ant_v, ant_w);
  if (res) goto EXIT;

  /* write the data back out again */
  n_read = fwrite(cc_data,size_cc,1,fpout_cc);
  if (n_read != 1) {
    fprintf(stderr,"%s: Failed to write cross correlations of size %d\n",__func__,(int)size_cc);
    res = 1;
  }

  if (debug) fprintf(fpd,"%s: processed scan %d\n",__func__,scan_count);

EXIT:
  if (cc_data != NULL) free(cc_data);
  return res;
}


/****************************
*****************************/
void parse_cmdline(const int argc,char * const argv[], const char *optstring) {
    int result=0;
    char arrayloc[80],*lon,*lat;

    arrayloc[0]='\0';

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'S': options.stationfilename = optarg;
            break;
          case 'o': options.outccfilename = optarg;
            break;
          case 'c': options.crosscor_filename = optarg;
            break;
          case 'd': debug = 1;
            fprintf(fpd,"Debugging on...\n");
            break;
          case 'I': options.configfilename = optarg;
            break;
          case 'H': options.header_filename = optarg;
            break;
          case 'F': options.flagfilename=optarg;
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
        options.arr_lat_rad = atof(lat);    /* convert string in degrees to float */
        options.arr_lon_rad = atof(lon);
        fprintf(fpd,"User specified array lon,lat: %g, %g (degs)\n",options.arr_lon_rad,options.arr_lat_rad);
        options.arr_lat_rad *= (M_PI/180.0); // convert to radian
        options.arr_lon_rad *= (M_PI/180.0);
    }

    /* auto flagging requires autocorrelations */
    if (options.do_flag) {
        fprintf(stderr,"WARNING: no autoflag this.\n");
        exit(1);
    }
}


/***************************
 ***************************/
void printusage(const char *progname) {
  fprintf(stderr,"Usage: %s [options]\n\n",progname);
  fprintf(stderr,"options are:\n");
  fprintf(stderr,"-c filename\tThe name of the input cross-correlation data file. Default: stdin.\n");
  fprintf(stderr,"-o filename\tThe name of the output cross correlation file. Default: stdout.\n");
  fprintf(stderr,"-S filename\tThe name of the file containing antenna name and local x,y,z. Default: %s\n",
                    options.stationfilename);
  fprintf(stderr,"-I filename\tThe name of the file containing instrument config. Default: %s\n",options.configfilename);
  fprintf(stderr,"-H filename\tThe name of the file containing observing metadata. Default: %s\n",options.header_filename);
  fprintf(stderr,"-A lon,lat \tSpecify East Lon and Lat of array center (degrees). Comma separated, no spaces. Default: MWA\n");
  fprintf(stderr,"-l         \tLock the phase center to the initial HA/DEC\n");
  fprintf(stderr,"-F filename\tOptionally apply global flags as specified in filename.\n");
  fprintf(stderr,"-d         \tturn debugging on.\n");
  fprintf(stderr,"-v         \treturn revision number and exit.\n");
  exit(1);
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

  mjd = data->date[0] - 2400000.5;  // get Modified Julian date of scan.
  lmst = palRanorm(palGmst(mjd) + options.arr_lon_rad);  // local mean sidereal time, given array location

  /* if no RA was specified in the header,then calculate the RA based on lmst and array location
     and update the ra */
  if (header->ra_hrs < -98.0 ) {
    // set the RA to be for the middle of the scan
//    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start + header->n_scans*header->integration_time*1.00274/(3600.0*2);   // include 1/2 scan offset
    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start;  // match existing code. RA defined at start of scan
    if (debug) fprintf(fpd,"Calculated RA_hrs: %g of field centre based on HA_hrs: %g and lmst_hrs: %g\n",
                        header->ra_hrs,header->ha_hrs_start,lmst*(12.0/M_PI));
  }

  return 0;
}


/**************************
***************************/
void checkInputs(Header *header,uvdata *data,InpConfig *inputs) {
    int total_ants;

    if(inputs->n_inputs != header->n_inputs) {
        fprintf(stderr,"%s ERROR: mismatch between the number of inputs in %s (%d) and header (%d)\n",__func__,
                options.configfilename,inputs->n_inputs,header->n_inputs);
        exit(1);
    }

    total_ants = countPresentAntennas(inputs);

    if (total_ants > data->array->n_ant) {
        fprintf(stderr,"ERROR: mismatch between the number of antennas in %s (%d) and %s (%d)\n",
                options.stationfilename,data->array->n_ant,options.configfilename,total_ants);
    }

    if (header->corr_type != 'C') {
        // silently override the header type since we only care about the cross correlations
        if (debug) fprintf(fpd,"INFO: Overriding input type to 'C'\n");
        header->corr_type = 'C';
    }

    if (options.crosscor_filename==NULL && debug!=0 ) {
        fprintf(fpd,"INFO: reading from stdin\n");
    }
    if (options.outccfilename==NULL && debug!=0 ) {
        fprintf(fpd,"INFO: writing to stdout\n");
    }
    if (header->ha_hrs_start == -99.0 && lock_pointing) {
        fprintf(stderr,"ERROR: HA must be specified in header if -l flag is used.\n");
        exit(1);
    }
    
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
  strcpy(data->array->name,"UNKNOWN");
}


