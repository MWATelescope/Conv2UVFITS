/* structure definitions and utilities for reading & writing
 * UV FITS files in C. Randall Wayth. Sep, 2006
 * 2013: Update for "iterator" single time step operation to reduce memory footprint
 * 2014: Update to re-sync this version and the one in the RTS
 */

#ifndef UVFITS_H
#define UVFITS_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define SIZE_ANT_NAME 8
#define SIZE_SOURCE_NAME 16
#define SIZE_POL_TYPE 4

/* constants for WGS84 Geoid model for the Earth */
#define EARTH_RAD_WGS84 6378137.0 /* meters */
#define E_SQUARED 6.69437999014e-3

typedef struct  _ant_table {
    char name[SIZE_ANT_NAME+1];
    double xyz_pos[3];			/* position relative to array centre  */
    float  xyz_deriv[3];
    /* skip orbital params for now */
    int  station_num;
    int  mount_type;
    double axis_offset[3];
    float pol_angleA;
    float pol_calA;
    float pol_angleB;
    float pol_calB;
    char pol_typeA[SIZE_POL_TYPE];
    char pol_typeB[SIZE_POL_TYPE];
} ant_table;

typedef struct _source_table {
    char name[SIZE_SOURCE_NAME+1];
    int  id;
    int  qual;
    char calcode[4];
    int  freq_id;
    double ra;    /* decimal hours */
    double dec;   /* decimal degrees */
} source_table;

typedef struct _array_table {
    int   n_ant;
    double xyz_pos[3];    /* the X,Y,Z coord of the array in conventional radio astronomy units */
    char  name[16];       /* name of the telescope (turns into TELESCOP keyword in fits file */
    char  instrument[16]; /* name of the instrument (turns into the INSTRUME keyword in fits file */
    ant_table *antennas;  /* a pointer to an array of ant_tables the size of the number of antennas */
    double arr_lon_rad;   /* array centre lon in radian */
    double arr_lat_rad;   /* array centre lat in radian */
} array_data;

typedef struct _fq_table {
    double freq; // The reference frequency (is a double in FQ table)
    float chbw; // The individual channel bandwidth (correlator channel size)
    float bandwidth; // total bandwidth for this IF
    int sideband; // -1 = lower, +1 = upper
} fq_table;


/* see notes on iterator mode below for how this data structure is used in that case */
/* this data structure reflects the old way the code was originally written where all data
   was read or written in a single chunk- hence was all in memory at the same time.
   The new iterator methods use this same structure, but only with 1 time step of data.
   Hence, if this all looks more complicated than it need to be, that is because it stays
   backwards compatible even though some of these arrays don't need to be arrays any more, 
   just single numbers.
*/
typedef struct _uvdata {
    int  n_pol;
    int  pol_type;        /* index to signal what kind of data: 1: Stokes, -1: Circ, -5: Linear */
    int  n_IF;            /* number of IFs in the data */
    int  n_freq;          /* number of freq channels if there is only 1 IF */
    int  n_vis;           /* number of sets of visibilities, one for each time instant */
    float cent_freq;      /* Hz */
    float freq_delta;
    double *date;         /* Julian date. array the size of n_vis */
    int  *n_baselines;    /* array the size of n_vis. number of baselines for each scan. */
    source_table *source; /* a pointer to a source table */
    array_data *array;    /* a pointer to an array struct */
    float **visdata;      /* array the size of n_vis whose elements point to arrays of visibiliites
                             the size of n_pol*n_freq*n_baselines complex floats. The data are ordered so
                             that pol changes most quickly, then freq, and baseline most slowly. Index with:
                             visdata[baseline*(n_pol*n_freq) + freq*n_pol + pol] as a float complex */
    float **weightdata;   /* weight data for visibiliites. Same data ordering as above, just float,
                             not complex */
    float **baseline;     /* same ordering again. encoded baseline using Miriad encoding convention */
    double **u;           /* arry the size of n_vis whose elements point to arrays of uvw
                             data the size of n_baselines */
    double **v;
    double **w;
    fq_table *fq;         /* pointer to array of FQ table for multi-IF data */
} uvdata;

/* create a mapping for the index of the various items that are in the visibilities. Since there are optional
   items, the index of an item in the group can change so this maps the item to the index in the group.
*/
typedef struct {
    int u;
    int v;
    int w;
    int date;
    int date0;  // try to support the non-standard CASA double-date format where this is the base date
    int bl;
    int su;
    int fq;
} ptype_mapping;


/* context data structure for iterative reading of UVFITS files */
typedef struct {
    int pcount,gcount,grp_index;
    float pscal[3],pzero[3];  // for u,v,w
    double base_date;         // date in PZERO5 which is the JD of the obs.
    float *grp_par, *grp_row;
    void *fptr;
    ptype_mapping ptype_map;
} uvReadContext;

/* context data structure for iterative writing of UVFITS files */
typedef struct {
    void *fptr;
    long ngroups,ntimes, n_grp_params;
    double jd_day_trunc;
} uvWriteContext;

/* public function prototypes */
int writeUVFITS(char *fname, uvdata *data); // old method - deprecated.
/* iterator mode API:
  In this mode, writeUVFITSiterator is called to set up the iterative method of writing 
  the data which uses much less memory. writeUVinstant is then use to write a single time
  instant of data. Finally, writeUVFITSfinalise is called to write the antenna table and close up
*/
int writeUVFITSiterator(char *filename, uvdata *data, uvWriteContext **iter);
int writeUVinstant(uvWriteContext *iter, uvdata *data, double jd_frac, int i);
int writeUVFITSfinalise(uvWriteContext *iter, uvdata *data);

/* reading API */
int readUVFITS(char *fname, uvdata **data); // old method - deprecated

/* iterator mode API:
   In this mode, the uvfits file is opened once and read many times, once for each time chunk
   in the file. This conserves system memory and allows processing of large uvfits files limited by
   disk read speed, not memory.

   To use this mode, don't use readUVFITS(). Instead use readUVFITSInitIterator() once to open the file
   and read the array table, then do a while loop on readUVFITSnextIter() to get a single time chunk of
   data from the file. See test_readuvfits.c for a simple example.

   In this case, the uvdata data structure only contains a single time instant (n_vis=1), so all the arrays of arrays
   have only 1 element in them and they all should be accessed via array index 0.
*/
int readUVFITSInitIterator(char *filename, uvdata **data, uvReadContext **iterator);
int readUVFITSnextIter(uvdata *obj, uvReadContext *iter);
int readUVFITSCloseIter(uvReadContext *iter);

/* utilities */
void printUVData(uvdata *data, FILE *fp);
void JD_to_Cal(double jd, int *year, int *month, int *day);
void JD_get_GSTIA0(double jd, double *GSTIA0);
void Cal_to_JD(int year, int month, int day, double *jd);
void uvfitsSetDebugLevel(int in_debug);
void freeUVFITSdata(uvdata *data);
void printAntennaData(array_data *array,FILE *fp);
void EncodeBaseline(int b1, int b2, float *result);
void DecodeBaseline(float blcode, int *b1, int *b2);
void Geodetic2XYZ(double lat_rad, double lon_rad, double height_meters, double *X, double *Y, double *Z);
void ENH2XYZ_absolute(double E,double N, double H, double lat_rad, double lon_rad, double *X, double *Y, double *Z);
void ENH2XYZ_local(double E,double N, double H, double lat, double *X, double *Y, double *Z);

#endif
