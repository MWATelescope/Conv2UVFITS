/************************
*************************/

#ifndef CONVUTIL_H
#define CONVUTIL_H

#include "uvfits.h"

#define MAX_ANT 256
#define MAX_LINE 1024
#define MWA_LAT -26.703319        // Array latitude. degrees North
#define MWA_LON 116.67081         // Array longitude. degrees East
#define MWA_HGT 377               // Array altitude. meters above sea level
#define EARTH_RAD 6378100.0       // meters
#define SIZ_PRODNAME 8            // size of char array with text type of pol products.
#define SIZ_TELNAME  8            // size of char with text for telescope and instrument names
#define VLIGHT 299792458.0        // speed of light. m/s
#define AUTOFLAG_N_NEIGHBOURS 3   // number of neighbours to use to form a median for channel flagging.
#define VEL_FACTOR  1.204         // the velocity factor of electic fields in RG-6 like coax.

#define CHAN_ALL_ANT_ALL_TIME      "CHAN_ALL_ANT_ALL_TIME"

typedef struct _header {
    int     n_inputs;               // the number of inputs to the corrleator
    int        n_scans;                // the number of time samples
    int        n_chans;                // the number of spectral channels
    int     corr_type;              // cross correlation, auto or both 'C', 'A', 'B'
    int     invert_freq;            // flag to indicate that freq decreases with increasing channel number.
    int     conjugate;              // conjugate the final vis to correct for any sign errors
    int     geom_correct;           // apply geometric phase correction
    float   integration_time;        // per time sample, in seconds
    double  cent_freq,bandwidth;    // observing central frequency and bandwidth (MHz)
    double  ra_hrs,dec_degs;        // ra,dec of phase centre.
    double  ha_hrs_start;           // the HA of the phase center at the start of the integration
    double  ref_el,ref_az;          // the el/az of the normal to the plane of the array (radian)
    int     year,month,day;            // date/time in UTC.
    int     ref_hour,ref_minute;
    float   ref_second;
    char    field_name[SIZE_SOURCE_NAME+1];
    char    pol_products[SIZ_PRODNAME+1];
    char    telescope[SIZ_TELNAME+1];
    char    instrument[SIZ_TELNAME+1];
} Header;

typedef struct {
    int n_inputs;
    int ant_index[MAX_ANT*2];
    float cable_len_delta[MAX_ANT*2];
    char pol_index[MAX_ANT*2];
    char inpFlag[MAX_ANT*2];
} InpConfig;


/* public function prototypes */
int createPolIndex(char *polprods, int *index);
int readArray(char *filename, double lat, double lon, array_data *array);
int readHeader(char *filename, Header *header);
int readInputConfig(char *filename, InpConfig *inp);
void calcUVW(double ha,double dec,double x,double y,double z,double *u,double *v,double *w);
int decodePolChar(int pol_char);
int decodePolIndex(int pol1, int pol2);
void azel2xyz(double az, double el, double *x, double *y, double *z);
int autoFlag(uvdata *uvdata, float sigma, int n_neighbours,int flag_type);
int flag_antenna(uvdata *uvdata, const int ant, const int pol, const int chan, const int t);
int flag_all_antennas(uvdata *uvdata, const int chan, const int t);
int applyFlagsFile(char *flagfilename,uvdata *uvdata);
int countPresentAntennas(InpConfig *inputs);
int checkAntennaPresent(InpConfig *inputs, int ant_index);
void precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
         double *xp, double *yp, double *zp, double lmst2000);
void rotate_radec(double rmat[3][3], double ra1, double dec1,
          double *ra2, double *dec2);
void aber_radec_rad(double eq, double mjd, double ra1, double dec1,
            double *ra2, double *dec2);
void stelaber(double eq, double mjd, double v1[3], double v2[3]);
void mat_transpose(double rmat1[3][3], double rmat2[3][3]);
void unitvecs_j2000(double rmat[3][3], double xhat[3], double yhat[3], double zhat[3]);
void ha_dec_j2000(double rmat[3][3], double lmst, double lat_rad, double ra2000,
                  double dec2000, double *newha, double *newlat, double *newlmst);
int makeBaselineLookup(InpConfig *inps, Header *header, array_data *array, int bl_ind_lookup[MAX_ANT][MAX_ANT]);
void setConvDebugLevel(int level);
int calcAntPhases(double mjd, Header *header, array_data *array,double ant_u[], double ant_v[], double ant_w[]);
int correctPhases(double mjd, Header *header, InpConfig *inps, array_data *array, int bl_ind_lookup[MAX_ANT][MAX_ANT], float *ac_data, float complex *cc_data,double ant_u[MAX_ANT],double ant_v[MAX_ANT],double ant_w[MAX_ANT]);

#endif
