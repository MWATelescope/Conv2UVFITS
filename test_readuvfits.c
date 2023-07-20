/*
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <complex.h>
#include "uvfits.h"

/* globals */
char *infilename=NULL;
int debug=0;
FILE *fpd;

void usage(char * const argv[]) {
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"%s <options> -i inputfilename\n",argv[0]);
    fprintf(stderr,"\t-d debug_level (0=no debugging)\n");
    exit(0);
}


void parse_cmdline(const int argc,char * const argv[]) {
    int result=0;
    const char optstring[] = "d:i:";

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'i': infilename = optarg;
            break;
          case 'd': debug = atoi(optarg);
            break;
          default:
              fprintf(stderr,"unknown option: %c\n",result);
              usage(argv);
        }
    }
    if (infilename==NULL) {
        usage(argv);
    }
}


int main(int argc, char *argv[]) {
    int res=0,chunk=0,i;
    uvdata *data_new;
    uvReadContext *iter;

    fpd = stderr;
    if (argc < 2) usage(argv);
    parse_cmdline(argc,argv);

    // set uvfits debugging output
    uvfitsSetDebugLevel(debug);

    if (debug) fprintf(fpd,"Reading file: %s\n",infilename);
    //res = readUVFITS(infilename,&data_old);
    if (res !=0) {
        fprintf(stderr,"readUVFITS failed with error %d\n",res);
    }

    /* now try the new way... */
/* */
    res = readUVFITSInitIterator(infilename, &data_new, &iter);
    if (res !=0) {
        fprintf(stderr,"readUVFITSInitIterator failed with error %d\n",res);
        exit(1);
    }
    if (debug) fprintf(fpd,"n_vis: %d, n_pol: %d, gcount: %d\n",data_new->n_vis,data_new->n_pol, iter->gcount);
    fflush(stdout);
    while ((res=readUVFITSnextIter(data_new,iter)) ==0) {
        fprintf(stdout,"Chunk %d. Time: %f. baselines: %d, grp_index: %d\n",chunk++,data_new->date[0],data_new->n_baselines[0],iter->grp_index);
    }
    if (res != 1) {
        fprintf(stderr,"readUVFITSnextIter returned %d\n",res);
        exit(1);
    }
    readUVFITSCloseIter(iter);
    if (res < 0) exit(res);

    fprintf(fpd,"N pols: %d, Pol type: %d\n",data_new->n_pol,data_new->pol_type);

    /* print a few vis values */
    for (i=0; i<10; i++) {
        float complex *visdata,temp;
        int f,p;
        visdata = (float complex *)data_new->visdata[0];
        fprintf(stdout,"u,v,w: %g,%g,%g. Baseline: %.6f.\n",data_new->u[0][i],data_new->v[0][i],data_new->w[0][i],data_new->baseline[0][i]);
        for (f=0; f<data_new->n_freq; f++) {
            fprintf(stdout,"Freq %g ",data_new->cent_freq +data_new->freq_delta*(f-data_new->n_freq/2));
            for(p=0; p<data_new->n_pol; p++) {
                temp=visdata[p+data_new->n_pol*f];
                fprintf(stdout,"%g,%g ",crealf(temp),cimagf(temp));
            }
            fprintf(stdout,". Weights: ");
            for(p=0; p<data_new->n_pol; p++) {
                fprintf(stdout,"%g ",data_new->weightdata[0][p+data_new->n_pol*f]);
            }
            fprintf(stdout,"\n");
        }
    }

    freeUVFITSdata(data_new);

    return 0;
}
