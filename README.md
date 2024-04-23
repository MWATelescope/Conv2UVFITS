# Conv2UVFITS
Convert raw correlator data to UVFITS format

This code was originally written for MWA 32T in 2008 by Randall Wayth, and expanded/upgraded a few times since. In particular the corrections for precession were added by Al Levine. The code originally used the SLALIB library, but was changed to use Starlink PAL, which has equivalent routines.

This repo also contains uvfits.c, which is a fairly self-contained library for reading and writing UVFITS files. The tool test_readuvfits is a simple standalone tool for reading a UVFITS file and dumping the output for testing.

## Requirements
cfitsio v3, Starlink PAL (Positional astronomy library).
Ubuntu packages: libcfitsio-dev libstarlink-pal-dev build-essential

## Compiling
make corr2uvfits

run the program with no command-line arguments for a usage summary.

## Metadata
There are 3 metadata input files required:

- a "header" file, which contains information about the time, frequency, phase centre etc of the observation.
- an "antenna_locations" file, which contains the names and coords of the antennas relative to the defined array centre
- an "instrument_config" file, which maps the logical ordering of signal "inputs" to antennas/polarisations, and hence correlation products. This can also specify flags for specific inputs and optionally a static delay correction factor that can be used to correct for differential delays in the signals (e.g. unequal cable lengths), applied as a frequency dependent phase.

Example files with the formatting etc are included in this repo. The actual files can be called anything and are specified on the command line, but default names are looked for also.

## Raw data format
The raw data are in "L-file" format (the origins of which are lost to history, but passed through Frank Briggs to the MWA project).  L-files are raw binary single-precision float files, with separate files for auto and cross correlations. Historically, these files have had suffixes ".LACSPC" and ".LCCSPC", but the files can be called anything. A simple and fast GPU-based correlator that generates this output format was also used in MWA 32T days, [found here](https://github.com/MWATelescope/corr_GPU).

The implicit dimensions of the binary data in the autocorrelation file is [time][input][channel], where channel changes most quickly in the data and time most slowly.
"input" here is a logical data input, so for dual polarisation antennas, there will be 2 inputs per antenna. A common ordering for "input" will thus be something like Ant1-X, Ant1-Y, Ant2-X, Ant2-Y etc.
Mapping of logical inputs to actual antennas is done in the metadata files as described aboove.

The implicit dimensions of the binary data in the cross correlation file is [time][product_index][channel], again where channel changes most quickly and time most slowly.
Cross correlations are single precision float complex numbers. The implicit ordering of inputs to products is inp1-inp2, inp1-inp3, ... inp2-inp3, inp2-inp4 etc.

Note that the code does not read ahead or look at the size of the input data, so named pipes can be used for the .LACSPC and .LCCSPC files to avoid creating intermediate files if desired in a streaming environment.

## Example usage

