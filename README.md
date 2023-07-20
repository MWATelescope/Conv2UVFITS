# Conv2UVFITS
Convert raw correlator data to UVFITS format

This code was originally written for MWA 32T in 2008 by Randall Wayth, and expanded/upgraded a few times since. In particular the corrections for precission were added by Al Levine.

This repo also contains uvfits.c, which is a fairly self-contained library for reading and writing UVFITS files. The tool test_readuvfits is a simple standalone tool for reading a UVFITS file and dumping the output for testing.

## Requirements
cfitsio v3, Starlink PAL (Positional astronomy library)
Ubuntu packages: libcfitsio3-dev libstarlink-pal-dev build-essential

## Compiling
make corr2uvfits
## Metadata

## Example usage

