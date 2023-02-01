#!/bin/bash
module load beagle-lib
module load X11 
module load Beast/1.10.5pre-GCCcore-8.3.0
beast -threads 4 -save_every 2000000 -beagle $1
