#!/bin/bash
# load newer version of java
#module load Java/1.8.0_181
module load beagle-lib

bash beast/bin/beast -threads 4 -beagle $1
