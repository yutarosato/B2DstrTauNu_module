#!/bin/bash

bsub -q e ./make_script_gMC_single.sh uds
bsub -q e ./make_script_gMC_single.sh charm
bsub -q e ./make_script_gMC_single.sh mixed
bsub -q e ./make_script_gMC_single.sh charged


