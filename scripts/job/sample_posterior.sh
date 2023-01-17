#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=12:00:00
#PJM -L rscgrp=gg57

# pjsub -X KEY=VALUE,KEY=VALUE -N ${JOB_NAME} -o ${LOG_FILE} -j --comment hoge job.sh

echo running $ARG
echo output $OUTPUT
module load python/3.7.3
/work/00/gg57/j29006/dbgphmm/target/release/sample_posterior $ARG > $OUTPUT
