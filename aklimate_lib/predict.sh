#!/bin/bash

COHORT=$1
DATATYPE=$2
SAMPLE_DATA_FILE=$3

if [[ -n $COHORT ]] && [[ -n $DATATYPE ]] && [[ -n $SAMPLE_DATA_FILE ]]; then
  echo "COHORT: $COHORT"
  echo "DATATYPE: $DATATYPE"
  echo "SAMPLE_DATA_FILE: $SAMPLE_DATA_FILE"
  Rscript --vanilla /aklimate/10_predict_samples_using_reduced_models.R $COHORT $DATATYPE $SAMPLE_DATA_FILE
else
  echo "There are three parameters:"
  echo "1) cohort code."
  echo "    - Allowed values are ACC, BLCA, BRCA, CESC, COADREAD, ESCC, GEA, HNSC, KIRCKICH, KIRP, LGGGBM, LIHCCHOL, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, SARC, SKCM, TGCT, THCA, THYM, UCEC, UVM."
  echo "2) datatype code."
  echo "    - Allowed values are GEXP, CNVR, METH, MULTI, TOP."
  echo "    - TOP will select the model with the highest prediction performance in 5-fold cross validation."
  echo "3) path to a sample data file where:"
  echo "    - first row gives feature IDs."
  echo "    - first column gives sample IDs."
  echo ""
  echo "If running from DOCKER, try something like:"
  echo "docker run --rm -v \`pwd\`:/data (IMAGE ID) (COHORT) (DATATYPE) /data/(SAMPLE_DATA_FILE)"
  echo "This command mounts the host current directory to /data in the Docker container so that files can be read and written."
fi
