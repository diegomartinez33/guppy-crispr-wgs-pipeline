#!/bin/bash
# 1. WGS individual
WGS_JOB=$(sbatch --parsable crispresso_wgs.sh)
echo "WGS job: $WGS_JOB"

# 2. Aggregate por grupo — después del WGS
AGG_JOB=$(sbatch --parsable \
  --dependency=afterok:${WGS_JOB} \
  crispresso_wgs_aggregate.sh)
echo "Aggregate job: $AGG_JOB"