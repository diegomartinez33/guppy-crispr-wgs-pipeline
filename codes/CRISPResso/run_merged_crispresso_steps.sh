#!/bin/bash
mkdir -p logs/

# 1. Merge BAMs
#MERGE_JOB=$(sbatch --parsable merge_bams.sh)
#echo "Merge job: $MERGE_JOB"

# 2. CRISPResso2 on-target — después del merge
CRISPR_JOB=$(sbatch --parsable \
  crispresso_ontarget_merged.sh)
echo "CRISPResso2 job: $CRISPR_JOB"

# 3. Compare — después de CRISPResso2
COMP_JOB=$(sbatch --parsable \
  --dependency=afterok:${CRISPR_JOB} \
  crispresso_compare_merged.sh)
echo "Compare job: $COMP_JOB"