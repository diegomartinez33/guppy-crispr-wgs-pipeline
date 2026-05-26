#!/bin/bash
# Después de MarkDuplicates:

# 1. On-target (en paralelo con GATK)
# OT_JOB=$(sbatch --parsable \
#   crispresso_ontarget.sh)

# # 2. Predicción off-targets (también en paralelo)
# OFF_JOB=$(sbatch --parsable \
#   casoffinder.sh)

# 3. Pooled por grupo — después de on-target
POOL_F=$(sbatch --parsable \
  crispresso_pooled_groups.sh)

# 4. Compare grupos pooled — después de pooled
COMP_F=$(sbatch --parsable \
  --dependency=afterok:${POOL_F} \
  crispresso_compare_groups.sh)

# # 4. WGS off-target (después de cas-offinder y markdup)
# sbatch --dependency=afterok:${OT_JOB}:${OFF_JOB} \
#   crispresso_wgs.sh