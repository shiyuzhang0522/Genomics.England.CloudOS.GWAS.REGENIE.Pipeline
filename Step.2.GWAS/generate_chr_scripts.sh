#!/usr/bin/env bash
set -euo pipefail


################################################################################
#
# Generate chromosome-specific REGENIE Step 2 launcher scripts
#
# Author:
#   Shelley
#
# Date:
#   2026-08-13
#
################################################################################


SCRIPT_DIR="$(pwd)"

COMMAND_DIR="${SCRIPT_DIR}/commands"


mkdir -p "${COMMAND_DIR}"


for CHR in {1..22}

do

cat > "${COMMAND_DIR}/run_chr${CHR}.sh" << EOF
#!/usr/bin/env bash

bash ${SCRIPT_DIR}/run_REGENIE_step2_chr.sh ${CHR}

EOF


chmod +x "${COMMAND_DIR}/run_chr${CHR}.sh"


done


echo
echo "Generated chromosome scripts:"
ls -lh "${COMMAND_DIR}"
