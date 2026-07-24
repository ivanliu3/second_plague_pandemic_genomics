#!/usr/bin/env bash

set -eo pipefail

DEMO_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="${DEMO_DIR}/../mahalanobisD/pcaDistCalc.mahalanobis.py"

cd "${DEMO_DIR}"

python3 "${SCRIPT}"  \
-tf ancient_pca.tsv  \
-tp 4:10 \
-ti 0 \
-rf reference_pca.tsv \
-rp 4:10 \
-rg 1 \
-ri 0

echo
echo "Demo completed successfully."
echo "Generated output files:"
echo "  data.mahadist.csv"
echo "  data.pvalue.csv"
