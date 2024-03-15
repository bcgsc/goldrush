#!/bin/bash

set -eu -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldRush"
goldrush run reads=test_reads G=1e6 t=4 p=goldrush_test shared_mem=./tmp -B

l50=$(abyss-fac goldrush_test_golden_path.goldpolish-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa |awk '{print $3}' |tail -n1)

if [ -e goldrush_test_golden_path.goldpolish-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa ] && [ ${l50} -eq 1 ]; then
  echo -e "\nTest successful!"
else
  echo -e "\nTest failed - please check your installation"
  exit 1
fi

exit 0
