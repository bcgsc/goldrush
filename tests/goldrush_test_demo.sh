#!/bin/bash

set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldRush"
goldrush run reads=test_reads G=1e6 t=4 p=goldrush_test -B

if [ -e goldrush_test_golden_path.goldrush-edit-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa ]; then
  echo "Test successful"
else
  echo "Final expected file not found - please check your installation"
  exit(1)
fi

exit(0)
