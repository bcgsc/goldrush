#!/bin/bash

set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldRush - this test may take 5-10min"
goldrush run reads=test_reads G=7e6 dev=True t=4 p=goldrush_test

if [ -e goldrush_test_golden_path.goldrush-edit-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa ]; then
  echo "Test successful"
else
  echo "Final expected file not found - please check your installation"
  exit(1)
fi

rm test_reads.fq

exit(0)
