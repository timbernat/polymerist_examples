#!/bin/bash

example_dirs=(
    '1-polymerization'
    '2-preparation'
    '3-workflows'
    '4-miscellaneous'
)
batch_status="passing"

for dirname in "${example_dirs[@]}"; do
    pushd "./$dirname"
    source "./notebook_tests.sh" || batch_status="failing"
    popd
done

if [ $batch_status = "failing" ]; then
    exit 1 # bubble up failure to CI if any of batches fail
fi