#!/bin/bash

example_dirs=(
    '1-polymerization'
    '2-preparation'
    '3-workflows'
    '4-miscellaneous'
)
for dirname in "${example_dirs[@]}"; do
    pushd "./$dirname"
    bash "./notebook_tests.sh"
    popd
done