#!/bin/sh

for filename in Xist/*_rep1.map; do
    python deltaSHAPE.py $filename "${filename/_rep1.map/_rep2.map}" --noplot -o "${filename/_rep1.map/.txt}"
done
