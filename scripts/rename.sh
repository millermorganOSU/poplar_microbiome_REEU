#!/bin/sh

for f in blast_*; do
f1=$(head -n 1 "$f" | cut -c8-12)
name="blastID_${f1}"
echo "$name"
cp -n "$f" "blasts/$name"
done

