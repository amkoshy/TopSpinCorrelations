#!/bin/sh

set -e

if [ "$1" = "" ]; then
    echo "Please provide a name - your analysis folders will be renamed"
    echo "For example: scripts/renameSelection testCut"
    echo " this would mv selectionRoot selectionRoot_testCut"
    echo " also the BTag, Plots, preunfolded, and Unfolding folders will be renamed"
    exit
else
    name="$1"
fi

for i in BTagEff Plots selectionRoot SVD UnfoldingResults preunfolded; do
    echo "renaming $i to ${i}_$name"
    mv "$i" "${i}_$name"
done

