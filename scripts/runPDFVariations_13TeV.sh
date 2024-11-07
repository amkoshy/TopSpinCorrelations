#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

for ch in emu mumu ee; do
    for no in `seq 0 100`; do # don't care if maximum > number of pdf sets! It will just do nothing - for NNPDF3.0
        w
        $LA --pdf $no -c $ch -f ttbarsignalplustau.root &
    done
done

wait
echo "Done!"

