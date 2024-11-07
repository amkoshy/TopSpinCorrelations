#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

for c in ee emu mumu; do
    w
    $LA -f ttbarsignalplustau.root -c $c &
    $LA -f ttbarsignalplustau.root -c $c --bgviatau &
    $LA -f dy -d 11 -c $c &
    $LA -f dy -d 13 -c $c &
    $LA -f dy -d 15 -c $c &
    $LA -f ttbarZ -c $c &
done

for c in ee emu mumu; do
    w
    $LA -f ${c}_run2016B -c $c &
    $LA -f ${c}_run2016C -c $c &
    $LA -f ${c}_run2016D -c $c &
    $LA -f ${c}_run2016E -c $c &
    #$LA -f ${c}_run2016F -c $c &
    $LA -f ${c}_run2016F1 -c $c &
    $LA -f ${c}_run2016F2 -c $c &
    $LA -f ${c}_run2016G -c $c &
    $LA -f ${c}_run2016H -c $c &
    for s in se smu; do
        $LA -f ${s}_run2016B -c $c &
        $LA -f ${s}_run2016C -c $c &
        $LA -f ${s}_run2016D -c $c &
        $LA -f ${s}_run2016E -c $c &
        #$LA -f ${s}_run2016F -c $c &
        $LA -f ${s}_run2016F1 -c $c &
        $LA -f ${s}_run2016F2 -c $c &
        $LA -f ${s}_run2016G -c $c &
        $LA -f ${s}_run2016H -c $c &
    done
done

for i in zz qcd tw.root ttbarbg.root wtol ww wz ttbarW ttgjets; do
    w
    if [ "$isNAF" = 1 ]; then
        $LA -f $i
    else
        $LA -f $i -c ee&
        $LA -f $i -c emu&
        $LA -f $i -c mumu&
    fi
done

wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analy"
else
    echo "Processing all nominal samples finished!"
fi


