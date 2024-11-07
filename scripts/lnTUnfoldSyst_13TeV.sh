#!/bin/zsh

BASEDIR=`pwd`

cd UnfoldingHistos
echo "Creating links for alternative ttbar MC samples ... " 
ln -s AMCATNLOFXFX POWHEG
ln -s POWHEGV2HERWIG MCATNLO
ln -s MADGRAPHMLM POWHEGHERWIG
cd $BASEDIR

foreach era (2017UL 2018UL 2016ULpreVFP 2016ULpostVFP)

    foreach channel (ee emu mumu)
    
        foreach syst (BFRAG_DOWN \
            BFRAG_UP \
            BSEMILEP_DOWN \
            BSEMILEP_UP \
            BTAG_DOWN \
            BTAG_ETA_DOWN \
            BTAG_ETA_UP \
            BTAG_LJET_DOWN \
            BTAG_LJET_ETA_DOWN \
            BTAG_LJET_ETA_UP \
            BTAG_LJET_PT_DOWN \
            BTAG_LJET_PT_UP \
            BTAG_LJET_UP \
            BTAG_PT_DOWN \
            BTAG_PT_UP \
            BTAG_UP \
            ERDON \
            ERDONRETUNE \
            GLUONMOVETUNE \
            JER_DOWN \
            JER_UP \
            JES_DOWN \
            JES_UP \
            KIN_DOWN \
            KIN_UP \
            LEPT_DOWN \
            LEPT_UP \
            MASS_DOWN \
            MASS_UP \
            MATCH_DOWN \
            MATCH_UP \
            MCATNLO \
            MEFACSCALE_DOWN \
            MEFACSCALE_UP \
            MERENSCALE_DOWN \
            MERENSCALE_UP \
            MESCALE_DOWN \
            MESCALE_UP \
            PDF_ALPHAS_DOWN \
            PDF_ALPHAS_UP \
            POWHEG \
            POWHEGHERWIG \
            PSFSRSCALE_DOWN \
            PSFSRSCALE_UP \
            PSISRSCALE_DOWN \
            PSISRSCALE_UP \
            PU_DOWN \
            PU_UP \
            TRIG_DOWN \
            TRIG_ETA_DOWN \
            TRIG_ETA_UP \
            TRIG_UP \
            UETUNE_DOWN \
            UETUNE_UP \
            UNCLUSTERED_DOWN \
            UNCLUSTERED_UP
        )
            if [ -d UnfoldingHistos/$era/$syst/$channel ] ; then
                cd UnfoldingHistos/$era/$syst/$channel/
                echo
                echo "Creating Links in ... " 
                pwd
                echo 
                ln -s ../../../$era/Nominal/$channel/*.root .
                cd $BASEDIR
            else
                mkdir -p UnfoldingHistos/$era/$syst/$channel
                cd UnfoldingHistos/$era/$syst/$channel/
                echo
                echo "Creating Links in ... " 
                pwd
                echo 
                ln -s ../../../$era/Nominal/$channel/*.root .
                cd $BASEDIR
            fi
        end

    end
end
