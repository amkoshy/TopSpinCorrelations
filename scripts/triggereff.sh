#!/bin/sh

# 2016preVFP UL
nohup ./install/bin/triggereff -y 2016preVFP_UL -t Data -s run2016_ULB &> nohuplogs/2016preVFP_UL_Data_B.out 
nohup ./install/bin/triggereff -y 2016preVFP_UL -t Data -s run2016_ULC &> nohuplogs/2016preVFP_UL_Data_C.out 
nohup ./install/bin/triggereff -y 2016preVFP_UL -t Data -s run2016_ULD &> nohuplogs/2016preVFP_UL_Data_D.out 
nohup ./install/bin/triggereff -y 2016preVFP_UL -t Data -s run2016_ULE &> nohuplogs/2016preVFP_UL_Data_E.out 
nohup ./install/bin/triggereff -y 2016preVFP_UL -t Data -s run2016_ULF1 &> nohuplogs/2016preVFP_UL_Data_F1.out 

#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 0 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_0.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 1 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_1.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 2 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_2.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 3 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_3.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 4 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_4.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 5 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_5.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 6 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_6.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 7 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_7.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 8 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_8.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 9 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_9.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 10 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_10.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 11 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_11.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 12 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_12.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 13 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_13.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 14 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_14.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 15 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_15.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 16 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_16.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 17 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_17.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 18 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_18.out
#nohup ./install/bin/triggereff -y 2016preVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 19 &> nohuplogs/2016preVFP_UL_MC_ttbarsignalplustau_fromDilepton_19.out

# 2016postVFP UL
nohup ./install/bin/triggereff -y 2016postVFP_UL -t Data -s run2016_ULF2 &> nohuplogs/2016postVFP_UL_Data_F2.out 
nohup ./install/bin/triggereff -y 2016postVFP_UL -t Data -s run2016_ULG &> nohuplogs/2016postVFP_UL_Data_G.out 
nohup ./install/bin/triggereff -y 2016postVFP_UL -t Data -s run2016_ULH &> nohuplogs/2016postVFP_UL_Data_H.out 

#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 0 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_0.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 1 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_1.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 2 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_2.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 3 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_3.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 4 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_4.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 5 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_5.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 6 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_6.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 7 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_7.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 8 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_8.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 9 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_9.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 10 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_10.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 11 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_11.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 12 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_12.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 13 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_13.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 14 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_14.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 15 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_15.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 16 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_16.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 17 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_17.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 18 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_18.out
#nohup ./install/bin/triggereff -y 2016postVFP_UL -t MC -s ttbarsignalplustau_fromDilepton -n 19 &> nohuplogs/2016postVFP_UL_MC_ttbarsignalplustau_fromDilepton_19.out

# 2017 UL
nohup ./install/bin/triggereff -y 2017_UL -t Data -s run2017_ULB &> nohuplogs/2017_UL_Data_B.out 
nohup ./install/bin/triggereff -y 2017_UL -t Data -s run2017_ULC &> nohuplogs/2017_UL_Data_C.out 
nohup ./install/bin/triggereff -y 2017_UL -t Data -s run2017_ULD &> nohuplogs/2017_UL_Data_D.out 
nohup ./install/bin/triggereff -y 2017_UL -t Data -s run2017_ULE &> nohuplogs/2017_UL_Data_E.out 
nohup ./install/bin/triggereff -y 2017_UL -t Data -s run2017_ULF &> nohuplogs/2017_UL_Data_F.out 

#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 0 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_0.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 1 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_1.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 2 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_2.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 3 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_3.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 4 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_4.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 5 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_5.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 6 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_6.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 7 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_7.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 8 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_8.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 9 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_9.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 10 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_10.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 11 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_11.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 12 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_12.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 13 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_13.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 14 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_14.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 15 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_15.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 16 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_16.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 17 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_17.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 18 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_18.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 19 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_19.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 20 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_20.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 21 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_21.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 22 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_22.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 23 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_23.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 24 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_24.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 25 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_25.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 26 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_26.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 27 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_27.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 28 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_28.out
#nohup ./install/bin/triggereff -y 2017_UL -t MC -s ttbarsignalplustau_fromDilepton -n 29 &> nohuplogs/2017_UL_MC_ttbarsignalplustau_fromDilepton_29.out

# 2018 UL
nohup ./install/bin/triggereff -y 2018_UL -t Data -s run2018_ULA &> nohuplogs/2018_UL_Data_A.out 
nohup ./install/bin/triggereff -y 2018_UL -t Data -s run2018_ULB &> nohuplogs/2018_UL_Data_B.out 
nohup ./install/bin/triggereff -y 2018_UL -t Data -s run2018_ULC &> nohuplogs/2018_UL_Data_C.out 
nohup ./install/bin/triggereff -y 2018_UL -t Data -s run2018_ULD &> nohuplogs/2018_UL_Data_D.out 

#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 0 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_0.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 1 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_1.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 2 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_2.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 3 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_3.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 4 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_4.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 5 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_5.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 6 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_6.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 7 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_7.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 8 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_8.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 9 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_9.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 10 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_10.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 11 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_11.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 12 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_12.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 13 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_13.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 14 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_14.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 15 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_15.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 16 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_16.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 17 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_17.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 18 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_18.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 19 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_19.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 20 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_20.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 21 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_21.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 22 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_22.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 23 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_23.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 24 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_24.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 25 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_25.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 26 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_26.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 27 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_27.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 28 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_28.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 29 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_29.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 30 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_30.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 31 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_31.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 32 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_32.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 33 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_33.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 34 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_34.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 35 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_35.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 36 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_36.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 37 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_37.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 38 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_38.out
#nohup ./install/bin/triggereff -y 2018_UL -t MC -s ttbarsignalplustau_fromDilepton -n 39 &> nohuplogs/2018_UL_MC_ttbarsignalplustau_fromDilepton_39.out

#nohup ./install/bin/triggereff -y 2017 -t Data -s run2017B &> nohuplogs/2017_Data_B.out 
#nohup ./install/bin/triggereff -y 2017 -t Data -s run2017C &> nohuplogs/2017_Data_C.out 
#nohup ./install/bin/triggereff -y 2017 -t Data -s run2017D &> nohuplogs/2017_Data_D.out 
#nohup ./install/bin/triggereff -y 2017 -t Data -s run2017E &> nohuplogs/2017_Data_E.out 
#nohup ./install/bin/triggereff -y 2017 -t Data -s run2017F &> nohuplogs/2017_Data_F.out 
#nohup ./install/bin/triggereff -y 2017 -t MC -s ttbarsignalplustau_fromDilepton &> nohuplogs/2017_MC_ttbarsignalplustau_fromDilepton.out
