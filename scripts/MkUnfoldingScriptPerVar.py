import sys
import argparse

def main () :

    parser = argparse.ArgumentParser()
    parser.add_argument("mode")
    parser.add_argument("-m", "--mode", action="store_true", help="Which mode to run in :: Nominal, Bootstrap, AllSysts, Plotting")

    args = parser.parse_args()
    mode = args.mode

    with open("TUnfoldVariablesList.txt", "r") as ipfile :
        varnames = ipfile.readlines()

    eras = [
    "2016ULpreVFP",
    "2016ULpostVFP",
    "2016UL",
    "2017UL",
    "2018UL",
    "fullRun2UL",
    ]

    for era in eras:
        unffile  = open("scripts/Unfold_vars_" + era + "_" + mode +".sh", "w")
        
        for var in varnames :
            var = var.strip('\n')

            if   mode == "Bootstrap" :
                unffile.write("nohup ./install/bin/MatrixUnfControl " + era + " " + var + " Nominal combined 1 1 1 1 0 0 0 0 ")
                unffile.write("&> nohuplogs/Unfold_vars_" + era + "_" + var + "_" + mode + ".out")
                unffile.write("\n")

            elif mode == "Nominal"  :
                unffile.write("nohup ./install/bin/MatrixUnfControl " + era + " " + var + " Nominal combined 0 1 1 0 0 0 1 0 ")
                unffile.write("&> nohuplogs/Unfold_vars_" + era + "_" + var + "_" + mode +".out")
                unffile.write("\n")
            
            elif mode == "AllSysts" :
                unffile.write("nohup ./install/bin/MatrixUnfControl " + era + " " + var + " Nominal combined 0 1 1 0 0 2 1 0 ")
                unffile.write("&> nohuplogs/Unfold_vars_" + era + "_" + var + "_" + mode +".out")
                unffile.write("\n")

            elif mode == "Plotting" :
                unffile.write("nohup ./install/bin/MatrixUnfControl " + era + " " + var + " Nominal combined 0 1 1 0 0 2 1 2 ")
                unffile.write("&> nohuplogs/Unfold_vars_" + era + "_" + var + "_" + mode +".out")
                unffile.write("\n")

            else :
                print("Mode not supported, exiting now")
                sys.exit()
    
        unffile.close()

if __name__ == "__main__" :
    main()
