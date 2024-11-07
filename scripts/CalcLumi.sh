mkdir -p RunLumis
brilcalc lumi -c web --begin 273150 --end 278807 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016preVFP.txt
brilcalc lumi -c web --begin 273150 --end 275376 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016B.txt
brilcalc lumi -c web --begin 275656 --end 276283 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016C.txt
brilcalc lumi -c web --begin 276315 --end 276811 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016D.txt
brilcalc lumi -c web --begin 276831 --end 277420 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016E.txt
brilcalc lumi -c web --begin 277932 --end 278807 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016F1.txt

brilcalc lumi -c web --begin 278769 --end 284044 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016postVFP.txt
brilcalc lumi -c web --begin 278769 --end 278808 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016F2.txt
brilcalc lumi -c web --begin 278820 --end 280385 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016G.txt
brilcalc lumi -c web --begin 281613 --end 284044 -u /fb  --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i ../../python/Data/Run2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt > RunLumis/2016H.txt