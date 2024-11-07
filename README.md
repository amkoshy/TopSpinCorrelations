# TopSpinCorrelations
The final version of spin correlation analysis framework developed during my PhD. I inherited an older version of the framework and I made changes necessary for my analysis; that is measuring spin correlation between top quark and anti-quark as a function of extra-jet multiplicity.

The analysis was run at computer clusters at Purdue and CERN, often on many 100s of CPUs at a time. The analysis runs only in CMSSW environment, which is a custom environment developed by CERN scientists to analyse the large amount of physics data produced by CERN. 

The analysis code is written primarily in C++. Jobs are run in the clusters by using shell scripts, and SLURM scripts. Python is used in creating many of these scripts, and for automating various tasks.
