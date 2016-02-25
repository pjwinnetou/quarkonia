echo compiling....
root -l -b -q test1_doFitUpsilon.C++

echo entering to the loop...

for i in  0 2    # 0=kPPDATA, 2=kAADATA
do
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,6.5,0,1.2,0,200)'  #Should be no spaces between arguments
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,6.5,1.2,2.4,0,200)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',6.5,10,0,1.2,0,200)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',6.5,10,1.2,2.4,0,200)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,30,0,1.2,0,200)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,30,1.2,2.4,0,200)'
done

# Centrlaity dependence 
for i in 2    # 0=kPPDATA, 2=kAADATA
do
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,10,0,1.2,0,60)'  #Should be no spaces between arguments
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,10,1.2,2.4,0,60)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,10,0,1.2,60,200)'  #Should be no spaces between arguments
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',3,10,1.2,2.4,0,200)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,50,0,1.2,0,60)'  #Should be no spaces between arguments
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,50,1.2,2.4,0,60)'
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,50,0,1.2,60,200)'  #Should be no spaces between arguments
    ./condor_root.sh fitResults 'test1_doFitUpsilon.C+('$i',10,50,1.2,2.4,0,200)'
done 
