mkdir eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance
cp onia2ySkim.C doFitUpsilonOneCB.C PsetCollection.h cutsAndBin.h eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance

root -l -b <<EOF
.L doFitUpsilonOneCB.C++
.q
EOF

echo entering to the loop...

for pt in '0,5' '5,100'
do
    for rap in '0,1.2' '1.2,2.4'
    do
	for i in 0 2   # 0=kPPDATA, 2=kAADATA
	do
	    ./condor_root.sh eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance 'doFitUpsilonOneCB.C+('$i','$pt','$rap',0,200,0,0.5,0)'  
	done
# Centrlaity dependence 
	for cent in '0,60' '60,200'
	do
	    ./condor_root.sh eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0,0.5,0)'    # PbPb only
	done 
# Event Plane dependence
	for cent in '20,120' 
	do
	    ./condor_root.sh eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0,0.167,0)'    # PbPb only
	    ./condor_root.sh eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0.167,0.333,0)'    # PbPb only
	    ./condor_root.sh eventPlane_Jan26_0-5-100GeV_0-120cent_3EventPlaneBin_looseAcceptance 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0.333,0.5,0)'    # PbPb only
	done 

    done
done

