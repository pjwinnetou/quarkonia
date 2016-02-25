# int collId = kPPDATA,                   float ptLow=5, float ptHigh=100,
# float yLow=0, float yHigh=1.2,                   int cLow=0, int cHigh=20,
# float dphiEp2Low=0, float dphiEp2High=0.5,   // In unit of PI!!
# float muPtCut=4.0,                   bool fixParameters=false


muPtCut='4'
fixParam='1'

outputDir='fitResults/fixParam'${fixParam}'MuPt'${muPtCut}'_Feb02-test'
mkdir -p $outputDir
cp onia2ySkim.C doFitUpsilonOneCB.C PsetCollection.h cutsAndBin.h $outputDir

root -l -b <<EOF
.L doFitUpsilonOneCB.C++
.q
EOF

echo entering to the loop...


for pt in '0,100'
do
    for rap in '0,2.4'
    do
	for i in 0 2  # 0=kPPDATA, 2=kAADATA 6=kAADATAPeri 7=kAADATACentL3
	do
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+('$i','$pt','$rap',0,200,0,0.5,'$muPtCut','$fixParam')'
	done
# Centrlaity dependence 
	for cent in '0,60' '60,200'
	do
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0,0.5,'$muPtCut','$fixParam')'    # PbPb only
	done 
# Event Plane dependence
	for cent in '20,120' 
	do
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0,0.167,'$muPtCut','$fixParam')'    # PbPb only
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0.167,0.333,'$muPtCut','$fixParam')'    # PbPb only
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0.333,0.5,'$muPtCut','$fixParam')'    # PbPb only
	done 
    done
done

