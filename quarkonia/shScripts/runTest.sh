# int collId = kPPDATA,                   float ptLow=5, float ptHigh=100,
# float yLow=0, float yHigh=1.2,                   int cLow=0, int cHigh=20,
# float dphiEp2Low=0, float dphiEp2High=0.5,   // In unit of PI!!
# float muPtCut=4.0,                   bool fixParameters=false


root -l -b <<EOF
.L doFitUpsilonOneCB.C++
.q
EOF


muPtCut='4'
fixParam='1'

outputDir='fitResults/fixParam'${fixParam}'MuPt'${muPtCut}'GeV_2'
mkdir -p $outputDir
cp runUpsilon.sh onia2ySkim.C doFitUpsilonOneCB.C PsetCollection.h cutsAndBin.h $outputDir


echo entering to the loop...


for pt in '0,5' '5,100' 
do
    for rap in '0,1.2'
    do
	for cent in '0,200'
	do
	    ./condor_root.sh $outputDir 'doFitUpsilonOneCB.C+(2,'$pt','$rap','$cent',0,0.5,'$muPtCut','$fixParam')'
	done 
    done
done

