# void onia2ySkim( int nevt = 2000,         int fileID = kPPDATA,
#                  int trigId=kL1DoubleMu0, int epSelection = kEPOppositeHF, 
#                  bool saveTracks=false 

#int kPPDATA = 0 ;
#int kPADATA = 1 ;
#int kAADATA = 2 ; // L1 doubleMu 0
#int kPPMC = 3 ;
#int kPAMC = 4 ;
#int kAAMC = 5 ;
#int kAADATAPeri = 6 ;
#int kAADATACentL3 = 7 ;


#int kNoSelection      =0;
#int kL1DoubleMu0      =1;
#int kL3JpsiCentral    =2 ;
#int kL3UpsilonCentral =3;
#int kL1DoubleMu0Peripheral=4;
#int kL1DoubleMu10     =5;


#int kEPl2HF = 0;
#int kEPOppositeHF = 1;
#int kEPSameSideHF = 2;

root -l -b <<EOF
.L onia2ySkim.C++
.q
EOF

# pp 
./condor_root.sh . 'onia2ySkim.C+(-1,0,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,0,0)'

# PbPb L1_DoubleMu 0
./condor_root.sh . 'onia2ySkim.C+(-1,2,1,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,2,2,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,2,3,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,2,4,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,2,0,1)'

# PbPb Peri and Cent triggers..
#./condor_root.sh . 'onia2ySkim.C+(-1,6,1,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,6,4,1)'
#./condor_root.sh . 'onia2ySkim.C+(-1,7,1,1)'
./condor_root.sh . 'onia2ySkim.C+(-1,7,3,1)'

