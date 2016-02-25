'''
int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;



int kNoTrigSel       =0;
int kL1DoubleMu0      =1;
int kL3JpsiCentral    =2;
int kL3UpsilonCentral =3;
int kL1DoubleMu0Peripheral=4;
int kL1DoubleMu10     =5;


int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;
'''

gitVer=$(git show | head -1 | awk '{print $2}')
#void onia2ySkim( int nevt = 2000, int fileID = kPPDATA, int trigId=kL1DoubleMu0, int epSelection = kEPOppositeHF, bool saveTracks=false, TString skimVersion="unIdentified"  ) {
root -l -q -b 'onia2ySkim.C++(-1, 6, 1,  1,0,"'$gitVer'")'  # Peripheral PD + L1_DoubleM0 HLT
root -l -q -b 'onia2ySkim.C+(-1, 6, 4,  1,0,"'$gitVer'")'   # Peripheral PD + peripheral HLT 
root -l -q -b 'onia2ySkim.C+(-1, 0, 1,  1,0,"'$gitVer'")'   # pp data
root -l -q -b 'onia2ySkim.C+(-1, 2, 1,  1,0,"'$gitVer'")'   # AA data
root -l -q -b 'onia2ySkim.C+(-1, 3, 1,  1,0,"'$gitVer'")'   # pp mc
root -l -q -b 'onia2ySkim.C+(-1, 5, 1,  1,0,"'$gitVer'")'   # pbpb mc



