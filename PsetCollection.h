#ifndef PsetCollection_h
#define PsetCollection_h

#include "cutsAndBin.h"

PSet3SingleCB getUpsilonPsets( int collId = kPPDATA,
			       float ptLow=5, float ptHigh=100,
			       float yLow=1.2, float yHigh=2.4,
             int cLow=0, int cHigh=160,
			       float muPtCut=4
			       )
{ 
  PSet3SingleCB ret;
  cout << " pt = " << ptLow<<" - " << ptHigh << endl;
  cout << " y = " << yLow<<" - " << yHigh << endl;
  cout << " Single muon pt > " << muPtCut << endl;

  ret.setNAlphaSigma(  -1,-1,-1,-1,-1,-1,-1,-1,-1 ); // set n, alpha, sigma for each 1S 2S 3S -> 3*3 = 9 parameter set
  ret.setNAlphaSigmaCB2(  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 ); // set two (n, alpha, sigma) for each 1S 2S 3S -> 2*3*3 = 18 parameter set

  ret.setSig1sF21NBkg(0,0,0);  // same as setSig1sF321NBkg excluding 3S/1S ratio
  ret.setSig1sF321NBkg(0,0,0,0);  // both for Toy MC (nSig1S, f21, f31, nBkg) --> lead to nSig1S, nSig2S, nSig3S, nBkg for Toy MC generation

  ret.setParMC( -1, -1, -1, -1, -1, -1); // set parameter from MC : n, alpha, sigma, m0, f, x
  ret.setParBkg( -1, -1, -1); // set bkg parameter : mu, sigma, lambda
  ret.setParBkgRes( -1, -1, -1, -1); // set bkg parameter from nominal result : mu, sigma, lambda, m0

  if (muPtCut == (float)4 )  
  {
    //***************************************************************
    //***PP & (PbPb : 0-100% + integrated // pT // rap binning)******
    //***************************************************************
 
    //Integrated Bin // No centrlaity dependence 
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {
      //Parameter from Chad // Apr 18th
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {
        ret.setParMC( 3.65, 1.52, 0.0675,    
          9.4577, 0.556, 1.94);
        ret.setParBkg(8.565,1.055,8.22);
      }
      else 
      {
        ret.setParMC( 3.74, 1.48, 0.0693,
		       9.4577, 0.573, 1.90 );
        ret.setParBkg(7.946,1.103,6.06);
      }

      ret.setNAlphaSigmaCB2( 8.83664, 6.98899, 8.13599, 1.10179, 7.00382, 8.75107, 
          1.21009, 1.12278, 1.15277, 1.13759, 7.54582, 7.52483, 
          0.103963, 0.0500101, 0.121012, 0.0720748, 0.139509, 0.0763136 );
      ret.setNAlphaSigma( 9.61274, 1.1, 9.15923,
          1.13981, 1.15, 7.83501,
          0.0843173, 0.0952217, 0.0934913 );
    } 
    
    //Full Rap & pT bin
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {
      //Parameter from Chad // Apr 18th
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {
        ret.setParMC( 3.56, 1.55,  0.0654,
          9.4577,  0.536, 1.99 );
        ret.setParBkg(8.98,0.97,6.29);
      }
      else
      {ret.setParMC(  2.48, 1.60, 0.0683, 
          9.4577, 0.573, 1.91 );
        ret.setParBkg(8.50,0.798,4.49);
      }
    
   //   ret.setParMC_PP( 2.95734, 1.58758, 0.0716493, 
   //     9.4577, 0.671119, 2.12786 );
    } 
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {
      //Parameter from Chad
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {ret.setParMC( 3.76, 1.49, 0.0690, 
          9.4577, 0.582, 1.93 );
        ret.setParBkg(0,0,21.3);
      }
      else
      {ret.setParMC( 6.67, 1.35, 0.0677,
          9.4577, 0.563, 1.92);
        ret.setParBkg(0,0,8.52);
      }
      
      //ret.setParMC_PP( 2.12817, 1.64774, 0.0735665, 
      //  9.45566, 0.678357, 2.03577 );
    } 
    else if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {
      //Parameter from Chad April 18th
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {ret.setParMC( 3.39, 1.55, 0.0694,
          9.4577, 0.556, 1.85 );
        ret.setParBkg(0,0,12.7);
      }
      else
      {ret.setParMC( 3.61, 1.52, 0.0750, 
          9.4577, 0.657, 1.905);
        ret.setParBkg(0,0,18.9);
      }
      
    }

    //Full pT & y bin 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )
    {
      //Parameter from Chad April 18th
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {ret.setParMC( 1.72, 1.82, 0.0612,
          9.4577, 0.641, 1.72) ;
        ret.setParBkg(8.54,1.48,11.2);
      }
      else
      {ret.setParMC( 1.49, 1.85, 0.0591,
          9.4577, 0.548, 1.69);
        ret.setParBkg(7.82,1.05,8.23);
      }
    } 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )
    {
      //Parameter from Chad April 18th
      if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)
      {ret.setParMC( 1.89, 1.91, 0.105, 
          9.4577, 0.563, 1.55) ;
        ret.setParBkg(8.64,1.10,4.50);
      }
      else
      {ret.setParMC( 2.11, 1.905, 0.110,
          9.4577, 0.618, 1.55);
        ret.setParBkg(8.18,1.20,2.97);
      }

    } 
    
 

   //old bin
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)2.5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 3.15106, 1.57896, 0.146825, 
        9.45898, 0.354881, 0.4833 );} 
    else if ( ( ptLow == (float)2.5 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 2.14371, 1.69434, 0.0722432, 
        9.45696, 0.654263, 2.04668 );} 
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)8 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 1.86268, 1.70799, 0.144922, 
        9.45598, 0.355838, 0.492459 );} 
    else if ( ( ptLow == (float)8 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 1.88988, 1.68954, 0.0718718, 
        9.45565, 0.636682, 1.98271 );} 
/*    else if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 2.0111, 1.66968, 0.0769498, 
        9.45554, 0.744286, 2.078 );} */
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)0.4 ) )
    {ret.setParMC( 1.29673, 1.88882, 0.0566235, 
        9.45916, 0.934458, 2.20436 );} 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0.4 ) && (yHigh == (float)0.8 ) )
    {ret.setParMC( 1.47658, 1.8233, 0.0685656, 
        9.45896, 0.902513, 1.94238 );} 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0.8 ) && (yHigh == (float)1.2 ) )
    {ret.setParMC( 1.31014, 1.98302, 0.147403, 
        9.45604, 0.14006, 0.600046 );} 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)1.6 ) )
    {ret.setParMC( 1.0219, 2.13453, 0.153924, 
        9.45063, 0.277966, 0.665491 );} 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.6 ) && (yHigh == (float)2 ) )
    {ret.setParMC( 1.29582, 2.00072, 0.178302, 
        9.44855, 0.340889, 0.652577 );} 
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)2 ) && (yHigh == (float)2.4 ) )
    {ret.setParMC( 1, 2.20011, 0.22998, 
        9.44562, 0.269098, 0.651062 );} 

    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
    {ret.setNAlphaSigma( 5.15366, 1.10287, 1.10287, 
        1.22781, 1.15052, 1.15052, 
        0.0850295, 0.0932266, 0.0932266 );} 





    /*
  ret.setNAlphaSigma(  -1,-1,-1,-1,-1,-1,-1,-1,-1 ); 
  
  if (muPtCut == (float)4 )  {  // set for 3rd (a0,a1,a2) background function
 
if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
 {ret.setNAlphaSigma( 5.15366, 1.10287, 1.10287, 
1.22781, 1.15052, 1.15052, 
0.0850295, 0.0932266, 0.0932266 );} 

    // Feb. 22nd .  pt in  '0,100'  '0,5' '5,100' '0,7' '7,100'    rap in '0,2.4' '0,1.6' '1.6,2.4'
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.6 ) )
      {ret.setNAlphaSigma( 10, 9.88409, 9.88409,
			   1.30449, 9.84178, 9.84178,
			   0.0800495, 0.0901807, 0.0901807 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)1.6 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 9.99951, 2.47988, 2.47988,
			   1.53688, 7.53012, 7.53012,
			   0.14, 0.14, 0.14 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)1.6 ) )
	    {ret.setNAlphaSigma( 9.99997, 9.998, 9.998,
				 1.29165, 7.59924, 7.59924,
				 0.0782797, 0.0895556, 0.0895556 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 9.99998, 4.83963, 4.83963,
			     1.2, 7.27032, 7.27032,
			     0.0875047, 0.10848, 0.10848 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)1.6 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 9.99887, 2.00035, 2.00035,
			     1.54407, 4.86455, 4.86455,
			     0.14, 0.14, 0.14 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)7 ) && (yLow == (float)0 ) && (yHigh == (float)1.6 ) )
	{ret.setNAlphaSigma( 2.00006, 9.99263, 9.99263,
			     1.20039, 4.08075, 4.08075,
			     0.139996, 0.117633, 0.117633 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)7 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 2.00026, 3.20551, 3.20551,
			     1.56116, 9.70192, 9.70192,
			     0.0904289, 0.108403, 0.108403 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)7 ) && (yLow == (float)1.6 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 2.12861, 2.0002, 2.0002,
			     1.8452, 9.99156, 9.99156,
			     0.14, 0.14, 0.14 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.6 ) )
	{ret.setNAlphaSigma( 7.27208, 9.99999, 9.99999,
			     1.34492, 6.65137, 6.65137,
			     0.0823384, 0.0914062, 0.0914062 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 2.00215, 10, 10,
			   1.60226, 6.00508, 6.00508,
			   0.0923425, 0.106126, 0.106126 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)1.6 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 2.02033, 2.53265, 2.53265,
			     2.04584, 9.02749, 9.02749,
			     0.14, 0.14, 0.14 );}
    else if ( ( ptLow == (float)7 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.6 ) )
	{ret.setNAlphaSigma( 2.00008, 2.00133, 2.00133,
			     8.26198, 1.20568, 1.20568,
			     0.0850329, 0.0845603, 0.0845603 );}
    else if ( ( ptLow == (float)7 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 2.00025, 10, 10,
			     1.65326, 5.91693, 5.91693,
			     0.0929326, 0.106411, 0.106411 );}
    else if ( ( ptLow == (float)7 ) && (ptHigh == (float)100 ) && (yLow == (float)1.6 ) && (yHigh == (float)2.4 ) )
	{ret.setNAlphaSigma( 2.00001, 2.00109, 2.00109,
			     1.93639, 8.17018, 8.17018,
			     0.14, 0.14, 0.14 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 2.08724, 9.94783, 9.94783,
			   1.58759, 8.98771, 8.98771,
			   0.0911935, 0.107455, 0.107455 );}

    // (Feb. 18th)
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)7 ) && (yLow == (float)0 ) && (yHigh == (float)0.4 ) )
      {ret.setNAlphaSigma( 2.00548, 2.00002, 2.00002,
			   1.80805, 1.91796, 1.91796,
			   0.058962, 0.0684566, 0.0684566 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)7 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 2.23089, 9.32961, 9.32961,
			   1.48476, 6.54795, 6.54795,
			  0.0906494, 0.103696, 0.103696 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)7 ) && (yLow == (float)0.4 ) && (yHigh == (float)0.8 ) )
      {ret.setNAlphaSigma( 2.00614, 2, 2,
			   1.49066, 1.63401, 1.63401,
			  0.0694319, 0.0693002, 0.0693002 );}
    // Rapidity < 1.93  ( Feb. 16th)
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.93 ) )
      {ret.setNAlphaSigma( 9.99977, 9.67821, 9.67821,
			   1.28714, 9.32963, 9.32963,
			   0.0849484, 0.0951778, 0.0951778 );}
    else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.93 ) )
      {ret.setNAlphaSigma( 8.86768, 9.95155, 9.95155,
			   1.30559, 9.82849, 9.82849,
			   0.0866596, 0.0950049, 0.0950049 );}
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)1.93 ) )
      {ret.setNAlphaSigma( 9.99983, 9.32882, 9.32882,
			   1.25728, 6.58999, 6.58999,
			   0.0838298, 0.101416, 0.101416 );}
	
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )  
      {
	ret.setNAlphaSigma(   8.1, 2.41, 4.5,
			      1.288, 1.21, 7.2,
			      0.0703, 0.078, 0.0849); }
    
    else if ( (ptLow== (float)5 ) && (ptHigh== (float)100 ) && (yLow== (float)0 ) && (yHigh== (float)1.2 ) )
      {
	ret.setNAlphaSigma(   5.3, 7.8, 2,
			      1.469, 10, 1.2,
			      0.0759, 0.078, 0.0772); }
    
    else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )  
      {
	ret.setNAlphaSigma(   5.3, 9.5, 2,
			      1.484, 8.5, 1.2,
			      0.126, 0.14, 0.131);}
    
    else if ( (ptLow== (float)5 ) && (ptHigh== (float)100 ) && (yLow== (float)1.2 ) && (yHigh== (float)2.4 ) )
      {
	ret.setNAlphaSigma( 10, 9.5, 2.0, 
			    1.367, 6.7, 1.24, 
			    0.1228, 0.14, 0.1213);}
  */
  
    // background systematics 
    if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   { 

      
      // from real data, May 2
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkgPol3( 0.153043, -0.238099, 0.0814439);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( 0.0470781, -0.253079, 0.10792);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( -0.118002, -0.292458, 0.139035);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( 0.200494, -0.433492, 0.157105);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( -0.342874, 0.0490931, 0.04333);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( -0.162706, 0.0252931, 0.0212731);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkgPol4( 0.153406, -0.236557, 0.0814901, 0.00335914);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol4( 0.0449198, -0.261034, 0.109096, -0.0170839);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol4( -0.124454, -0.313708, 0.149625, -0.0434754);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol4( 0.199545, -0.434946, 0.156449, -0.00440958);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol4( -0.345239, 0.0271187, 0.0568589, -0.0341808);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol4( -0.164412, 0.00479173, 0.0299365, -0.0342342);


      //result from the nominal fit

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkgRes(8.49532,1.17973,11.8666,9.45209 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgRes( 8.46524,1.055,8.22475,9.45023 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgRes( 8.58147,1.00163,4.88473,9.4415);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgRes( 6.64382,0.8062,6.64382,9.45084);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgRes( 0, 0,8.03611,9.4493);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgRes( 0, 0,17.0224,9.44918);


      /* SS MC
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkgPol3( 0.15057, -0.221047, 0.0584276);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( 0.0486762, -0.249424, 0.0923441);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( -0.150275, -0.305512, 0.152585);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( 0.172673, -0.425264, 0.146032);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgPol3( -0.214283, 0.0331719, 0.0316277);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
      ret.setParBkgPol3( -0.137293, 0.0259804, 0.00347688);    
      */
      
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkg( 8.53959, 1.47378, 11.2341);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg( 8.5165, 1.27541, 7.95591);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg( 8.64007, 1.09979, 4.49538);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg( 8.97432, 0.969171, 6.2889);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg( 0, 0, 12.7024);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg( 0, 0, 21.3048);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkg2ErrExp( 9.09684, 1.20723, 4.55871,
			      7.93726, 2.66578, 45.6162, 0.545745);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg2ErrExp( 7.92787, 1.1246, 6.32217,
			      8.97375, 1.13974, 9.64848, 0.55323);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg2ErrExp( 8.81906, 1.27244, 4.01575,
			      8.43383, 0.718196, 5.47444, 0.276265);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg2ErrExp( 5.59077, 2.32819, 1.68291,
			      8.97434, 0.969083, 6.28955, 1);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg2ErrExp( 0, 0, 12.8255,
			      0, 0, 12.6111, 0.587102);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkg2ErrExp( 0, 0, 49.9993,
			      0, 0, 2.84152, 0.0876578);

    }
    
    else {  // PbPb
      /*if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 7.8194, 1.04522, 8.22826);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkg( 7.85668, 0.927897, 6.22555);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 7.85887, 1.0159, 6.08148);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkg( 8.13523, 1.10166, 7.76745);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkg( 7.7818, 1.13994, 5.70851);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkg( 9.67679, 1.95593, 2.99213);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)160 )  )
	ret.setParBkg( 8.88746, 1.01721, 4.65605);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkg( 7.84871, 1.00154, 5.82828);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkg( 7.8803, 1.21752, 6.11419);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkg( 7.96869, 1.03374, 6.08945);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
	ret.setParBkg( 8.1093, 0.496274, 8.58243);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 8.1816, 1.19908, 2.97442);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 8.49738, 0.797861, 4.4877);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 0, 0, 18.8578);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkg( 0, 0, 8.51591);
*/
      // from real data fit, May 2
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.139622, -0.140123, 0.0800708);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkgPol3( -0.236093, -0.120547, 0.097607);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.237187, -0.138513, 0.0934057);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkgPol3( -0.0987132, -0.267165, 0.137858);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkgPol3( -0.277973, -0.117348, 0.0836331);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkgPol3( 0.149861, -0.24264, 0.194676);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( 0.247142, -0.408111, 0.292681);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkgPol3( -0.263367, -0.139223, 0.0890883);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkgPol3( -0.228942, -0.135845, 0.0780673);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkgPol3( -0.206513, -0.171631, 0.101607);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
	ret.setParBkgPol3( -0.112646, -0.20485, 0.137838);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.502781, -0.137075, 0.12317);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.118782, -0.340394, 0.202157);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.207017, 0.00926542, 0.0235402);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol3( -0.342847, 0.020315, 0.00752558);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.140653, -0.150509, 0.0833484, -0.0187622);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkgPol4( -0.237267, -0.135218, 0.104524, -0.0255536);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.238673, -0.15568, 0.101417, -0.030147);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkgPol4( -0.104161, -0.304131, 0.154502, -0.0722361);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkgPol4( -0.278845, -0.128008, 0.0887048, -0.0184208);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkgPol4( 0.138267, -0.277842, 0.197145, -0.065954);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( 0.265569, -0.383514, 0.299676, 0.0636796);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkgPol4( -0.265001, -0.165885, 0.102533, -0.046685);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkgPol4( -0.230305, -0.157211, 0.0876892, -0.0370089);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkgPol4( -0.207043, -0.176545, 0.103677, -0.00875852);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
	ret.setParBkgPol4( -0.111812, -0.201952, 0.137158, 0.00503467);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.504231, -0.172404, 0.15257, -0.0569553);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.12644, -0.36678, 0.212134, -0.0561921);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.209296, -0.0145804, 0.0319277, -0.0398372);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgPol4( -0.342828, 0.017376, 0.00898467, -0.0047585);
     

      //result from the nominal 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 7.88537, 1.24244, 7.74079,9.4523);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkgRes( 7.8167, 1.103, 6.46753,9.4579);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 7.94553,1.103,6.05766,9.4503);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkgRes( 8.30924,1.103,6.31127,9.4580);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkgRes( 7.85177,1.103,5.83277,9.4533);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkgRes( 8.36736,1.103,9.03956,9.4471);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 9.04258,1.103,5.3108,9.457);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkgRes( 7.95566,1.103,5.728,9.4474);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkgRes( 7.89851,1.103,6.39386,9.4498);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkgRes( 8.12265,1.103,5.71561,9.4492);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
	ret.setParBkgRes( 8.17478,1.103,7.14846,9.4419);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 8.21679,1.01966,3.36472,9.4412);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 8.59429,0.819975,4.68143,9.4510);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 0,0,13.743,9.4550);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgRes( 0,0,8.54481,9.4473);

      /* MC SS 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.15974, -0.124924, 0.0795904);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkgPol3( -0.278328, -0.111362, 0.0817727);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.270729, -0.117976, 0.0905687);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkgPol3( -0.0698811, -0.178187, 0.1256);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkgPol3( -0.295367, -0.102098, 0.0957806);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkgPol3( -0.196448, -0.256386, 0.0844397);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.0207115, -0.408187, 0.152843);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkgPol3( -0.298478, -0.114252, 0.0848217);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkgPol3( -0.229667, -0.128292, 0.10228);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkgPol3( -0.245873, -0.149978, 0.0864309);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
      ret.setParBkgPol3( -0.139084, -0.143909, 0.116901);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.59852, -0.0971782, 0.119311);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.196155, -0.309863, 0.202749);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.181371, -0.0298397, -0.0346417);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
	ret.setParBkgPol3( -0.337143, 0.0230273, 0.0158519);
      */

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setParBkg2ErrExp( 8.11858, 0.632277, 46.608,
			    9.99696, 2.58143, 2.26691, 0.445348);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
      ret.setParBkg2ErrExp( 7.81141, 2.13067, 6.29515,
			    8.65364, 0.504939, 1.3533, 0.120596);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setParBkg2ErrExp( 7.52226, 2.99916, 7.19687,
			    8.09777, 0.842954, 4.61882, 0.612696);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
      ret.setParBkg2ErrExp( 8.30565, 0.686399, 1.51663,
			    9.51109, 0.562063, 47.138, 0.589143);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
      ret.setParBkg2ErrExp( 9.45175, 0.764216, 7.97158,
			    8.38944, 0.710573, 1.00212, 0.354892);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
      ret.setParBkg2ErrExp( 9.57507, 1.49902, 3.54996,
			    9.68228, 2.4717, 1.19673, 0.140101);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
      ret.setParBkg2ErrExp( 8.91417, 0.982643, 4.58231,
			    8.7967, 1.19317, 4.9076, 0.219837);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
      ret.setParBkg2ErrExp( 8.0072, 0.790722, 4.03563,
			    9.88213, 2.80972, 4.46411, 0.394269);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
      ret.setParBkg2ErrExp( 9.99952, 0.901286, 30.5598,
			    8.62586, 1.06864, 1.62388, 0.578181);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
      ret.setParBkg2ErrExp( 9.37942, 2.90797, 3.83128,
			    8.33842, 0.400136, 4.51689, 0.291002);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
      ret.setParBkg2ErrExp( 8.61136, 0.400398, 1.08194,
			    9.99588, 2.40842, 4.93665, 0.796494);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setParBkg2ErrExp( 9.17895, 1.61718, 2.84503,
			    8.04599, 0.612087, 2.15591, 0.405453);

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setParBkg2ErrExp( 8.54773, 0.626949, 4.38699,
			    9.66643, 2.02949, 2.82492, 0.334572);
    
    if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setParBkg2ErrExp( 0, 0, 18.8776,
			    0, 0, 18.8776, 0.260287);
    
    if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setParBkg2ErrExp( 0, 0, 8.48803,
			    0, 0, 8.55283, 0.438691);
    } // PbPb

    //  Set 2nd background systematics // May 2nd
    if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   {  // pp 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setParBkgErrExpExp( 8.71604, 0.990639, 13.0202, 11, 0.182202);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgErrExpExp( 8.49627, 1.03955, 8.17772,18.9281, 0.0376224);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgErrExpExp( 8.67837, 0.906286, 4.84363, 6.25284, 0.110066);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgErrExpExp( 8.91188, 0.801246, 6.65364, 19.9997, 0.01);

      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgErrExpExp( 0, 0, 19.9994, 3.94923, 0.373256);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setParBkgErrExpExp( 0, 0, 21.3048, 10.9989, 0.260729);
    } // end of pp 
    else {
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 8.12285, 1.01632, 5.72878, 8.8278, 0.161383);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
	ret.setParBkgErrExpExp( 8.33566, 0.932331, 7.50459,2.18136, 0.157304);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
	ret.setParBkgErrExpExp( 7.95585, 1.22543, 5.38293,9.95731, 0.111084);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
	ret.setParBkgErrExpExp( 8.49392, 0.648108, 6.24239, 5.63387, 0.435291);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
	ret.setParBkgErrExpExp( 8.25398, 0.855212, 6.83401, 5.86311, 0.274226);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
	ret.setParBkgErrExpExp( 8.64092, 0.840586, 5.79928, 5.96887, 0.394377);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
	ret.setParBkgErrExpExp( 8.5581, 0.957217, 6.94003, 2.17766, 0.106471);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
	ret.setParBkgErrExpExp( 8.49186, 0.842133, 5.73034,10.9843, 0.100001);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
	ret.setParBkgErrExpExp( 11.0307, 2.29225, 2.99213,10.9999, 0.100001);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 8.94408, 0.632479, 6.66266,10.9979, 0.1);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 8.00472, 1.18449, 7.5388,10.7105, 0.102206);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 8.38809, 1.02568, 2.97442, 8.36066, 0.136021);

      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 8.66849, 0.765402, 4.41278,11, 0.1);
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 0, 0, 16.8826, 8.54001, 1);
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
	ret.setParBkgErrExpExp( 0, 0, 19.9817, 9.60377, 0.426122);
    } // pbpb
    
    
    
    //  Set signal numbers // Updated from Chad's result On Aug 1st
    if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   {  // pp 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setSig1sF21NBkg( 20656.7, 0.32375, 77583.6 ) ;
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 14266.4, 0.317267, 43852.3 ) ;
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 16847.1, 0.302235, 78762.4 );
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 13249.9, 0.311204, 36589.6 );
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 4357.89, 0.380202, 7108.71);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg(34453, 0.319887, 122179 );

    } // end of pp 
    else { // pbpb 
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 2634.67, 0.0966339, 49689.6) ; 
    if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 2383.06, 0.0756197, 48647.6 ) ; 
    if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 807.607, 0.142482, 5899.82) ; 

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF321NBkg( 5793,  0.0994625, 0.00943264, 103811 ) ;

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 3632.06,  0.110272 , 76327.3 ) ;
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 2232.06, 0.0629488, 27911.1 );

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
      ret.setSig1sF21NBkg( 1076.58, 0.0335189, 24044 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
      ret.setSig1sF21NBkg( 871.894, 0.120531, 20196.5 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
      ret.setSig1sF21NBkg( 1331.63, 0.0959053, 28706.2 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
      ret.setSig1sF21NBkg( 1009.57, 0.127288, 15984.4 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
      ret.setSig1sF21NBkg( 849.816, 0.123842, 10377.1 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
      ret.setSig1sF21NBkg( 502.413, 0.13744, 4762.33 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
      ret.setSig1sF21NBkg( 293.383, 0.110187,  2235);
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
      ret.setSig1sF21NBkg( 198.671, 0.164113 , 953.692 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 65.8021, 0.416928, 480.779 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)200 )  )
      ret.setSig1sF21NBkg( 102.348, 0.248146, 1019.19);


    } // end of pbpb 
  /*  
    //  Set signal numbers // Updated from Chad's result On Apr 27th
    if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   {  // pp 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 )  )
	ret.setSig1sF21NBkg( 20657, 0.323742, 77582.9 ) ;
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 14265.8, 0.31731, 43854.5 ) ;
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg(16847, 0.30224, 78762.2 );
      if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg(13250, 0.311233, 36589.9 );
      if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg( 4357, 0.380275, 7108.31);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
	ret.setSig1sF21NBkg(34450.7, 0.319878, 122184 );

    } // end of pp 
    else { // pbpb 
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 2611, 0.0998969, 49288) ; 
    if ( ( ptLow == (float)5 ) && (ptHigh == (float)12 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 2379, 0.0746074, 48633 ) ; 
    if ( ( ptLow == (float)12 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 805,  0.142994, 5897.8) ; 

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF321NBkg( 5793,  0.0994625, 0.00943264, 103811 ) ;

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 3611.45,  0.110555 , 76054.1 ) ;
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg(2225.16, 0.0641333, 27773 );

    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)10 )  )
      ret.setSig1sF21NBkg( 1075.78, 0.0322151, 24045.4 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)10 ) && (cHigh == (int)20 )  )
      ret.setSig1sF21NBkg( 893.958, 0.120287, 20177.4 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)40 )  )
      ret.setSig1sF21NBkg( 1308.78, 0.0978241, 28729.4 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)40 ) && (cHigh == (int)60 )  )
      ret.setSig1sF21NBkg( 1001.53, 0.128813, 15991.6 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)80 )  )
      ret.setSig1sF21NBkg( 864.603, 0.126703, 10359.5 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)80 ) && (cHigh == (int)100 )  )
      ret.setSig1sF21NBkg( 507.686, 0.139436, 4755.28 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)120 )  )
      ret.setSig1sF21NBkg( 290.273, 0.109301,  2238.18);
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)120 ) && (cHigh == (int)140 )  )
      ret.setSig1sF21NBkg( 194.736, 0.18231 , 950 );
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)140 ) && (cHigh == (int)160 )  )
      ret.setSig1sF21NBkg( 63.8036, 0.408185, 484.195 );


    } // end of pbpb 
*/
    //  Set signal numbers for 3S upper limit setSig1sF321NBkg --> Nsig1S, Nsig2S (f21), Nsig3S (f31), Nbkg for 0-10%, 10-30%, 30-50%, 50-80% and pp
    /*if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   
    {  // pp 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
        ret.setSig1sF321NBkg(34450.7, 0.319878, 0.168, 122184 );
    } // end of pp 
    else 
    { // pbpb 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)20 )  )
        ret.setSig1sF321NBkg( 1964.61, 0.0724332, 0.0331024, 44229.9 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)60 )  )
        ret.setSig1sF321NBkg( 2309.96, 0.110981, 4.79439*TMath::Power(10,-8), 44721.4 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)100 )  )
        ret.setSig1sF321NBkg( 1370.84, 0.132736, 0.0388728, 15112.6);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)160 )  )
        ret.setSig1sF321NBkg( 540.209, 0.161275, 0.0275389, 3696.81);
    } // end of pbpb 
    */


    //Aug 1st
    if(collId == kPPDATA || collId == kPPMC || collId == kPPMCUps1S || collId == kPPMCUps2S || collId == kPPMCUps3S)   
    {  // pp 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 )  )
        ret.setSig1sF321NBkg(34450.7, 0.319878, 0.168, 122184 );
    } // end of pp 
    else 
    { // pbpb 
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)0 ) && (cHigh == (int)20 )  )
        ret.setSig1sF321NBkg( 1948.46, 0.0736699, 0.0456728, 44236.9 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)20 ) && (cHigh == (int)60 )  )
        ret.setSig1sF321NBkg( 2341.4, 0.109522, 1.35731*TMath::Power(10,-8), 44689.9 );
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)60 ) && (cHigh == (int)100 )  )
        ret.setSig1sF321NBkg( 1351.18, 0.1303, 0.040567, 15136.6);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)160 )  )
        ret.setSig1sF321NBkg( 554.471, 0.161745, 0.0222137, 3682.45);
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) && ( cLow == (int)100 ) && (cHigh == (int)200 )  )
        ret.setSig1sF321NBkg( 590, 0.149433, 0.0262152, 4218.25);
    } // end of pbpb 

  } // 4GeV muon cut
  
  return ret;
}

#endif
