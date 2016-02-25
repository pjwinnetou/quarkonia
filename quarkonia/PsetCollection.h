#ifndef PsetCollection_h
#define PsetCollection_h

#include "cutsAndBin.h"

PSet3SingleCB getUpsilonPsets( int collId = kPPDATA,
			       float ptLow=5, float ptHigh=100,
			       float yLow=1.2, float yHigh=2.4,
			       float muPtCut=4
			       )
{ 
  PSet3SingleCB ret;
  cout << " pt = " << ptLow<<" - " << ptHigh << endl;
  cout << " y = " << yLow<<" - " << yHigh << endl;
  cout << " Single muon pt > " << muPtCut << endl;

  ret.setNAlphaSigma(  -1,-1,-1,-1,-1,-1,-1,-1,-1 ); 
  
  if (muPtCut == (float)4 )  {  // set for 3rd (a0,a1,a2) background function
   
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
  }
  
  else if (muPtCut == (float)3.5 )  {
    
    if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )
      {ret.setNAlphaSigma( 7.70205, 4.37304, 4.37304,
			   1.30535, 7.86261, 7.86261,
			     0.0739032, 0.0809132, 0.0809132 );}
     else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 9.85872, 8.94052, 8.94052,
			   7.54106, 4.81017, 4.81017,
			   0.0945045, 0.0985856, 0.0985856 );}
     else if ( ( ptLow == (float)0 ) && (ptHigh == (float)100 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 8.15202, 8.38118, 8.38118,
			   1.38085, 8.56071, 8.56071,
			   0.126096, 0.14, 0.14 );}
     else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )
      {ret.setNAlphaSigma( 9.9677, 5.43065, 5.43065,
			   1.25754, 9.17669, 9.17669,
			   0.0717877, 0.0759169, 0.0759169 );}
     else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 6.07335, 2, 2,
			   1.2, 1.2, 1.2,
			   0.0871961, 0.0989585, 0.0989585 );}
     else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 5.27611, 9.99995, 9.99995,
			   8.83383, 4.49077, 4.49077,
			   0.130046, 0.14, 0.14 );}
     else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )
      {ret.setNAlphaSigma( 2.76482, 2, 2,
			   1.43454, 1.30664, 1.30664,
			   0.075975, 0.0820177, 0.0820177 );}
     else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)0 ) && (yHigh == (float)2.4 ) )
      {ret.setNAlphaSigma( 9.99996, 2.96964, 2.96964,
			   1.20002, 9.99432, 9.99432,
			   0.0902858, 0.105487, 0.105487 );}
     else if ( ( ptLow == (float)5 ) && (ptHigh == (float)100 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )
       {ret.setNAlphaSigma( 8.8612, 4.48035, 4.48035,
			    9.99973, 2.94075, 2.94075,
			    0.12575, 0.136888, 0.136888 );}
  }
  /*  Legacy parameters for pt1 >4GeV, pt2 >3.5GeV
      if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)0 ) && (yHigh == (float)1.2 ) )
      {
      ret.setNAlphaSigma(   4.99,  2.7,  0.5, 
      1.361, 5.5, 8.7,       
      0.0720, 0.0776, 0.0800); }
      else if ( (ptLow== (float)5 ) && (ptHigh== (float)100 ) && (yLow== (float)0 ) && (yHigh== (float)1.2 ) )
      {
      ret.setNAlphaSigma(   4.71,       0.5,        3.6,
      1.384,   1.625,     7.3,
      0.07443,   0.0824,   0.0835);   }
      
      else if ( ( ptLow == (float)0 ) && (ptHigh == (float)5 ) && (yLow == (float)1.2 ) && (yHigh == (float)2.4 ) )
      {
      ret.setNAlphaSigma(   5, 3.3, 2,
      1.388, 8.2, 5.6,
      0.1272, 0.15, 0.114 ); }
      
      else if ( (ptLow== (float)5 ) && (ptHigh== (float)100 ) && (yLow== (float)1.2 ) && (yHigh== (float)2.4 ) )
      {
      ret.setNAlphaSigma(   4.9999, 4.99, 1.5,
      1.56, 5.0, 1.35,
      0.1229, 0.1406, 0.1130);  }
			      }
  */
  
  
  return ret;
}

#endif
