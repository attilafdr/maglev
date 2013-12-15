clc;

% /www/usr/fisher/helpers/mkfilter -Ch -9.0000000000e-01 -Lp -o 10 -a 7.0000000000e-02 0.0000000000e+00 -Z 1.0000000000e-01
% /www/usr/fisher/helpers/mkfilter -Ch -5.0000000000e-01 -Hp -o 10 -a 1.0000000000e-01 0.0000000000e+00


HPZ = [1 1 1 1 1 1 1 1 1 1]; 

HPP = [0.7972873800 + 1i * -0.5764317455;
	   0.7319194498 + 1i * -0.5977786312;
	   0.5925428607 + 1i * -0.6488758499;
	   0.2962312300 + 1i * -0.6707891898;
	  -0.1783936189 + 1i * -0.3737007802;
	  -0.1783936189 + 1i *  0.3737007802;
	   0.2962312300 + 1i *  0.6707891898;
	   0.5925428607 + 1i *  0.6488758499;
	   0.7319194498 + 1i *  0.5977786312;
	   0.7972873800 + 1i *  0.5764317455]';
   
 LPZ = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1] %0.8090169944+1i*0.587785252 0.8090169944+1i*-0.5877852523];
 
LPP = [0.8961729025+1i*0.4210355198;
	   0.8957002115+1i*0.3760297025;
	   0.9076061028+1i*0.2978263464;
	   0.9234138270+1i*0.1916570037;
	   0.9342053995+1i*0.0662045396;
	   0.9342053995+1i*-0.0662045396;
	   0.9234138270+1i*-0.1916570037;
	   0.9076061028+1i*-0.2978263464;
	   0.8957002115+1i*-0.3760297025;
	   0.8961729025+1i*-0.4210355198;
       %0;
       %0
       ]';
   
    
     HPfilter = zpk(HPZ,HPP,1/1.668454352e+01,10000)
     LPfilter = zpk(LPZ,LPP,1/3.785114626e+08,10000)

    [HPnum, HPden] = tfdata(HPfilter)
    [LPnum, LPden] = tfdata(LPfilter)







