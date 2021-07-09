function [Y]  = fsigmaT_free(PinA,Temp,FullTemp) %values in brackets get imported


Sf = PinA(1);
Su = PinA(2);
H  = PinA(3);
Tm = PinA(4);


Y = Su + (  (Sf - Su)./(1+ exp( ((H*1000)/-8.314)*((1./Temp)-(1/(Tm*100)))) ) );


%Y = Y/Ynorm;

