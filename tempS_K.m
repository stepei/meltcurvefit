% script to plot extended curve points of peptide A and B into a surfboard

SfA=PinA(1);
SuA=PinA(2);

SfB=PinB(1);
SuB=PinB(2);


    for i=1:1:length(FullTemp)
        XdataL(i)=(SuA-YAstretch(i))*(YBstretch(i)-SfB);
        YdataL(i)=(YAstretch(i)-SfA)*(SuB-YBstretch(i)); 
    end
   
  N = min(numel(shiftcalcA,shiftcalcB)); %if datasets A and B are not equal in length
  shiftcalcKA = shiftcalcA(1:N)
  shiftcalcKB = shiftcalcB(1:N)
    
   for i=1:1:N
        XdataS(i)=(SuA-shiftcalcKA(i))*(shiftcalcKB(i)-SfB);
        YdataS(i)=(shiftcalcKA(i)-SfA)*(SuB-shiftcalcKB(i)); 
    end


[Ydata_maxV, Ydata_maxP]=max(YdataL);
figure(101), clf
plot(XdataL, YdataL, '-r', XdataS, YdataS, 'xb') % , XdataL(Ydata_maxP:length(XdataL)), YdataL(Ydata_maxP:length(YdataL)), 'x-r')

line([min(XdataL),max(XdataL)],[min(YdataL),max(YdataL)])
k = (max(YdataL))/(max(XdataL));
title(['y = ' num2str(k) '*x'])
 xlabel('(sUA-sA)*(sB-sFB)'), ylabel('(sA-sFA)*(sUB-sB)')
 
 


