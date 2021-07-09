close all;
echo off;
warning ('off')

Sf = [2.5, 2.8]; 
Su = [-1, -1];
H = [22, 30]; % *10.000
Tm = [2.5, 2.5]; % *100

PinA = [Sf(1) Su(1) H(1) Tm(1)] %makes 1x4 double with starting suggestions for PA
PinB = [Sf(2) Su(2) H(2) Tm(2)] %makes 1x4 double with starting suggestions for PB
VTdata; %imports VT data in a ..x3 double (=ans)

Temp=A(:,1);
shiftPA=A(:,2);
TempB=B(:,1);
shiftPB=B(:,2);

Fxn_name = 'fsigmaT_free'

lbA = [-10,-10,-inf,0]; %lower bounds for ALL variables, [Sf,Su,H,Tm] for PA
ubA = [10,10,inf,6]; %upper bounds for ALL variables, [Sf,Su,H,Tm] for PA
%lbB = [0,0,20000,300]; %lower bounds for ALL variables, [Sf,Su,H,Tm] for PB
%ubB = [10,10,30000,360]; %upper bounds for ALL variables, [Sf,Su,H,Tm] for PB

Iterations = 1000;
OPTIONS=optimset('Display','off','MaxFunEvals',100000,'TolFun',1e-1,'MaxIter',10000,'Algorithm','levenberg-marquardt','DiffMaxChange',0.01);

% PEPTIDE A (HALOGEN BONDING)
for j=1:Iterations

    PinA = lsqcurvefit(Fxn_name,PinA,Temp,shiftPA,lbA,ubA,OPTIONS)
    
    
end

YA = PinA
shiftcalcA = feval(Fxn_name,YA,Temp)
error = shiftPA - shiftcalcA

figure(1), clf
hold on
line (Temp,shiftcalcA,'Color','blue')       %calc fxn as line
plot (Temp,shiftPA,'.k')    %real data points % k=black
title(['Pep A: Sf= ' num2str(YA(1)) ' Su= ' num2str(YA(2)) ' Hm= ' num2str(YA(3)*1000) ' Tm= ' num2str(YA(4)*100)])
xlabel('°T(K)'), ylabel('shift (ppm)');
hold off

%plot(Temp,error,'o'), grid
%xlabel('°T(K)'), ylabel('d-sigma');
%title('Error Peptide A');

% PEPTIDE B (REFERENCE)
Temp = TempB;

for j=1:Iterations

    PinB = lsqcurvefit(Fxn_name,PinB,Temp,shiftPB,lbA,ubA,OPTIONS)
    
end

YB = PinB
shiftcalcB = feval(Fxn_name,YB,Temp)
error = shiftPB - shiftcalcB

figure(2), clf
hold on
line (Temp,shiftcalcB,'Color','red')       %calc fxn as line
plot (Temp,shiftPB,'.k')    %real data points 
title(['Pep B: Sf= ' num2str(YB(1)) ' Su= ' num2str(YB(2)) ' Hm= ' num2str(YB(3)*1000) ' Tm= ' num2str(YB(4)*100)])
xlabel('°T(K)'), ylabel('shift (ppm)');
hold off

%plot(Temp,error,'o'), grid
%xlabel('°T(K)'), ylabel('d-sigma');
%title('Error Peptide B');

FullTemp = linspace(101,600,200);
YAstretch = Stefan_stretch(YA,FullTemp);
YA = YB;
YBstretch = Stefan_stretch(YA,FullTemp);


figure(3), clf
hold on
line (FullTemp,YAstretch,'Color','blue')      
line (FullTemp,YBstretch,'Color','red')    
title('Overlap: Peptide A: blue; Peptide B: red')
xlabel('°T(K)'), ylabel('shift (ppm)');
hold off

tempS_K;