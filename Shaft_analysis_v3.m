%% MEE 342 - Shaft Analysis v3     ~ Eduardo Alvarez

close all ; clear ; clc ;

%% Inputs 

ns = 2 ;% static F.S. 
nf = 2 ;% fatigue F.S.
Lt_i = input('Enter the total length of the input shaft (must be greater than 7.5 inches). ') ;
dy = -260 ; % lb 
dx = 50.25 ; % lb
Td = 463.2 ; % in-lb ; 38.6 ft-lb

%% Shaft Lengths 

La = (Lt_i - 4) * (.25/11.5) ; % Gap to First Bearing
Lb = (Lt_i -4) * (.5) ; % Length to Gear Loads
Lc = (Lt_i - 4) - La ; 

%% Solving Rxns

Rb = input('Enter the radius of the gear (inches). ') ;
cx = - dx ; 
Tb = - Td ; 
bz = Tb/Rb ; 
B = bz / cos(20*pi/180) ;
by = B * sin(20*pi/180) ;
cy = ( (dy*(Lt_i-La))-(by*(Lb-La)) ) / (Lc) ; 
ay = dy - by - cy ;
cz = (-bz*(Lb-La)) / (Lc) ; 
az = -bz - cz ;

%% Moment Equations

% My Diagran
x1 = La:Lt_i/10000:Lb ;
x2 = Lb:Lt_i/10000:Lc ;
Myb = ay*(Lb-La) ;
y1 = (((Myb)./(Lb - La)).*x1) - (((Myb)./(Lb - La)).*La) ; %piece 1
y2 = -(((Myb)./(Lc - Lb)).*x2) + (((Myb)./(Lc - Lb)).*Lc) ; %piece 2 
figure(1) ; 
plot(x1,y1,'k',x2,y2,'k') ; 
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (My)') ;

% Mz Diagram
Mzb = -az*(Lb-La) ; 
Mzc = dy*(4+La) ;
y3 = (((Mzb)./(Lb - La)).*x1) - (((Mzb)./(Lb - La)).*La) ; %piece 1
y4 = (((Mzc-Mzb)/(Lc-Lb)).*x2) + Mzb - (((Mzc-Mzb)/(Lc-Lb))*Lb) ; %piece 2
x3 = Lc:Lt_i/1000:Lt_i ;
y5 = ((-Mzc/(Lt_i-Lc)).*x3) + ((Mzc/(Lt_i-Lc)).* Lt_i) ; %piece 3
figure(2) ; 
plot(x1,y3,'k',x2,y4,'k') ;
hold on ;
plot(x3,y5,'k') ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (Mz)') ;

% Mr Diagram

yr1 = sqrt((y1).^2 + (y3).^2) ;
yr2 = sqrt((y2).^2 + (y4).^2) ;
yr3 = y5 ;
figure(3) ;
plot(x1,yr1,'k',x2,yr2,'k') ; 
hold on ;
plot(x3,yr3,'k') ;
hold on ;
plot([x2(length(x2)) x2(length(x2))], [yr2(length(yr2)) yr3(1)],'k')
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Resultant Moment Diagram (Mr)') ;

% Torque Diagram 

Ti = -Tb ; % Internal Torque
figure(4) ;
plot([0 Lb],[0 0], 'k') ;
hold on ;
plot([Lb Lb],[0 Ti], 'k') ;
plot([Lb Lt_i],[Ti Ti],'k') ;
xlabel('X Distance [in]') ; ylabel('Torque [lb-in]') ; title('Torque Diagram (Mx)') ;

%% Choosing Moments

Mr1 = yr1(round(length(yr1)/2)) ; 
Mr2 = yr1(length(yr1)) ; 
Mr3 = yr3(round(length(yr3)/2)) ; 
save('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;


