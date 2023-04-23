
%% MEE 342 - Countershaft  

clear ; clc ;

load('InputShaft.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf', 'Lb', 'x1', 'x2', 'x3', 'Lt_i','Tb') ;

%% Inputs

Lt_c = input('Enter the total length of the counter shaft (must be greater than 5 inches). ') ;
fprintf(" \nSecond Gear must be smaller than first gear \nFirst gear was selected to have a radius of %5.3f inches \n",Rb)
Rg = input('Enter the radius of the first gear (inches). ') ;

%% Solving Rxns

Le = (.75/11.5)*Lt_c ;  % Ratios from Example in Book
Lf = (2.75/11.5)*Lt_c ;
Lg = (8.5/11.5)*Lt_c ;
Lh = (10.75/11.5)*Lt_c ;
Rf = Rb ; 

Gz = -Tb/Rg ; 
G = Gz / cos(20*180/pi) ;
Gy = G * sin(20*180/pi) ;
Tg = Gz * Rg ;

Tf = -Tg ;
Fz = Tf/Rf ;
F = Fz / cos(20*180/pi) ;
Fy = F * sin(20*180/pi) ;

Hy = ( (-Fy*(Lf-Le)) - (Gy*(Lg-Le)) ) / (Lh - Le) ; % Moment from E
Hz = ( (-Fz*(Lf-Le)) - (Gz*(Lg-Le)) ) / (Lh - Le) ;

Ey = -Fy - Gy - Hy ;
Ez = -Fz - Gz - Hz ;

%% Shear Diagrams

figure(3) ;
subplot(3,1,1) ;
plot([Le Le],[0 Ey],'k', [Le Lf], [Ey Ey], 'k') ;
hold on ;
plot([Lf Lf],[Ey (Ey+Fy)],'k', [Lf Lg], [(Ey+Fy) (Ey+Fy)], 'k') ;
hold on ;
plot([Lg Lg],[(Ey+Fy) (Ey+Fy+Gy)],'k', [Lg Lh], [(Ey+Fy+Gy) (Ey+Fy+Gy)], 'k') ;
hold on ;
plot([Lh Lh],[(Ey+Fy+Gy) (Ey+Fy+Gy+Hy)],'k') ;
xlabel('X Distance [in]') ; ylabel('Shear Force [Fy]') ; title('Shear Force along the Y-axis') ;
hold off ;

subplot(3,1,2) ;
plot([Le Le],[0 Ez],'k', [Le Lf], [Ez Ez], 'k') ;
hold on ;
plot([Lf Lf],[Ez (Ez+Fz)],'k', [Lf Lg], [(Ez+Fz) (Ez+Fz)], 'k') ;
hold on ;
plot([Lg Lg],[(Ez+Fz) (Ez+Fz+Gz)],'k', [Lg Lh], [(Ez+Fz+Gz) (Ez+Fz+Gz)], 'k') ;
hold on ;
plot([Lh Lh],[(Ez+Fz+Gz) (Ez+Fz+Gz+Hz)],'k') ;
xlabel('X Distance [in]') ; ylabel('Shear Force [Fz]') ; title('Shear Force along Z-axis') ;
hold off ;

subplot(3,1,3) ;
plot([Lf Lg],[-Tg -Tg],'k', [Lf Lf],[0 -Tg],'k') ;
hold on ;
plot([Lg Lg],[0 -Tg],'k') ;
xlabel('X Distance [in]') ; ylabel('Torque [lb-in]') ; title('Torque Diagram (Mx)') ;
hold off

%% Moment Equations 

x4 = Le:Lt_c/10000:Lf ;
x5 = Lf:Lt_c/10000:Lg ;
x6 = Lg:Lt_c/10000:Lh ;
x_tot = 0:Lt_c/10000:Lt_c ;

Myf = -Ez*(Lf-Le) ;
Mzf = Ey*(Lf-Le) ;

Myg = -Hz*(Lh-Lg) ;
Mzg = Hy*(Lh-Lg) ;

y4 = (Myf/(Lf-Le)).*x4 - ( (Myf*Le) / (Lf-Le)) ;
y4b = (Mzf/(Lf-Le)).*x4 - ( (Mzf*Le) / (Lf-Le)) ; % b is z plane
yr4 = sqrt(y4.^2 + y4b.^2) ;

y5 = ((Myg-Myf)/(Lg-Lf)).*x5 + Myf - (((Myg-Myf)/(Lg-Lf))*Lf) ;
y5b = ((Mzg-Mzf)/(Lg-Lf)).*x5 + Mzf - (((Mzg-Mzf)/(Lg-Lf))*Lf) ;
yr5 = sqrt(y5.^2 + y5b.^2) ;

y6 = (-Myg/(Lh-Lg)).*x6 + ((Myg/(Lh-Lg))*Lh) ;
y6b = (-Mzg/(Lh-Lg)).*x6 + ((Mzg/(Lh-Lg))*Lh) ;
yr6 = sqrt(y6.^2 + y6b.^2) ;

figure(4) ; 
subplot(3,1,1) ;
plot(x4,y4,'k',x5,y5,'k') ;
hold on ;
plot(x6,y6,'k') ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (My)') ;

subplot(3,1,2) ;
plot(x4,y4b,'k',x5,y5b,'k') ;
hold on ;
plot(x6,y6b,'k') ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (Mz)') ;

subplot(3,1,3) ;
plot(x4,yr4,'k',x5,yr5,'k') ;
hold on ;
plot(x6,yr6,'k') ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (Mr)') ;


%% Diameter Change Locations

% D1 from x = 0 to x = (1.25/11.5)*Lt_c
% D2 from x = (1.25/11.5)*Lt_c   to x = (3.5/11.5)*Lt_c 
% D3 from x = (3.5/11.5)*Lt_c  to   x = (7.5/11.5)*Lt_c
% D4 from x = (7.5/11.5)*Lt_c  to   x = (10.25/11.5)*Lt_c 
% D5 from x = (10.25/11.5)*Lt_c   to   x = Lt_c

%% Choosing Moments

% Choose 4 Moments for Each Diameter Change

Mr4 = yr4(round((.5/2)*length(yr4))) ; % Mr @ x = (1.25/11.5)*Lt_c
Mr5 = yr5(round((.75/(8.5-2.75))*length(yr5))) ; % Mr @ x = (3.5/11.5)*Lt_c
Mr6 = yr5(round(((7.5-2.75)/(8.5-2.75))*length(yr5))) ; % Mr @ x = (7.5/11.5)*Lt_c
Mr7 = yr6(round(((10.25-8.5)/(10.75-8.5))*length(yr6))) ; % Mr @ x = (10.25/11.5)*Lt_c

%% Output Variables
Ti_c = -Tg ;
save("CountershaftStaticLoad.mat")

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% MEE 342 - Static Diameter Loop v2 

clear ; clc ; 
load("CountershaftStaticLoad.mat") ;
load('variables2.mat','Sut','Sy') ;


%% Algorithm ------ Concentration 1

D4 = (ns*sqrt((32*Mr4)^2 + (3*(16*0)^2)) / (Sy* pi) )^(1/3)      % initial guesses based off Von Mises @ Geometry Change
                                                                 % No T on D4
D5 = (ns*sqrt((32*Mr5)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3)   % D5 based on Mr5 and T! Not Mr4
D5b = (ns*sqrt((32*Mr5)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3) ;

if Mr5 > Mr6                                                     % Optimizes for most conservative guess
    D6 = (ns*sqrt((32*Mr5)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3)
    D6b = (ns*sqrt((32*Mr5)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3) ;
else
    D6 = (ns*sqrt((32*Mr6)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3)
    D6b = (ns*sqrt((32*Mr6)^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3) ;
end

D7 = (ns*sqrt((32*max(yr5))^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3) % use max moment to prevent hr = 0
D7b = (ns*sqrt((32*max(yr5))^2 + (3*(16*Ti_c)^2)) / (Sy* pi) )^(1/3) ;

D8 =(ns*sqrt((32*Mr7)^2 + (3*(16*0)^2)) / (Sy* pi) )^(1/3) 
                                                                 % No T on D8
                                                                 
  save("OPS.mat","D4","D5","D5b","D6","D6b","D7","D7b","D8");                                                               
%%

h1 = D5-D4 ; 
fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.25)
r2 = input('Enter Fillet Radius in Inches for concentration 1: ') ;
hr1 = h1/r2 ;
FS_verify = 0 ;

while FS_verify < 2 % Diameters 4 & 5
    
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r2)^(1/2)) - (0.131 * h1/r2);
         c2 = 0.022 - (3.405 * (h1/r2)^(1/2)) + (0.915 * h1/r2);
         c3 = 0.869 + (1.777 * (h1/r2)^(1/2)) - (0.555 * h1/r2);
         c4 = -.810 + (.422 * (h1/r2)^(1/2)) - (0.260 * h1/r2);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r2)^(1/2)) - (0.008 * h1/r2);
         c2 = -3.813 + (.968 * (h1/r2)^(1/2)) - (0.260 * h1/r2);
         c3 = 7.423 - (4.868 * (h1/r2)^(1/2)) + (0.869 * h1/r2);
         c4 = -3.839 + (3.070 * (h1/r2)^(1/2)) - (0.6 * h1/r2);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r2)^(1/2)) - (0.075 * h1/r2);
         c6 = -0.437 - (1.969 * (h1/r2)^(1/2)) + (0.553 * h1/r2);
         c7 = 1.557 + (1.073 * (h1/r2)^(1/2)) - (0.578 * h1/r2);
         c8 = -1.061 + (.171 * (h1/r2)^(1/2)) + (0.086 * h1/r2);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D5)) + (c3 * (2*h1/D5)^2) + (c4 * (2*h1/D5)^3) ;
    kts1 = c5 + (c6 * (2*h1/D5)) + (c7 * (2*h1/D5)^2) + (c8 * (2*h1/D5)^3) ;
    FS_verify = Sy / sqrt((32*Mr4*kt1 / (pi * D4^3))^2 + (3*(16*0*kts1 / (pi * D4^3))^2)) ; % Von Mises
    
    if FS_verify < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr1 < 3.9 
            D5 = D5 + .01 ; %in
            h1 = D5-D4 ;
            hr1 = h1/r2 ;
        else
            D4 = D4 + .01 ; %in
            h1 = D5-D4 ;
            hr1 = h1/r2 ;
        end
        
    else
        continue 
    end
end

% FS_verify
% D4
% D5

%% --------------------Concentration 2-----------------------------
clc;

load("OPS.mat","D4","D5","D5b","D6","D6b","D7","D7b","D8"); 
load("CountershaftStaticLoad.mat") ;
load('variables2.mat','Sut','Sy') ;

FS_verify2 = 0 ;
h2 = D6 - D5b ;
fprintf("\nChoose a fillet radius value between %5.3f and %5.3f \n (not including the bounds)\n", h2/4, h2/.25)
r3 = input('Enter Fillet Radius in Inches for concentration 2: ') ;
hr2 = h2/r3 ;

while FS_verify2 < 2 % Diameters 5b & 6
    
    % Bending kt values
    if hr2 >= .1 && hr2 <= 2 
         c1 = 0.947 + (1.206 * (h2/r3)^(1/2)) - (0.131 * h2/r3);
         c2 = 0.022 - (3.405 * (h2/r3)^(1/2)) + (0.915 * h2/r3);
         c3 = 0.869 + (1.777 * (h2/r3)^(1/2)) - (0.555 * h2/r3);
         c4 = -.810 + (.422 * (h2/r3)^(1/2)) - (0.260 * h2/r3);
    
    elseif hr2 >= 2 && hr2 <= 20
         c1 = 1.232 + (.832 * (h2/r3)^(1/2)) - (0.008 * h2/r3);
         c2 = -3.813 + (.968 * (h2/r3)^(1/2)) - (0.260 * h2/r3);
         c3 = 7.423 - (4.868 * (h2/r3)^(1/2)) + (0.869 * h2/r3);
         c4 = -3.839 + (3.070 * (h2/r3)^(1/2)) - (0.6 * h2/r3);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr2 >= .25 && hr2 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r3)^(1/2)) - (0.075 * h2/r3);
         c6 = -0.437 - (1.969 * (h2/r3)^(1/2)) + (0.553 * h2/r3);
         c7 = 1.557 + (1.073 * (h2/r3)^(1/2)) - (0.578 * h2/r3);
         c8 = -1.061 + (.171 * (h2/r3)^(1/2)) + (0.086 * h2/r3);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clc ;
    end

    kt2 = c1 + (c2 * (2*h2/D6)) + (c3 * (2*h2/D6)^2) + (c4 * (2*h2/D6)^3) ;
    kts2 = c5 + (c6 * (2*h2/D6)) + (c7 * (2*h2/D6)^2) + (c8 * (2*h2/D6)^3) ;
    FS_verify2 = Sy / sqrt((32*Mr5*kt2 / (pi * D5b^3))^2 + (3*(16*Ti_c*kts2 / (pi * D5b^3))^2)) ;  % Von Mises on D3
    
    if FS_verify2 < 2 
        if hr2 < 3.9 
            D5b = D5b + .01 ; %in
            h2 = D6 - D5b ;
            hr2 = h2/r3 ;
        else 
            D6 = D6 + .01 ;
            h2 = D6 - D5b ;
            hr2 = h2/r3 ;
        end
    else
        continue 
    end
end

% FS_verify2
% D5b
% D6

%% --------------------Concentration 3-----------------------------
clc;
FS_verify3 = 0 ;
h3 = D7 - D6b ;
fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h3/4, h3/.25)
r4 = input('Enter Fillet Radius in Inches for concentration 3: ') ;
hr3 = h3/r4 ;
while FS_verify3 < 2 % Diameters 6b & 7
    
    % Bending kt values
    if hr3 >= .1 && hr3 <= 2 
         c1 = 0.947 + (1.206 * (h3/r4)^(1/2)) - (0.131 * h3/r4);
         c2 = 0.022 - (3.405 * (h3/r4)^(1/2)) + (0.915 * h3/r4);
         c3 = 0.869 + (1.777 * (h3/r4)^(1/2)) - (0.555 * h3/r4);
         c4 = -.810 + (.422 * (h3/r4)^(1/2)) - (0.260 * h3/r4);
    
    elseif hr3 >= 2 && hr3 <= 20
         c1 = 1.232 + (.832 * (h3/r4)^(1/2)) - (0.008 * h3/r4);
         c2 = -3.813 + (.968 * (h3/r4)^(1/2)) - (0.260 * h3/r4);
         c3 = 7.423 - (4.868 * (h3/r4)^(1/2)) + (0.869 * h3/r4);
         c4 = -3.839 + (3.070 * (h3/r4)^(1/2)) - (0.6 * h3/r4);
    else 
         fprintf('Your fillet radius is invalid. Hr3 ratio is %4.3f \n',hr3)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr3 >= .25 && hr3 <= 4.0 
         c5 = 0.905 + (.783 * (h3/r4)^(1/2)) - (0.075 * h3/r4);
         c6 = -0.437 - (1.969 * (h3/r4)^(1/2)) + (0.553 * h3/r4);
         c7 = 1.557 + (1.073 * (h3/r4)^(1/2)) - (0.578 * h3/r4);
         c8 = -1.061 + (.171 * (h3/r4)^(1/2)) + (0.086 * h3/r4);
    else 
         fprintf('Your fillet radius is invalid. Hr3 ratio is %4.3f \n',hr3)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h3/D7)) + (c3 * (2*h3/D7)^2) + (c4 * (2*h3/D7)^3) ;
    kts2 = c5 + (c6 * (2*h3/D7)) + (c7 * (2*h3/D7)^2) + (c8 * (2*h3/D7)^3) ;
    FS_verify3 = Sy / sqrt((32*Mr6*kt2 / (pi * D6b^3))^2 + (3*(16*Ti_c*kts2 / (pi * D6b^3))^2)) ;  % Von Mises on D3
    
    if FS_verify3 < 2 
        if hr3 < 3.9 
            D7 = D7 + .01 ;
            h3 = D7 - D6b ;
            hr3 = h3/r4 ;
        else 
            D6b = D6b + .01 ; %in
            h3 = D7 - D6b ;
            hr3 = h3/r4 ;
        end
    else
        continue 
    end
end

% FS_verify3
% D6b
% D7

%% --------------------------Concentration 4------------------------------------
clc;
h4 = D7b-D8 ; 
fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h4/4, h4/.25)
r5 = input('Enter Fillet Radius in Inches for concentration 4: ') ;
hr4 = h4/r5 ;
FS_verify4 = 0 ;

while FS_verify4 < 2 % Diameters 7b & 8
    
    % Bending kt values
    if hr4 >= .1 && hr4 <= 2 
         c1 = 0.947 + (1.206 * (h4/r5)^(1/2)) - (0.131 * h4/r5);
         c2 = 0.022 - (3.405 * (h4/r5)^(1/2)) + (0.915 * h4/r5);
         c3 = 0.869 + (1.777 * (h4/r5)^(1/2)) - (0.555 * h4/r5);
         c4 = -.810 + (.422 * (h4/r5)^(1/2)) - (0.260 * h4/r5);
    
    elseif hr4 >= 2 && hr4 <= 20
         c1 = 1.232 + (.832 * (h4/r5)^(1/2)) - (0.008 * h4/r5);
         c2 = -3.813 + (.968 * (h4/r5)^(1/2)) - (0.260 * h4/r5);
         c3 = 7.423 - (4.868 * (h4/r5)^(1/2)) + (0.869 * h4/r5);
         c4 = -3.839 + (3.070 * (h4/r5)^(1/2)) - (0.6 * h4/r5);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr4)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr4 >= .25 && hr4 <= 4.0 
         c5 = 0.905 + (.783 * (h4/r5)^(1/2)) - (0.075 * h4/r5);
         c6 = -0.437 - (1.969 * (h4/r5)^(1/2)) + (0.553 * h4/r5);
         c7 = 1.557 + (1.073 * (h4/r5)^(1/2)) - (0.578 * h4/r5);
         c8 = -1.061 + (.171 * (h4/r5)^(1/2)) + (0.086 * h4/r5);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr4)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h4/D7b)) + (c3 * (2*h4/D7b)^2) + (c4 * (2*h4/D7b)^3) ;
    kts1 = c5 + (c6 * (2*h4/D7b)) + (c7 * (2*h4/D7b)^2) + (c8 * (2*h4/D7b)^3) ;
    FS_verify4 = Sy / sqrt((32*Mr7*kt1 / (pi * D8^3))^2 + (3*(16*0*kts1 / (pi * D8^3))^2)) ; % Von Mises
    
    if FS_verify4 < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr4 < 3.9 
            D7b = D7b + .01 ; %in
            h4 = D7b-D8 ;
            hr4 = h4/r5 ;
        else
            D8 = D8 + .01 ; %in
            h4 = D7b-D8 ;
            hr4 = h4/r5 ;
        end
        
    else
        continue 
    end
end

% FS_verify4
% D7b
% D8

%% Keyhole Concentration 

kt_hole = 2.14 ; % Approximation from Table 7-1 
kts_hole = 3 ;
D5_hole = ((ns/(Sy*pi)) * sqrt( (32*yr4(length(yr4))*kt_hole)^2 + (3*((16*Ti_c*kts_hole)^2)) ))^(1/3) ;
D7_hole = ((ns/(Sy*pi)) * sqrt( (32*yr5(length(yr5))*kt_hole)^2 + (3*((16*Ti_c*kts_hole)^2)) ))^(1/3) ;

%% Display 

% D4 
% %
% if D5 > D5b && D5 > D5_hole
%     D5
% elseif D5b > D5 && D5b > D5_hole
%     D5b 
% elseif D5_hole > D5 && D5_hole > D5b
%     D5_hole
% end
% %
% if D6 > D6b
%     D6
% else 
%     D6b
% end
% %
% if D7 > D7b && D7 > D7_hole
%     D7
% elseif D7b > D7 && D7b > D7_hole
%     D7b 
% elseif D7_hole > D7 && D7_hole > D7b
%     D7_hole
% end
% 
% D8
%     
% FS_verify
% FS_verify2
% FS_verify3
% FS_verify4

%% Output Variables

save('contershaft_analysis.mat') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% MEE 342 - Fatigue Diameter Loop v1     

clear ; clc ; 

load('contershaft_analysis.mat') ;

%% Inputs

fprintf("Type in one of the following values for the corresponding finish \nGround = 1 \nMachined/Cold-drawn = 2 \nHot-rolled = 3 \nAs-forged = 4 \n")
A = input("Select from list above ") ;
fprintf("\n Input operating temperature in degrees fahrenheit \n")
Tf = input("Temperature = ") ;

%% Algorithm ------ Concentration 1

% h1 = D5-D4 ; 
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.24)
% r2 = input('Enter Fillet Radius in Inches: ') 
% hr1 = h1/r2 ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 1: ") ;
FS_verify_f = 0 ;

while FS_verify_f < 2 % Diameters 4 & 5
    
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r2)^(1/2)) - (0.131 * h1/r2);
         c2 = 0.022 - (3.405 * (h1/r2)^(1/2)) + (0.915 * h1/r2);
         c3 = 0.869 + (1.777 * (h1/r2)^(1/2)) - (0.555 * h1/r2);
         c4 = -.810 + (.422 * (h1/r2)^(1/2)) - (0.260 * h1/r2);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r2)^(1/2)) - (0.008 * h1/r2);
         c2 = -3.813 + (.968 * (h1/r2)^(1/2)) - (0.260 * h1/r2);
         c3 = 7.423 - (4.868 * (h1/r2)^(1/2)) + (0.869 * h1/r2);
         c4 = -3.839 + (3.070 * (h1/r2)^(1/2)) - (0.6 * h1/r2);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r2)^(1/2)) - (0.075 * h1/r2);
         c6 = -0.437 - (1.969 * (h1/r2)^(1/2)) + (0.553 * h1/r2);
         c7 = 1.557 + (1.073 * (h1/r2)^(1/2)) - (0.578 * h1/r2);
         c8 = -1.061 + (.171 * (h1/r2)^(1/2)) + (0.086 * h1/r2);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D5)) + (c3 * (2*h1/D5)^2) + (c4 * (2*h1/D5)^3) ;
    kts1 = c5 + (c6 * (2*h1/D5)) + (c7 * (2*h1/D5)^2) + (c8 * (2*h1/D5)^3) ;    
    
    % Fatigue Concentration Factors----------------------------------------
    
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r2) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r2) ) ) ;

    
    % Endurance Limit------------------------------------------------------
    
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
       
% for ka
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;

% for kb

    if D4 >= .1 && D4 <= 2
        %.1 <= D1 <= 2
        kb = 0.879*D4^-0.107 ; 
    else
        kb = 0.91*D4^-0.157 ;
    end
    
% for kc

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% ke,kf_coeff

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
    %----------------------------------------------------------------------
%Goodman 

% for Stress Amplitude

    Sa_D4 = sqrt((32*Mr4*kf_1 / (pi * D4^3))^2) ;

% for Mean Stress

    Sm_D4 = sqrt((3*(16*0*kfs_1 / (pi * D4^3))^2)) ; %No Torque at first concentration
    
% Factor of Safety via Modified Goodman
    FS_verify_f = ( (Sa_D4/(Se*10^3)) + (Sm_D4/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr1 < 3.9 
            D5 = D5 + .01 ; %in
            h1 = D5-D4 ;
            hr1 = h1/r2 ;
        else
            D4 = D4 + .01 ; %in
            h1 = D5-D4 ;
            hr1 = h1/r2 ;
        end
        
    else
        continue 
    end
end

% FS_verify
% D4
% D5

%% --------------------Concentration 2-----------------------------

% h2 = D6 - D5b ;
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h2/4, h2/.24)
% r3 = input('Enter Fillet Radius in Inches: ') 
% hr2 = h2/r3 ;

fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 2: ") ;
FS_verify_f2 = 0 ;
while FS_verify_f2 < 2 % Diameters 5b & 6
    
    % Bending kt values
    if hr2 >= .1 && hr2 <= 2 
         c1 = 0.947 + (1.206 * (h2/r3)^(1/2)) - (0.131 * h2/r3);
         c2 = 0.022 - (3.405 * (h2/r3)^(1/2)) + (0.915 * h2/r3);
         c3 = 0.869 + (1.777 * (h2/r3)^(1/2)) - (0.555 * h2/r3);
         c4 = -.810 + (.422 * (h2/r3)^(1/2)) - (0.260 * h2/r3);
    
    elseif hr2 >= 2 && hr2 <= 20
         c1 = 1.232 + (.832 * (h2/r3)^(1/2)) - (0.008 * h2/r3);
         c2 = -3.813 + (.968 * (h2/r3)^(1/2)) - (0.260 * h2/r3);
         c3 = 7.423 - (4.868 * (h2/r3)^(1/2)) + (0.869 * h2/r3);
         c4 = -3.839 + (3.070 * (h2/r3)^(1/2)) - (0.6 * h2/r3);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr2 >= .25 && hr2 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r3)^(1/2)) - (0.075 * h2/r3);
         c6 = -0.437 - (1.969 * (h2/r3)^(1/2)) + (0.553 * h2/r3);
         c7 = 1.557 + (1.073 * (h2/r3)^(1/2)) - (0.578 * h2/r3);
         c8 = -1.061 + (.171 * (h2/r3)^(1/2)) + (0.086 * h2/r3);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h2/D6)) + (c3 * (2*h2/D6)^2) + (c4 * (2*h2/D6)^3) ;
    kts2 = c5 + (c6 * (2*h2/D6)) + (c7 * (2*h2/D6)^2) + (c8 * (2*h2/D6)^3) ;

        % Fatigue Concentration Factors----------------------------------------
    
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r3) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r3) ) ) ;

    
    % Endurance Limit------------------------------------------------------
    
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
       
% for ka
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;

% for kb

    if D5b >= .1 && D5b <= 2
        kb = 0.879*D5b^-0.107 ; 
    else
        kb = 0.91*D5b^-0.157 ;
    end
    
% for kc

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% ke,kf_coeff

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
    %----------------------------------------------------------------------
%Goodman 

% for Stress Amplitude

    Sa_D5b = sqrt((32*Mr5*kf_1 / (pi * D5b^3))^2) ;

% for Mean Stress

    Sm_D5b = sqrt((3*(16*Ti_c*kfs_1 / (pi * D5b^3))^2)) ; 
    
% Factor of Safety via Modified Goodman
    FS_verify_f2 = ( (Sa_D5b/(Se*10^3)) + (Sm_D5b/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f2 < 2 
        if hr2 < 3.9 
            D5b = D5b + .01 ; %in
            h2 = D6 - D5b ;
            hr2 = h2/r3 ;
        else 
            D6 = D6 + .01 ;
            h2 = D6 - D5b ;
            hr2 = h2/r3 ;
        end
    else
        continue 
    end
end

% FS_verify2
% D5b
% D6

%% --------------------Concentration 3-----------------------------

fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 3: ") ;
FS_verify_f3 = 0 ;
h3 = D7 - D6b ;
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h3/4, h3/.24)
% r4 = input('Enter Fillet Radius in Inches: ') 
hr3 = h3/r4 ; 
while FS_verify_f3 < 2 % Diameters 6b & 7
    
    % Bending kt values
    if hr3 >= .1 && hr3 <= 2 
         c1 = 0.947 + (1.206 * (h3/r4)^(1/2)) - (0.131 * h3/r4);
         c2 = 0.022 - (3.405 * (h3/r4)^(1/2)) + (0.915 * h3/r4);
         c3 = 0.869 + (1.777 * (h3/r4)^(1/2)) - (0.555 * h3/r4);
         c4 = -.810 + (.422 * (h3/r4)^(1/2)) - (0.260 * h3/r4);
    
    elseif hr3 >= 2 && hr3 <= 20
         c1 = 1.232 + (.832 * (h3/r4)^(1/2)) - (0.008 * h3/r4);
         c2 = -3.813 + (.968 * (h3/r4)^(1/2)) - (0.260 * h3/r4);
         c3 = 7.423 - (4.868 * (h3/r4)^(1/2)) + (0.869 * h3/r4);
         c4 = -3.839 + (3.070 * (h3/r4)^(1/2)) - (0.6 * h3/r4);
    else 
         fprintf('Your fillet radius is invalid. Hr3 ratio is %4.3f \n',hr3)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr3 >= .25 && hr3 <= 4.0 
         c5 = 0.905 + (.783 * (h3/r4)^(1/2)) - (0.075 * h3/r4);
         c6 = -0.437 - (1.969 * (h3/r4)^(1/2)) + (0.553 * h3/r4);
         c7 = 1.557 + (1.073 * (h3/r4)^(1/2)) - (0.578 * h3/r4);
         c8 = -1.061 + (.171 * (h3/r4)^(1/2)) + (0.086 * h3/r4);
    else 
         fprintf('Your fillet radius is invalid. Hr3 ratio is %4.3f \n',hr3)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h3/D7)) + (c3 * (2*h3/D7)^2) + (c4 * (2*h3/D7)^3) ;
    kts2 = c5 + (c6 * (2*h3/D7)) + (c7 * (2*h3/D7)^2) + (c8 * (2*h3/D7)^3) ;

    % Fatigue Concentration Factors----------------------------------------
    
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r4) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r4) ) ) ;

    
    % Endurance Limit------------------------------------------------------
    
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
       
% for ka
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;

% for kb

    if D6 >= .1 && D6 <= 2
        kb = 0.879*D6^-0.107 ; 
    else
        kb = 0.91*D6^-0.157 ;
    end
    
% for kc

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% ke,kf_coeff

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
    %----------------------------------------------------------------------
%Goodman 

% for Stress Amplitude

    Sa_D6b = sqrt((32*Mr6*kf_1 / (pi * D6b^3))^2) ;

% for Mean Stress

    Sm_D6b = sqrt((3*(16*Ti_c*kfs_1 / (pi * D6b^3))^2)) ; 
    
% Factor of Safety via Modified Goodman
    FS_verify_f3 = ( (Sa_D6b/(Se*10^3)) + (Sm_D6b/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f3 < 2 
        if hr3 < 3.9 
            D7 = D7 + .01 ;
            h3 = D7 - D6b ;
            hr3 = h3/r4 ;
        else 
            D6b = D6b + .01 ; %in
            h3 = D7 - D6b ;
            hr3 = h3/r4 ;
        end
    else
        continue 
    end
end

% FS_verify3
% D6b
% D7

%% --------------------------Concentration 4------------------------------------

fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 4: ") ;
h4 = D7b-D8 ; 
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h4/4, h4/.24)
% r5 = input('Enter Fillet Radius in Inches: ') 
hr4 = h4/r5 ;
FS_verify_f4 = 0 ;

while FS_verify_f4 < 2 % Diameters 7b & 8
    
    % Bending kt values
    if hr4 >= .1 && hr4 <= 2 
         c1 = 0.947 + (1.206 * (h4/r5)^(1/2)) - (0.131 * h4/r5);
         c2 = 0.022 - (3.405 * (h4/r5)^(1/2)) + (0.915 * h4/r5);
         c3 = 0.869 + (1.777 * (h4/r5)^(1/2)) - (0.555 * h4/r5);
         c4 = -.810 + (.422 * (h4/r5)^(1/2)) - (0.260 * h4/r5);
    
    elseif hr4 >= 2 && hr4 <= 20
         c1 = 1.232 + (.832 * (h4/r5)^(1/2)) - (0.008 * h4/r5);
         c2 = -3.813 + (.968 * (h4/r5)^(1/2)) - (0.260 * h4/r5);
         c3 = 7.423 - (4.868 * (h4/r5)^(1/2)) + (0.869 * h4/r5);
         c4 = -3.839 + (3.070 * (h4/r5)^(1/2)) - (0.6 * h4/r5);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr4)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr4 >= .25 && hr4 <= 4.0 
         c5 = 0.905 + (.783 * (h4/r5)^(1/2)) - (0.075 * h4/r5);
         c6 = -0.437 - (1.969 * (h4/r5)^(1/2)) + (0.553 * h4/r5);
         c7 = 1.557 + (1.073 * (h4/r5)^(1/2)) - (0.578 * h4/r5);
         c8 = -1.061 + (.171 * (h4/r5)^(1/2)) + (0.086 * h4/r5);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr4)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h4/D7b)) + (c3 * (2*h4/D7b)^2) + (c4 * (2*h4/D7b)^3) ;
    kts1 = c5 + (c6 * (2*h4/D7b)) + (c7 * (2*h4/D7b)^2) + (c8 * (2*h4/D7b)^3) ;

        % Endurance Limit------------------------------------------------------
    
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
       
% for ka
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;

% for kb

    if D8 >= .1 && D8 <= 2
        kb = 0.879*D8^-0.107 ; 
    else
        kb = 0.91*D8^-0.157 ;
    end
    
% for kc

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% ke,kf_coeff

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
    %----------------------------------------------------------------------
%Goodman 

% for Stress Amplitude

    Sa_D8 = sqrt((32*Mr7*kf_1 / (pi * D8^3))^2) ;

% for Mean Stress

    Sm_D8 = sqrt((3*(16*0*kfs_1 / (pi * D8^3))^2)) ; 
    
% Factor of Safety via Modified Goodman
    FS_verify_f4 = ( (Sa_D8/(Se*10^3)) + (Sm_D8/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f4 < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr4 < 3.9 
            D7b = D7b + .01 ; %in
            h4 = D7b-D8 ;
            hr4 = h4/r5 ;
        else
            D8 = D8 + .01 ; %in
            h4 = D7b-D8 ;
            hr4 = h4/r5 ;
        end
        
    else
        continue 
    end
end

% FS_verify4
% D7b
% D8

%% Output Variables

save('countershaft_fatigue.mat') ;


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% MEE 342 - Key Hole Concentration     

clear ; clc ;
load('countershaft_fatigue.mat') ;

%% Fatigue Algorithm at Key Hole F

FS_verify_f_holeF = 0 ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for key 1: ") ;

while FS_verify_f_holeF < 2 
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    r_key = 0.2 * D5_hole ;
    kf_hole = 1 + ( (kt_hole - 1) / (1 + a_b/sqrt(r_key) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_hole = 1 + ( (kts_hole - 1) / (1 + a_t/sqrt(r_key) ) ) ;
    
    Sa_D5_hole = sqrt((32*yr4(length(yr4))*kf_hole / (pi * D5_hole^3))^2) ;
    Sm_D5_hole = sqrt((3*(16*Ti_c*kfs_hole / (pi * D5_hole^3))^2)) ;
    
    % ------------------------------ Se ----------------------------------
    
    % ka --------------------
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;
    
    % kb ---------------------------
    if D5_hole >= .1 && D5_hole <= 2
        kb = 0.879*D5_hole^-0.107 ; 
    else
        kb = 0.91*D5_hole^-0.157 ;
    end  
    
    % for kc ---------------------------

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
% End of Se Loop
    
    FS_verify_f_holeF = ( (Sa_D5_hole/(Se*10^3)) + (Sm_D5_hole/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f_holeF < 2
        D5_hole = D5_hole + .01 ;
    else 
        continue 
    end
end 

% FS_verify_f_holeF
% D5_hole

%% Fatigue Algorithm at Key Hole G

FS_verify_f_holeG = 0 ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for key 2: ") ;

while FS_verify_f_holeG < 2 
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    r_key = 0.2 * D7_hole ;
    kf_hole = 1 + ( (kt_hole - 1) / (1 + a_b/sqrt(r_key) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_hole = 1 + ( (kts_hole - 1) / (1 + a_t/sqrt(r_key) ) ) ;
    
    Sa_D7_hole = sqrt((32*yr5(length(yr5))*kf_hole / (pi * D7_hole^3))^2) ;
    Sm_D7_hole = sqrt((3*(16*Ti_c*kfs_hole / (pi * D7_hole^3))^2)) ;
    
    % ------------------------------ Se ----------------------------------
    
    % ka --------------------
    if A == 1
        a = 1.34 ; b = -0.085 ; %Ground
    elseif A== 2
        a = 2.7 ; b = -0.265 ; % Machined or Cold drawn
    elseif A == 3
        a = 14.4 ; b = -0.718 ; % hot-rolled
    elseif A == 4
        a = 219.9 ; b = -0.995 ; % as-forged
    end
    
    ka = (a*Sut^b) ;
    
    % kb ---------------------------
    if D7_hole >= .1 && D7_hole <= 2
        kb = 0.879*D7_hole^-0.107 ; 
    else
        kb = 0.91*D7_hole^-0.157 ;
    end  
    
    % for kc ---------------------------

    if B == 1
        kc = .59 ;
    elseif B == 2
        kc = 1 ;
    elseif B == 3
        kc = .85 ;
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % for kd

    kd = 0.975 + 0.432*(10^-3)*Tf - 0.115*(10^-5)*Tf^2 + 0.104*(10^-8)*Tf^3 - 0.595*(10^-12)*Tf^4 ;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ke = 1 ;
    kf_coeff = 1 ; 

    Se = Se1*ka*kb*kc*kd*ke*kf_coeff ; % in Ksi
    
% End of Se Loop
    
    FS_verify_f_holeG = ( (Sa_D7_hole/(Se*10^3)) + (Sm_D7_hole/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f_holeG < 2
        D7_hole = D7_hole + .01 ;
    else 
        continue 
    end
end 

% FS_verify_f_holeG
% D7_hole

%% Output Variables

save('countershaft_key_fatigue.mat') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


%% MEE 342 - Diameter Outputs Countershaft    

clear ; clc ;

load('countershaft_key_fatigue.mat') ;

%% 

if D5 > D5b && D5 > D5_hole
    D5_final = D5 ;
elseif D5b > D5 && D5b > D5_hole
    D5_final = D5b ;
elseif D5_hole > D5 && D5_hole > D5b
    D5_final = D5_hole ;
end

if D6 >= D6b 
    D6_final = D6 ;
elseif D6b > D6 
    D6_final = D6b ;
end

if D7 > D7b && D7 > D7_hole
    D7_final = D7 ;
elseif D7b > D7 && D7b > D7_hole
    D7_final = D7b ;
elseif D7_hole > D7 && D7_hole > D7b
    D7_final = D7_hole ;
end

%%

fprintf('The countershaft is stepped with five diameters. \nAll diameters have been optimized to minimize size and weight yet meet the safety factor requirements.\n')
fprintf('Diameter 4 from x = 0 to x = %5.4f is %5.4f inches. \n', x4(round((.5/2)*length(yr4))), ceil(D4 * 40) / 40) ;
fprintf('Diameter 5 from x = %5.4f to %5.4f is %5.4f inches. \n', x4(round((.5/2)*length(yr4))), x5(round((.75/(8.5-2.75))*length(yr5))), ceil(D5_final * 40) / 40) ;
fprintf('Diameter 6 from x = %5.4f to %5.4f is %5.4f inches. \n', x5(round((.75/(8.5-2.75))*length(yr5))), x5(round(((7.5-2.75)/(8.5-2.75))*length(yr5))), ceil(D6_final * 40) / 40) ;
fprintf('Diameter 7 from x = %5.4f to %5.4f is %5.4f inches. \n', x5(round(((7.5-2.75)/(8.5-2.75))*length(yr5))), x6(round(((10.25-8.5)/(10.75-8.5))*length(yr6))), ceil(D7_final * 40) / 40) ;
fprintf('Diameter 8 from x = %5.4f to %5.4f is %5.4f inches. \n', x6(round(((10.25-8.5)/(10.75-8.5))*length(yr6))), Lt_c, ceil(D8 * 40) / 40) ;


