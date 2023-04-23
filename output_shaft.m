
clc;
clear;


%% Static Analysis Output Shaft

clc;
clear ;

%% Inputs
load("InputShaft.mat","Rb","ns","nf")
load("CountershaftStaticLoad.mat","Gz","Gy","Ti_c")

a = 6 ;
L_OS = 0;
while L_OS <= a 

fprintf("Please give a shaft length larger than %5.3f inches \n",a)
L_OS = input ("Output Shaft Length = ") ;

if L_OS <= a 
   fprintf("\nThat input value will not work \n") 
end
end
fprintf("\n");

R_G5 = Rb;


%% Distance Between Forces

L1 = (L_OS - 4) * (.25/11.5) ; % Gap to First Bearing
L2 = (L_OS -4) * (.5) ; % Length to Gear Loads
L3 = (L_OS - 4) - L1 ; 

%% Reaction Forces

F_gy = -Gy ;
F_gz = -Gz ;

R_b2z =  -(F_gz*(L2-L1))/(L3-L1) ;
R_b1z = -(F_gz + R_b2z) ;
 
R_b2y = [F_gy*(L2-L1)]/(L3-L1) ; 
R_b1y = -F_gy - R_b2y ;

%% Bending Moments

% My Diagran
x1 = L1:L_OS/10000:L2 ;
x2 = L2:L_OS/10000:L3 ;
Myb = R_b1y*(L2-L1) ;
y1 = (((Myb)./(L2 - L1)).*x1) - (((Myb)./(L2 - L1)).*L1) ; %piece 1
y2 = -(((Myb)./(L3 - L2)).*x2) + (((Myb)./(L3 - L2)).*L3) ; %piece 2 
figure(1) ; 
subplot(2,1,1) ;
plot(x1,y1,'k',x2,y2,'k') ;
xlim([0 L_OS])
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (My)') ;

% Mz Diagram
Mzb = -R_b1z*(L2-L1) ; 
y3 = (((Mzb)./(L2 - L1)).*x1) - (((Mzb)./(L2 - L1)).*L1) ; %piece 1
y4 = -(((Mzb)./(L3 - L2)).*x2) + (((Mzb)./(L3 - L2)).*L3) ; %piece 2
x3 = L3:L_OS/1000:L_OS ;
subplot(2,1,2) ;
plot(x1,y3,'k',x2,y4,'k') ;
xlim([0 L_OS]) ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (Mz)') ;

% Mr Diagram

yr1 = sqrt((y1).^2 + (y3).^2) ;
yr2 = sqrt((y2).^2 + (y4).^2) ;
figure(2)
subplot(2,1,1) ;
plot(x1,yr1,'k',x2,yr2,'k') ; 
hold on
plot([x2(length(x2)) L_OS], [ 0 0 ], 'k' )
hold off
xlim([0 L_OS]) ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Resultant Moment Diagram (Mr)') ;


% Torque Diagram 

T_G = -Ti_c ; % Internal Torque
subplot(2,1,2) ;
plot([0 0],[0 T_G], 'k') ;
hold on ;
plot([0 L2],[T_G T_G], 'k') ;
plot([L2 L2],[T_G 0],'k') ;
xlim([0 L_OS]) ;
xlabel('X Distance [in]') ; ylabel('Torque [lb-in]') ; title('Torque Diagram (Mx)') ;

% Shear in Y Diagram

figure(3)
subplot(2,1,1) ;
hold on
plot([ 0 L1 ],[ 0 0 ], 'k')
plot([ L1 L1 ],[ 0 R_b1y ], 'k')
plot([ L1 (L2) ],[ R_b1y R_b1y ], 'k')
plot([ (L2) (L2) ],[ R_b1y (R_b1y + F_gy) ], 'k')
plot([ (L2) (L3) ],[ (R_b1y + F_gy) (R_b1y + F_gy) ], 'k')
plot([ (L3) (L3) ],[ (R_b1y + F_gy) (R_b1y + F_gy + R_b2y) ], 'k')
plot([ (L3) (L_OS) ],[ (R_b1y + F_gy + R_b2y) (R_b1y + F_gy + R_b2y) ], 'k')
hold off
xlabel('X Distance [in]') ; ylabel('Shear Force[lb]') ; title('Shear Force along Y-axis (Fy)') ;

% Shear in Z Diagram

subplot(2,1,2) ;
hold on
plot([ 0 L1 ],[ 0 0 ], 'k')
plot([ L1 L1 ],[ 0 R_b1z ], 'k')
plot([ L1 (L2) ],[ R_b1z R_b1z ], 'k')
plot([ (L2) (L2) ],[ R_b1z (R_b1z + F_gz) ], 'k')
plot([ (L2) (L3) ],[ (R_b1z + F_gz) (R_b1z + F_gz) ], 'k')
plot([ (L3) (L3) ],[ (R_b1z + F_gz) (R_b1z + F_gz + R_b2z) ], 'k')
plot([ (L3) (L_OS) ],[ (R_b1z + F_gz + R_b2z) (R_b1z + F_gz + R_b2z) ], 'k')
hold off
xlabel('X Distance [in]') ; ylabel('Shear Force[lb]') ; title('Shear Force along Z-axis (Fz)') ;


%% Choosing Moments

Mr1 = yr1(round(length(yr1)/2)) ; 
Mr2 = yr1(length(yr1)) ; 
Mr3 = yr2(round(length(yr2)/2)) ;

save("OutputShaftStatic.mat","Mr1","Mr2","Mr3","T_G","ns","nf","Rb","yr1","yr2","x1","x2","L_OS")

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


clc;
clear;

%% Output Shaft --- Static FS Diameters 

load("OutputShaftStatic.mat","Mr1","Mr2","Mr3","T_G","ns","nf","Rb","yr1","yr2","x1","x2","L_OS") ;
Ti = T_G ;
n = 2 ;
Sy = 50000 ; % psi
Sut = 68000 ; % psi
Mr8 = Mr1 ;
Mr9 = Mr3 ;

%% Algorithm 

D9 = (n*sqrt((32*Mr8)^2 + (3*(16*0)^2)) / (Sy* pi) )^(1/3) % initial guesses based off Von Mises
                                                           % No T on D9
D10 = (n*sqrt((32*Mr2)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3)
D10b = D10 ; % duplicate for preservation of Static Diameter ;
D11 = (n*sqrt((32*Mr9)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3)
h1 = D10-D9 ; 

fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.25)
r6 = input('Enter Fillet Radius in Inches: ') 
hr1 = h1/r6 ;
FS_verify_8 = 0 ;

while FS_verify_8 < 2 % Diameters 1 & 2
    
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r6)^(1/2)) - (0.131 * h1/r6);
         c2 = 0.022 - (3.405 * (h1/r6)^(1/2)) + (0.915 * h1/r6);
         c3 = 0.869 + (1.777 * (h1/r6)^(1/2)) - (0.555 * h1/r6);
         c4 = -.810 + (.422 * (h1/r6)^(1/2)) - (0.260 * h1/r6);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r6)^(1/2)) - (0.008 * h1/r6);
         c2 = -3.813 + (.968 * (h1/r6)^(1/2)) - (0.260 * h1/r6);
         c3 = 7.423 - (4.868 * (h1/r6)^(1/2)) + (0.869 * h1/r6);
         c4 = -3.839 + (3.070 * (h1/r6)^(1/2)) - (0.6 * h1/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r6)^(1/2)) - (0.075 * h1/r6);
         c6 = -0.437 - (1.969 * (h1/r6)^(1/2)) + (0.553 * h1/r6);
         c7 = 1.557 + (1.073 * (h1/r6)^(1/2)) - (0.578 * h1/r6);
         c8 = -1.061 + (.171 * (h1/r6)^(1/2)) + (0.086 * h1/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D10)) + (c3 * (2*h1/D10)^2) + (c4 * (2*h1/D10)^3) ;
    kts1 = c5 + (c6 * (2*h1/D10)) + (c7 * (2*h1/D10)^2) + (c8 * (2*h1/D10)^3) ;
    FS_verify_8 = Sy / sqrt((32*Mr8*kt1 / (pi * D9^3))^2 + (3*(16*0*kts1 / (pi * D9^3))^2)) ; % Von Mises
    
    if FS_verify_8 < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr1 < 3.9 
            D10 = D10 + .01 ; %in
            h1 = D10-D9 ;
            hr1 = h1/r6 ;
        else
            D9 = D9 + .01 ; %in
            h1 = D10-D9 ;
            hr1 = h1/r6 ;
        end
        
    else
        continue 
    end
end

% FS_verify
% D9
% D10

%% ----------------------Concentration 2---------------------------
FS_verify_9 = 0 ;
h2 = D10b - D11 ;
hr2 = h2/r6 ;
while FS_verify_9 < 2 % Diameters 2 & 3
    
    % Bending kt values
    if hr2 >= .1 && hr2 <= 2 
         c1 = 0.947 + (1.206 * (h2/r6)^(1/2)) - (0.131 * h2/r6);
         c2 = 0.022 - (3.405 * (h2/r6)^(1/2)) + (0.915 * h2/r6);
         c3 = 0.869 + (1.777 * (h2/r6)^(1/2)) - (0.555 * h2/r6);
         c4 = -.810 + (.422 * (h2/r6)^(1/2)) - (0.260 * h2/r6);
    
    elseif hr2 >= 2 && hr2 <= 20
         c1 = 1.232 + (.832 * (h2/r6)^(1/2)) - (0.008 * h2/r6);
         c2 = -3.813 + (.968 * (h2/r6)^(1/2)) - (0.260 * h2/r6);
         c3 = 7.423 - (4.868 * (h2/r6)^(1/2)) + (0.869 * h2/r6);
         c4 = -3.839 + (3.070 * (h2/r6)^(1/2)) - (0.6 * h2/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr2 >= .25 && hr2 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r6)^(1/2)) - (0.075 * h2/r6);
         c6 = -0.437 - (1.969 * (h2/r6)^(1/2)) + (0.553 * h2/r6);
         c7 = 1.557 + (1.073 * (h2/r6)^(1/2)) - (0.578 * h2/r6);
         c8 = -1.061 + (.171 * (h2/r6)^(1/2)) + (0.086 * h2/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h2/D10b)) + (c3 * (2*h2/D10b)^2) + (c4 * (2*h2/D10b)^3) ;
    kts2 = c5 + (c6 * (2*h2/D10b)) + (c7 * (2*h2/D10b)^2) + (c8 * (2*h2/D10b)^3) ;
    FS_verify_9 = Sy / sqrt((32*Mr9*kt2 / (pi * D11^3))^2 + (3*(16*Ti*kts2 / (pi * D11^3))^2)) ;  % Von Mises on D11
    
    if FS_verify_9 < 2 
        if hr2 < 3.9 
            D10b = D10b + .01 ; %in
            h2 = D10b-D11 ;
            hr2 = h2/r6 ;
        else 
            D11 = D11 + .01 ;
            h2 = D10b-D11 ;
            hr2 = h2/r6 ;
        end
    else
        continue 
    end
end

% FS_verify2
% D10b
% D11

%% Keyhole Concentration 

kt_hole = 2.14 ; % Approximation from Table 7-1 
kts_hole = 3 ;
D10_hole = ((n/(Sy*pi)) * sqrt( (32*Mr2*kt_hole)^2 + (3*((16*Ti*kts_hole)^2)) ))^(1/3) ;

%% Display 

D9 
if D10 > D10b 
    D10
elseif D10b > D10 
    D10b 
elseif D10_hole
    D10_hole
end
D11
FS_verify_8
FS_verify_9


save("OutputShaftDiameter_Guess.mat",'yr1','yr2','Ti','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r6', 'D9','D10','D10b', 'D10_hole', 'D11', 'Sut','Mr8','Mr9','Mr2','x1','x2','L_OS') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% MEE 342 - Output Shaft Fatigue

clear ; clc ;
load("OutputShaftDiameter_Guess.mat") ;

%% Inputs

fprintf("Type in one of the following values for the corresponding finish \nGround = 1 \nMachined/Cold-drawn = 2 \nHot-rolled = 3 \nAs-forged = 4 \n")
A = input("Select from list above ") ;
fprintf("\n Input operating temperature in degrees fahrenheit \n")
Tf = input("Temperature = ") ;

%% Algorithm ------ Concentration 1

h1 = D10-D9 ; 
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.24)
% r6 = input('Enter Fillet Radius in Inches: ') 
hr1 = h1/r6 ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 1: ") ;
FS_verify_f8 = 0 ;

while FS_verify_f8 < 2 % Diameters 4 & 5
    
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r6)^(1/2)) - (0.131 * h1/r6);
         c2 = 0.022 - (3.405 * (h1/r6)^(1/2)) + (0.915 * h1/r6);
         c3 = 0.869 + (1.777 * (h1/r6)^(1/2)) - (0.555 * h1/r6);
         c4 = -.810 + (.422 * (h1/r6)^(1/2)) - (0.260 * h1/r6);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r6)^(1/2)) - (0.008 * h1/r6);
         c2 = -3.813 + (.968 * (h1/r6)^(1/2)) - (0.260 * h1/r6);
         c3 = 7.423 - (4.868 * (h1/r6)^(1/2)) + (0.869 * h1/r6);
         c4 = -3.839 + (3.070 * (h1/r6)^(1/2)) - (0.6 * h1/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r6)^(1/2)) - (0.075 * h1/r6);
         c6 = -0.437 - (1.969 * (h1/r6)^(1/2)) + (0.553 * h1/r6);
         c7 = 1.557 + (1.073 * (h1/r6)^(1/2)) - (0.578 * h1/r6);
         c8 = -1.061 + (.171 * (h1/r6)^(1/2)) + (0.086 * h1/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D10)) + (c3 * (2*h1/D10)^2) + (c4 * (2*h1/D10)^3) ;
    kts1 = c5 + (c6 * (2*h1/D10)) + (c7 * (2*h1/D10)^2) + (c8 * (2*h1/D10)^3) ;    
    
    % Fatigue Concentration Factors----------------------------------------
    
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r6) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r6) ) ) ;

    
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

    if D9 >= .1 && D9 <= 2
        %.1 <= D1 <= 2
        kb = 0.879*D9^-0.107 ; 
    else
        kb = 0.91*D9^-0.157 ;
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

    Sa_D9 = sqrt((32*Mr8*kf_1 / (pi * D9^3))^2) ;

% for Mean Stress

    Sm_D9 = sqrt((3*(16*Ti*kfs_1 / (pi * D9^3))^2)) ; 
    
% Factor of Safety via Modified Goodman
    FS_verify_f8 = ( (Sa_D9/(Se*10^3)) + (Sm_D9/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f8 < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr1 < 3.9 
            D10 = D10 + .01 ; %in
            h1 = D10-D9 ;
            hr1 = h1/r6 ;
        else
            D9 = D9 + .01 ; %in
            h1 = D10-D9 ;
            hr1 = h1/r6 ;
        end
        
    else
        continue 
    end
end

% FS_verify
% D9
% D10

%% --------------------Concentration 2-----------------------------

h2 = D10b - D11 ;
% fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h2/4, h2/.24)
% r6 = input('Enter Fillet Radius in Inches: ') 
hr6 = h2/r6 ;

fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 2: ") ;
FS_verify_f9 = 0 ;
while FS_verify_f9 < 2 % Diameters 5b & 6
    
    % Bending kt values
    if hr6 >= .1 && hr6 <= 2 
         c1 = 0.947 + (1.206 * (h2/r6)^(1/2)) - (0.131 * h2/r6);
         c2 = 0.022 - (3.405 * (h2/r6)^(1/2)) + (0.915 * h2/r6);
         c3 = 0.869 + (1.777 * (h2/r6)^(1/2)) - (0.555 * h2/r6);
         c4 = -.810 + (.422 * (h2/r6)^(1/2)) - (0.260 * h2/r6);
    
    elseif hr6 >= 2 && hr6 <= 20
         c1 = 1.232 + (.832 * (h2/r6)^(1/2)) - (0.008 * h2/r6);
         c2 = -3.813 + (.968 * (h2/r6)^(1/2)) - (0.260 * h2/r6);
         c3 = 7.423 - (4.868 * (h2/r6)^(1/2)) + (0.869 * h2/r6);
         c4 = -3.839 + (3.070 * (h2/r6)^(1/2)) - (0.6 * h2/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr6 ratio is %4.3f \n',hr6)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr6 >= .25 && hr6 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r6)^(1/2)) - (0.075 * h2/r6);
         c6 = -0.437 - (1.969 * (h2/r6)^(1/2)) + (0.553 * h2/r6);
         c7 = 1.557 + (1.073 * (h2/r6)^(1/2)) - (0.578 * h2/r6);
         c8 = -1.061 + (.171 * (h2/r6)^(1/2)) + (0.086 * h2/r6);
    else 
         fprintf('Your fillet radius is invalid. Hr6 ratio is %4.3f \n',hr6)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h2/D10b)) + (c3 * (2*h2/D10b)^2) + (c4 * (2*h2/D10b)^3) ;
    kts2 = c5 + (c6 * (2*h2/D10b)) + (c7 * (2*h2/D10b)^2) + (c8 * (2*h2/D10b)^3) ;

        % Fatigue Concentration Factors----------------------------------------
    
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r6) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r6) ) ) ;

    
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

    if D10b >= .1 && D10b <= 2
        kb = 0.879*D10b^-0.107 ; 
    else
        kb = 0.91*D10b^-0.157 ;
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

    Sa_D10b = sqrt((32*Mr9*kf_1 / (pi * D10b^3))^2) ;

% for Mean Stress

    Sm_D10b = sqrt((3*(16*0*kfs_1 / (pi * D10b^3))^2)) ; % No torque on second concentration
    
% Factor of Safety via Modified Goodman
    FS_verify_f9 = ( (Sa_D10b/(Se*10^3)) + (Sm_D10b/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f9 < 2 
        if hr6 < 3.9 
            D10b = D10b + .01 ; %in
            h2 = D10b - D11 ;
            hr6 = h2/r6 ;
        else 
            D10b = D10b + .01 ;
            h2 = D10b - D11 ;
            hr6 = h2/r6 ;
        end
    else
        continue 
    end
end

% FS_verify2
% D10b
% D11

%% Keyhole Fatigue Concentration

%% Fatigue Algorithm at Key Hole F

FS_verify_f_hole_out = 0 ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for key 1: ") ;

while FS_verify_f_hole_out < 2 
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    r_key = 0.2 * D10_hole ;
    kf_hole = 1 + ( (kt_hole - 1) / (1 + a_b/sqrt(r_key) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_hole = 1 + ( (kts_hole - 1) / (1 + a_t/sqrt(r_key) ) ) ;
    
    Sa_D10_hole = sqrt((32*Mr2*kf_hole / (pi * D10_hole^3))^2) ;
    Sm_D10_hole = sqrt((3*(16*Ti*kfs_hole / (pi * D10_hole^3))^2)) ;
    
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
    if D10_hole >= .1 && D10_hole <= 2
        kb = 0.879*D10_hole^-0.107 ; 
    else
        kb = 0.91*D10_hole^-0.157 ;
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
    
    FS_verify_f_hole_out = ( (Sa_D10_hole/(Se*10^3)) + (Sm_D10_hole/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f_hole_out < 2
        D10_hole = D10_hole + .01 ;
    else 
        continue 
    end
end 

% FS_verify_f_hole_out
% D10_hole

%% Final Output

if D10 > D10b && D10 > D10_hole
    D10_final = D10 ;
elseif D10b > D10 && D10b > D10_hole 
    D10_final = D10b ;
elseif D10_hole > D10 && D10_hole > D10b
    D10_final = D10_hole ;
end


fprintf('The output shaft is stepped with three diameters. \nAll diameters have been optimized to minimize size and weight yet meet the safety factor requirements.\n')
fprintf('Diameter 9 from x = 0 to x = %5.4f is %5.4f inches. \n', x1(round(length(yr1)/2)), ceil(D9 * 40) / 40) ;
fprintf('Diameter 10 from x = %5.4f to %5.4f is %5.4f inches. \n', x1(round(length(yr1)/2)), x2((round(length(yr2)/2))), ceil(D10_final * 40) / 40) ;
fprintf('Diameter 11 from x = %5.4f to %5.4f is %5.4f inches. \n', x2(round(length(yr2)/2)), L_OS, ceil(D11 * 40) / 40) ;

