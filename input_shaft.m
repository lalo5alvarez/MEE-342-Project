
%% MEE 342 - Shaft Analysis v3     

close all ; clear ; clc ;

%Input Shaft---------------------------------------------------------------

%% Inputs 

ns = 2 ;% static F.S. 
nf = 2 ;% fatigue F.S.

a = 7.5 ;
L_OS = 0;
while L_OS <= a 

fprintf("Please give a shaft length larger than %5.3f inches \n",a)
L_OS = input ("Input Shaft Length = ") ;

if L_OS <= a 
   fprintf("\nThat input value will not work \n") 
end
end
fprintf("\n");

Lt_i = L_OS ;
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
B = bz / cos(20*pi/180) ; % 
by = B * sin(20*pi/180) ;
cy = ( (dy*(Lt_i-La))-(by*(Lb-La)) ) / (Lc-La) ; 
ay = -dy - by - cy ;
cz = (-bz*(Lb-La)) / (Lc) ; 
az = -bz - cz ;

%% Shear & Torque Diagrams

figure(1) ;
subplot(3,1,1) ;
plot([La La],[0 ay],'k', [La Lb], [ay ay], 'k') ;
hold on ;
plot([Lb Lb],[ay (ay+by)],'k', [Lb Lc], [(ay+by) (ay+by)], 'k') ;
hold on ;
plot([Lc Lc],[(ay+by) (ay+by+cy)],'k', [Lc Lt_i], [(ay+by+cy) (ay+by+cy)], 'k') ;
hold on ;
plot([Lt_i Lt_i],[(ay+by+cy) (ay+by+cy+dy)],'k') ;
xlabel('X Distance [in]') ; ylabel('Shear Force [lb]') ; title('Shear Force along the Y-axis [Fy]') ;
hold off ;

subplot(3,1,2) ;
plot([La La],[0 az],'k', [La Lb], [az az], 'k') ;
hold on ;
plot([Lb Lb],[az (az+bz)],'k', [Lb Lc], [(az+bz) (az+bz)], 'k') ;
hold on ;
plot([Lc Lc],[(az+bz) (az+bz+cz)],'k', [Lc Lt_i], [(az+bz+cz) (az+bz+cz)], 'k') ;
hold on ;
plot([Lt_i Lt_i],[(az+bz+cz) (az+bz+cz+0)],'k') ;
xlabel('X Distance [in]') ; ylabel('Shear Force [lb]') ; title('Shear Force along the Z-axis [Fz]') ;
hold off ;

% Torque Diagram 

Ti = -Tb ; % Internal Torque
subplot(3,1,3) ;
plot([0 Lb],[0 0], 'k') ;
hold on ;
plot([Lb Lb],[0 Ti], 'k') ;
plot([Lb Lt_i],[Ti Ti],'k') ;
xlabel('X Distance [in]') ; ylabel('Torque [lb-in]') ; title('Torque Diagram (Mx)') ;

%% Moment Equations

% My Diagran
x1 = La:Lt_i/10000:Lb ;
x2 = Lb:Lt_i/10000:Lc ;
Myb = ay*(Lb-La) ;
y1 = (((Myb)./(Lb - La)).*x1) - (((Myb)./(Lb - La)).*La) ; %piece 1
y2 = -(((Myb)./(Lc - Lb)).*x2) + (((Myb)./(Lc - Lb)).*Lc) ; %piece 2  
figure(2) ;
subplot(3,1,1) ;
plot(x1,y1,'k',x2,y2,'k') ;
xlim([0 Lt_i])
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (My)') ;

% Mz Diagram
Mzb = -az*(Lb-La) ; 
Mzc = dy*(4+La) ;
y3 = (((Mzb)./(Lb - La)).*x1) - (((Mzb)./(Lb - La)).*La) ; %piece 1
y4 = (((Mzc-Mzb)/(Lc-Lb)).*x2) + Mzb - (((Mzc-Mzb)/(Lc-Lb))*Lb) ; %piece 2
x3 = Lc:Lt_i/1000:Lt_i ;
y5 = ((-Mzc/(Lt_i-Lc)).*x3) + ((Mzc/(Lt_i-Lc)).* Lt_i) ; %piece 3
subplot(3,1,2) ;
plot(x1,y3,'k',x2,y4,'k') ;
hold on ;
plot(x3,y5,'k') ;
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Moment Diagram (Mz)') ;

% Mr Diagram

yr1 = sqrt((y1).^2 + (y3).^2) ;
yr2 = sqrt((y2).^2 + (y4).^2) ;
yr3 = sqrt(y5.^2) ;
subplot(3,1,3) ;
plot(x1,yr1,'k',x2,yr2,'k') ; 
hold on ;
plot(x3,yr3,'k') ;
hold on ;
plot([x2(length(x2)) x2(length(x2))], [yr2(length(yr2)) yr3(1)],'k')
xlabel('X Distance [in]') ; ylabel('Bending Moment [lb-in]') ; title('Resultant Moment Diagram (Mr)') ;

%% Choosing Moments

Mr1 = yr1(round(length(yr1)/2)) ; 
Mr2 = yr1(length(yr1)) ; 
Mr3 = yr3(round(length(yr3)/2)) ; 

%% Output Variables
save('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf', 'Lb', 'x1', 'x2', 'x3', 'Lt_i','Tb') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% MEE 342 - Static Diameter Loop v2 

clear ; clc ;
load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;

%% Inputs 

n = 2 ;
Sy = 50000 ; % psi
Sut = 68000 ; % psi

%% Algorithm 

D1 = (n*sqrt((32*Mr1)^2 + (3*(16*0)^2)) / (Sy* pi) )^(1/3) ; % initial guesses based off Von Mises
                                                           % No T on D1
D2 = (n*sqrt((32*Mr2)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3) ;
D2_static = D2 ; % duplicate for preservation of Static Diameter ;
D3 = (n*sqrt((32*Mr3)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3) ;
h1 = D2-D1 ; 

fprintf("The following diamaters are inital estimates \nin order to run an iterative process \nD1 = %5.3f inches\nD2 = %5.3f inches\nD3 = %5.3f inches \n",D1,D2,D3)

fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.25)
r1 = input('Enter Fillet Radius in Inches: ') ;
hr1 = h1/r1 ;
FS_verify = 0 ;

while FS_verify < 2 % Diameters 1 & 2
    
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r1)^(1/2)) - (0.131 * h1/r1);
         c2 = 0.022 - (3.405 * (h1/r1)^(1/2)) + (0.915 * h1/r1);
         c3 = 0.869 + (1.777 * (h1/r1)^(1/2)) - (0.555 * h1/r1);
         c4 = -.810 + (.422 * (h1/r1)^(1/2)) - (0.260 * h1/r1);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r1)^(1/2)) - (0.008 * h1/r1);
         c2 = -3.813 + (.968 * (h1/r1)^(1/2)) - (0.260 * h1/r1);
         c3 = 7.423 - (4.868 * (h1/r1)^(1/2)) + (0.869 * h1/r1);
         c4 = -3.839 + (3.070 * (h1/r1)^(1/2)) - (0.6 * h1/r1);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
         %close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r1)^(1/2)) - (0.075 * h1/r1);
         c6 = -0.437 - (1.969 * (h1/r1)^(1/2)) + (0.553 * h1/r1);
         c7 = 1.557 + (1.073 * (h1/r1)^(1/2)) - (0.578 * h1/r1);
         c8 = -1.061 + (.171 * (h1/r1)^(1/2)) + (0.086 * h1/r1);
    else 
         fprintf('Your fillet radius is invalid. Hr1 ratio is %4.3f \n',hr1)
%          close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D2)) + (c3 * (2*h1/D2)^2) + (c4 * (2*h1/D2)^3) ;
    kts1 = c5 + (c6 * (2*h1/D2)) + (c7 * (2*h1/D2)^2) + (c8 * (2*h1/D2)^3) ;
    FS_verify = Sy / sqrt((32*Mr1*kt1 / (pi * D1^3))^2 + (3*(16*0*kts1 / (pi * D1^3))^2)) ; % Von Mises
    
    if FS_verify < 2 % Iterating for new diameter
            % If loop keeps hr1 ratio within equation parameters.
        if hr1 < 3.9 
            D2 = D2 + .01 ; %in
            h1 = D2-D1 ;
            hr1 = h1/r1 ;
        else
            D1 = D1 + .01 ; %in
            h1 = D2-D1 ;
            hr1 = h1/r1 ;
        end
        
    else
        continue 
    end
end

% FS_verify
% D1
% D2

%% -------------------------------------------------
D2b = D2 ;
FS_verify2 = 0 ;
h2 = D2b - D3 ;
hr2 = h2/r1 ;
while FS_verify2 < 2 % Diameters 2 & 3
    
    % Bending kt values
    if hr2 >= .1 && hr2 <= 2 
         c1 = 0.947 + (1.206 * (h2/r1)^(1/2)) - (0.131 * h2/r1);
         c2 = 0.022 - (3.405 * (h2/r1)^(1/2)) + (0.915 * h2/r1);
         c3 = 0.869 + (1.777 * (h2/r1)^(1/2)) - (0.555 * h2/r1);
         c4 = -.810 + (.422 * (h2/r1)^(1/2)) - (0.260 * h2/r1);
    
    elseif hr2 >= 2 && hr2 <= 20
         c1 = 1.232 + (.832 * (h2/r1)^(1/2)) - (0.008 * h2/r1);
         c2 = -3.813 + (.968 * (h2/r1)^(1/2)) - (0.260 * h2/r1);
         c3 = 7.423 - (4.868 * (h2/r1)^(1/2)) + (0.869 * h2/r1);
         c4 = -3.839 + (3.070 * (h2/r1)^(1/2)) - (0.6 * h2/r1);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end
    
        % Torsion kts values
    if hr2 >= .25 && hr2 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r1)^(1/2)) - (0.075 * h2/r1);
         c6 = -0.437 - (1.969 * (h2/r1)^(1/2)) + (0.553 * h2/r1);
         c7 = 1.557 + (1.073 * (h2/r1)^(1/2)) - (0.578 * h2/r1);
         c8 = -1.061 + (.171 * (h2/r1)^(1/2)) + (0.086 * h2/r1);
    else 
         fprintf('Your fillet radius is invalid. Hr2 ratio is %4.3f \n',hr2)
         %close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h2/D2b)) + (c3 * (2*h2/D2b)^2) + (c4 * (2*h2/D2b)^3) ;
    kts2 = c5 + (c6 * (2*h2/D2b)) + (c7 * (2*h2/D2b)^2) + (c8 * (2*h2/D2b)^3) ;
    FS_verify2 = Sy / sqrt((32*Mr3*kt2 / (pi * D3^3))^2 + (3*(16*Ti*kts2 / (pi * D3^3))^2)) ;  % Von Mises on D3
    
    if FS_verify2 < 2 
        if hr2 < 3.9 
            D2b = D2b + .01 ; %in
            h2 = D2b-D3 ;
            hr2 = h2/r1 ;
        else 
            D3 = D3 + .01 ;
            h2 = D2b-D3 ;
            hr2 = h2/r1 ;
        end
    else
        continue 
    end
end

% FS_verify2
% D2b
% D3

%% Keyhole Concentration 

kt_hole = 2.14 ; % Approximation from Table 7-1 
kts_hole = 3 ;
D2_hole = ((n/(Sy*pi)) * sqrt( (32*Mr2*kt_hole)^2 + (3*((16*Ti*kts_hole)^2)) ))^(1/3) ;

%% Display 

% D1 
% if D2 > D2b 
%     D2
% elseif D2b > D2 
%     D2b 
% elseif D2_hole
%     D2_hole
% end
% D3
% FS_verify
% FS_verify2

%% Output Variables

save('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut','Sy')

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% MEE 342 - Fatigue Diameter Loop v1   

clear ; clc ;
load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;
load('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut')
fprintf("The following inputs will correspond to the first stress concentration area\n")
fprintf("Type in one of the following values for the corresponding finish \nGround = 1 \nMachined/Cold-drawn = 2 \nHot-rolled = 3 \nAs-forged = 4 \n")
A = input("Select from list above: ") ;
fprintf("\n Input operating temperature in degrees fahrenheit \n")
Tf = input("Temperature = ") ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentation 1: ") ;

%% Fatigue Algorithm Check at Concentration 1
FS_verify_f = 0 ;
h1 = D2-D1 ;
hr1 = h1/r1 ;

while FS_verify_f < 2 
        
    % Bending kt values
    if hr1 >= .1 && hr1 <= 2 
         c1 = 0.947 + (1.206 * (h1/r1)^(1/2)) - (0.131 * h1/r1);
         c2 = 0.022 - (3.405 * (h1/r1)^(1/2)) + (0.915 * h1/r1);
         c3 = 0.869 + (1.777 * (h1/r1)^(1/2)) - (0.555 * h1/r1);
         c4 = -.810 + (.422 * (h1/r1)^(1/2)) - (0.260 * h1/r1);
    
    elseif hr1 >= 2 && hr1 <= 20
         c1 = 1.232 + (.832 * (h1/r1)^(1/2)) - (0.008 * h1/r1);
         c2 = -3.813 + (.968 * (h1/r1)^(1/2)) - (0.260 * h1/r1);
         c3 = 7.423 - (4.868 * (h1/r1)^(1/2)) + (0.869 * h1/r1);
         c4 = -3.839 + (3.070 * (h1/r1)^(1/2)) - (0.6 * h1/r1);
    else 
         fprintf('Your fillet radius is invalid')
         close all ; clear ;
    end
    
    % Torsion kts values
    if hr1 >= .25 && hr1 <= 4.0 
         c5 = 0.905 + (.783 * (h1/r1)^(1/2)) - (0.075 * h1/r1);
         c6 = -0.437 - (1.969 * (h1/r1)^(1/2)) + (0.553 * h1/r1);
         c7 = 1.557 + (1.073 * (h1/r1)^(1/2)) - (0.578 * h1/r1);
         c8 = -1.061 + (.171 * (h1/r1)^(1/2)) + (0.086 * h1/r1);
    else 
         fprintf('Your fillet radius is invalid')
         close all ; clear ;
    end

    kt1 = c1 + (c2 * (2*h1/D2)) + (c3 * (2*h1/D2)^2) + (c4 * (2*h1/D2)^3) ;
    kts1 = c5 + (c6 * (2*h1/D2)) + (c7 * (2*h1/D2)^2) + (c8 * (2*h1/D2)^3) ;
    
    % Fatigue Concentration Factors
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_1 = 1 + ( (kt1 - 1) / (1 + a_b/sqrt(r1) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_1 = 1 + ( (kts1 - 1) / (1 + a_t/sqrt(r1) ) ) ;

    % Endurance Limit
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for kb

    if D1 >= .1 && D1 <= 2
        %.1 <= D1 <= 2
        kb = 0.879*D1^-0.107 ; 
    else
        kb = 0.91*D1^-0.157 ;
    end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for kc

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
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%Goodman 

% for Stress Amplitude

    Sa_D1 = sqrt((32*Mr1*kf_1 / (pi * D1^3))^2) ;

% for Mean Stress

    Sm_D1 = sqrt((3*(16*0*kfs_1 / (pi * D1^3))^2)) ; %No Torque at first concentration
    
% Factor of Safety via Modified Goodman
    FS_verify_f = ( (Sa_D1/(Se*10^3)) + (Sm_D1/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f < 2 % Iterating for new diameter
        if hr1 < 3.9
            D2 = D2 + .01 ; %in
            h1 = D2-D1 ;
            hr1 = h1/r1 ;
        else
            D1 = D1 + .01 ; %in
            h1 = D2-D1 ;
            hr1 = h1/r1 ;
        end
    else
        continue 
    end
end

% FS_verify_f
% hr1
% D2
% D1

%% Output Variables

save('variables3.mat','FS_verify_f', 'hr1', 'D2', 'D1') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% MEE 342 - Fatigue Diameter Loop v1    

clear ; clc ;
load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;
load('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut') ; 
fprintf("The following inputs will correspond to the second stress concentration area\n")
fprintf("Type in one of the following values for the corresponding finish \nGround = 1 \nMachined/Cold-drawn = 2 \nHot-rolled = 3 \nAs-forged = 4 \n")
A = input("Select from list above: ") ;
fprintf("\n Input operating temperature in degrees fahrenheit \n")
Tf = input("Temperature = ") ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above for concentration 2: ") ;

%% Fatigue Algorithm Check at Concentration 2

FS_verify_f2 = 0 ;
h2 = D2b-D3 ;
hr2 = h2/r1 ;

while FS_verify_f2 < 2 
        
    % Bending kt values
    if hr2 >= .1 && hr2 <= 2 
         c1 = 0.947 + (1.206 * (h2/r1)^(1/2)) - (0.131 * h2/r1);
         c2 = 0.022 - (3.405 * (h2/r1)^(1/2)) + (0.915 * h2/r1);
         c3 = 0.869 + (1.777 * (h2/r1)^(1/2)) - (0.555 * h2/r1);
         c4 = -.810 + (.422 * (h2/r1)^(1/2)) - (0.260 * h2/r1);
    
    elseif hr2 >= 2 && hr2 <= 20
         c1 = 1.232 + (.832 * (h2/r1)^(1/2)) - (0.008 * h2/r1);
         c2 = -3.813 + (.968 * (h2/r1)^(1/2)) - (0.260 * h2/r1);
         c3 = 7.423 - (4.868 * (h2/r1)^(1/2)) + (0.869 * h2/r1);
         c4 = -3.839 + (3.070 * (h2/r1)^(1/2)) - (0.6 * h2/r1);
    else 
         fprintf('Your fillet radius is invalid')
         close all ; clear ;
    end
    
    % Torsion kts values
    if hr2 >= .25 && hr2 <= 4.0 
         c5 = 0.905 + (.783 * (h2/r1)^(1/2)) - (0.075 * h2/r1);
         c6 = -0.437 - (1.969 * (h2/r1)^(1/2)) + (0.553 * h2/r1);
         c7 = 1.557 + (1.073 * (h2/r1)^(1/2)) - (0.578 * h2/r1);
         c8 = -1.061 + (.171 * (h2/r1)^(1/2)) + (0.086 * h2/r1);
    else 
         fprintf('Your fillet radius is invalid')
         close all ; clear ;
    end

    kt2 = c1 + (c2 * (2*h2/D2b)) + (c3 * (2*h2/D2b)^2) + (c4 * (2*h2/D2b)^3) ;
    kts2 = c5 + (c6 * (2*h2/D2b)) + (c7 * (2*h2/D2b)^2) + (c8 * (2*h2/D2b)^3) ;
    
    % Fatigue Concentration Factors
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    kf_2 = 1 + ( (kt2 - 1) / (1 + a_b/sqrt(r1) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_2 = 1 + ( (kts2 - 1) / (1 + a_t/sqrt(r1) ) ) ;

    % Endurance Limit
    if Sut <= 200  % in kpsi
        Se1 = .5*Sut ;
    else
        Se1 = 100 ;
    end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for kb

    if D3 >= .1 && D3 <= 2
        %.1 <= D1 <= 2
        kb = 0.879*D3^-0.107 ; 
    else
        kb = 0.91*D3^-0.157 ;
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% for kc

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
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%Goodman 

% for Stress Amplitude

    Sa_D3 = sqrt((32*Mr3*kf_2 / (pi * D3^3))^2) ;

% for Mean Stress

    Sm_D3 = sqrt((sqrt(3))*(16*Ti*kfs_2 / (pi * D3^3))) ;
    
% Factor of Safety via Modified Goodman
    FS_verify_f2 = ( (Sa_D3/(Se*10^3)) + (Sm_D3/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f2 < 2 % Iterating for new diameter
        if hr2 < 3.9
            D2b = D2b + .01 ; %in
            h2 = D2b-D3 ;
            hr2 = h2/r1 ;
        else
            D3 = D3 + .01 ; %in
            h2 = D2b-D3 ;
            hr2 = h2/r1 ;
        end
    else
        continue 
    end
end
% hr2
% D2b 
% D3
% FS_verify_f2

%% Output Variables

save('variables4.mat','FS_verify_f2', 'hr2', 'D2b', 'D3', 'A','B','Se1','Tf') ;

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% MEE 342 - Key Hole Concentration   

clear ; clc ;

load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;
load('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut') ;
load('variables3.mat','FS_verify_f', 'hr1', 'D2', 'D1') ;
load('variables4.mat','FS_verify_f2', 'hr2', 'D2b', 'D3', 'A','B','Se1','Tf') ; 

%% Fatigue Algorithm at Key Hole

FS_verify_f_hole = 0 ;

while FS_verify_f_hole < 2 
    a_b = 0.246 - 3.08*(10^-3)*Sut + 1.51*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % axial/bending notch sensitivity 
    r_key = 0.2 * D2_hole ;
    kf_hole = 1 + ( (kt_hole - 1) / (1 + a_b/sqrt(r_key) ) ) ;
    a_t = .19 - 2.5*(10^-3)*Sut + 1.35*(10^-5)*Sut^2 - 2.67*(10^-8)*Sut^3 ; % torsion notch sensitivity
    kfs_hole = 1 + ( (kts_hole - 1) / (1 + a_t/sqrt(r_key) ) ) ;
    
    Sa_D2_hole = sqrt((32*Mr2*kf_hole / (pi * D2_hole^3))^2) ;
    Sm_D2_hole = sqrt((3*(16*Ti*kfs_hole / (pi * D2_hole^3))^2)) ;
    
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
    if D2_hole >= .1 && D2_hole <= 2
        kb = 0.879*D2_hole^-0.107 ; 
    else
        kb = 0.91*D2_hole^-0.157 ;
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
    
    FS_verify_f_hole = ( (Sa_D2_hole/(Se*10^3)) + (Sm_D2_hole/(Sut*10^3)) )^-1 ; %converted to psi
    
    if FS_verify_f_hole < 2
        D2_hole = D2_hole + .01 ;
    else 
        continue 
    end
end 

% FS_verify_f_hole 
% D2_hole

fprintf('Calculation complete. Please Proceed. \n') ;

%% Output Variables

save('variables5.mat','FS_verify_f_hole','D2_hole')

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% MEE 342 - Diameter Outputs  

clear ; clc ;

load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf', 'Lb', 'x1', 'x2', 'x3', 'Lt_i','Tb') ;
load('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut') ;
load('variables3.mat','FS_verify_f', 'hr1', 'D2', 'D1') ;
load('variables4.mat','FS_verify_f2', 'hr2', 'D2b', 'D3', 'A','B','Se1','Tf') ; 
load('variables5.mat','FS_verify_f_hole','D2_hole') ;

if D2 > D2b && D2 > D2_hole
    D2_final = D2 ;
elseif D2b > D2 && D2b > D2_hole
    D2_final = D2b ;
elseif D2_hole > D2 && D2_hole > D2b
    D2_final = D2_hole ;
end

fprintf('The input shaft is stepped with three diameters. \nAll diameters have been optimized to minimize size and weight yet meet the safety factor requirements.\n')
fprintf('Diameter 1 from x = 0 to x = %5.4f is %5.4f inches. \n', x1(round(length(x1)/2)), ceil(D1 * 40) / 40) ;
fprintf('Diameter 2 from x = %5.4f to %5.4f is %5.4f inches. \n', x1(round(length(x1)/2)), x3(round(length(x3)/2)), ceil(D2_final * 40) / 40) ;
fprintf('Diameter 3 from x = %5.4f to %5.4f is %5.4f inches. \n', x3(round(length(x3)/2)), Lt_i, ceil(D3 * 40) / 40) ;

%% Output Variables

save('InputShaft.mat') ;


