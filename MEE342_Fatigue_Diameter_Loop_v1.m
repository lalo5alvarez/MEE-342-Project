%% MEE 342 - Fatigue Diameter Loop v1     ~ Eduardo Alvarez

close all ; clear ; clc ;
load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;
load('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'r1' , 'D1','D2','D2b','D3')
load('variables3.mat','Sut') ; 
fprintf("Type in one of the following values for the corresponding finish \nGround = 1 \nMachined/Cold-drawn = 2 \nHot-rolled = 3 \nAs-forged = 4 \n")
A = input("Select from list above ") ;
fprintf("\n Input operating temperature in degrees fahrenheit \n")
Tf = input("Temperature = ") ;
fprintf("\nTorsion = 1 \nBending = 2 \nAxial = 3 \n")
B = input("Select from list above ") ;

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


%     if .1 <= D2 <= 2
%         kb2 = 0.879*D2^-0.107 ; 
%     else
%         kb2 = 0.91*D2^-0.157 ;
%     end
% 
% 
%     if .1 <= D2b <= 2
%         kb2b = 0.879*D2b^-0.107 ; 
%     else
%         kb2b = 0.91*D2b^-0.157 ;
%     end
% 
% 
%     if .1 <= D3 <= 2
%         kb3 = 0.879*D3^-0.107 ; 
%     else
%         kb3 = 0.91*D3^-0.157 ;
%     end
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

FS_verify_f
hr1
D2
D1
