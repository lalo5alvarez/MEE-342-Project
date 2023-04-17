%% MEE 342 - Key Hole Concentration     ~ Eduardo Alvarez

clear ; close all ; clc ;

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

%% Output Variables

save('variables5.mat','FS_verify_f_hole','D2_hole')

