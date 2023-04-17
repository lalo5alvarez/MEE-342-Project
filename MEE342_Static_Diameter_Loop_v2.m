%% MEE 342 - Static Diameter Loop v2 ~ Eduardo Alvarez 

close all ; clear ; clc ;
load('variables.mat','Mr1','Mr2', 'Mr3','Ti', 'Rb','ns','nf') ;

%% Inputs 

n = 2 ;
Sy = 50000 ; % psi
Sut = 68000 ; % psi

%% Algorithm 

D1 = (n*sqrt((32*Mr1)^2 + (3*(16*0)^2)) / (Sy* pi) )^(1/3) % initial guesses based off Von Mises
                                                           % No T on D1
D2 = (n*sqrt((32*Mr2)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3)
D2_static = D2 ; % duplicate for preservation of Static Diameter ;
D3 = (n*sqrt((32*Mr3)^2 + (3*(16*Ti)^2)) / (Sy* pi) )^(1/3)
h1 = D2-D1 ; 

fprintf("Choose a fillet radius value between %5.3f and %5.3f \n (not including the bounds) \n", h1/4, h1/.24)
r1 = input('Enter Fillet Radius in Inches: ') 
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

D1 
if D2 > D2b 
    D2
elseif D2b > D2 
    D2b 
elseif D2_hole
    D2_hole
end
D3
FS_verify
FS_verify2

%% Output Variables

save('variables2.mat','kt1', 'kts1' ,'kt2','kts2', 'kt_hole','kts_hole','r1', 'D1','D2','D2b', 'D2_hole', 'D3', 'Sut')
