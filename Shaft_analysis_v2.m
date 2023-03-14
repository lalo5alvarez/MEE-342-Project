%% MEE 342 - Shaft Analysis v2     ~ Eduardo Alvarez

close all ;
clear ;
clc ;

%% Inputs

E = input("Enter your material's Young Modulus.(Enter purely the numerical value in base units) ") ; 
Sy = input("Enter your material's Yield Stress. (Enter purely the numerical value in base units) ") ;
Se = input("Enter your material's Endurance Limit. (Enter purely the numerical value in base units) ") ;
Sut = input("Enter your material's Ultimate Strength. (Enter purely the numerical value in base units) ") ;

n = input('Enter the desired factor of safety ') ;
% fail_theory = input('Which failure theory would you like to use? Enter 1 for Distortion Energy Theory. Enter 2 for Maximum Shear Stress Theory ') ;

L1 = input('Enter the total length of the input shaft. (Enter purely the numerical value in base units) ') ;
L2 = input('Enter the total length of the countershaft. (Enter purely the numerical value in base units) ') ;
L3 = input('Enter the total length of the output shaft. (Enter purely the numerical value in base units) ') ;

% wg1 = input('Enter the width of gears 1 and 2 (Enter purely the numerical value in base units) ') ;
% wg2 = input('Enter the width of gears 3 and 4 (Enter purely the numerical value in base units) ') ;
% wb = input('Enter the width of the bearings (Enter purely the numerical value in base units) ') ;

R1 = input('Enter the radius of the input Shaft''s Gear (Enter purely the numerical value in base units) ') ; 
F1_r = input('Enter the radial force component on the first gear. (Enter purely the numerical value in base units) ') ;
F1_t = input('Enter the tangential force component on the first gear. (Enter purely the numerical value in base units) ') ;

R2 = input('Enter the radius of the countershaft''s first Gear (Enter purely the numerical value in base units) ') ;
F2_r = -F1_r ;
F2_t = -F1_t ;

R3 = input('Enter the radius of the countershaft''s second Gear (Enter purely the numerical value in base units) ') ;
F3_r = input('Enter the radial force component on the third gear. (Enter purely the numerical value in base units) ') ;
F3_t = -F2_t * R2 / R3 ;

R4 = input('Enter the radius of the output Shaft''s Gear (Enter purely the numerical value in base units) ') ;
F4_r = -F3_r ;
F4_t = -F3_t ;

%% Calculations for Input Shaft

Lf_i = L1/2 ; % Location of Forces 
Lb1_i = L1 * (.25/11.5) ; % Location of First Bearing
Lb2_i = L1 - Lb1_i ; % Location of Second Bearing

R_ay = -F1_r / 2 ; % Reaction Force of first bearing in y direction
R_az = -F1_t / 2 ; 
R_by = R_ay ;
R_bz = R_az ;

My = -.5 * R_az * (L1 - (2*Lb1_i)) ; % Moment about y axis at gear's location
Mz = .5* R_ay * (L1 - (2*Lb1_i)) ; % Moment about z axis at gear's location
Mx = F1_t * R1 ; % Torque on Shaft

%% Input shaft Moment Diagrams

x1 = Lb1_i : L1/1000 : L1/2 ;
x2 = L1/2 : L1/1000 : Lb2_i ;

y1_y = (My / ( (L1/2) - Lb1_i)).*x1 - ( (Lb1_i .* My) / (( (L1/2) - Lb1_i))) ;
y2_y = -(My / ( Lb2_i - (L1/2) )).*x2 + ( (My * Lb2_i) / (Lb2_i - (L1/2)));

y1_z = (Mz / ( (L1/2) - Lb1_i)).*x1 - ( (Lb1_i .* Mz) / (( (L1/2) - Lb1_i))) ; 
y2_z = -(Mz / ( Lb2_i - (L1/2) )).*x2 + ( (Mz * Lb2_i) / (Lb2_i - (L1/2)));

y1_x = zeros(1,length(2 * x1)) ;
y1_x(1:length(x1)) = Mx ;
y1_x(1+length(x1) : (2*length(x1))) = 0 ;

y1_r = sqrt((y1_y.^2) + (y1_z.^2)) ;
y2_r = sqrt((y2_y.^2) + (y2_z.^2)) ;
yr = [y1_r y2_r] ;

%% Plot

figure(1) ;
subplot(4,1,1) ;
plot(x1,y1_y,'r',x2,y2_y,'r') ;
xlabel('Location on Input Shaft') ; ylabel('Moment about Y-axis') ; title('Bending Moment Diagram: Input Shaft') ;
subplot(4,1,2) ;
plot(x1,y1_z,'b',x2,y2_z,'b') ;
xlabel('Location on Input Shaft') ; ylabel('Moment about Z-axis') ; title('Bending Moment Diagram: Input Shaft') ;
subplot(4,1,3) ;
plot(x1,y1_r,'g',x2,y2_r,'g') ;
xlabel('Location on Input Shaft') ; ylabel('Resultant Moment') ; title('Bending Moment Diagram: Input Shaft') ;
subplot(4,1,4) ;
plot(x1,-y1_x(1:length(x1)),'c') ;
hold on ;
plot([L1/2 L1/2], [0 -Mx], 'c') ;
xlabel('Location on Input Shaft') ; ylabel('Internal Torque') ; title('Torque Diagram: Input Shaft') ; xlim([0 L1]) ;

%% Calculating Diameters Based on Static Loading
% In a future version of the code, this section will consider static
% loading conditions. It will inlcude an iterative process for determining
% the static concentration factors while simultaneously determining the
% diameters needed for each section. 

%% Calculating Diameters considering Fatigue & Stress Concentration
% Note: Version 2 will include equations to have the program solve the
% concentration factors and diameters simultaneously and iteratively.

k_key = input('Enter the concentration factor for the key ') ;
kf_1 = input('Enter the concentration factor for the first shoulder ') ;
kf_2 = input('Enter the concentration factor for the second shoulder ') ;
kfs_1 = input('Enter the torsional concentration factor for the first shoulder ') ;
kfs_2 = input('Enter the torsional concentration factor for the second shoulder ') ;

% Using Goodman Equation (Before First Shoulder)

A = 2*kf_1*yr( round(length(x1) / 2)) ;
B = sqrt(3) * kfs_1 * sqrt(Mx^2) ;
d1 = (16*n/pi*((A/Se) + (B/Sut)))^(1/3) ;

% b/w shoulders kf = kfs = 1

A2 = 2 * yr(length(x1)) ; 
B2 = sqrt(3) * sqrt(Mx^2) ; 
d2 = (16*n/pi*((A2/Se) + (B2/Sut)))^(1/3) ;

% After Second Shoulder 

A3 = 2*kf_2*yr( round(length(2*x1) * .75)) ;
B3 = sqrt(3) * kfs_2 * sqrt(Mx^2) ;
d3 = (16*n/pi*((A3/Se) + (B3/Sut)))^(1/3) ;

fprintf("The stepped diameters on the input shaft are as follows:\nBefore the first shoulder: %0.2f \nBetween the shoulders: %0.2f \nAfter the second shoulder: %0.2f", d1,d2,d3)

%% Diameter Output
% In a future version of the code, this section will compare the required
% diameters obtained from both static and fatigue analyses. The most
% conservative diameter will be elected.

%% Deflection Analysis
% In a future version of the code, a section will be added here analyzing
% the deflection of the shaft. It will take into account the deflection
% constraints. Deflection will be calculated via an FEA. 
%% ----------------------------------------------------- Output Shaft--------------------------------------------------------------------------------------------

%% Calculations for Output Shaft

Lf_o = L3/2 ; % Location of Forces 
Lb1_o = L3 * (.25/11.5) ; % Location of First Bearing
Lb2_o = L3 - Lb1_o ; % Location of Second Bearing

R_ay = -F4_r / 2 ; % Reaction Force of first bearing in y direction
R_az = -F4_t / 2 ; 
R_by = R_ay ;
R_bz = R_az ;

My = -.5 * R_az * (L2 - (2*Lb1_o)) ; % Moment about y axis at gear's location
Mz = .5* R_ay * (L2 - (2*Lb1_o)) ; % Moment about z axis at gear's location
Mx = F4_t * R4 ; % Torque on Shaft

%% Output Shaft Moment Diagrams

x1 = Lb1_o : L3/1000 : L3/2 ;
x2 = L3/2 : L3/1000 : Lb2_o ;

y1_y = (My / ( (L3/2) - Lb1_o)).*x1 - ( (Lb1_o .* My) / (( (L3/2) - Lb1_o))) ;
y2_y = -(My / ( Lb2_o - (L3/2) )).*x2 + ( (My * Lb2_o) / (Lb2_o - (L3/2)));

y1_z = (Mz / ( (L3/2) - Lb1_o)).*x1 - ( (Lb1_o .* Mz) / (( (L3/2) - Lb1_o))) ; 
y2_z = -(Mz / ( Lb2_o - (L3/2) )).*x2 + ( (Mz * Lb2_o) / (Lb2_o - (L3/2)));

y1_x = zeros(1,length(2 * x1)) ;
y1_x(1:length(x1)) = Mx ;
y1_x(1+length(x1) : (2*length(x1))) = 0 ;

y1_r = sqrt((y1_y.^2) + (y1_z.^2)) ;
y2_r = sqrt((y2_y.^2) + (y2_z.^2)) ;
yr = [y1_r y2_r] ;

%% Plot

figure(2) ;
subplot(4,1,1) ;
plot(x1,y1_y,'r',x2,y2_y,'r') ;
xlabel('Location on Output Shaft') ; ylabel('Moment about Y-axis') ; title('Bending Moment Diagram: Output Shaft') ;
subplot(4,1,2) ;
plot(x1,y1_z,'b',x2,y2_z,'b') ;
xlabel('Location on Output Shaft') ; ylabel('Moment about Z-axis') ; title('Bending Moment Diagram: Output Shaft') ;
subplot(4,1,3) ;
plot(x1,y1_r,'g',x2,y2_r,'g') ;
xlabel('Location on Output Shaft') ; ylabel('Resultant Moment') ; title('Bending Moment Diagram: Output Shaft') ;
subplot(4,1,4) ;
plot(x2,-y1_x(1:length(x2)),'c') ;
hold on ;
plot([L3/2 L3/2], [0 -Mx], 'c') ;
xlabel('Location on Output Shaft') ; ylabel('Internal Torque') ; title('Torque Diagram: Output Shaft') ; xlim([0 L3]) ;

%% Calculating Diameters Based on Static Loading
% In a future version of the code, this section will consider static
% loading conditions. It will inlcude an iterative process for determining
% the static concentration factors while simultaneously determining the
% diameters needed for each section.

%% Calculating Diameters considering Fatigue & Stress Concentration - 
% Note: Version 2 will include equations to have the program solve the
% concentration factors and diameters simultaneously and iteratively.

k_key = input('Enter the concentration factor for the key ') ;
kf_1 = input('Enter the concentration factor for the first shoulder ') ;
kf_2 = input('Enter the concentration factor for the second shoulder ') ;
kfs_1 = input('Enter the torsional concentration factor for the first shoulder ') ;
kfs_2 = input('Enter the torsional concentration factor for the second shoulder ') ;

% Using Goodman Equation (Before First Shoulder)

A = 2*kf_1*yr( round(length(x1) / 2)) ;
B = sqrt(3) * kfs_1 * sqrt(Mx^2) ;
d1 = (16*n/pi*((A/Se) + (B/Sut)))^(1/3) ;

% b/w shoulders kf = kfs = 1

A2 = 2 * yr(length(x1)) ; 
B2 = sqrt(3) * sqrt(Mx^2) ; 
d2 = (16*n/pi*((A2/Se) + (B2/Sut)))^(1/3) ;

% After Second Shoulder 

A3 = 2*kf_2*yr( round(length(2*x1) * .75)) ;
B3 = sqrt(3) * kfs_2 * sqrt(Mx^2) ;
d3 = (16*n/pi*((A3/Se) + (B3/Sut)))^(1/3) ;

fprintf("\nThe stepped diameters on the output shaft are as follows:\nBefore the first shoulder: %0.2f \nBetween the shoulders: %0.2f \nAfter the second shoulder: %0.2f", d1,d2,d3)

%% Diameter Output
% In a future version of the code, this section will compare the required
% diameters obtained from both static and fatigue analyses. The most
% conservative diameter will be elected.

%% Deflection Analysis
% In a future version of the code, a section will be added here analyzing
% the deflection of the shaft. It will take into account the deflection
% constraints. Deflection will be calculated via an FEA.

%% -------------------------------------------------------------Countershaft----------------------------------------------------------------------------------------
% 1.) The first step of the analysis will consist of solving the reaction
% forces on the countershaft
% 2.) The moments at both gears will be calculated in both planes
% 3.) The resultant moment will be calculated
% 4.) The Torsion will be calculated 
% 5.) All moment and torque diagrams will be plotted
% 6.) Diameters will iteratively be calculated based on static loading conditions and
% concentration factors 
% 7.) Diameters will itervatively be calculated based on fatigue and
% fatigue concentration factors
% 8.) The most conservative diameters will be elected and outputted 
% 9.) A deflection analysis will be performed via FEA
