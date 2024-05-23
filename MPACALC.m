%MATLAB Microstrip Patch Antenna Calculator
%Formulas are from Antenna Theory Analysis and Design by C. Balanis

%Constant Parameters
clear all;
clc;
c = 299792458; % Speed of Light
Zf = 50; % Characteristic Impedance
fr = input("Input Frequency in (GHz):"); %Frequency
fr_2 = fr*10^9; %Resonant Frequency
th_ck = input("Input Thickness in mm:"); %Substrate Thickness
th_ck2 = th_ck*10^-3;
e_r = input("Input Dielectric Constant:"); %Dielectric Constant

W = (c/(2*fr_2))*sqrt(2/(e_r+1)); %Patch Width
Ereff=((e_r+1)/2)+((e_r-1)/2)*(1/sqrt(1+(12*(th_ck2/W)))); %Effective Dielectric Constant
L_d = th_ck2*0.412*(( (Ereff + 0.3)*((W/th_ck2)+0.264))) /( (Ereff - 0.258)*((W/th_ck2)+0.8)); %Delta_L
L = (c/(2 * fr_2 * sqrt(Ereff))) - (2*(L_d)); %Actual Length
Wg = (6*10^-3)+ W; %Width of Ground and Substrate
Lg = (6*10^-3) + L; %Length of Ground and Substrate
lambda_0 = c/fr_2; %Wavelength 
k_0 = (2*pi)/lambda_0;
x_0 = k_0 * W;
i_1=(-2)+cos(x_0)+(x_0 *sinint(x_0))+(sin(x_0)/x_0);
    
G_1=i_1/(120*pi^2); %Self conductance
fun_G_12 = @(th) ((sin(((k_0.*W)./2).*cos(th))./cos(th)).^2) .* (besselj(0,(k_0 .* L .*sin(th))).*sin(th).^3);
G_12 = (1/(120*pi^2))*integral(fun_G_12,0,pi); %Mutual Conductance
R_in = 1/(2*(G_1+G_12)); %Input Impedance 
y_0 = (L/pi)*acos(sqrt(Zf/R_in)); %Inset Depth 


syms x
eqn = (120*pi)/(sqrt(Ereff)*((x/th_ck2)+1.393+0.667*log((x/th_ck2)+1.444))) == 50;
sol = solve(eqn, x); %Feed Width 1
B = (377*pi)/(2*Zf*sqrt(e_r));
w_f = ((2*th_ck2)/pi)*( B-1-log(2*B-1)+(((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r))));

fr_22 = ['Frequency: ', num2str(fr_2), ' MHz'];
th_ck22 = ['Substrate Thickness: ', num2str(th_ck2*10^3), ' mm'];
e_r2 = ['Dielectric Constant: ', num2str(e_r)];
W_2 = ['Patch Width: ', num2str(W*10^3), ' mm'];
L_2 = ['Patch Length: ', num2str(L*10^3), ' mm'];
Wg_2 = ['Width of Ground and Substrate: ' num2str(Wg*10^3), ' mm'];
Lg_2 = ['Length of Ground and Substrate: ' num2str(Lg*10^3), ' mm'];
G_11 = ['Self Conductance: ', num2str(G_1)];
G_122 = ['Mutual Conductance: ', num2str(G_12)];
R_in1 = ['Input Impedance: ', num2str(R_in), ' ohms'];
y_01 = ['Inset Depth: ' num2str(y_0*10^3), ' mm'];
sol_2 = ['Feed Width: ', num2str(double(sol)*10^3), ' mm, Balanis'];
w_f2 = ['Feed Width2: ', num2str(w_f*10^3), ' mm, Pozar'];
disp(fr_22);
disp(th_ck22);
disp(e_r2);
disp(W_2);
disp(L_2);
disp(Wg_2);
disp(Lg_2);
disp(G_11);
disp(G_122);
disp(R_in1);
disp(y_01);
disp(sol_2);
disp(w_f2);