function [n_b] = Sellemeir (lambda,T,String); 
%INPUT
% lambda is wavelength in nanometers
% T is temperature in celsius
%-----------------
if (String == 'o')
    %ordinary Parameters:
    A1 = 4.582;
    A2 = 9.921e4;
    A3 = 2.109e2;
    A4 = 2.192e-8;
    B1 = 5.2716e-2;
    B2 = -4.9143e-5;
    B3 = 2.2971e-7;
   %----------------------
else
    %extraordinary Parameters:
    A1 = 4.9048;
    A2 = 1.1775e5;
    A3 = 2.1802e2;
    A4 = 2.7153e-8;
    B1 = 2.2314e-2;
    B2 = -2.9671e-5;
    B3 = 2.1429e-8;
    %----------------------
end
F = (T-24.5)*(T+570.5);
n_sq = A1 + (A2+B1*F)/(lambda^2 - A3+B2*F) + B3*F - A4*lambda^2;
n_b = sqrt(n_sq);
end
