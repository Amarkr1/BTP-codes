function  [neff] = calc_neff(w,h,lambda_0,String)
%INPUTS
% w,h are in micrometers
%lambda_0 is the wavelenth in nm
syms x y
%x,y give alpha_y, alpha_z
assume(y,'real')
assume(y>=0)
assume(x,'real')
assume(x>=0)

T = 25;
% T is temperature in celsius

F = (T-24.5)*(T+570.5);
if(String =='o')
    [A1,A2,A3,A4,B1,B2,B3] = sellemeir_coefficients('o');
else
    [A1,A2,A3,A4,B1,B2,B3] = sellemeir_coefficients('e');
end

n_bsq = A1 + (A2+B1*F)/(lambda_0^2 - A3+B2*F) + B3*F - A4*lambda_0^2;
n_b = sqrt(n_bsq);

delta_n = 0.0038 ;
k_0 = (2*pi)*1e9 / lambda_0 ;

%Equations
A = n_b^2;
B = - 1 / ( (k_0)^2 * w^2 ) ;
C = + 3 / ( (k_0)^2 * h^2 );
D = 8 * n_b * delta_n ;
f(x,y) = A + B * x^2 + C * y^2 + D * x * y^3 * (2*y^2 + 1)^(-3/2) * (2*x^2 + 1)^(-1/2);

S_1 = 4*k_0^2*w^2*n_b*delta_n;
S_2 = 4*k_0^2*h^2*n_b*delta_n;
eqn_1 = (S_1*S_2^2*y^2/(2*y^2+1)^5) - (S_1^0.5*S_2^1.5/(2*y^2+1)^3) - 2*S_1 ==0+0i;
alpha_z = solve(eqn_1,y);
size(alpha_z,1);

for i = 1:size(alpha_z,1)
eqn_2 = x*(2*x^2+1)^1.5 + 0i == S_1*alpha_z(i)^3/(2*alpha_z(i)^2+1)^1.5;
% alpha_z(i)
alpha_y(i) = solve(eqn_2,x);
% [solv, solu] = solve([eqn_1, eqn_2], [x, y])
end
% alpha_y
% neff_sq = zeros(1,size(alpha_z,1));
for i = 1:size(alpha_z,1)
    neff_sq(i)=f(alpha_z(i),alpha_y(i));
end
[max_neffsq,Idx] = max(neff_sq);
max_neffsq;
alpha_y(Idx);
alpha_z(Idx);
neff = sqrt(max_neffsq);

% psi(y,z) = sqrt((16*alpha_y(Idx)*alpha_z(Idx))/(pi*w*h))*alpha_z(Idx)*(z/h)*exp(-alpha_y(Idx)^2*y^2/w^2)*exp(-alpha_z(Idx)^2*z^2/h^2)
% ezsurf(psi)
% shading interp
end
