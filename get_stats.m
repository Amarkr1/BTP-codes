function [ n , alpha_y, alpha_z ] = give_stats(lambda,w,h,T,String)%w,h in um and lambda in nm
syms x y z
assume(y,'real')
assume(y>=0)
assume(x,'real')
assume(x>=0)

%INPUTS
w = w*1e-6 ;
h = h*1e-6;
lambda_0 = lambda; 
n_b = Sellemeir(lambda_0,T,String);  % Sellemeir (lambda,T,String)
delta_n = 0.0038 ;
k_0 = (2*pi)*1e9 / lambda_0 ;


%Writing function f ,i.e., (n_eff)^2
A = n_b^2;
B = - 1 / ( (k_0)^2 * w^2 ) ;
C = + 3 / ( (k_0)^2 * h^2 );
D = 8 * n_b * delta_n ;
f(x,y) = A + B * x^2 + C * y^2 + D * x * y^3 * (2*y^2 + 1)^(-3/2) * (2*x^2 + 1)^(-1/2);

%Calculating maximum of f
S_1 = 4*k_0^2*w^2*n_b*delta_n;
S_2 = 4*k_0^2*h^2*n_b*delta_n;
eqn_1 = (S_1*S_2^2*y^2/(2*y^2+1)^5) - (S_1^0.5*S_2^1.5/(2*y^2+1)^3) - 2*S_1 ==0+0i;
alpha_z = solve(eqn_1,y);

for i = 1:size(alpha_z,1)
    eqn_2 = x*(2*x^2+1)^1.5 + 0i == S_1*alpha_z(i)^3/(2*alpha_z(i)^2+1)^1.5;
    alpha_y(i) = solve(eqn_2,x);
  
end


for i = 1:size(alpha_z,1)
    neff_sq(i)=f(alpha_z(i),alpha_y(i));
end

[max_neffsq,Idx] = max(neff_sq);
n = sqrt( max_neffsq);
alpha_y = alpha_y(Idx);
alpha_z = alpha_z(Idx);

end
