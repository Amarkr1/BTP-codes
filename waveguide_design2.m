function [w,h,difference,plt_neff_p,plt_neff_s,plt_neff_i]  = waveguide_design2 (lambda_p,sigma,lambda_s,L)
close all
syms wave_pump wave_signal

assume(wave_pump,'real')
assume(wave_pump>0)
assume(wave_signal,'real')
assume(wave_signal>0)
%INPUT
% lambda_p = pump wavelength
% lambda_s = desired signal wavelength
% signma is pump bandwidth

lhs = 1/(sigma*sigma);
lambda_i = (lambda_p*lambda_s)/(lambda_s - lambda_p);
num = 4;
% w = linspace(3e-6,6e-6,num);
% h = linspace(4e-6,6e-6,num);
% ll = 6;
% ul = 16;
% w = linspace(ll*1e-6,ul*1e-6,num);
% h = linspace(ll*1e-6,ul*1e-6,num);
w = linspace(6.5*1e-6,7.5*1e-6,num);
h = linspace(7.5*1e-6,8.5*1e-6,num);
width_vector = size(w);
% width_vector(2)
height_vector = size(h);
gamma = 0.04822;
c = 3e8;
plt_neff_p = zeros(width_vector(2),height_vector(2));
plt_neff_s = zeros(width_vector(2),height_vector(2));
plt_neff_i = zeros(width_vector(2),height_vector(2));
difference = zeros(width_vector(2),height_vector(2));
for i=1:width_vector(2)
    for j= 1:height_vector(2)
       % i,j
        %pump -> ordinary ; signal -> extraordinary;   idler->ordinary
        neff_p = calc_neff(w(i),h(j),lambda_p,'o');
        neff_s = calc_neff(w(i),h(j),lambda_s,'e');
        neff_i = calc_neff(w(i),h(j),lambda_i,'o');
            
        %these variables are to plot
        plt_neff_p(i,j) = neff_p;
        plt_neff_s(i,j) = neff_s;
        plt_neff_i(i,j) = neff_i;        
%         plt_neff_s(i,j) = calc_neff(w(i),h(j),lambda_s,'e');
%         plt_neff_i(i,j) = calc_neff(w(i),h(j),lambda_i,'o');
        
        
        first_term_rhs = (neff_p-(lambda_p*num_differentiation(w(i),h(j),lambda_p,'o')))-(neff_s-(lambda_s*num_differentiation(w(i),h(j),lambda_s,'e')));
        second_term_rhs = (neff_p-(lambda_p*num_differentiation(w(i),h(j),lambda_p,'o')))-(neff_i-(lambda_i*num_differentiation(w(i),h(j),lambda_i,'o')));
        rhs = double ( (-gamma*(L^2)/c)*first_term_rhs*second_term_rhs);
        difference(i,j) = abs(lhs-rhs);
    end
%     i

end
% surf(w,h,plt_neff_p)
% difference
% w
% h
% surf(difference)
