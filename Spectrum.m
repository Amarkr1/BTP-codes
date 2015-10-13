function [] = Spectrum(I_p, pump ,main_signal ,w ,h , T )  %I_p is intensity of pump in W ; pump in nm ; L is the interaction length;  w, h in um ; T in deg. C
syms sgnl idlr
% assume(sgnl,'real')
assume(sgnl>=0)
% assume(idlr,'real')
assume(idlr>=0)
%Constants
d_24 = 5e-12; % xe-12 m/V
c =  3e8 ; 
hbar  = 1.0545718e-34 ;
epsilon_0 = 8.8541e-12 ;

%Parameters
L = 1e-2; % Interaction Length in m

%Idler wavelength as a functin of signal using energy coservation
idler(sgnl) = solve(1/pump == (1/sgnl) + (1/idlr) , idlr );

main_signal;
main_idler = double(solve(1/pump == 1/main_signal + 1/idlr , idlr ));
[i_eo , n_po , n_se, n_io] = I_eo(pump, main_signal , main_idler , w , h, T); %Considering narrow bandwidth


% Lambda = Lambda_QPM (pump , main_signal , w ,h ,T );

Lambda = double((main_signal *pump)/ ((n_po*main_signal) - (n_se*pump) - (n_io*main_signal) + (n_io*pump)));
K = double(2*pi/Lambda);

n_po = double(n_po);
n_se = double(n_se);
n_io = double(n_io);

delta_k (sgnl) = (-2*pi*n_po/pump) +(2*pi*n_se/sgnl) + (2*pi*n_io/idler) +(2*pi/Lambda);
m_sgnl = double(solve(delta_k == 0 , sgnl))
snc(sgnl)= sinc(delta_k(sgnl)*L/(2*pi));
% snc(sgnl)  = sin(delta_k*L/2)/(delta_k*L/2);
spctrm(sgnl) = 16 * pi * d_24^2 * hbar * I_p * L^2 *i_eo^2 * (sinc(delta_k(sgnl)*L/(2*pi)))^2 / (n_se * n_io * sgnl^4 * idler * epsilon_0);

delta = 2;
Low = main_signal - delta;
High = main_signal + delta;
wavelength = linspace (Low,High,1000);
y = zeros(1,length(wavelength));
spectrum = zeros(1,length(wavelength));
id = zeros(1,length(wavelength));
dk = zeros(1,length(wavelength));

% value = double(snc(2000))
% value2 = double(solve(snc^2==0.2,sgnl))
% check = double(delta_k(780))


for i = 1:length(wavelength)
    id(i) = pump*wavelength(i)/(wavelength(i)-pump);
    dk(i) = ((-2*pi*n_po/pump) +(2*pi*n_se/wavelength(i)) + (2*pi*n_io/id(i)) +(2*pi/Lambda))*1e9;
    y(i) = (sinc(dk(i)*L/(2)))^2;
    spectrum(i) = 16 * (pi^3)*c * d_24^2 * hbar * I_p * L^2 *i_eo^2 *y(i) / (n_se * n_io * n_po * (wavelength(i)*1e-9)^4 * (id(i)*1e-9) * epsilon_0);
end
% id
% dk
% y
% wavelength
close all
figure
subplot(121)
% ezplot(snc^2);
plot(wavelength,y);
title('sinc^2(\Delta kL/2)');
xlabel('Signal wavelength(nm)');
ylabel('sinc^2(\Delta kL/2)');
subplot(122)
plot(wavelength,spectrum);
title('Spectrum');
xlabel('Signal wavelength(nm)');
end
