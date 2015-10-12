function [] = Spectrum(I_p, pump ,main_signal ,w ,h , T )  %I_p is intensity of pump ; pump in nm ; L is the interaction length;  w, h in um ; T in deg. C
syms sgnl idlr

%Constants
d_24 = 5e-12; % xe-12 m/V
c =  3e8 ; 
hbar  = 1.0545718e-34 ;
epsilon_0 = 8.8541e-12 ;

%Parameters
L = 1; % Interaction Length

%Idler wavelength as a functin of signal using energy coservation
idler(sgnl) = solve(1/pump == 1/sgnl + 1/idlr , idlr ); 

%QPM
Lambda = Lambda_QPM (pump , main_signal , w ,h ,T );
K = 2*pi/Lambda;


main_idler = solve(1/pump == 1/main_signal + 1/idlr , idlr ); 

%I_eo
[i_eo , n_po , n_se, n_io] = I_eo(pump, main_signal , main_idler , w , h, T); %Considering narrow bandwidth

delta_k (sgnl) = 2*pi*n_po/pump -(2*pi*n_se/sgnl + 2*pi*n_io/idler +K) ;

snc = sin(delta_k*L/2)/(delta_k*L/2);
spctrm = 16 * pi * d_24^2 * hbar * I_p * L^2 *i_eo^2 * (snc)^2 / (n_se * n_io * sgnl^4 * idler * epsilon_0);


%Plotting
% power = -9; %-9 for nm
% scale = 1.5 ;
% lim = scale*Lambda ;
% X = main_signal * 10^(power) ;
% ezplot(spctrm,[X-lim,X+lim])
figure
ezplot(snc^2);
figure
ezplot(spctrm);

end
