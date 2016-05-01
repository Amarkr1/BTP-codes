pump = 519;
signal = 1050;
c = 3e8;
L = 1e3;%mm
[w,h,diff] = waveguide_design(pump,50,signal,L);
[M,I] = min(diff(:));
[I_row, I_col] = ind2sub(size(diff),I)
width_min = w(I_row)
height_min = h(I_col)

%pump
lambda_p = linspace(pump-4,pump+4,5); %value of lambda is in nm
vector_length_p = size(lambda_p);
neff = zeros(1,vector_length_p(2));
for i=1:vector_length_p(2)
    neff(i) = calc_neff(width_min,height_min,lambda_p(i),'o');
end
pwg_p= polyfit(lambda_p,neff,1);
x1 = lambda_p;
y1 = polyval(pwg_p,x1);
slope_pump = pwg_p(1);
neff_central_pump = calc_neff(width_min,height_min,pump,'o');

%signal
lambda_s = linspace(signal-4,signal+4,5); %value of lambda is in nm
vector_length_s = size(lambda_s);
neff = zeros(1,vector_length_s(2));
for i=1:vector_length_s(2)
    neff(i) = calc_neff(width_min,height_min,lambda_s(i),'e');
end
pwg_s= polyfit(lambda_s,neff,1);
x1 = lambda_s;
y1 = polyval(pwg_s,x1);
slope_signal = pwg_s(1);
neff_central_signal = calc_neff(width_min,height_min,signal,'e');

%idler
idler = signal - pump;
% lambda_i = lambda_s - lambda_p; %value of lambda is in nm
lambda_i = linspace(idler-4,idler+4,5); %value of lambda is in nm
vector_length_p = size(lambda_i);
neff = zeros(1,vector_length_p(2));
for i=1:vector_length_p(2)
    neff(i) = calc_neff(width_min,height_min,lambda_i(i),'o');
end
pwg_i= polyfit(lambda_i,neff,1);
x1 = lambda_p;
y1 = polyval(pwg_i,x1);
slope_idler = pwg_i(1);
neff_central_idler = calc_neff(width_min,height_min,idler,'o');


%inverse of group velocities
k_dash_pump = (neff_central_pump/c1) - ((slope_pump)/(c1*((1/pump) + (slope_pump/neff_central_pump))));
k_dash_signal = (neff_central_signal/c1) - ((slope_signal)/(c1*((1/signal) + (slope_signal/neff_central_signal))));
k_dash_idler = (neff_central_idler/c1) - ((slope_idler)/(c1*((1/idler) + (slope_idler/neff_central_idler))));


Y=0.04822; %value of gamma can be controlled to adjust the pump bandwidth

pump_wave = pump/1000 ;%um
signal_wave = signal/1000;%um

np = neff_central_pump;
ni = neff_central_idler;
ns = neff_central_idler;

velocity_pump = c/np;
k1p  = k_dash_pump;
k1i  = k_dash_idler;
k1s  = k_dash_signal;

vgp = 1/k1p
vgi = 1/k1i
vgs = 1/k1s


s0 = signal_wave; %um
ws0 = 2*pi*c/(s0);
i0 = signal_wave;
wi0 = 2*pi*c/(i0);

central = s0;
shift = 0.15;
low = s0 - shift;
high = s0 + shift;
signal = linspace(low,high,100);
idler = linspace(low,high,100);

% alpha = zeros(length(signal),length(idler));

sigma = 1e3/sqrt(-Y*L^2*(k1p-k1s)*(k1p-k1i))

for j = 1: length(signal)
    for k = 1: length(idler)
        kp = 2*pi*np*( 1/idler(k) + 1/signal(j));
        ks = 2*pi*ns/signal(j);
        ki = 2*pi*ni/idler(k);
        delta_k = kp-ks-ki;
        phi(j,k) = exp(-Y*1e-6*L^2*((2*pi*c/signal(j)-ws0)*(k1p-k1s)+((2*pi*c/idler(k)-wi0)*(k1p-k1i)))^2)*exp(-i*(L/2)*(delta_k));
        
%         [nso,nse] = Sellemeir_BBO (signal(j));
%         [nio,nie] = Sellemeir_BBO (idler(k));
        omega_s = 2*pi*c/(signal(j))-ws0;
        omega_i = 2*pi*c/(idler(k))-wi0;
        alpha(j,k) = exp(-((omega_s+omega_i)/sigma)^2);
    end
end
colormap gray
subplot 131
imagesc(signal,idler,alpha)
tt=title('Pump Envelope')
xx = xlabel('Signal (\mu m)');
yy = ylabel('Idler (\mu m)');
set(xx, 'FontSize', 14);
set(tt, 'FontSize', 14);
set(yy, 'FontSize', 14);
set(gca,'XDir','reverse')
set(gca, 'FontSize', 14);

subplot 132
Phase = abs(phi);
imagesc(signal,idler,Phase)
tt=title('Phase Matching')
xx=xlabel('Signal (\mu m)');
yy=ylabel('Idler (\mu m)');
set(gca,'XDir','reverse')
set(gca, 'FontSize', 14);
set(xx, 'FontSize', 14);
set(tt, 'FontSize', 14);
set(yy, 'FontSize', 14);

subplot 133
F = Phase.*alpha;
sum1 = sum(F);
sum2 = sum(sum1);
F = F/sum2;
S = F.^2;
imagesc(signal,idler,S)
tt=title('Spectral Density')
xx=xlabel('Signal (\mu m)');
yy=ylabel('Idler (\mu m)');
set(gca,'XDir','reverse')
set(gca, 'FontSize', 14);
set(xx, 'FontSize', 14);
set(tt, 'FontSize', 14);
set(yy, 'FontSize', 14);
%Schmidt decomposition
[U,S,V] = svd(F);
diag(S);
