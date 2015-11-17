function [wf] = wave_fn(lambda,w,h,T,String)
syms y z
[n,alpha_y,alpha_z]= give_stats(lambda,w,h,T,String);
 w = w*1e-6;
 h = h*1e-6;
wf(y,z) = -1*sqrt((16*alpha_y*alpha_z)/(pi*w*h))*alpha_z*(z/h)*exp(-alpha_y^2*y^2/w^2)*exp(-alpha_z^2*z^2/h^2);

%Plotting
scale = 1.5;
w_lim = scale*w;
h_lim= -scale*h;
ezsurf(wf,[-w_lim,w_lim],[h_lim,0])
set(gca, 'FontSize', 12);

xx = xlabel('y (m)');
set(xx, 'FontSize', 14);
zz = ylabel('z (m)');
set(zz, 'FontSize', 14);
shading interp;

end
