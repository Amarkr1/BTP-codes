function [I_eo , n_po , n_se, n_io] = I_eo(pump, signal , idler , w , h, T)%wavelengths , T temperature in deg C
[n_po, alpha_poy, alpha_poz] = give_stats(pump,w, h, T,'o');
[n_se, alpha_sey, alpha_sez] = give_stats(signal,w, h, T,'e');
[n_io, alpha_ioy, alpha_ioz] = give_stats(idler,w, h, T,'o');
I_eo = (32 * (alpha_poz*alpha_sez*alpha_ioz)^1.5*(alpha_poy*alpha_sey*alpha_ioy)^0.5 )/( pi * sqrt(w*h) * (alpha_poz^2+alpha_sez^2+alpha_ioz^2)^2*(alpha_poy^2+alpha_sey^2+alpha_ioy^2)^0.5 );
end
