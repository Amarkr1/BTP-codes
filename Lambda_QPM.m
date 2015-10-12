function [Lambda] = Lambda_QPM (pump , signal , w ,h ,T ) % p,s and i in nm ; w and h in um ; T in deg. C

syms idlr Lambda_0
[n_po, alpha_poy, alpha_poz] = give_stats(pump,w, h, T,'o');
[n_se, alpha_sey, alpha_sez] = give_stats(signal,w, h, T,'e');

idler = solve (1/pump == 1/signal + 1/idlr , idlr); %Conservation of Energy

[n_io, alpha_ioy, alpha_ioz] = give_stats(idler,w, h, T,'o');

Lambda = solve( 1/Lambda_0 == n_po/pump + n_se/signal +n_io/idler ,Lambda_0 );
end
%Result in nm
