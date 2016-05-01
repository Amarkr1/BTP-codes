function [slope] = slope_wg(LAMBDA,w,h,String)

lambda = linspace(LAMBDA-5,LAMBDA+5,5); %value of lambda is in nm
vector_length = size(lambda);
neff = zeros(1,vector_length(2));
for i=1:vector_length(2)
    neff(i) = calc_neff(w,h,lambda(i),String);
end
p= polyfit(lambda,neff,1);
slope = p(1);
% x1 = lambda;
% y1 = polyval(p,x1);
% figure
% plot(lambda,neff_p,'x')
% hold on
% plot(x1,y1)
% plot(lambda,neff,'-ro',x1,y1,'-.b')
% hleg1 = legend('data','fit');
end
