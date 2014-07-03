function [ P,R ] = plot_fit_line( aircraft, behr )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[P, R] = polyfit_R2(aircraft,behr,1);
line(0:1e16:1e17,polyval(P,0:1e16:1e17),'color','k','linestyle','--','linewidth',2);
line(0:1e16:1e17,1:1e16:1e17,'color','r','linestyle',':','linewidth',2);
legend('All points',sprintf('Fit: %.4fx + %.2g \nR^2 = %.4f',P(1),P(2),R),'1:1')
end

