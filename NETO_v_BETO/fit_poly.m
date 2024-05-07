function [y, R, P]  = fit_poly(t,v,ord)

P = polyfit(t(v~=0),v(v~=0),ord);
y = polyval(P,t(v~=0));
R = v(v~=0)-y;