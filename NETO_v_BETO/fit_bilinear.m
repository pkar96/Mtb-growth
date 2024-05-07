% calculates bilinear fit and returns predicted y, residuals and fit
% parameters
function [y, R, beta]  = fit_bilinear(t,v)

    [beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(t(v~=0),v(v~=0),@piecewise,[0, 5, 0.02]);
    y=piecewise(beta,t(v~=0));
end
