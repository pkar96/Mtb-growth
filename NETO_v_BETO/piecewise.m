function y=piecewise(coef,x)
% return piecewise linear fit given [b1,m1,c,m2] as coefficient array and vector x
% c is breakpoint, b1 is intercept of first section, m1,m2 are two segment slopes
%     b1=coef(1); m1=coef(2);  % reduce parentheses clutter....
%     c=coef(3);  m2=coef(4);
%     ix=x>=c;                  % location > breakpoint
%     y(ix)=[b1+c*(m1-m2)+m2*x(ix)];
%     y(~ix)=[b1+m1*x(~ix)];
%     y=reshape(y,size(x));

    b1=coef(1);  % reduce parentheses clutter....
    c=coef(2);  m1=coef(3);
    ix=x>=c;                  % location > breakpoint
    y(~ix)=[b1];
    y(ix)=[b1+m1*(x(ix)-c)];
    y=reshape(y,size(x));

end