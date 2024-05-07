function [es, time, grow, lin_o_bilin] = find_es_time_auto(td, coeff, P, option, aic, bic)

es = nan; time = nan; grow = nan; lin_o_bilin = nan;

if option ==0
    if aic>0 && bic>0
        lin_o_bilin = 1;
        es = coeff(3);
        time = coeff(2);
        grow = coeff(3)*(td-coeff(2));
    else
        lin_o_bilin = 0;
        es = P(1);
        time = 0;
        grow = es*td;
    end
elseif option~=0
    if aic>0 && bic>0
        lin_o_bilin = 1;
        es = coeff(3);
        time = coeff(2);
        grow = coeff(3)*(td- time);
    else
        if -P(2)/P(1)>=0
            lin_o_bilin = -1;
            es = P(1);
            time = -P(2)/P(1);
            grow = es*(td- time);
        else
            lin_o_bilin = 0;
            es = P(1);
            time = 0;
            grow = es*td;
        end
    end
end
