function v_pred = predict_v(time, time_trans, es_rate, constant, intercept)

v_pred = nan(1, length(time));

if time_trans>time(1)
    for t=1:length(time)
        if time(t)<=time_trans
            v_pred(t) = constant;
        else
            v_pred(t) = constant + es_rate*(time(t)-time_trans);
        end
    end
else
    for t=1:length(time)
            v_pred(t) = intercept + es_rate*time(t);
    end 
end


