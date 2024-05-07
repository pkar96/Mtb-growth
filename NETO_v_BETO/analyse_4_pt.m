% We fit length vs time traces of only old pole, and only new pole growth 
% for each individual cell. We use 2 different models- linear, bilinear. We
% calculate various characteristics of old pole and new pole growth such as
% when they start growing and how much do they grow from each pole. We
% consider only those cells which have no-HADA label region at both poles 
% from start. For these cells, col 1 ~= col 2 ~= col 3 ~= col 4 throughout.
% Select the growth medium to analyze. file_call = 1 for pH 5.9, 2 for pH 
% 7. Output is cell_write.

clear

file_call = 1; % choose 1 for pH 5.9 , 2 for pH 7 here

file_name = ["pH_59.xlsx", "pH_7.xlsx"];
[num,txt,raw] = xlsread(file_name(file_call));

indt1 = find(num(:,1)==1); % start points of each cell (time t =1)
indte = [indt1(2:end)-1; length(num)]; % end time of each cell
TotsingleCell = num(:,2:5);

% store growth for both poles in an array
keep_cells =[]; % cells which go from unipolar -> bipolar growth
vs_old_pole = NaN(150, length(indt1)); % stores the growth from old pole
vs_new_pole = NaN(150, length(indt1)); % stores the growth from new pole
vs = NaN(150, length(indt1)); % stores the total length
ts = NaN(150, length(indt1)); % stores the time from birth
td = []; % stores the generation time
cnt = 1; % counter

for i=1:length(indt1)

    t = num(indt1(i):indte(i),1);
    ls = TotsingleCell(indt1(i):indte(i),:);

    lbe = ls(1,:);
    if ~(lbe(1)~=lbe(2) && lbe(2)~=lbe(3) && lbe(3)~=lbe(4)) % only bipolar cells
        continue
    end

    % record lengths
    keep_cells = [keep_cells; i];
    td = [td length(ls(:,1))];
    vs_old_pole(1:td(end),cnt) = ls(1:td(end),4)-ls(1:td(end),3)-(ls(1,4)-ls(1,3));
    vs_new_pole(1:td(end),cnt) = ls(1:td(end),2)-ls(1:td(end),1)-(ls(1,2)-ls(1,1));
    vs(1:td(end),cnt) = ls(1:td(end),4);
    ts(1:td(end),cnt) = 1:td(end);    

    cnt=cnt+1;
    
end

vs_old_pole(:,cnt:end)=[];
vs_new_pole(:,cnt:end)=[];
vs(:,cnt:end)=[];
ts(:,cnt:end)=[];

%%

% fitting and calculating residuals
[~,n] = size(vs_old_pole);
cnt = 1;

% Output
cell_write= {'No','Cell number in original file', 'Td', 'Old-Time at growth (h)', 'Old-Elongation speed(um/h)','Old-Amt of growth(um)', 'New-Time at growth (h)', 'New-Elongation speed(um/h)','New-Amt of growth(um)', 'BEITO(0)/NETO(1)/OETO(2)'};

for i=1:n 

    v_old = vs_old_pole(~isnan(vs(:,i)),i);
    v_new = vs_new_pole(~isnan(vs(:,i)),i);
    t = ts(~isnan(vs(:,i)),i); 

    % bilinear fit
    [fitobj_bilin_old , resid_bilin_old, beta_old] = fit_bilinear(t, v_old);
    [fitobj_bilin_new , resid_bilin_new, beta_new] = fit_bilinear(t, v_new);
    % linear fit
    [fitobj_lin_old, resid_lin_old, P_old] = fit_poly(t,v_old,1);
    [fitobj_lin_new, resid_lin_new, P_new] = fit_poly(t,v_new,1);
    % Mean squared residual adjusted
    MSR_bilin_old = sum(resid_bilin_old.^2)/(length(find(v_old~=0))); MSR_lin_old = sum(resid_lin_old.^2)/(length(find(v_old~=0)));
    MSR_bilin_new = sum(resid_bilin_new.^2)/(length(find(v_new~=0))); MSR_lin_new = sum(resid_lin_new.^2)/(length(find(v_new~=0)));
    % AIC and BIC calculation
    Delta_AIC_old = 2*2+12/(length(find(v_old~=0))-3)+length(find(v_old~=0))*log(MSR_lin_old)-(2*3+24/(length(find(v_old~=0))-4)+length(find(v_old~=0))*log(MSR_bilin_old)); % AIC_lin-AIC_bilin
    Delta_AIC_new = 2*2+12/(length(find(v_new~=0))-3)+length(find(v_new~=0))*log(MSR_lin_new)-(2*3+24/(length(find(v_new~=0))-4)+length(find(v_new~=0))*log(MSR_bilin_new)); % AIC_lin-AIC_bilin
    Delta_BIC_old = log(length(find(v_old~=0)))*2+length(find(v_old~=0))*log(MSR_lin_old)-(log(length(find(v_old~=0)))*3+length(find(v_old~=0))*log(MSR_bilin_old)); % BIC_lin-BIC_bilin
    Delta_BIC_new = log(length(find(v_new~=0)))*2+length(find(v_new~=0))*log(MSR_lin_new)-(log(length(find(v_new~=0)))*3+length(find(v_new~=0))*log(MSR_bilin_new)); % BIC_lin-BIC_bilin

    if Delta_BIC_old*Delta_AIC_old<0
        'BIC and AIC signs exchanged for old'
    end
    if Delta_BIC_new*Delta_AIC_new<0
        'BIC and AIC signs exchanged for new'
    end

    [es_old, time_old, grow_old, lin_o_bil_old] = find_es_time_auto(td(i), beta_old, P_old, 0, Delta_AIC_old, Delta_BIC_old);
    [es_new, time_new, grow_new, lin_o_bil_new] = find_es_time_auto(td(i), beta_new, P_new, 0, Delta_AIC_new, Delta_BIC_new);
    
    % find if cell is BEITO, NETO or OETO
    if time_new<1 && time_old<1
        ty = 0;
    elseif time_old<time_new
        ty= 1;
    else
        ty=2;
    end
    
    cell_write(end+1, :)= {cnt,strcat('Cell ',int2str(keep_cells(i))), td(i), time_old, es_old, grow_old, time_new, es_new, grow_new, ty};
    cnt = cnt+1;

    % plot old pole and new pole growth
    % Get values for plotting
    v_pred_old = predict_v(t(v_old~=0), time_old, es_old, beta_old(1), P_old(2));
    v_pred_new = predict_v(t(v_new~=0), time_new, es_new, beta_new(1), P_new(2));

%     figure
%     scatter(t,v_old, 40, 'filled', MarkerFaceColor='b')
%     hold on
%     pl=plot(t(v_old~=0), v_pred_old, 'Color','b', 'LineWidth',2);
%     scatter(t,v_new, 40, 'fill', MarkerFaceColor='r')    
%     pl=plot(t(v_new~=0), v_pred_new, 'Color','r', 'LineWidth',2);
%     set(gca, 'FontSize', 30)
%     legend off

end

% xlswrite(strcat("HADA_anl_4pt_",file_name(file_call)),cell_write)