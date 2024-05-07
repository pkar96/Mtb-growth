% We fit length vs time traces of only old pole, and only new pole growth 
% for each individual cell. We use 2 different models- linear, bilinear. We
% calculate various characteristics of old pole and new pole growth such as
% when they start growing and how much do they grow from each pole. We
% consider only those cells which go from no-HADA label region at one pole 
% to both. For these cells, col 1 ~= col 2 = col 3 ~= col 4 in the 
% beginning while at the end all columns are unequal. Select the growth 
% medium to analyze. file_call = 1 for pH 5.9, 2 for pH 7. Output is 
% cell_write.

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
time_to_both = []; % stores the earliest time when col 1 ~= col 2 ~= col 3 ~= col 4
cnt = 1; % counter

for i=1:length(indt1)

    t = num(indt1(i):indte(i),1);
    ls = TotsingleCell(indt1(i):indte(i),:);

    lbe = ls(1,:);
    if ~(lbe(1)~=lbe(2) && lbe(2)==lbe(3) && lbe(3)~=lbe(4)) % only unipolar to begin
        continue
    end
    
    lind= ls(end,:);
    if ~(lind(1)~=lind(2) && lind(2)~=lind(3) && lind(3)~=lind(4)) % only bipolar to end
        continue
    end

    % first time point recorded
    vs_old_pole(1,cnt) = 0;
    vs_new_pole(1,cnt) = 0;
    vs(1,cnt) = ls(1,4);
    ts(1,cnt) = 1;

    for j = 2:length(ls(:,1))
        lbeg = ls(j,:);

        % check when the cell starts growing from both poles
        if lbeg(1)~=lbeg(2) && lbeg(2)~=lbeg(3) && lbeg(3)~=lbeg(4)
            time_to_both =[time_to_both; j];
            keep_cells = [keep_cells; i];
            break
        end

        % record lengths before both poles start growing
        vs_old_pole(j,cnt) = ls(j,4)-ls(j,3)-(ls(1,4)-ls(1,3));
        vs_new_pole(j,cnt) = 0;
        vs(j,cnt) = ls(j,4);
        ts(j,cnt) = j;
    end

    % record lengths after both poles start growing
    j = time_to_both(end);

    td = [td length(ls(:,1))];

    vs_old_pole(j:td(end),cnt) = ls(j:td(end),4)-ls(j:td(end),3)-(ls(1,4)-ls(1,3));
    vs_new_pole(j:td(end),cnt) = ls(j:td(end),2)-ls(j:td(end),1);
    vs(j:td(end),cnt) = ls(j:td(end),4);
    ts(j:td(end),cnt) = j:td(end);    

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
    t = ts(~isnan(vs(:,i)),i); t_old = t; t_new = t;

    % some data points are ignored or in some cases whole thing is ignored
    if file_call == 1
        if i==2 || i==3 || i == 16 || i == 84
            continue
        elseif i == 38
            v_old([44]) = []; t_old([44]) = []; v_new([44]) = []; t_new([44]) = [];
        elseif i == 39
            v_old([46]) = []; t_old([46]) = []; v_new([46]) = []; t_new([46]) = [];
        elseif i == 44
            v_old([60]) = []; t_old([60]) = []; v_new([60]) = []; t_new([60]) = [];
        end
    elseif file_call==2
        if i==32 || i== 78 || i==86 || i==87
            continue
        elseif i == 18 
            v_new([3,4]) = []; t_new([3,4]) = [];
        elseif i == 31
            v_old([90]) = []; t_old([90]) = []; v_new([3,90]) = []; t_new([3,90]) = [];
        elseif i == 37
            v_old([7]) = []; t_old([7]) = [];
        elseif i == 52
            v_old([28]) = []; t_old([28]) = [];
        end
    end

    % bilinear fit
    [fitobj_bilin_old , resid_bilin_old, beta_old] = fit_bilinear(t_old, v_old);
    [fitobj_bilin_new , resid_bilin_new, beta_new] = fit_bilinear(t_new, v_new);
    % linear fit
    [fitobj_lin_old, resid_lin_old, P_old] = fit_poly(t_old,v_old,1);
    [fitobj_lin_new, resid_lin_new, P_new] = fit_poly(t_new,v_new,1);
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
    [es_new, time_new, grow_new, lin_o_bil_new] = find_es_time_auto(td(i), beta_new, P_new, length(find(v_new==0)), Delta_AIC_new, Delta_BIC_new);
    
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
    v_pred_old = predict_v(t_old(v_old~=0), time_old, es_old, beta_old(1), P_old(2));
    v_pred_new = predict_v(t_new(v_new~=0), time_new, es_new, beta_new(1), P_new(2));

%     figure
%     scatter(t_old,v_old, 40, 'filled', MarkerFaceColor='b')
%     hold on
%     pl=plot(t_old(v_old~=0), v_pred_old, 'Color','b', 'LineWidth',2);
%     scatter(t_new,v_new, 40, 'fill', MarkerFaceColor='r')    
%     pl=plot(t_new(v_new~=0), v_pred_new, 'Color','r', 'LineWidth',2);
%     set(gca, 'FontSize', 30)
%     legend off

end

% xlswrite(strcat("HADA_anl_3pt_",file_name(file_call)),cell_write)