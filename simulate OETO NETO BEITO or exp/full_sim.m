% The simulation strategy is mentioned in Simulations section of the paper.
% Size control strategy depends on experiments. Poles have
% different growth rates. Old and new pole can start growing after certain 
% time. Growth rate distributions are sampled from experiments. vs_old and
% vs_new also have experimental errors added to them. Choose num_in 1 to
% simulate pH 5.9 condition, 2 for pH 7

clear
% close

num_in = 1;

[num,txt,raw] = xlsread('params.xlsx');
sheet_name =["pH_59";"pH_7"];
[num1,txt1,raw1] = xlsread('sample_time_from_here.xlsx',sheet_name(num_in));

% Calculate cell cycle characteristics
ngen = num(num_in,1); % number of generation
n3pt = num(num_in,2); % new pole have no HADA label in the beginning
n4pt = num(num_in,3); % new pole have HADA label in the beginning
alpha= num(num_in, 4); % cell size regulation parameter
bbd= num(num_in, 5); % average length at birth
cve = num(num_in, 6); % standard deviation (std dev) of meausrement error
cvr= num(num_in, 7); % std dev of division noise
cvt2= sqrt(((num(num_in, 9)*num(num_in, 10))^2-cve^2)*alpha*(2-alpha)-4*cvr^2*num(num_in,9)^2)*2; % std dev of size additive division noise
tau=num(num_in,8); % generation time in hr
gr_lin_old = num(num_in, 11); % mean linear growth rate for old pole 
gr_lin_new = num(num_in, 12); % mean linear growth rate for new pole
cvl_old = num(num_in, 13); % CV of linear growth rate for old pole
cvl_new = num(num_in, 14); % CV of linear growth rate for new pole
prob_t = num(num_in, 17); % probability of time series growing from beginning

% Output 
vs_old = NaN(floor(tau*10), ngen); % records vol for old pole at each tStep time t from birth
vs_new = NaN(floor(tau*10), ngen); % records vol for new pole at each tStep time t from birth
ts = NaN(floor(tau*10), ngen); % all time slots
vs = NaN(floor(tau*10), ngen); % records the total measured volume
npgrs=[]; % record growth from new pole
opgrs=[]; % record growth from old pole
oldpol_ts =[]; % record time for old pole growth
newpol_ts =[]; % record time for new pole growth

%-----Advance in time for %gens generations
tStep = 0.1; % in units of hrs
delt=1; % hrs after which length is recorded
gen_cnt=1; % counts the number of generation passed

while gen_cnt <= ngen

    time_cnt=1;
    next_rec=delt;
    
    t=1;
    vb = num(num_in,9) + randn()*num(num_in,9)*num(num_in,10);
    v=vb;
    opgr= 0;
    npgr= 0;
    vd = (2*(1-alpha)*vb + 2*alpha*bbd)+randn()*cvt2;
    while vd<vb
        vd = (2*(1-alpha)*vb + 2*alpha*bbd)+randn()*cvt2;
    end
    rate_op = max(gr_lin_old + randn()*cvl_old*gr_lin_old, gr_lin_old*0.05);
    rate_np = max(gr_lin_new + randn()*cvl_new*gr_lin_new, gr_lin_new*0.05);

    if rand()<prob_t
        old_pole_gr_t = 0;
        new_pole_gr_t = 0;
    else
        indt = randsample(size(num1,1),1);
        old_pole_gr_t = num1(indt,1);
        new_pole_gr_t = num1(indt,2);
    end

    old_pole_gr_t_c = old_pole_gr_t; new_pole_gr_t_c = new_pole_gr_t;

    while v-vd < 0

        %-----Record event and output
        if next_rec>=delt
            ts(time_cnt, gen_cnt)= (t);
            vs_old(time_cnt, gen_cnt) = opgr;
            vs_new(time_cnt, gen_cnt) = npgr;
            time_cnt=time_cnt+1;
            next_rec=0;
        end

        %-----Step
        next_rec=next_rec+ tStep;
        t = t+tStep;
        if old_pole_gr_t<=0
            v = v + tStep*rate_op; % grow linear old pole
            opgr=opgr+tStep*rate_op;
        end
        if new_pole_gr_t<=0
            v = v + tStep*rate_np; % grow linear new pole
            npgr=npgr+tStep*rate_np;
        end
        old_pole_gr_t = old_pole_gr_t-tStep; new_pole_gr_t = new_pole_gr_t-tStep;

    end
    
    if old_pole_gr_t>0 || new_pole_gr_t>0 % if one of the poles has not started growing, ignore that cell

        vs_old(:,gen_cnt) = nan*ones(length(vs_old(:,gen_cnt)),1);
        vs_new(:,gen_cnt) = nan*ones(length(vs_new(:,gen_cnt)),1);
        ts(:,gen_cnt) = nan*ones(length(ts(:,gen_cnt)),1);

    else
        % record length measurement
        % record n3pt prop of cells where the new pole's HADA label is very small to be registered for a while 
        if rand()<n3pt/ngen
            if new_pole_gr_t_c<=1
                [n,tl,r] = xlsread('sample_new_pl_zeros.xlsx'...
                    ,strcat("lin_",sheet_name(num_in))); % column 1 stores time, column 2 stores a bias
                indt = randsample(size(n,1),1); % get a random time from experimental data for which vs_new = 0
                tnnan=ts(~isnan(ts(:,gen_cnt)),gen_cnt);
                while n(indt,1)>=tnnan(end) % if the time generated is greater than the doubling time, generate again
                    indt = randsample(size(n,1),1);
                end
                vs_new(:,gen_cnt) = vs_new(:,gen_cnt)+n(indt,2);
                vs_new(1:n(indt,1), gen_cnt) = 0;
            else
                [n,tl,r] = xlsread('sample_new_pl_zeros.xlsx'...
                    ,strcat("bi_cos_0_",sheet_name(num_in))); % same thing repeated for bilinear cases
                indt = randsample(size(n,1),1);
                tnnan=ts(~isnan(ts(:,gen_cnt)),gen_cnt);
                while new_pole_gr_t_c+n(indt,1)>=tnnan(end)  %n(indt,1)<new_pole_gr_t_c && n(indt,1)>=tnnan(end)
                    indt = randsample(size(n,1),1);
                end
                vs_new(1:round(n(indt,1)+new_pole_gr_t_c), gen_cnt) = 0;
            end
        end

        vs(:,gen_cnt) = vb+vs_old(:,gen_cnt)+vs_new(:,gen_cnt)+randn(length(vs(:,gen_cnt)),1)*cve; % total length with measurement error added
        % adding measurement errors to vs_old and vs_new
        tnnan=ts(~isnan(ts(:,gen_cnt)),gen_cnt);
        indnz = find(vs_new(1:round(tnnan(end)),gen_cnt)~=0);
        vs_old(2:round(tnnan(end)),gen_cnt) = vs_old(2:round(tnnan(end)),gen_cnt)+randn(round(tnnan(end))-1,1)*cve;
        vs_new(indnz,gen_cnt) = vs_new(indnz,gen_cnt)+randn(length(indnz),1)*cve;
        % recording growth amount and time of start of growth
        npgrs=[npgrs npgr];
        opgrs=[opgrs opgr];
        oldpol_ts = [oldpol_ts old_pole_gr_t_c];
        newpol_ts = [newpol_ts new_pole_gr_t_c];
        gen_cnt=gen_cnt+1;
    end

end

[m,n] = size(vs);
indrm=[];
for i=1:n
    v=vs(~isnan(vs(:,i)),i);
    if length(v)<=5
        indrm=[indrm i];
    end
end
vs(:,indrm)=[];
ts(:,indrm)=[];

%%
% analysis on collected data
% getting values of lb, ld, td values
[m,n] = size(vs);
ld=[];lb=[]; td=[];
for i=1:n
    v=vs(~isnan(vs(:,i)),i);
    t=ts(1:length(v),i);
    ld=[ld; v(end)];lb=[lb; v(1)]; td=[td; t(end)];
    
end

v_inp=mean(lb);
gr_lin = mean((ld-lb)./td);
cvgrlin= std((ld-lb)./td)/gr_lin;

%%

% calculate growth rate and elongation speed vs age.

npt = 100; % number of points in the cell cycle to be binned first
tot_bin=30; % total number of bins
edge_bin = 0:1/tot_bin:1; 
[m,n] = size(vs);
a_f = NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % age
el_rate_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % elongation speed of all cells
el_rate_l_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % growth rate of all cells 

cnt=1;

for ind=1:n
    v=vs(~isnan(vs(:,ind)),ind);
    t=ts(1:length(v),ind);
    
    ad=length(v);
    el_rate_l_i= 1./v(1:end-1).*(v(2:end)-v(1:end-1))./(t(2:end)-t(1:end-1));
    el_rate_i=(v(2:end)-v(1:end-1))./(t(2:end)-t(1:end-1));
    
    xv= t(1):(t(end-1)-t(1))/npt:t(end-1);
    vq_l = interp1(t(1:end-1),el_rate_l_i,xv);
    vq = interp1(t(1:end-1),el_rate_i,xv);
    a_i= xv/t(end);
    el_rate_l_i= vq_l;
    el_rate_i= vq;
    
    [N, edges]= histcounts(a_i, edge_bin);
    bin=zeros(1,length(edges)-1);
    da=zeros(1,length(bin));
    for i=1:length(bin)
        da(i)=mean(el_rate_l_i(edges(i)<=a_i & edges(i+1)>=a_i));
        bin(i)=(edges(i)+edges(i+1))/2;
    end
    
    el_rate_l_f(cnt:cnt+(edge_bin(end)-edge_bin(1))*tot_bin-1)= da;
    a_f(cnt:cnt+(edge_bin(end)-edge_bin(1))*tot_bin-1) = bin;
    
    [N, edges]= histcounts(a_i, edge_bin);
    bin=zeros(1,length(edges)-1);
    da=zeros(1,length(bin));
    for i=1:length(bin)
        da(i)=mean(el_rate_i(edges(i)<=a_i & edges(i+1)>=a_i));
        bin(i)=(edges(i)+edges(i+1))/2;
    end
    
    el_rate_f(cnt:cnt+(edge_bin(end)-edge_bin(1))*tot_bin-1)= da;
    cnt=cnt+(edge_bin(end)-edge_bin(1))*tot_bin;
    
end

%%
% plot growth rate vs age
figure
[bin,da,yfit, P, err]=binning_with_error_1(el_rate_l_f,a_f, edge_bin); 
hold on
ind= find(bin<2);
s=scatter(bin(ind),da(ind), 100, 'filled');
s.MarkerEdgeColor='k';
s.MarkerFaceColor='k'; %[0.85, 0.33, 0.1];
eb = errorbar(bin(ind),da(ind),err(ind), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', 'k');
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('\lambda (hr^{-1})');
xlabel('Age');
xlim([0.1, 1])
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
ax = gca; % axes handle
ax.YAxis.Exponent = -2;
box on

%%
% plot dL/dt vs age
figure
[bin,da,yfit, P, err]=binning_with_error_1(el_rate_f,a_f, edge_bin); 
hold on
ind= find(bin<2);
s=scatter(bin(ind),da(ind), 100, 'filled');
s.MarkerEdgeColor='k';
s.MarkerFaceColor='k'; %[0.85, 0.33, 0.1];
eb = errorbar(bin(ind),da(ind),err(ind), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', 'k');
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('dL/dt (\mum hr^{-1})');
xlabel('Age');
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
ax = gca; % axes handle
ax.YAxis.Exponent = -1;
box on
xlim([0.1, 1])

%%
% plot growth asymmetry
figure
histogram(opgrs./(npgrs+opgrs))
xlabel('Growth asymmetry')
ylabel('Count')
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
box on

mean(opgrs./(npgrs+opgrs));
std(opgrs./(npgrs+opgrs));

%%
% plot time of start of old pole growth and new pole growth in non-BEITO cases 
figure
histogram(oldpol_ts(oldpol_ts>1|newpol_ts>1))
hold on
histogram(newpol_ts(oldpol_ts>1|newpol_ts>1))
xlabel('Time of old pole growth')
ylabel('Count')
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
box on
