% The code finds the mean squared residuals between the binned data from 
% experiments and simulations of growth rate vs age and elongation speed vs
% age plots. The mean squared residuals are calculated for ntimes (=200) 
% runs of the simulations. The minimum value is shown for each of the two 
% plots. The simulation strategy is mentioned in Simulations section of the
% paper. Size control strategy depends on experiments. Poles have
% different growth rates. Old and new pole can start growing after certain 
% time. Growth rate distributions are sampled from experiments. vs_old and
% vs_new also have experimental errors added to them. Choose num_in 1 to
% simulate pH 5.9 condition, 2 for pH 7.

clear
% close

[num,txt,raw] = xlsread('params.xlsx');
num_in=1;
sheet_name =["pH_59";"pH_7"];
[num1,txt1,raw1] = xlsread('sample_time_from_here.xlsx',sheet_name(num_in));

[numl,txtl,rawl] = xlsread('binned_data_exp.xlsx',strcat(sheet_name(num_in),'_gr_v_age'));
[numel,txtel,rawel] = xlsread('binned_data_exp.xlsx',strcat(sheet_name(num_in),'_dLdt_v_age'));
el_rate_l_exp = numl(1:30,2);
el_rate_el_exp = numel(1:30,2);

ntimes = 200;

% Outputs
ssr_l = NaN(1,ntimes);
ssr_el = NaN(1,ntimes);

for sim_cnt = 1:ntimes

% Calculate cell cycle characteristics
ngen = num(num_in,1); % number of generation
n3pt = num(num_in,2); % new pole have no HADA label in the beginning
n4pt = num(num_in,3); % new pole have HADA label in the beginning
alpha= num(num_in,4); % cell size regulation parameter
bbd= num(num_in, 5); % average length at birth
cve = num(num_in,6); % standard deviation (std dev) of meausrement error
cvr= num(num_in, 7); % std dev of division noise
cvt2= sqrt(((num(num_in, 9)*num(num_in, 10))^2-cve^2)*alpha*(2-alpha)-4*cvr^2*num(num_in,9)^2)*2; % std dev of size additive division noise
cve = 0;
tau=num(num_in,8); % generation time in hr
gr_lin_old = num(num_in, 11); % mean linear growth rate for old pole 
gr_lin_new = num(num_in, 12); % mean linear growth rate for new pole
cvl_old = num(num_in, 13); % CV of linear growth rate for old pole
cvl_new = num(num_in, 14); % CV of linear growth rate for new pole
prob_t = num(num_in, 17); % probability of time series growing from beginning

% Output 
vs_old = NaN(floor(tau*30), ngen); % records vol for old pole at each tStep time t from birth
vs_new = NaN(floor(tau*30), ngen); % records vol for new pole at each tStep time from birth
ts = NaN(floor(tau*30), ngen); % all time slots
vs = NaN(floor(tau*30), ngen); % records the measured volume after changes are done
npgrs=[]; % record growth from new pole
opgrs=[]; % record growth from old pole
oldpol_ts =[]; % record time for old pole growth
newpol_ts =[]; % record time for new pole growth

%-----Advance in time for %gens generations
tStep = 0.05; %in units of hrs
delt=0.9; % hrs after which length is recorded
gen_cnt=1; % countsthe number of generation passed

while gen_cnt <= ngen

    time_cnt=1;
    next_rec=delt;
    
    t=0;
    vb = num(num_in,9) + randn()*num(num_in,9)*num(num_in,10);
    while vb<0
        vb = num(num_in,9) + randn()*num(num_in,9)*num(num_in,10);
    end
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
        vs(:,gen_cnt) = vb+vs_old(:,gen_cnt)+vs_new(:,gen_cnt);
        gen_cnt=gen_cnt+1;
    end

end

ts(vs==0)= NaN;
vs(vs==0)= NaN;

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

% calculate growth rate and elongation speed vs age according to Nordholt et al.

npt = 100; % number of points in the cell cycle to be binned first
tot_bin=30;
edge_bin = 0:1/tot_bin:1; 
[m,n] = size(vs);
a_f = NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1);
el_rate_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % elongation rate final for all cells
el_rate_l_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % 

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

[bin_l,da_l,yfit, P, err_l]=binning_with_error_1(el_rate_l_f,a_f, edge_bin);
bin_l = bin_l'; da_l = da_l'; err_l = err_l';
ssr_l(sim_cnt) = mean((da_l - el_rate_l_exp).^2,'omitnan');

[bin,da,yfit, P, err]=binning_with_error_1(el_rate_f,a_f, edge_bin); 
bin = bin'; da = da'; err = err';
ssr_el(sim_cnt) = mean((da - el_rate_el_exp).^2,'omitnan');

end

[min(ssr_l) min(ssr_el)]
ssr_el = ssr_el'; ssr_l = ssr_l';
