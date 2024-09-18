% Division centric model simulated for multiple generations over a single 
% lineage. There is no initiation data. Size control strategy depends on 
% experiments. Returns growth rate, elongation speed and age for a
% simulation with zero experimental error

% clear
% close

function [el_rate_l_f_bin, el_rate_f_bin, a_f_bin] = birth_to_div_gr_v_age_sims(num_in)

[num,txt,raw] = xlsread('params.xlsx');

% Input parameters
tau=num(num_in,1); % generation time in minutes
gr_lin = num(num_in, 12); % lin growth rate
ngen=500; % number of generation
alpha= num(num_in, 7); % Ld = 2(1-alpha) Lb + 2 alpha Delta
bbd= num(num_in, 8); % change in volume from birth to division

%-----Noise parameters
cvl = 0; % num(num_in, 13)*num(num_in, 12); % std linear growth rate
cvr= num(num_in, 17); % std dev of division noise
cve= num(num_in, 16);   % std dev of experimental noise
cvt2= 0; % sqrt(((num(num_in, 2)*num(num_in, 3))^2-cve^2)*alpha*(2-alpha)-4*cvr^2*num(num_in,2)^2)*2; % std dev of size additive division noise

%-----Outputs
% record volume of cell and time from birth(t=0 signifies birth)
vs = NaN(floor(tau*50), ngen); % records vol at each tStep time t from birth
ts = NaN(floor(tau*50), ngen); % all time slots

%-----Advance in time for %gens generations
tStep = 0.1; % in units of hrs
delt=1; % hrs after which length is recorded
gen_cnt=1; % countsthe number of generation passed

while gen_cnt <= ngen
    
    time_cnt=1;
    next_rec=delt;
    t=0;

    vb = num(num_in,2);% + randn()*num(num_in,2)*num(num_in,3);
    v=vb;
    vd = num(num_in,4); %(2*(1-alpha)*vb + 2*alpha*bbd)+randn()*cvt2;
    rate = max(gr_lin + randn()*cvl, gr_lin*0.05);

    while v-vd<0

        %-----Record event and output
        if next_rec>=delt
            ts(time_cnt, gen_cnt)= t;
            vs(time_cnt, gen_cnt)= (v+randn()*cve);
            time_cnt=time_cnt+1;
            next_rec=0;
        end
    
        %-----Step
        next_rec=next_rec+ tStep;
        t = t + tStep;
        v = v + tStep*rate; % grow linear
    end

    gen_cnt = gen_cnt+1;
    
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

% Calculation codes 
[m,n] = size(vs);
ts_tr=zeros(m*n,1); % time from birth
v_f=zeros(m*n,1); % time from birth
el_rate=zeros(m*n,1); % elongation rate or growth rate and each time t
el_rate_l=zeros(m*n,1);
rate=[];
ld=[];lb=[]; td=[];
cnt1=1;
for i=1:n
    v=vs(~isnan(vs(:,i)),i);
    t=ts(1:length(v),i);
    ld=[ld; v(end)];lb=[lb; v(1)]; td=[td; t(end)];
    
    ad=length(v);
    ts_tr(cnt1:cnt1+ad-2,1)=t(1:end-1);
    v_f(cnt1:cnt1+ad-2,1)=v(2:end);
    el_rate_l(cnt1:cnt1+ad-2,1)= 1./v(1:end-1).*(v(2:end)-v(1:end-1))./(t(2:end)-t(1:end-1));
    el_rate(cnt1:cnt1+ad-2,1)=(v(2:end)-v(1:end-1))./(t(2:end)-t(1:end-1));
    cnt1=cnt1+ad-1;
end
ts_tr=ts_tr(1:cnt1-1,1);
v_f=v_f(1:cnt1-1,1);
el_rate_l=el_rate_l(1:cnt1-1,1);
el_rate=el_rate(1:cnt1-1,1);
v_inp=mean(lb);

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

gr_lin= mean((ld-lb)./td);
cvgrlin= std((ld-lb)./td)/gr_lin;
gr=mean(rate);
cv=std(rate)/mean(rate);

[N, edges]= histcounts(a_f, edge_bin);
el_rate_f_bin=zeros(1,length(edges)-1);
el_rate_l_f_bin=zeros(1,length(edges)-1);
a_f_bin=zeros(1,length(el_rate_f_bin));
for i=1:length(a_f_bin)
    el_rate_f_bin(i)=mean(el_rate_f(edges(i)<=a_f & edges(i+1)>=a_f));
    el_rate_l_f_bin(i)=mean(el_rate_l_f(edges(i)<=a_f & edges(i+1)>=a_f));
    a_f_bin(i)=(edges(i)+edges(i+1))/2;
end

