% Exponential growth is simulated. Size control strategy depends on 
% experiments. Choose num_in 1 to simulate pH 5.9 condition, 2 for pH 7

clear
% close

num_in=1;

[num,txt,raw] = xlsread('params.xlsx');

% Calculate cell cycle characteristics
ngen = num(num_in,1); % number of generation
alpha= num(num_in, 4); % cell size regulation parameter
bbd= num(num_in, 5); % average length at birth
cve = num(num_in, 6); % standard deviation (std dev) of meausrement error
cvr= num(num_in, 7); % std dev of division noise
cvt2= sqrt(((num(num_in, 9)*num(num_in, 10))^2-cve^2)*alpha*(2-alpha)-4*cvr^2*num(num_in,9)^2)*2; % std dev of size additive division noise
tau=num(num_in,8); % generation time in hr
gr = num(num_in,15); % mean growth rate
cvl = num(num_in,16)*gr; % std dev of growth rate

% Output 
ts = NaN(floor(tau*20), ngen); % all time slots
vs = NaN(floor(tau*20), ngen); % records the total measured volume

%-----Advance in time for %gens generations
tStep = 0.01; %in units of hrs
delt=1; % hrs after which length is recorded
gen_cnt=1; % countsthe number of generation passed

while gen_cnt <= ngen

    time_cnt=1;
    next_rec=delt;
    
    t=1;
    vb = num(num_in,9) + randn()*num(num_in,9)*num(num_in,10);
    v=vb;
    vd = (2*(1-alpha)*vb + 2*alpha*bbd)+randn()*cvt2;
    while vd<vb
        vd = (2*(1-alpha)*vb + 2*alpha*bbd)+randn()*cvt2;
    end
    rate = max(gr+randn()*cvl, gr*0.1);

    while v-vd < 0

        %-----Record event and output
        if next_rec>=delt
            ts(time_cnt, gen_cnt)= (t);
            vs(time_cnt, gen_cnt) = v+randn()*cve;
            time_cnt=time_cnt+1;
            next_rec=0;
        end

        %-----Step
        next_rec=next_rec+ tStep;
        t = t+tStep;
        v = v*exp(tStep*rate); % grow exponentially
    end
    
    gen_cnt=gen_cnt+1;

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

