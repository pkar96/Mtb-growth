clear

% plot growth rate vs age info as shown in the main text of the paper for
% acidic (pH 5.9) and neutral medium (pH 7). Select the growth medium to 
% analyze. file_call = 1 for pH 5.9, 2 for pH 7

file_call = 1; % choose 1 for pH 5.9 , 2 for pH 7 here

file_name = ["pH_59.xlsx", "pH_7.xlsx"];
[num,txt,raw] = xlsread(file_name(file_call));

indt1 = find(num(:,1)==1); % start points of each cell (time t =1)
indte = [indt1(2:end)-1; length(num)]; % end time of each cell

% stores volume and times of all cells
vs=NaN(150,length(indt1));
ts=NaN(150,length(indt1));

TotsingleCell = num(:,2:5);

for i=1:length(indt1)
    vs(1:indte(i)-indt1(i)+1, i) = TotsingleCell(indt1(i):indte(i),4);
    ts(1:indte(i)-indt1(i)+1, i) = num(indt1(i):indte(i),1);
end

ts(vs==0)= NaN;
vs(vs==0)= NaN;

%%
% analysis codes

% estimating measurement noise. Method used in Appendix 3 of Messelink et
% al.

delt=1;
[m,n] = size(vs);
ld=[]; lb=[]; td=[];
for i=1:n
    v=vs(~isnan(vs(:,i)),i);
    t=ts(1:length(v),i);
    ld=[ld; v(end)]; lb=[lb; v(1)]; td=[td; t(end)];
end

min_td=min(td); std_gr=std((ld-lb)./td)^2;

lmt= vs(2:min_td/delt,:)-vs(1,:);
var_t=var(lmt,0,2);
x_t= 1:min_td/delt-1;

mdl = fittype(@(a,b,x) a + b*x.^2,'independent','x');
fittedmdl = fit(x_t',var_t,mdl,'start',[10^-4, 10^-2]);

figure
plot(fittedmdl)
hold on
plot(1:min_td/delt-1,var_t,'ro', 'MarkerFaceColor','r')
ylabel('Var(l_m(\Delta t)-l_m(0))');
xlabel('\Delta t (hr)');
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
ax= gca;
ax.YAxis.Exponent = -2;
cve = sqrt(fittedmdl.a/2);

%%
% calculate growth rate and elongation speed vs age according to Nordholt et al. 

npt = 100; % number of points in the cell cycle to be binned first
tot_bin=30;
edge_bin = 0:1/tot_bin:1; 
[m,n] = size(vs);
a_f = NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1);
el_rate_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % elongation rate final for all cells
el_rate_l_f= NaN(round((edge_bin(end)-edge_bin(1))*tot_bin*n),1); % growth rate final for all cells

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

v_inp = mean(ld-lb);
gr_lin= mean((ld-lb)./td); % plots linear growth rate
cvgrlin = std((ld-lb)./td)/gr_lin;

% collect simulation results to make plots for comparing against
% experiments

repeat = 150;
el_rate_f_sim = NaN(repeat,30); el_rate_l_f_sim =NaN(repeat,30); a_f_sim=NaN(repeat,30);
for i=1:repeat
    [el_rate_l_f_sim_e, el_rate_f_sim_e, a_f_sim_e] = birth_to_div_gr_v_age_sims(file_call+5);
    el_rate_l_f_sim(i,:)=el_rate_l_f_sim_e;
    el_rate_f_sim(i,:)=el_rate_f_sim_e;
    a_f_sim(i,:) = a_f_sim_e;
end

el_rate_f_sim_f = mean(el_rate_f_sim); el_rate_l_f_sim_f = mean(el_rate_l_f_sim); a_f_sim_f = mean(a_f_sim);

% Plots growth rate vs age
figure
[bin,da,yfit, P, err]=binning_with_error_1(el_rate_l_f,a_f, edge_bin); 
hold on
ind= find(bin<2);
s=scatter(bin(ind),da(ind), 100, 'filled');
s.MarkerEdgeColor='k';
s.MarkerFaceColor=[0.85, 0.33, 0.1];
s1=plot(a_f_sim_f,el_rate_l_f_sim_f, 'LineWidth',4, 'Color','k');
eb = errorbar(bin(ind),da(ind),err(ind), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.1]);
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('\lambda (hr^{-1})');
xlabel('Age');
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
ax = gca;
ax.YAxis.Exponent = -2;
box on

%%
% Plots elongation speed vs age
figure
[bin,da,yfit, P, err]=binning_with_error_1(el_rate_f,a_f, edge_bin);
hold on
ind= find(bin<2);
s=scatter(bin(ind),da(ind), 150, 'filled');
s.MarkerEdgeColor='k';
s.MarkerFaceColor=[0.85, 0.33, 0.1];
s1=plot(a_f_sim_f,el_rate_f_sim_f, 'LineWidth',4, 'Color','k');
eb = errorbar(bin(ind),da(ind),err(ind), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', [0.85, 0.33, 0.1]);
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('dL/dt (\mum/hr)');
xlabel('Age');
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
box on