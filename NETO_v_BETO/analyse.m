% generate histograms of growth asymmetry and elongation speed of each pole 
% in acidic and neutral medium. Select file_call = 1 for acidic pH and 2
% for neutral

clear

file_call = 1; % choose 1 for pH 5.9 , 2 for pH 7 here

file_name = ["pH_59", "pH_7"];

sheet_name = file_name(file_call);
[num,txt,raw] = xlsread('HADA_anl_comb.xlsx', sheet_name);

td = num(:,3);
time_old = num(:,4);
es_old = num(:,5);
grow_old = num(:,6);
time_new = num(:,7);
es_new = num(:,8);
grow_new = num(:,9);
ty = num(:,10);% BEITO = 0, NETO =1 , OETO =2

gr_asym = grow_old./(grow_new+grow_old);
figure
% h = histogram(gr_asym);
% h(1).FaceColor = 'c';
% h(1).EdgeColor = 'k';
hold on
f = histfit(gr_asym,14);
f(1).FaceColor = 'c';
box on
set(gca, 'FontSize', 28)
set(gcf, 'Position',[276,42,777,602])
xlabel('growth asymmetry')
ylabel('Relative frequency')

figure
h = histogram(es_old);
h(1).FaceColor = 'b';
h(1).EdgeColor = 'k';
hold on
h = histogram(es_new);
h(1).FaceColor = 'r'; 
h(1).EdgeColor = 'k'; 
box on
set(gca, 'FontSize', 28)
set(gcf, 'Position',[276,42,777,602])
legend('old pole', 'new pole')
xlabel('Elongation speed (\mum/h)')
ylabel('cell count')