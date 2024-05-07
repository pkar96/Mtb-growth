% clear

% gets data from File- ./processed data/(Christin
% 09.03.2021)CDCssbGFPmcherry.csv and returns length vs at time ts for all
% cells being analyzed. Also returns if the cell is accelerator or
% alternator

function [vs, ts, ty] = analysis_09_03()

vs=NaN(75,220);
ts=NaN(75,220);
ty=NaN(1,220);
cnt=1;

% Read in the data
num=[];
txt=[];
raw=[];
[num,txt,raw] = xlsread('./data/processed data/(Christin 09.03.2021)CDCssbGFPmcherry.csv');
txt= txt(2:end,1);

ind= 1; % index whose cell number is searched. Starts with cell at first row
while ind<= length(num)
    [row, column, dataValues] = find(num(:,1)==num(ind,1) & num(:,2)==num(ind,2) & strcmp(txt, txt{ind}));
    if num(ind,3)==3
        ty(cnt) = 2;
    else
        ty(cnt)=num(ind, 3);
    end
    vs(1:length(row), cnt) = sqrt((num(row,5)-num(row, 8)).^2 + (num(row,6)-num(row, 9)).^2 + (num(row,7)-num(row, 10)).^2);
    ts(1:length(row), cnt) = 0:length(row)-1;
    cnt=cnt+1;
    ind=row(end)+1;
end

% post processing
ty(cnt:end)=[];
vs(:,cnt:end)=[];
ts(:,cnt:end)=[];

ts(vs==0)= NaN;
vs(vs==0)= NaN;
