function makeomovie(movname,varhead,file,figname,islog,step)
% function MAKEMOVIE
% Make movie with a series of variables with name 'varhead' in 'file'.
% Usage: makemovie(movname,varhead,file,figname)
% Inputs:
%       movname ------ name of output file in format of '*.avi'
%       varhead ------ key charaters of variable name
%          file ------ Matlab mat-file to be load from (not with .mat)
%       figname ------ name of figure
% Outputs:
%       an avi-file with name movname.avi
%
% by skiloop@126.com
% 2011/12/19 09:40
%

% check if inputs are all chars
if ~ischar(varhead)||~ischar(file)||~ischar(movname)||~ischar(figname)
    error('All inputs are required to be char!');
end

% check if file exist
file = strcat(file,'.mat');
if ~exist(file,'file');
    error('File not exist!\nPlease check the file and rerun the program!');
end

% load variables
varheads=strcat(varhead,'*');
%load(file,varheads);
li = whos('-file',file,varheads);
if isempty(li)
    warning('No variable names with  varhead found in the file.');
    return;
end
load(file,li(1).name);
% initial the figure
h = figure('NumberTitle','OFF','Name',figname);
%l = round(dim*27):round(dim*175);
tvar= eval(li(1).name);
% if islog
%     tvar(tvar<1)=1;
%     ih= imagesc(log(tvar(l,l)));
% else
%     ih= imagesc(tvar(l,l));
% end
ih= imagesc(tvar);
set(gca,'YDir','Normal');
ax = gca;
colorbar;
title(figname);
mov(1)=getframe(h);
% yt = get(gca,'ytick');
% xt = get(gca,'xtick');
% start to make movie
len = length(li);
if len>1000
    len=1000;
end

for iter = 1:len
    ffname=strcat(varhead,int2str(iter*step));
    load(file,ffname);
    tvar= eval(ffname);
    if islog
        tvar(tvar<1)=1;
        tvar = log(tvar);
    end
    %set(ih,'CData',tvar(l,l));
    set(ih,'CData',tvar);
    xlabel(ax,['Time : ',datestr(clock)]);
    title(ax,[figname,' step = ',int2str(iter*step)]);
    drawnow;
    clear(li(iter).name);
    mov(iter) = getframe(h);
end
movie2avi(mov, movname, 'compression', 'None');



