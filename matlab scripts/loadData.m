len=100000;
% for i=1:len
%     ffname=strcat('cdat/E_Field_',int2str(i*10),'.txt');
%     load(ffname);
% end
for i=40000:len
    ffname=strcat('cdat/E_Field_',int2str(i*10),'.txt');
    load(ffname);
end