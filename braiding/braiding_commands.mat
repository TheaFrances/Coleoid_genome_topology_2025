clear

addpath '~/matlab/braidlab-release-3.2.5'
import braidlab.*


pre='esc-obi.10m'

maxhom=50;

fn=strcat(pre,'.braid.map');

map = readtable(fn, "FileType","text",'Delimiter', '\t');

ids=unique(map(:,1));

rows=size(ids,1);
 res=[];

for i=1:rows
 id=ids(i,1).Var1;
 submap=map(map.Var1==id,:);
 

mapa0=cat(2,submap.Var10,submap.Var11,submap.Var12);
mapa=mapa0;
if length(mapa)>maxhom;
 ratio=length(mapa)/maxhom;
 mapa=mapa0(1:floor(ratio):end,:);
end

bb=[];
bb(1,:,:)=transpose(mapa(:,2));
bb(2,:,:)=transpose(mapa(:,3));
bb(1,2,:)=repelem(1,length(mapa));
bb(2,2,:)=repelem(1,length(mapa));

r=braid(bb,0);
rc=r.compact;
r.entropy;
rc.entropy;

plot(r,'LineWidth',1);
lenr=length(r);
if lenr==0;
 lenr=length(mapa);
end
 
set(gca,'DataAspectRatio', [1 0.5*lenr/length(mapa) 1]);
saveas(gcf,strcat('./braidplots/',fn,'.',int2str(i),'.pdf'));


l=loopcoords(r);

disp(strcat(num2str(i),'::',fn,': ',num2str(r.entropy)));
res=[res;[i r.entropy]];

end

dlmwrite(strcat('./braidplots/',fn,'.',int2str(i),'.tsv'),res,'\t')

