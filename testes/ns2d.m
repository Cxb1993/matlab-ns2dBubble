path('/home/gustavo/projects/matlab/ns2dBubble',path)
workDir='/home/gustavo/projects/';
workDirVtk='/home/gustavo/projects/';

Re=10;
Sc=4;
We=0.25;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% utilizacao:                                               %             
% test: test (step,cavity,couette)                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

m1=Model2d();
%m1=test(m1,29,15,'cavityBubble');
%m1=testCirc(m1,8,21,'cavityBubble');
m1=testCirc(m1,8,31,'cavityBubble');
m1=vtkMeshOut(m1,workDirVtk,'oscillatingBubble');

show(m1);


s1=Simulator2d(m1);
s1=init(s1);


cfl=1;
dt=cfl*sqrt((max(m1.Y)-min(m1.Y))*(max(m1.X)-min(m1.X))/s1.nvert)/max([s1.us;1])*10;

i=1
s1=stepInterface(s1,dt,true,'uncoupled',Re,Sc,We);
s1=step(s1,dt,true,'uncoupled',Re,Sc);
%saveDump(s1,workDir,'sim',i)
%saveSol(s1,workDir,'sim',i)

for i=2:200
    i
    s1=stepInterface(s1,dt,true,'uncoupled',Re,Sc,We);
    s1=step(s1,dt,true,'uncoupled',Re,Sc);
    if(mod(i,2)==0) 
        s1=update(s1);
    end;
    %saveSol(s1,workDir,'sim',i)
    %showScalar(s1,s1.kappa)
    %vtkCompleteOut(s1,workDirVtk,'field',i)
    %vort(s1);
    show(s1)
end;