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
m1=testCirc(m1,12,41,'cavityBubble');
m1=vtkMeshOut(m1,workDir,'teste')

