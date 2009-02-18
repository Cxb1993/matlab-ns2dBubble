function s = stepInterface(s,dt,comp,steptype,Re,Sc,We)

%modelo de documentacao a partir de:
%http://www.engin.umd.umich.edu/CIS/course.des/cis400/matlab/oop.html

%SIMULATOR simulator class constructor.
%   s = Simulator(m) creates a simulator object from the mesh object

%Name: step
%Location: <path>/@Simulator
%Purpose: this is the main program,

% modificado em 13/03/2006
% revisado   em 09/04/2007

IEN = getIEN(s.m);
X= getX(s.m);
Y=getY(s.m);
Z=getZ(s.m);

nelem=size(IEN,1);
nvert=size(X,1)-nelem;
nnodes=size(X,1);

velu=s.us;
velv=s.vs;
velc=s.cs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  alpha = 0   -> explicito                                     %
%  alpha = 0.5 -> crank-nicholson                               %
%  alpha = 1   -> implicito                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

alpha=1;

%% coeficiente do laplaciano - viscosidade e concentracao
k=1/Re;
kc=0.5;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% metodo acoplado                                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'coupled'))

    Mclump=diag(sparse(sum(s.Mc,2)));
    matc=1/dt*Mclump+alpha*kc*s.Kc;
    Ac=matc;
    b1c=vc;
    [Atc b1c] = setCoupledCBC(s,Ac,b1c);
    c=Atc\b1c;
    cs=c;

end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Metodo da Projecao discreto baseado em decomposicao LU        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if(strcmp(steptype,'uncoupled'))

    if(comp)
        Mclump=diag(sparse(sum(s.Mc,2)));
        matc=1/dt*Mclump+alpha*kc*s.Kc;
        [Atc b1c ipc] = setUncoupledCBC(s,matc,velc*0);

        rc = symrcm(-Atc);
        Atc=Atc(rc,rc);
        Rc = cholinc(-Atc,1e-4);

        s.Atc=Atc;
        s.b1c=b1c;
        s.ipc=ipc;
        s.Rc=Rc;
        s.rc=rc;

    else
        Atc=s.Atc;
        b1c=s.b1c;
        ipc=s.ipc;
        Rc=s.Rc;
        rc=s.rc;
    end;
    
    
    
	% retorna os indices com s.cs=0.5
    surface=find(s.cs==0.5);
	% procura
	% X(surface) sao os valores de X para os vertices em surface, ou seja,
    % sao os valores de X dos vertices onde s.cs=0.5
	%  - vetor com mesma dimensao de surface
	% Y(surface) sao os valores de Y para os vertices em surface, ou seja,
    % sao os valores de Y dos vertices onde s.cs=0.5
    
	% closer - vetor com mesma dimensao de X e Y
    % closer - para cada vertice da malha procura-se o vertice da bolha
    % mais proximo retornando o numero do vertice
    closer=surface(dsearchn([X(surface) Y(surface)],[X Y]));
    distance=sqrt(((X-X(closer)).*(X-X(closer)))+(Y-Y(closer)).*(Y-Y(closer)));
    distance=distance(1:nvert).*(s.cs-0.5)*2;
    

    cs=s.cs;
    Mclump=diag(sparse(sum(s.Mc,2)));
    for i=1:1
        vc=((1/dt)*Mclump-(1-alpha)*kc*s.Kc)*[distance];
        b1c=b1c+vc.*ipc;
        cs(rc,1) = pcg(Atc,b1c(rc),1e-6,20,Rc',Rc);
    end;

end;

% invM=diag(sparse(1./sum(s.M,2)));
% g = invM*s.G*cs;
% %g=s.M\(s.G*cs);
% modG=sqrt(g(1:nnodes).*g(1:nnodes)+g(1+nnodes:2*nnodes).*g(1+nnodes:2*nnodes));
% gn(1:nnodes,1)=g(1:nnodes)./(modG+10e-2);
% gn(nnodes+1:2*nnodes,1)=g(nnodes+1:2*nnodes)./(modG+10e-2);
 invMclump=diag(sparse(1./sum(s.Mc,2)));
% kappa1=invMclump*s.D*gn;
% 
% 
% g=s.M\(s.G*cs);
% modG=sqrt(g(1:nnodes).*g(1:nnodes)+g(1+nnodes:2*nnodes).*g(1+nnodes:2*nnodes));
% gn(1:nnodes,1)=g(1:nnodes)./(modG+10e-2);
% gn(nnodes+1:2*nnodes,1)=g(nnodes+1:2*nnodes)./(modG+10e-2);
% kappa=s.Mc\(s.D*gn);

%kappa3=-(s.Kc\cs)./(modG(1:nvert));
%kappa3=s.Mc\s.Kc*cs;
kappa=invMclump*(s.Kc*cs);

% kappa1=invMclump*(s.Kc*distance);
% kappa2=s.Mc\(s.Kc*distance);
% 
% teste1=s.Mc\(s.Kc*X(1:nvert));
% teste2=s.Mc\(s.Kc*ones(nvert,1));
% teste3=s.M(1:nnodes,1:nnodes)\(s.K(1:nnodes,1:nnodes)*X);
% teste4=s.M(1:nnodes,1:nnodes)\(s.K(1:nnodes,1:nnodes)*Y);

%kappa=0.5*(kappa+kappa1);




% 
%  vc=((1/dt)*Mclump-(1-alpha)*kc*s.Kc)*kappa;
%         b1c=b1c+vc.*ipc;
%         kappa(rc,1) = pcg(Atc,b1c(rc),1e-6,20,Rc',Rc);

%surface=find(s.cs==0.5);
%dsearch(X(surface),Y(surface),X,Y);
% plot3(X(surface),Y(surface),kappa1(surface),X(surface),Y(surface),kappa(surface)...
%    , X(surface),Y(surface),0.5*(kappa1(surface)+kappa(surface)));

plot3(X(surface),Y(surface),kappa(surface));


kappas=kappa(surface(dsearchn([X(surface) Y(surface)],[X Y])));

%kappa = setCentroid(s,[kappa;zeros(nnodes-nvert,1);kappa;zeros(nnodes-nvert,1)]);

% 
%  vc=((1/dt)*Mclump-(1-alpha)*kc*s.Kc)*kappas;
%         b1c=b1c+vc.*ipc;
%         kappas(rc,1) = pcg(Atc,b1c(rc),1e-6,20,Rc',Rc);

s.kappa=kappas*We;


