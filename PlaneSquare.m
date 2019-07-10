%This is an example FE model of a plane square structure
%composed of triangle elements

addpath(genpath('Functions'));

%% Load and set parameters
% load('nodeData')
% Node=NodeBeforeRecenter;


Node=[0 0; 1 0; 1 1;0 1];
DT = delaunayTriangulation(Node);
Element=DT.ConnectivityList;
NodeOnPressureFace=logical([0;1; 1; 0]);
NodeOnfixedBoundary=logical([1;0; 0; 1]);


E=10;   %Modulus
u=0.3;  %Poisson's ratio
pressure=-0.06; %Horizontal pressure
Em=1/E*[1, -u ,0; -u 1 0; 0 0 2+2*u];

%% Construct stiffness matrix
numE=size(Element,1);
DoF=6;  %Degree of freedom per element



% A=5;

Ke=zeros(DoF,DoF,numE);
for n=1:numE
    Ke(:,:,n)=PlaneTri2Stiffness( [Node(Element(n,1),:),Node(Element(n,2),:),Node(Element(n,3),:)]', Em);
end



eft=zeros(numE,DoF);
for i=1:numE
   eft(i,:)=[2*Element(i,1)-1,2*Element(i,1),2*Element(i,2)-1,2*Element(i,2),2*Element(i,3)-1,2*Element(i,3)]; %3 node per element
end

Kglob=zeros(size(Node(:),1));

for n=1:numE
   for i=1:DoF
       for j=1:DoF
           Kglob(eft(n,i),eft(n,j))= Kglob(eft(n,i),eft(n,j)) + Ke(i,j,n);
       end
   end   
end





%% Construct BCs

stiffMsize=size(Kglob,1);
%zero DBC nodes
temp=[2*find(NodeOnfixedBoundary)-1,2*find(NodeOnfixedBoundary)]';
DBCnodeIndex=(temp(:));

%Construct force BC
allEdges=[Element(:,1),Element(:,2); %Big edge matrix
          Element(:,2),Element(:,3);
          Element(:,3),Element(:,1);];
   Pedge=[] ;   %Edges with pressure
for i=1:size(allEdges,1)
   if NodeOnPressureFace(allEdges(i,1)) & NodeOnPressureFace(allEdges(i,2))
    Pedge=[Pedge;allEdges(i,:)];
   end
end



applyF=zeros(stiffMsize,1);

for i=1:size(Pedge,1)
edgeNode=Pedge(i,:);
n1=Node(edgeNode(1),:);
n2=Node(edgeNode(2),:);
edgeL=sqrt(sum((n2-n1).^2));
%Test simple pressure condition
applyF(2*edgeNode(1)-1)=1/2*edgeL*pressure ;       %Node 1 fx
applyF(2*edgeNode(1))=0    ;        %Node 1 fy

applyF(2*edgeNode(2)-1)= 1/2*edgeL*pressure ;     %Node 2 fx
applyF(2*edgeNode(2))= 0   ;      %Node 2 fy

end

%Apply DBC
%Apply 0 to ux1 uy1 and uy2

% zeroDBC=false(stiffMsize,1);
% zeroDBC(DBCnodeIndex)=true;

I=eye(stiffMsize);
for i=DBCnodeIndex
Kglob(i,:)=I(i,:);
Kglob(:,i)=I(:,i);
end

% SolvedDisp=(applyF/Kglob);
%output Format 1
 SolvedDisp=((Kglob)\applyF);
DisplacementList=[Node,SolvedDisp(1:2:end),SolvedDisp(2:2:end)];



% HolderX=dispMap(:,:,1);
% HolderY=dispMap(:,:,2);
% dispMagnitude=sqrt(HolderX.^2+HolderY.^2);
% surf(HolderY)
% title('displacement magnitude')



% tri = delaunay(DisplacementList(:,1),DisplacementList(:,2));
% z = DisplacementList(:,4);
% trisurf(tri,DisplacementList(:,1),DisplacementList(:,2),DisplacementList(:,3))

%%
figure;plot(Node(NodeOnfixedBoundary,1),Node(NodeOnfixedBoundary,2),'x','MarkerSize',15,'LineWidth',3)
hold on
plot(Node(NodeOnPressureFace,1),Node(NodeOnPressureFace,2),'o','MarkerSize',15,'LineWidth',3)

triplot(Element,Node(:,1),Node(:,2),'k') 
legend('Fixed Node','Forced Node','Bar')
hold off
xlim([-0.5,2])
ylim([-0.5,2])
title('Before')


NodeNew=Node+DisplacementList(:,3:4);
figure;plot(NodeNew(NodeOnfixedBoundary,1),NodeNew(NodeOnfixedBoundary,2),'x','MarkerSize',15,'LineWidth',3)
hold on
plot(NodeNew(NodeOnPressureFace,1),NodeNew(NodeOnPressureFace,2),'o','MarkerSize',15,'LineWidth',3)

triplot(Element,NodeNew(:,1),NodeNew(:,2),'k') 
legend('Fixed Node','Forced Node','Bar')
hold off

xlim([-0.5,2])
ylim([-0.5,2])
title('After')
