%----------------------------------------------------------------
% PURPOSE
%	Analysis of a plate.
% unit mm N MPa
%----------------------------------------------------------------
clear all;clc;close all;
%----- Topology -------------------------------------------------
%square plate with n element each side
n_length=100;%number of element in long side
m_width=100;%number of element in wide side
uni=1;%element size
num=0;
for j=1:m_width
    for i=1:n_length
        num=num+1;
        Enode(num,:)=[num (j-1)*(n_length+1)+i (j-1)*(n_length+1)+i+1 (j)*(n_length+1)+i+1 (j)*(n_length+1)+i];
        ex(num,:)=[(i-1)*uni i*uni i*uni (i-1)*uni];%unit mm
        ey(num,:)=[(j-1)*uni (j-1)*uni j*uni j*uni];
    end
end
ndof=3;
Edof=caldof(Enode,ndof);
%----- Element stiffness ----------------------------------------
E=0.5; %2.1e5MPa
v=0.48;%Poison ratio
D=hooke(E,v);%stress=D*strain %% Mechanics of Plate and Shell 
t=1; %thickness of the plate
ep=[t];% thickness of element
nie=size(Enode,1); %number of element
for i=1:nie
Ke(:,:,i)=platre(ex(i,:),ey(i,:),ep,D); %%% Mechanics of Plate and Shell 
end

%----- Assemble Ke into K ---------------------------------------
K=zeros(max(max(Edof))); K=assem(Edof,K,Ke);

%----- load vector f and boundary conditions bc -----------------
f=zeros(max(max(Edof)),1);	
f(18240)=-0.01;% force
num=0;
%find boundary node
for i=0:m_width 
    num=num+1;
    boundary_node(num)=(i)*(n_length+1)+1;
    num=num+1;
    boundary_node(num)=(i+1)*(n_length+1);
end

%map boundary constraint to edof(element degree of freedom)
for i=1:length(boundary_node)
    bc((i-1)*3+1,:)=[3*(boundary_node(i)-1)+1 0]; %0 is displacement
    bc((i-1)*3+2,:)=[3*(boundary_node(i)-1)+2 0];
    bc((i-1)*3+3,:)=[3*(boundary_node(i)-1)+3 0];
%     bc(i,:)=[3*(boundary_node(i)-1)+1 0]; %simple supported
end
%----- Solve the system of equations and compute reactions ------
[a]=solveq(K,f,bc); %return node displacement
Ed=extract(Edof,a); %write disp in element form
%----- Postprocess ----------------------------------------------
%deformation
num=0;
for i=1:(m_width+1)
    for j=1:(n_length+1)
        num=num+1;
        x(i,j)=(j-1)*uni;
        y(i,j)=(i-1)*uni;
        z(i,j)=abs(a((num-1)*3+1));
        XX(num)=x(i,j);
        YY(num)=y(i,j);     
    end
end
figure
surf(x,y,z)
title('Displacement')


%internel force
for i=1:nie
    F_node(i,:)=Ke(:,:,i)*Ed(i,:)'; % element internel force
end
%plot internal force with patch(matlab commmand)
Shear=zeros((m_width+1)*(n_length+1),1);
Mx=zeros((m_width+1)*(n_length+1),1);
My=zeros((m_width+1)*(n_length+1),1);
right_line_node=[(n_length+1):(n_length+1):m_width*(n_length+1)];
bottom_line_node=[(m_width*(n_length+1)+1):((m_width+1)*(n_length+1)-1)];
last_node=(m_width+1)*(n_length+1);
for i=1:nie
        Shear(Enode(i,2))=F_node(i,1);
        Mx(Enode(i,2))=F_node(i,2);
        My(Enode(i,2))=F_node(i,3);
end
for i=1:length(right_line_node)
    Shear(right_line_node(i))=-F_node(n_length+(i-1)*n_length,4);
    Mx(right_line_node(i))=F_node(n_length+(i-1)*n_length,5);
    My(right_line_node(i))=F_node(n_length+(i-1)*n_length,6);
end
for i=1:length(bottom_line_node)
    Shear(bottom_line_node(i))=-F_node(n_length+(m_width-2)*n_length+i,10);
    Mx(bottom_line_node(i))=F_node(n_length+(m_width-2)*n_length+i,11);
    My(bottom_line_node(i))=F_node(n_length+(m_width-2)*n_length+i,12);
end
Shear(last_node)=F_node(nie,7);
Mx(last_node)=F_node(nie,8);
My(last_node)=F_node(nie,9);
    
        
figure
title('Shear')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',Shear,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;

figure
title('Mx')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',Mx,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;
figure
title('My')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',My,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;
% %------------------------ end -----------------------------------

dlmwrite('f.csv',z,'delimiter',',','precision',6)










