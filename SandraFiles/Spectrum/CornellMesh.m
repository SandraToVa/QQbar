function [mesh,x] = CornellMesh(~)
% the endpoints of the integration interval:
system.a=0.001;   
system.b=22; 
% parameters of the boundary conditions:
system.A1= eye(3);
system.A2= zeros(3,3);
system.B1= eye(3);
system.B2= zeros(3,3);
% function handle to the function returning the potential matrix
system.V=@potentialMatrix;

tol=1.5e-8; %for sd->p/p1->s hybrid x=173 (bottom)
%construct a mesh corresponding to the user input tolerance tol
meshData = constructMesh(system,tol);

mesh=meshData;

x=[system.a system.a+cumsum(mesh.h)];
disp(x);

fid = fopen('x.txt','w');
fprintf(fid, '%.16e\n', x);
fclose(fid);

end

function r=potentialMatrix(x) % returns the potential matrix evaluated in x

r = zeros(3,3,4);

for i=1:4 
  % CAS QUARKONIUM S
  r(1,1,i)=R(x(i)); 
  r(2,2,i)=R(x(i)); 
  r(3,3,i)=R(x(i)); 
  r(1,2,i)=0;  
  r(2,1,i)=0;
  
end

end

function f1=Vg(x)
  %Eo=2.6984; %charm
  Eo=9.5325; %bottom   
  %%%%
  k_g=0.489;   %Full potential
  %k_g=0;   %Long distances
  %%%%
  f1=-k_g/x+0.187*x+Eo; 
end

function M2=R(x)
  %m=1.47; %charm
  m=4.88; %bottom
  M2=m*Vg(x);
end