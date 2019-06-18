function Z = Interface_FK(V,d,lam,r,step,Tend,IC,rho)
%numerical solution of the Fisher-Kolmogorov equation

global N dx D n L R RHO 
N=length(IC);
dx=1;
D=d;
n=step;
L=lam;
R=r;
RHO=rho;
options =odeset('RelTol',1e-4);
sol = ode23s(@fisher,[0 Tend], IC,options);

Y=sol.y;
Y(1,:)=4/3*Y(2,:)-1/3*Y(3,:);
Y(N,:)=4/3*Y(N-1,:)-1/3*Y(N-2,:);
Z=Y(:,end);
    
end


function dy = fisher(~,y)
global N dx D n L R RHO 
dy=zeros(N,1);
y(1)=4/3*y(2)-1/3*y(3);
y(N)=4/3*y(N-1)-1/3*y(N-2);
Uxx = @(y,i) 1/dx^2 * (y(i+1)-2*y(i)+y(i-1));
Ux = @(y,i) 1/(2*dx) * (y(i+1)-y(i-1));
for i=2:N-1
    dy(i) = D* Uxx(y,i) +RHO*Ux(y,i)*(2*y(i)-1)+ L*y(i).*(1-y(i)).^n-R*y(i);
end
end