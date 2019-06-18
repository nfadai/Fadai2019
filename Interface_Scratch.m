function Interface_Scratch(V,rm,rp,rd,rho,tf)

%% To compute 2D scratch assay simulations
% Inputs: V (compartment size, defined as m in the paper), rm (motility
% rate), rp (proliferation rate), rd (death rate), rho (movement bias,
% shown in SI; rho=0 retrieves standard scratch assay shown in Section
% 3.2), tf (final dimensional time)

n=240; %lattice nodes in x-direction
m=24; %lattice nodes in y-direction
vox=n/V;
vox2=m/V;
P1=rm; %motility rate
P2=rp; %proliferation rate
P3=rd; %death rate
step = 1; %how far do daughter cells go
alpha = 1/V; %probability of proliferation/migration ending up in neighbour voxel
Tend=tf; %dimensional time

Y=100; %number of realisations

CC=zeros(m,n);
CC(:,51:70)=1; % Initial condition of scratch assay

Q0=sum(sum(CC));


C00=CC;
if min(n/vox,m/vox2)>1
    Z=zeros(m/vox2,n/vox,vox*vox2);
    k=1;
    for i=1:vox
        for j=1:vox2
            Z(:,:,k)=CC(1+m/vox2*(j-1):m/vox2*j,1+n/vox*(i-1):n/vox*i);
            k=k+1;
        end
    end
    
    CC=(reshape(sum(sum(Z)),vox2,vox,1));
end


parfor p=1:Y
    

    T=0;
    j=1;
    tau=0;
    Q=Q0;
    C=CC;
    
    
    JJ=randi([1,vox*vox2],1,10000);
    W=rand(1,10000);
    W2=rand(1,10000);
    W3=rand(1,10000);
    y=1;
    Qend=Q0;
    while T<Tend && Qend<n*m && Qend>0
        
        tau(j+1)=tau(j)+log(1/W2(y))/((P1+P2+P3)*Q(j));
        R=(P1+P2+P3)*Q(j)*rand;
        Q(j+1)=Q(j);
        % find a random occupied voxel, weighted to maximum capacity
        % to find a uniformly random particle 
        while  W3(y) < ( 1 - C(JJ(y))/V^2 )
            y=y+1;
            if y==10001
                JJ=randi([1,vox*vox2],1,10000);
                W=rand(1,10000);
                W2=rand(1,10000);
                W3=rand(1,10000);
                y=1;
            end
        end
        J=JJ(y);
        
        
        if R<P1*Q(j)
            % cell movement
            I=W(y)*4;
            if I<=1  &&rand<alpha&& mod(J,vox2)~=0  && C(J+1)<(n/vox)*(m/vox2) && rand<(1-C(J+1)*vox*vox2/(n*m))
                C(J+1)= C(J+1)+1;
                C(J)=C(J)-1;
            elseif I<=2  &&rand<alpha&& I>1 && mod(J,vox2)~=1  && C(J-1)<(n/vox)*(m/vox2) && rand<(1-C(J-1)*vox*vox2/(n*m))
                C(J-1)= C(J-1)+1;
                C(J)=C(J)-1;
            elseif I<=3-rho  &&rand<alpha&& I>2 && J-vox2>=1 && C(J-vox2)<(n/vox)*(m/vox2) && rand<(1-C(J-vox2)*vox*vox2/(n*m))
                C(J-vox2)=C(J-vox2)+1;
                C(J)=C(J)-1;
            elseif I>3-rho  &&rand<alpha&& J+vox2<=vox*vox2 && C(J+vox2)<(n/vox)*(m/vox2) && rand<(1-C(J+vox2)*vox*vox2/(n*m))
                C(J+vox2)=C(J+vox2)+1;
                C(J)=C(J)-1;
            else %movement does not happen
            end
            
        elseif R<(P1+P2)*Q(j)
            % cell proliferation
            
    
                if rand<alpha
                    I=ceil(W(y)*4);
                    
                    
                    if I==1
                        
                        if J+vox2*step<=vox*vox2 && max(C(J+vox2:vox2:J+vox2*step))<(n/vox)*(m/vox2) &&...
                                rand<(1-C(J+vox2*step)*vox*vox2/(n*m))
                            C(J+vox2*step)=C(J+vox2*step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                        
                    elseif I==2
                        
                        if ceil(J/vox2)==ceil((J+step)/vox2) && max(C(J+1:J+step))<(n/vox)*(m/vox2) &&...
                                rand<(1-C(J+step)*vox*vox2/(n*m))
                            C(J+step)=C(J+step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                        
                    elseif I==3
                        
                        if  ceil(J/vox2)==ceil((J-step)/vox2) && max(C(J-step:J-1))<(n/vox)*(m/vox2)&&...
                                rand<(1-C(J-step)*vox*vox2/(n*m))
                            C(J-step)=C(J-step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                        
                    elseif I==4
                        
                        if J-vox2*step>=1 && max(C(J-vox2*step:vox2:J-vox2))<(n/vox)*(m/vox2) &&...
                                rand<(1-C(J-vox2*step)*vox*vox2/(n*m))
                            C(J-vox2*step)=C(J-vox2*step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    end
                else
                    if  C(J)<(n/vox)*(m/vox2) &&...
                            rand<(1-C(J)*vox*vox2/(n*m))
                        C(J)=C(J)+1;
                        Q(j+1)=Q(j+1)+1;
                    end
                end
                
            
        else
            %cell death
            C(J)=C(J)-1;
            Q(j+1)=Q(j+1)-1;
        end
        
        T=tau(j+1);
        Qend=Q(j+1);
        
        y=y+1;
        if y==10001
            JJ=randi([1,vox2*vox],1,10000);
            W=rand(1,10000);
            W2=rand(1,10000);
            W3=rand(1,10000);
            y=1;
        end
        
        j=j+1;
        
    end
    
   
    X=linspace(0,Tend,1000);
    [tau,U]=unique(tau,'first');
    Q=Q(U);
  
    Fi(:,:,p)=C;
    
    
    G2=griddedInterpolant(tau,Q/(n*m),'nearest');
    Qi(p,:)=G2(X);
end
QQ=mean(Qi,1);
F2=mean(Fi,3);
FF=mean(F2,1)'/V^2;

figure(515)
W=(1+V)/2:V:n+1-(1+V)/2;
plot(W,FF,'Color',[0 2/3 1],'LineWidth',2)
hold on

Z=Interface_FK(V,P1*V^2*alpha/(4),P2,P3,step,Tend,...
    mean(C00,1),(rho*alpha)*P1*V/(2)); %numerical solution of continuum limit

plot(1:240,Z,'k--','LineWidth',1)
set(gca,'FontSize',18)
xlabel('x')
ylabel('C(x,T), <Cm(x,T)>')
axis([0 240 0 1])
xticks(0:60:240);
yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})


end

