function Interface_Prolif(V,rm,rp,rd)


%% To compute 2D proliferation assay simulations
% Inputs: V (compartment size, defined as m in the paper), rm (motility
% rate), rp (proliferation rate), rd (death rate)



n=120; %number of lattice points in each direction
vox=n/V;
P1=rm; %motility rate
P2=rp; %proliferation rate
P3=rd; %death rate
step = 1; %how far do daughter cells go


Tend=10; %NB NONdimensional time
k=1;
D=reshape(1:n^2,n,n);
B=zeros(n/vox,n/vox,vox^2);
for i=1:vox
    for j=1:vox
        B(:,:,k)=D(1+n/vox*(j-1):n/vox*j,1+n/vox*(i-1):n/vox*i);
        k=k+1;
    end
end
Y=100; %number of realisations
m=0.05; % initial agent density

% Gillespie loop
parfor p=1:Y
    
    C=double(rand(n,n)<m);
    
    
    Q0=sum(sum(C));
    while Q0==0
        C=double(rand(n,n)<m);
        Q0=sum(sum(C));
    end
    Q=Q0;
    
    
    if n/vox>1
        Z=zeros(n/vox,n/vox,vox^2);
        k=1;
        
        for i=1:vox
            for j=1:vox
                Z(:,:,k)=C(1+n/vox*(j-1):n/vox*j,1+n/vox*(i-1):n/vox*i);
                k=k+1;
            end
        end
        
        C=(reshape(sum(sum(Z)),vox,vox,1));
    end
    

    T=0;
    j=1;
    tau=0;
    pp=1;
    
    
    JJ=randi([1,vox^2],1,10000);
    W=rand(1,10000);
    W2=rand(1,10000);
    W3=rand(1,10000);
    y=1;
    Qend=Q0;
    while T<Tend && Qend<n^2 && Qend>0
        
        tau(j+1)=tau(j)+(P2-P3)*log(1/W2(y))/((P1+P2+P3)*Q(j));
        R=(P1+P2+P3)*Q(j)*rand;
        Q(j+1)=Q(j);
        %find a random occupied site
        while  W3(y) < ( 1 - C(JJ(y))/V^2 )
            y=y+1;
            if y==10001
                JJ=randi([1,vox^2],1,10000);
                W=rand(1,10000);
                W2=rand(1,10000);
                W3=rand(1,10000);
                y=1;
            end
        end
        J=JJ(y);
        
        
        if R<P1*Q(j)
            % cell movement
            I=ceil(W(y)*4);
            
            if I==1 && J+vox<=vox^2 && C(J+vox)<(n/vox)^2 && rand<(1-C(J+vox)*vox^2/n^2)
                C(J+vox)=C(J+vox)+1;
                C(J)=C(J)-1;
            elseif I==2  && mod(J,vox)~=0  && C(J+1)<(n/vox)^2 && rand<(1-C(J+1)*vox^2/n^2)
                C(J+1)= C(J+1)+1;
                C(J)=C(J)-1;
            elseif I==3 && mod(J,vox)~=1  && C(J-1)<(n/vox)^2 && rand<(1-C(J-1)*vox^2/n^2)
                C(J-1)= C(J-1)+1;
                C(J)=C(J)-1;
            elseif I==4 && J-vox>=1 && C(J-vox)<(n/vox)^2 && rand<(1-C(J-vox)*vox^2/n^2)
                C(J-vox)=C(J-vox)+1;
                C(J)=C(J)-1;
            else %movement does not happen
            end
            
        elseif R<(P1+P2)*Q(j)
            % cell proliferation
            
          
                I=ceil(W(y)*4);
                
                
                if I==1
                    if rand<(vox/n)^pp
                        if J+vox*step<=vox^2 && rand<(1-C(J+vox*step)*vox^2/n^2)
                            
                            C(J+vox*step)=C(J+vox*step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    else
                        if  rand<(1-C(J)*vox^2/n^2)
                            C(J)=C(J)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    end
                    
                elseif I==2
                    if rand<(vox/n)^pp
                        if ceil(J/vox)==ceil((J+step)/vox) && rand<(1-C(J+step)*vox^2/n^2)
                            
                            C(J+step)=C(J+step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    else
                        if  rand<(1-C(J)*vox^2/n^2)
                            C(J)=C(J)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    end
                    
                elseif I==3
                    if rand<(vox/n)^pp
                        if  ceil(J/vox)==ceil((J-step)/vox) && rand<(1-C(J-step)*vox^2/n^2)
                            
                            C(J-step)=C(J-step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    else
                        if  rand<(1-C(J)*vox^2/n^2)
                            C(J)=C(J)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    end
                elseif I==4
                    if rand<(vox/n)^pp
                        if J-vox*step>=1 && rand<(1-C(J-vox*step)*vox^2/n^2)
                            
                            C(J-vox*step)=C(J-vox*step)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
                    else
                        if  rand<(1-C(J)*vox^2/n^2)
                            C(J)=C(J)+1;
                            Q(j+1)=Q(j+1)+1;
                        end
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
            JJ=randi([1,vox.^2],1,10000);
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
    G2=griddedInterpolant(tau,Q/n^2,'nearest');
    Qi(p,:)=G2(X);
end
QQ=mean(Qi,1);

X=linspace(0,Tend,1000);
Fish=@(x,y) fisher(x,y,P2,P3);
sol = ode23s(@(x,y) Fish(x,y),[0 Tend], m);
G=griddedInterpolant(sol.x,sol.y,'linear');

figure(55)
plot(X,G(X),'k-.','LineWidth',2)
hold on
plot(X,QQ,'g-.','LineWidth',2)
set(gca,'FontSize',18)
axis([0 Tend 0 1])
xlabel('T')
ylabel('C')


end
function dy = fisher(~,y,P2,P3)
%Logistic equation (Equation 5 in Section 3.1 of paper)
dy = y.*((P2/(P2-P3)*(1-y))-P3/(P2-P3));
end

