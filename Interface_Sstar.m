function Fmax=Interface_Sstar(rm,rp,rd,epsilon,q)

% Computes the threshold correlation radius s*(T), defined in Section 3.1,
% Equation 8. It takes inputs rm (motility rate), rp (proliferation rate),
% rd (death rate), epsilon (correlation tolerance; see Section 3.l), and q
% (initial cell density). It outputs the maximum of s*(T), denoted as
% sigma* in Equation 9.


Y=100;    %number of realisations
    P1=rm; %motility rate
    P2=rp; %proliferation rate
    P3=rd; %death rate
    step = 1; %how far do daughter cells go
    Tend=8; %NB NONdimensional time
parfor p=1:Y
    disp(p)
    m=q; % initial seeding
    n=120; %lattice nodes in each direction
    C=double(rand(n,n)<m);
    Q0=sum(sum(C));
    while Q0==0
        C=double(rand(n,n)<m);
        Q0=sum(sum(C));
    end
    Q=Q0;
  
    T=0;
    j=1;
    tau=0;
    
     Qend=Q0;
    
    JJ=randi([1,n^2],1,10000);
    W=rand(1,10000);
    y=1;
    FAA=1;
    while T<Tend && Qend<n^2 && Qend>0
        
        tau(j+1)=tau(j)+(P2-P3)*log(1/rand)/((P1+P2+P3)*Q(j));
        R=(P1+P2+P3)*Q(j)*rand;
        Q(j+1)=Q(j);
        %find a random occupied site
        while C(JJ(y))==0
            y=y+1;
            if y==10001
                JJ=randi([1,n.^2],1,10000);
                W=rand(1,10000);
                y=1;
            end
        end
        J=JJ(y);
        
        
        if R<P1*Q(j)
            % cell movement
            I=ceil(W(y)*4);
            
            C(J)=0;
            if I==1 && J+n<=n^2 && C(J+n)==0
                C(J+n)=1;
            elseif I==2  && mod(J,n)~=0  && C(J+1)==0
                C(J+1)=1;
            elseif I==3 && mod(J,n)~=1  && C(J-1)==0
                C(J-1)=1;
            elseif I==4 && J-n>=1 && C(J-n)==0
                C(J-n)=1;
            else
                C(J)=1; %movement does not happen
            end
            
        elseif R<(P1+P2)*Q(j)
            % cell proliferation
           
                I=ceil(W(y)*4);
                
                
                if I==1 && J+n*step<=n^2 && max(C(J+n:n:J+n*step))==0
                    C(J+n*step)=1;
                    Q(j+1)=Q(j+1)+1;
                elseif I==2 && ceil(J/n)==ceil((J+step)/n) && max(C(J+1:J+step))==0
                    C(J+step)=1;
                    Q(j+1)=Q(j+1)+1;
                elseif I==3 && ceil(J/n)==ceil((J-step)/n) && max(C(J-step:J-1))==0
                    C(J-step)=1;
                    Q(j+1)=Q(j+1)+1;
                elseif I==4 && J-n*step>=1 && max(C(J-n*step:n:J-n))==0
                    C(J-n*step)=1;
                    Q(j+1)=Q(j+1)+1;
                end

        else
            %cell death
            C(J)=0;
            Q(j+1)=Q(j+1)-1;
        end
        
        T=tau(j+1);
        Qend=Q(j+1);
        
        y=y+1;
        if y==10001
            JJ=randi([1,n.^2],1,10000);
            W=rand(1,10000);
            y=1;
        end
        %compute correlation function
        
        S=1;
        E=epsilon;
        %use shift functions to calculate FAA
        Q2=C.*circshift(C,[S 0]) + C.*circshift(C,[0 S]) + ...
            C.*circshift(C,[-S 0])+ C.*circshift(C,[0 -S]);
        while sum(sum(Q2(1+S:end-S,1+S:end-S)))/Qend^2*n^4/(4*n*(n-S)) > 1+E
            S=S+1;
            Q2=C.*circshift(C,[S 0]) + C.*circshift(C,[0 S]) + ...
                C.*circshift(C,[-S 0])+ C.*circshift(C,[0 -S]);
            if S>n/2
                break
            end
        end
        FAA(j+1)=S;
        j=j+1;
    end
   
    X=linspace(0,Tend,1000);
    [tau,U]=unique(tau,'first');
    Q=Q(U);
    G2=griddedInterpolant(tau,Q/n^2);
    Qi(p,:)=G2(X);
    G=griddedInterpolant(tau,FAA);
    F(p,:)=G(X);
end
X=linspace(0,Tend,1000);
FF=mean(F,1);
Fmax=max(FF);
figure(50)
plot(X,FF,'m','LineWidth',2)
hold on
 axis([0 8 0 5])
 xlabel('Nondimensional Time T')
 ylabel('Maximum Correlation Distance s*(T)')
end


