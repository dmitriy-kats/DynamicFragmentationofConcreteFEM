
clear; clc; close all;

eC=[10 20 25 30 50 60 80 100 150 250 500 1500 1800 2500 3500 5000 7500 10000];
NC=eC.*0;
NfC=eC.*0;

for ee=1:length(eC)
elem=eC(ee);

Initialization_concrete_meshconv


for n=1:nsteps   
    %Performs the Newmark Method
    uu1=uu+dt*vv+.5*dt*dt*aa;
    uuc=uuc+dt*vvc+.5*dt*dt*aac;
    uu1(1)=-strainrate*lengthB*n*dt/2; %LHS BC   
    uu1(end)=strainrate*lengthB*n*dt/2; %RHS BC  
    delta=uuc-uu1;    
    uu=uu1;

    for ii=1:elem  %Force Calculation
        if cohTracker(ii)==true
            if delta(ii)>=deltamax(ii) && delta(ii)<=deltaC(ii) && DD(ii)<1 && delta(ii)>=0
                deltamax(ii)=delta(ii);
                TT(ii)=SigmaC(ii)*(1-delta(ii)/deltaC(ii));
                TTmax(ii)=TT(ii);
                DD(ii)=min(deltamax(ii)/deltaC(ii),1);
                
            elseif delta(ii)>=0 && delta(ii)<deltamax(ii) && DD(ii)<1
                TT(ii)=TTmax(ii)*delta(ii)/deltamax(ii);       % (eq. 4b: IJSS 2015 Vocialta-Molinari)
                DD(ii)=min(deltamax(ii)/deltaC(ii),1);
                
            elseif delta(ii)<=0  || DD(ii)==1
                TT(ii)=0;
            elseif delta(ii)>=deltaC(ii)
                TT(ii)=0;
                TTmax(ii)=SigmaC(ii); 
                deltamax(ii)=deltaC(ii);
                DD(ii)=min(deltamax(ii)/deltaC(ii),1);
            else
                delta(ii)
                deltamax(ii)
                deltaC(ii)
                DD(ii)
                
                sprintf('Error?!')
            end
            
            NodeForceCoh(ii)=-TT(ii)*A+EE(ii)*A*(uu(ii+1)-uuc(ii))/h;
        end
        
        if ii==1 
            NodeForce(1)=EE(ii)*A*(uu(2)-uu(1))/h; 
            if cohTracker(elem)==false || n==1
                NodeForce(end)=-EE(elem)*A*(uu(end)-uu(end-1))/h;
            elseif cohTracker(elem)==true
                NodeForce(end)=-EE(elem)*A*(uu(end)-uuc(end-1))/h;
            end
        elseif cohTracker(ii-1)==false && cohTracker(ii)==false
            NodeForce(ii)=-EE(ii-1)*A*(uu(ii)-uu(ii-1))/h+EE(ii)*A*(uu(ii+1)-uu(ii))/h;
            AvgForce(ii)=(EE(ii-1)*A*(uu(ii)-uu(ii-1))/h+EE(ii)*A*(uu(ii+1)-uu(ii))/h)/2;
        elseif cohTracker(ii-1)==false &&cohTracker(ii)==true
            NodeForce(ii)=-EE(ii-1)*A*(uu(ii)-uu(ii-1))/h+TT(ii)*A;
            AvgForce(ii)=(EE(ii-1)*A*(uu(ii)-uu(ii-1))/h+TT(ii)*A)/2;
        elseif cohTracker(ii-1)==true && cohTracker(ii)==false  
            NodeForce(ii)=-EE(ii-1)*A*(uu(ii)-uuc(ii-1))/h+EE(ii)*A*(uu(ii+1)-uu(ii))/h;
            AvgForce(ii)=(EE(ii-1)*A*(uu(ii)-uuc(ii-1))/h+EE(ii)*A*(uu(ii+1)-uu(ii))/h)/2;
        elseif cohTracker(ii-1)==true &&cohTracker(ii)==true
            NodeForce(ii)=-EE(ii-1)*A*(uu(ii)-uuc(ii-1))/h+TT(ii)*A;
            AvgForce(ii)=(EE(ii-1)*A*(uu(ii)-uuc(ii-1))/h+TT(ii)*A)/2;
        else
            sprintf('Error?')
        end
        
    end
    

    %Newmark and BoundaryConditions
    aa1=NodeForce./mm;
    aac1=NodeForceCoh./mmc;
    aa1(1)=0; aac1(1)=0;         % extreme nodes have a=0 (v cte)
    aa1(end)=0; aac1(end)=0;     % extreme nodes have a=0
    vv1=vv+0.5*dt*(aa+aa1);
    vvc=vvc+0.5*dt*(aac+aac1);
    vv1(1)=-strainrate*lengthB/2; vvc1(1)=-strainrate*lengthB/2;         %  extreme nodes have v = cte = v0
    vv1(end)=strainrate*lengthB/2; vvc1(end)=strainrate*lengthB/2;       %  extreme nodes have v = cte = v0
    vv=vv1;
    aa=aa1;
    aac=aac1;
    
    t(n+1)=dt*n;
    KE(n+1)=1/2*sum(mm.*vv.^2);

    ExtWork(n+1)=ExtWork(n)+NodeForce(end)*vv1(end)*dt + NodeForce(1)*(vv1(1))*dt;   
       
    for ii=1:elem
        PE(n+1)=PE(n+1)+1/2*A*EE(ii)*(uu(ii+1)-uu(ii))^2/h;
        if cohTracker(ii)==true 
                E_diss(n+1)=E_diss(n+1)+0.5*SigmaC(ii)*deltamax(ii)*A;
            if DD(ii)<1 && delta(ii)<=deltamax(ii) && delta(ii)>=0 
                E_rev(n+1)=E_rev(n+1)+0.5*TT(ii)*delta(ii)*A;
            end
            PE(n+1)=PE(n+1)-1/2*A*EE(ii)*(uu(ii+1)-uu(ii))^2/h+1/2*A*EE(ii)*(uu(ii+1)-uuc(ii))^2/h; 
            KE(n+1)=KE(n+1)+1/2*mmc(ii)*vvc(ii)^2;
        elseif cohTracker(ii)==false && AvgForce(ii)>=(SigmaC(ii)*A) && ii>1
            mm(ii)=mm(ii)/2;
            mmc(ii)=mm(ii);
            uuc(ii)=uu(ii);
            vvc(ii)=vv(ii);
            aac(ii)=aa(ii);
            cohTracker(ii)=true;
            delta(ii)=0;
        end
    end
    
    numFrag(n+1)=length(find(DD==1))+1;
    
    if n>nstep_tolerance &&  (numFrag(n+1) == numFrag(n-nstep_tolerance) ) && numFrag(n+1)>1
            numFrag=numFrag(1:n+1);
             break
    end


end


breakpnt=coords(find(DD==1)); %Position of Break Points
if length(breakpnt)>1
    FragLength=[breakpnt(1) breakpnt(2:end)-breakpnt(1:end-1) coords(end)-breakpnt(end)];
elseif length(breakpnt)==1
    FragLength=[breakpnt(1) coords(end)-breakpnt(end)];
else
    FragLength=[];
end
NC(ee)=elem;
NfC(ee)=length(FragLength);
NfC(ee)
end

loglog(NC,NfC,'*')


%% Plotting after multiple runs of the above
errorbar(N,F,neg,pos)
set(gca,'yscale','log')
set(gca,'xscale','log')
errorbar(N,F,neg,pos)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Number of Elements')
ylabel('Number of Fragments')
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',14)