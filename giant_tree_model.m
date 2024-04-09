% Simple approaches to calculating the population of giant trees based on observed populations and demographic rates
%
% T. Pugh
% 09.04.24

%% Simplest model
% Assume constant mortality rates with tree size

Ngiant_init=20; %Stems per ha
Mgiant=0.01; %Fraction of stems per year
Rgiant=1; %New stems per ha per year

%--- Numerical solution ---

nyear=1000; %Number of years to simulate

Ngiant=NaN(nyear,1);
Ngiant(1)=Ngiant_init;
for nn=2:nyear
    Ngiant(nn) = Ngiant(nn-1) + Rgiant - (Ngiant(nn-1)*Mgiant);
end
clear nn

figure
plot(Ngiant)

% If we run for long enough to have a clear asymptote then can find the time at which the 95th percentile of tree number
% is reached
Ngiant_95th=Ngiant(nyear)*0.95;
diff_95th=Ngiant-Ngiant_95th;
time_to_95th=find(abs(diff_95th)==min(abs(diff_95th))); %Find the index (equivalent to year) for which the Ngiant value corresponds most closely to the 95th percentile of Ngiant


%--- Analytical solution ---

% dNgiant/dt = Rgiant - (Ngiant*Mgiant) = 0 (Eq.1)
% Ngiant = Rgiant/Mgiant (Eq. 2)
% Integrate Eq. 1 and solve for time based on e.g. 95% of Ngiant value from Eq. 2


%% Consider changing mortality rate with size
% Assume basic area of 1 ha.

ntree=20;
Drange=0.01:0.001:5; %Diameter range considered (m)
Dgiant_init=rand(1,ntree); %Some random initialisation of D per tree - this should be read in from the data
Ggiant=1./(1.01-exp(-(Drange)))*0.01 + 0.001; %A crappy function to change relative growth rates according to size (D) - should be replaced by one based on data
Mgiant=1./(1.01-exp(-(Drange)))*0.002 + 0.0005; %A crappy function to change mortality rates according to size (D) - should be replaced by one based on data
Rgiant=1.1; %New stems per ha per year
Rthres=0.1; %Recruitment diameter threshold for new trees (m)

nyear=1000; %Number of years to simulate

Dgiant=NaN(nyear,ntree);
Dgiant(1,:)=Dgiant_init;
for nn=2:nyear
    for tt=1:ntree
        if isnan(Dgiant(nn-1,tt))
            continue %Tree was not alive at previous timepoint so skip calculations for this tree
        end
        %Find size dependent rates
        size_diff=Drange-Dgiant(nn-1,tt);
        ind=find(abs(size_diff)==min(abs(size_diff)));
        Grate_tree=Ggiant(ind);
        Mrate_tree=Mgiant(ind);
        clear size_diff ind
        
        %First growth
        Dgiant(nn,tt) = Dgiant(nn-1,tt) + Dgiant(nn-1,tt)*Grate_tree;

        %Then death (stochastically defined, based on a uniform distribution between 0 and 1)
        %n.b. to make stochastic parts repeatable, then define a fixed initial seed for the simulations
        if rand(1)<Mrate_tree
            Dgiant(nn,tt)=NaN; % If tree is dead, set diameter to NaN
        end
    end
    
    %Finally recruitment
    newtrees=0;
    newtrees=newtrees+floor(Rgiant); %Increment tree number by the number of whole stems recruited each year
    if rand(1)<=rem(Rgiant,floor(Rgiant)) %Deal with the probability specified by the decimal part of Rgiant
        newtrees=newtrees+1;
    end
    Dgiant=cat(2,Dgiant,NaN(nyear,newtrees)); %Append new tree rows to array
    Dgiant(nn,ntree+1:ntree+newtrees)=Rthres;
    ntree=ntree+newtrees;

    clear Grate_tree Mrate_tree newtrees

    if mod(nn,100)==0
        fprintf('Completed year %d\n',nn)
    end
end
clear nn tt

figure
subplot(1,2,1)
histogram(Dgiant(1,:))
xlabel('D (m)')
ylabel('N trees')
title('Year 1')
subplot(1,2,2)
histogram(Dgiant(nyear,:))
xlabel('D (m)')
ylabel('N trees')
title('Final year')

% Number of trees over 80 cm each year
ntree_over_80cm=NaN(nyear,1);
for nn=1:nyear
    bigtree_ind=find(~isnan(Dgiant(nn,:))); % Find indexes of all trees that are alive (i.e. not NaN)
    ntree_over_80cm(nn)=length(bigtree_ind);
end

figure
plot(ntree_over_80cm)
xlabel('Year')
ylabel('N trees')

