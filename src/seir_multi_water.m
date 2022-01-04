function dxdt = seir_multi_water(t,x,param)
% Michael, Savannah, Sophia, Kaiyue
% SEIR model, multi-water, only considering droplets transmission

%% Force small numbers to 0
    %for i=1:length(x)
    %    if x(i)  < 5/100 %original: 1/(N*100), can be changed into 1/100
    %        x(i) = 0; 
    %    end
    %end

    
%% Get parameter vectors and parameters
    nWater   = param(1);
    toComm   = param(2); %single transit rate from one water site to the community 
    fromComm = param(3); %overall transit rate from the community to transit
    ktw      = param(4); %single
    ktc      = param(5); %single
    betaC = param(6);
    betaW = param(7);
    sigma = param(8); 
    gamma = param(9);

    
%% Water Site Transmission
    dxdt = zeros(length(x), 1);
    SC = x(1);
    EC = x(2); 
    IC = x(3);
    RC = x(4);
    
    dxdt(1) = - betaC*SC*(IC);
    dxdt(2) = + betaC*SC*(IC) - sigma*EC;
    dxdt(3) =                 + sigma*EC - gamma*IC;
    dxdt(4) =                            + gamma*IC;
    
%% Community-Transit-WaterSite
    for i = 1:nWater

        STcw = x(5+(i-1)*12);
        ETcw = x(6+(i-1)*12);
        ITcw = x(7+(i-1)*12);
        RTcw = x(8+(i-1)*12);
        
        SW = x(9 +(i-1)*12);
        EW = x(10+(i-1)*12);
        IW = x(11+(i-1)*12);
        RW = x(12+(i-1)*12);
        
        STwc = x(13+(i-1)*12);
        ETwc = x(14+(i-1)*12);
        ITwc = x(15+(i-1)*12);
        RTwc = x(16+(i-1)*12);
        
        dSTcw = + (fromComm/nWater)*SC - ktw*STcw; 
        dETcw = + (fromComm/nWater)*EC - ktw*ETcw;
        dITcw = + (fromComm/nWater)*IC - ktw*ITcw;
        dRTcw = + (fromComm/nWater)*RC - ktw*RTcw;
        
        dSW = - betaW*SW*(IW) + ktw*STcw - toComm*SW; 
        dEW = + betaW*SW*(IW) + ktw*ETcw - toComm*EW;
        dIW =                 + ktw*ITcw - toComm*IW;
        dRW =                 + ktw*RTcw - toComm*RW;
        
        dSTwc = + toComm*SW - ktc*STwc; 
        dETwc = + toComm*EW - ktc*ETwc;
        dITwc = + toComm*IW - ktc*ITwc;
        dRTwc = + toComm*RW - ktc*RTwc;
        
        dxdt((5+(i-1)*12):(16+(i-1)*12)) = [dSTcw dETcw dITcw dRTcw...
            dSW dEW dIW dRW dSTwc dETwc dITwc dRTwc];
        
        dxdt(1) = dxdt(1) + ktc*STwc - (fromComm/nWater)*SC;
        dxdt(2) = dxdt(2) + ktc*ETwc - (fromComm/nWater)*EC;
        dxdt(3) = dxdt(3) + ktc*ITwc - (fromComm/nWater)*IC;
        dxdt(4) = dxdt(4) + ktc*RTwc - (fromComm/nWater)*RC;
        
    end
end    
