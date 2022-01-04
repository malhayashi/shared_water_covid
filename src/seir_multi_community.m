%% Multi-community Function 10/5/2020

function dxdt = seir_multi_community(t,x,param)
%% Set Parameter Values
    nComm = param(1);
    toWater   = param(2); %single transit rate from one community to transit
    fromWater = param(3); %overall transit rate from water site to transit
    ktw       = param(4); %single
    ktc       = param(5); %single
    betaC = param(6);
    betaW = param(7);
    sigma = param(8); 
    gamma = param(9);
    
%% Water Site Transmission
    endpoint =  12*nComm + 4; 

    dxdt = zeros(length(x), 1);
    SW = x(endpoint-3);
    EW = x(endpoint-2); 
    IW = x(endpoint-1);
    RW = x(endpoint-0);
    
    dxdt(endpoint-3) = - betaW*SW*(IW);
    dxdt(endpoint-2) = + betaW*SW*(IW);
    dxdt(endpoint-1) = 0              ;
    dxdt(endpoint-0) = 0              ;
    
%% Community-Transit-WaterSite
    for i = 1:nComm
        SC = x(1+(i-1)*12);
        EC = x(2+(i-1)*12);
        IC = x(3+(i-1)*12);
        RC = x(4+(i-1)*12);
        
        STcw = x(5+(i-1)*12);
        ETcw = x(6+(i-1)*12);
        ITcw = x(7+(i-1)*12);
        RTcw = x(8+(i-1)*12);
        
        STwc = x(9 +(i-1)*12);
        ETwc = x(10+(i-1)*12);
        ITwc = x(11+(i-1)*12);
        RTwc = x(12+(i-1)*12);
        
        dSC = - betaC*SC*(IC)                       - toWater*SC + ktc*STwc; 
        dEC = + betaC*SC*(IC) - sigma*EC            - toWater*EC + ktc*ETwc;
        dIC =                 + sigma*EC - gamma*IC - toWater*IC + ktc*ITwc;
        dRC =                            + gamma*IC - toWater*RC + ktc*RTwc;
        
        dSTcw = + toWater*SC - ktw*STcw; 
        dETcw = + toWater*EC - ktw*ETcw;
        dITcw = + toWater*IC - ktw*ITcw;
        dRTcw = + toWater*RC - ktw*RTcw;
        
        dSTwc = + (fromWater/nComm)*SW - ktc*STwc; 
        dETwc = + (fromWater/nComm)*EW - ktc*ETwc;
        dITwc = + (fromWater/nComm)*IW - ktc*ITwc;
        dRTwc = + (fromWater/nComm)*RW - ktc*RTwc;

        dxdt((1+(i-1)*12):(12+(i-1)*12)) = [dSC dEC dIC dRC...
            dSTcw dETcw dITcw dRTcw dSTwc dETwc dITwc dRTwc];
        
        dxdt(endpoint-3) = dxdt(endpoint-3) + ktw*STcw - (fromWater/nComm)*SW ;
        dxdt(endpoint-2) = dxdt(endpoint-2) + ktw*ETcw - (fromWater/nComm)*EW ;
        dxdt(endpoint-1) = dxdt(endpoint-1) + ktw*ITcw - (fromWater/nComm)*IW ;
        dxdt(endpoint-0) = dxdt(endpoint-0) + ktw*RTcw - (fromWater/nComm)*RW ;
        
    end

end