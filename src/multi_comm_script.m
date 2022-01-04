%% Model 10/29/2020
AC = 250/0.0036; %[0.48 0.1 0.0036], for each single community  
kC = 5; %[2 5]
%% Assign time parameters

% Time spent in different settings, unit = day(s)
time_T  =  0.5/24;   

% The length of infectious period and exposed period
time_infectious = 9; %9 days of being infectious in total  
time_exposed    = 3; %5-day incubation period, infectious 1-2 days before symptoms develop  


%% Assign number of communities and create blank matrix

nCommVector = 1:10;
tRange = 1:2000;

WaterTotalMatrix   =   zeros(length(tRange), length(nCommVector));
CommTotalMatrix    =   zeros(length(tRange), length(nCommVector));
WaterDensityVector =   zeros(1, length(nCommVector));
CommDensityVector  =   zeros(1, length(nCommVector));

TotalEIMatrix =    zeros(length(tRange), length(nCommVector));
OutbreakTime  =    zeros(1, length(nCommVector));

TotalCumuInfection    = zeros(length(tRange), length(nCommVector));
MaxTotalInfectionTime = zeros(1, length(nCommVector));
MaxTotalInfection     = zeros(1, length(nCommVector));


%% For loop
for i = 1:length(nCommVector) 

    nComm  = nCommVector(i); 
    time_W = (20 + 10*nComm)/(60*24);
    time_C = 1 - time_W - 2*time_T; 
    
%% Transition rate between settings

toWater   = 1/time_C; %separate 
fromWater = 1/time_W; %total 
ktw = 1/time_T;
ktc = 1/time_T; 

%% Beta term & delta

% Surface area A (unit = m^2)
% AC = 250/0.0036; %[0.48 0.1 0.0036], for each single community  
AW = 25; 

% Probability of infection p
p  = 1/10;
%kC = 5; 
kW = 10; 

% Beta 
betaC = (p*kC)/AC;
betaW = (p*kW)/AW;

% Sigma and gamma
sigma = 1/time_exposed; 
gamma = 1/time_infectious; 


%% Set initial conditions

% Different states
x0 = zeros((12*nComm + 4), 1); 

% Number of people in different states 
N = 250*nComm;       %250 people for each single community 
x0(3) = 1;           %1 infectious people in community 1 
x0(1) = 250 - x0(3); %249 susceptible people in community 1
if nComm > 1
    for k = 2:nComm
        x0(12*(k - 1) + 1) = 250; 
    end
end 

%% Set parameters vector  

param = [nComm toWater fromWater ktw ktc betaC betaW sigma gamma]; 


%% Run simulation 
[t,x] = ode45(@seir_multi_community,tRange,x0,[],param);

%Water site SEIR
endpoint =  12*nComm + 4; 
SW = x(:, endpoint-3);
EW = x(:, endpoint-2); 
IW = x(:, endpoint-1);
RW = x(:, endpoint  );

%Total SEIR
Susceptible = zeros(length(tRange), 1); 
Exposed     = zeros(length(tRange), 1); 
Infectious  = zeros(length(tRange), 1); 
Recovered   = zeros(length(tRange), 1); 
for j = 1:nComm
        Susceptible = Susceptible + x(:, 1+(j-1)*12) + x(:, 5+(j-1)*12) + x(:, 9+(j-1)*12);
        Exposed     = Exposed     + x(:, 2+(j-1)*12) + x(:, 6+(j-1)*12) + x(:, 10+(j-1)*12);
        Infectious  = Infectious  + x(:, 3+(j-1)*12) + x(:, 7+(j-1)*12) + x(:, 11+(j-1)*12);
        Recovered   = Recovered   + x(:, 4+(j-1)*12) + x(:, 8+(j-1)*12) + x(:, 12+(j-1)*12);
end
Susceptible = Susceptible + SW; 
Exposed     = Exposed     + EW; 
Infectious  = Infectious  + IW;
Recovered   = Recovered   + RW;

%Community SEIR 
SC  = zeros(length(tRange), 1); 
EC  = zeros(length(tRange), 1); 
IC  = zeros(length(tRange), 1); 
RC  = zeros(length(tRange), 1); 
for j = 1:nComm
        SC = SC + x(:, 1+(j-1)*12);
        EC = EC + x(:, 2+(j-1)*12);
        IC = IC + x(:, 3+(j-1)*12);
        RC = RC + x(:, 4+(j-1)*12);
end

%Sum it up
WaterTotalMatrix(:, i) = SW + EW + IW + RW; 
CommTotalMatrix(:, i)  = SC + EC + IC + RC; 

WaterDensityVector(i) = (max(SW + EW + IW + RW) / (AW)      ) ; 
CommDensityVector(i)  = (max(SC + EC + IC + RC) / (AC*nComm)) ; 
fprintf(append(('\n'), ('The max average density in water site is '), string( WaterDensityVector(i) ), (' person(s)/m^2.')));
fprintf(append(('\n'), ('The max average density in communities is '), string( CommDensityVector(i) ), (' person(s)/m^2.')));

% Find the time of the outbreak and maximum of...
OutbreakTime(1, i) = find((Exposed + Infectious) == max(Exposed + Infectious), 1); %I + E

TotalInfection = Exposed + Infectious + Recovered; 
if max(TotalInfection) <= 1 
    MaxTotalInfectionTime(i) = 0;
else 
    MaxTotalInfectionTime(i) = find((TotalInfection) >= floor(max(TotalInfection) - 1), 1);
end
MaxTotalInfection(i) = max(TotalInfection);

% For loops
TotalCumuInfection(:, i) = (Exposed + Infectious + Recovered); 
TotalEIMatrix(:, i) = (Exposed + Infectious); 

% Report progress
fprintf(append(('\n'), (string(i)), ('/10 completed.'))); 

end 



%% Plots 
ColorM1 = [0   0   0; 
          100 30  22;
          148 49  38;
          176 58  46;
          203 67  53;
          231 76  60;
          236 112 99;
          241 148 138;
          245 183 177;
          %245 200 200; %newly added  
          255 233 233]; 
      
ColorB1 = [0   0   0; 
          0 0 128;
          0 0 205;
          0 0 255;
          30 144 255;
          0 191 255;
          135 206 250;
          176 224 230;
          200 225 240;
          240 248 255]; 
      
% Reverse the color matrix  
ColorM = zeros(10, 3); 
for l = 1:10
    ColorM(l, :) = ColorM1(11 - l, :); 
end 
ColorB = zeros(10, 3); 
for k = 1:10
    ColorB(k, :) = ColorB1(11 - k, :); 
end 

%two-axis plot
figure(999)
hold on
for i = 1:length(nCommVector)
    ColorV = ColorB(i, :)/255; 
yyaxis left
line(nCommVector, MaxTotalInfection/N, 'color', [0 0.4470 0.7410], 'LineWidth', 1); 
plot(nCommVector(:, i), MaxTotalInfection(:, i)/N, 'o', 'LineWidth', 5, 'color', ColorV);
axis([1 10 0 1])
end 
xlabel('Number of Communities');
ylabel('Maximum Attack Rate');

for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
yyaxis right
axis([1 10 0 1000])
line(nCommVector, tRange(OutbreakTime), 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1); 
plot(nCommVector(:, i), tRange(OutbreakTime(:, i)), 'o', 'LineWidth', 5, 'color', ColorV);

end 
xlabel('Number of Communities');
ylabel('Outbreak Time (days)');
title('Maximum Attack Rate & Peak Time with Varying Number of Communities'); 
hold off





%% Other Plots
%Water sites population density
figure(10)
hold on 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(nCommVector(:, i), WaterDensityVector(:, i), 'o', 'LineWidth', 2, 'color', ColorV);
line(nCommVector, WaterDensityVector, 'color', [231 76 60]/255); 

end 
xlabel('Number of Communities');
ylabel('Population Density at Water Site (person/m^2)');
%axis([1 10 0 2.5]); 
title('Population Density at Water Site with Different Number of Communities');
hold off



%People at water site
figure(11)
hold on 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot( t, (WaterTotalMatrix(:, i)/i ), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('People at Water Site');
%axis(); 
title('People at Water Site with Different Number of Communities');
legend('1 Community',  '2 Communities',...
       '3 Communities', '4 Communities', '5 Communities',...
       '6 Communities', '7 Communities', '8 Communities',...
       '9 Communities', '10 Communities'); 
hold off



%Population density at water site
figure(12)
hold on 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(t, (WaterTotalMatrix(:, i) / (i*AW) ), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('Population density at Water Site');
%axis(); 
title('Population density at Water Site with Different Number of Communities');
legend('1 Communities',  '2 Communities',...
       '3 Communities',  '4 Communities', '5 Communities',...
       '6 Communities',  '7 Communities', '8 Communities',...
       '9 Communities',  '10 Communities'); 
hold off



%Overall E+I
figure(1)
hold on 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(t, TotalEIMatrix(:, i)/(i*250), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('(E+I)% in Population');
axis([0 1000 0 0.6]); 
title('3-Setting Overall (E+I)% in Population with Different Number of Communities');
legend('1 Water Communities',  '2 Water Communities',...
       '3 Water Communities',  '4 Water Communities', '5 Water Communities',...
       '6 Water Communities',  '7 Water Communities', '8 Water Communities',...
       '9 Water Communities',  '10 Water Communities'); 
hold off



%Outbreak time v.s Number of Water Place(s)
figure(2)
hold on
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(nCommVector(:, i), tRange(OutbreakTime(:, i)), 'o', 'LineWidth', 2, 'color', ColorV);
line(nCommVector, tRange(OutbreakTime), 'color', [231 76 60]/255); 

end 
xlabel('Number of Communities');
ylabel('Outbreak Time (days)');
%axis([1 10 80 240])
title('Outbreak Time of (I+E) with Varying Number of Communities'); 
hold off



%Outbreak max scale v.s Number of Water Place(s)
TotalIvector = zeros(1, length(nCommVector)); 
figure(3)
hold on
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(nCommVector(:, i), TotalEIMatrix(OutbreakTime(:, i), i)/(i*250), 'o', 'LineWidth', 2, 'color', ColorV);
TotalIvector(:, i) = TotalEIMatrix(OutbreakTime(:, i), i)/(i*250);

end 
line(nCommVector(1, :), TotalIvector(1, :), 'color', [231 76 60]/255); 
xlabel('Number of Communities');
ylabel('Peak level of outbreak');
%axis([1 10 0 600])
title('Maximum of (I+E)% in Population with Varying Number of Communities'); 
hold off
      


%Overall transmission--cumulative infection
figure(4)
hold on 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(t, TotalCumuInfection(:, i)/double(250*i), 'LineWidth', 1, 'color', ColorV);
    
end 
xlabel('Days');
ylabel('Attack Rate%');
%axis([0 1000 0 0.9])
title('3-Setting Attack Rate with Different Number of Communities (%)'); 
legend('1 Community',  '2 Communities',...
       '3 Communities', '4 Communities', '5 Communities',...
       '6 Communities',  '7 Communities', '8 Communities',...
       '9 Communities', '10 Communities'); 
hold off



%Maximum of Overall Infection Time 
figure(5)
hold on
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(nCommVector(:, i), tRange(MaxTotalInfectionTime(:, i)), 'o', 'LineWidth', 2, 'color', ColorV);
line(nCommVector, tRange(MaxTotalInfectionTime), 'color', [231 76 60]/255); 

end 
xlabel('Number of Communities');
ylabel('Time for Maximum of Attack Rate (days)');
%axis([1 10 150 500])
title('Time for Maximum of Attack Rate with Different Number of Communities');
hold off



%Maximum level of Overall Infection
figure(6)
hold on
thing = zeros(1, length(nCommVector)); 
for i = 1:length(nCommVector)
    ColorV = ColorM(i, :)/255; 
plot(nCommVector(:, i), MaxTotalInfection(:, i)/double(250*i), 'o', 'LineWidth', 2, 'color', ColorV);
thing(:, i) = MaxTotalInfection(:, i)/double(250*i); 

end 
line(nCommVector, thing, 'color', [231 76 60]/255); 
xlabel('Number of Communities');
ylabel('Maximum Attack Rate');
%axis([1 10 0.45 0.9])
title('Maximum Attack Rate by Number of Communities'); 
hold off


