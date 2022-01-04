%% Model 10/29/2020
AC = 2500/0.0036; %[0.48 0.1 0.0036] [Mumbai Nairobi Pretoria]
kC = 5; %[2 5]
%% Assign time parameters

% Time spent in different settings, unit = day(s)
time_T  =  0.5/24;   

% The length of infectious period and exposed period
time_infectious = 9; %9 days of being infectious in total  
time_exposed    = 3; %5-day incubation period, infectious 1-2 days before symptoms develop  


%% Assign number of communities and create blank matrix

nWaterVector = 1:10;
tRange = 1:2000;

WaterTotalMatrix   =   zeros(length(tRange), length(nWaterVector));
CommTotalMatrix    =   zeros(length(tRange), length(nWaterVector));
WaterDensityVector =   zeros(1, length(nWaterVector));
CommDensityVector  =   zeros(1, length(nWaterVector));

TotalEIMatrix =    zeros(length(tRange), length(nWaterVector));
OutbreakTime  =    zeros(1, length(nWaterVector));

TotalCumuInfection    = zeros(length(tRange), length(nWaterVector));
MaxTotalInfectionTime = zeros(1, length(nWaterVector));
MaxTotalInfection     = zeros(1, length(nWaterVector));


%% For loop
for i = 1:length(nWaterVector) 

    nWater = nWaterVector(i); 
    time_W = (60 - 5*(nWater - 1))/(60*24);
    time_C = 1 - time_W - 2*time_T; 
    
%% Transition rate between settings

fromComm = 1/time_C; %total
toComm   = 1/time_W; %seperate 
ktw = nWater*(1/time_T);
ktc = nWater*(1/time_T); %modified on 10/29/2020 

%% Beta term & delta

% Surface area A (unit = m^2)
% AC = 2500/0.48; %[0.48 0.1 0.0036] [Mumbai Nairobi Pretoria]
AW = 25; 

% Probability of infection p
p  = 1/10;
%kC = 5; %[2 5]
kW = 10; 

% Beta 
betaC = (p*kC)/AC;
betaW = (p*kW)/AW;

% Sigma and gamma
sigma = 1/time_exposed; 
gamma = 1/time_infectious; 


%% Set initial conditions

% Different states
x0 = zeros((4 + 12*nWater), 1); 

% Number of people in different states 
N = 2500;       
x0(3) = 1;           %1 infectious people in community 
x0(1) = N - x0(3);   %(2500 - 1) susceptible people in community


%% Set parameters vector  

param = [nWater toComm fromComm ktw ktc betaC betaW sigma gamma]; 


%% Run simulation 
[t,x] = ode45(@seir_multi_water,tRange,x0,[],param);

%Community SEIR
SC = x(:, 1);
EC = x(:, 2); 
IC = x(:, 3);
RC = x(:, 4);

%Total SEIR
Susceptible = zeros(length(tRange), 1); 
Exposed     = zeros(length(tRange), 1); 
Infectious  = zeros(length(tRange), 1); 
Recovered   = zeros(length(tRange), 1); 
for j = 1:nWater
        Susceptible = Susceptible + x(:, 13+(j-1)*12) + x(:, 5+(j-1)*12) + x(:, 9+(j-1)*12);
        Exposed     = Exposed     + x(:, 14+(j-1)*12) + x(:, 6+(j-1)*12) + x(:, 10+(j-1)*12);
        Infectious  = Infectious  + x(:, 15+(j-1)*12) + x(:, 7+(j-1)*12) + x(:, 11+(j-1)*12);
        Recovered   = Recovered   + x(:, 16+(j-1)*12) + x(:, 8+(j-1)*12) + x(:, 12+(j-1)*12);
end
Susceptible = Susceptible + SC; 
Exposed     = Exposed     + EC; 
Infectious  = Infectious  + IC;
Recovered   = Recovered   + RC;

%Water Site SEIR 
SW  = zeros(length(tRange), 1); 
EW  = zeros(length(tRange), 1); 
IW  = zeros(length(tRange), 1); 
RW  = zeros(length(tRange), 1); 
for j = 1:nWater
        SW = SW + x(:, 5+(j-1)*12);
        EW = EW + x(:, 6+(j-1)*12);
        IW = IW + x(:, 7+(j-1)*12);
        RW = RW + x(:, 8+(j-1)*12);
end

%Sum it up
WaterTotalMatrix(:, i) = SW + EW + IW + RW; 
CommTotalMatrix(:, i)  = SC + EC + IC + RC; 

WaterDensityVector(i) = (max(SW + EW + IW + RW) / (AW*nWater)) ; 
CommDensityVector(i)  = (max(SC + EC + IC + RC) / (AC       )) ; 
fprintf(append(('\n'), ('The max average density in water sites is '), string( WaterDensityVector(i) ), (' person(s)/m^2.')));
fprintf(append(('\n'), ('The max average density in community is '),   string( CommDensityVector(i) ), (' person(s)/m^2.')));

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
ColorM = [0   0   0; 
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
      
ColorB = [0   0   0; 
          0 0 128;
          0 0 205;
          0 0 255;
          30 144 255;
          0 191 255;
          135 206 250;
          176 224 230;
          200 225 240;
          240 248 255]; 

%two-axis plot
figure(999)
hold on
for i = 1:length(nWaterVector)
    ColorV = ColorB(i, :)/255; 
yyaxis left
line(nWaterVector, MaxTotalInfection/N, 'color', [0 0.4470 0.7410], 'LineWidth', 1); 
plot(nWaterVector(:, i), MaxTotalInfection(:, i)/N, 'o', 'LineWidth', 5, 'color', ColorV);
axis([1 10 0 1]);
end 
xlabel('Number of Water Place(s)');
ylabel('Maximum Attack Rate');

for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
yyaxis right
axis([1 10 0 1000]);
line(nWaterVector, tRange(OutbreakTime), 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1); 
plot(nWaterVector(:, i), tRange(OutbreakTime(:, i)), 'o', 'LineWidth', 5, 'color', ColorV);

end 
xlabel('Number of Water Place(s)');
ylabel('Outbreak Time (days)');
title('Maximum Attack Rate & Peak Time with Varying Number of Water Place(s)'); 
hold off



%% Other Plots 
%Population density at water site
figure(10)
hold on 
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(nWaterVector(:, i), WaterDensityVector(:, i), 'o', 'LineWidth', 2, 'color', ColorV);
line(nWaterVector, WaterDensityVector, 'color', [231 76 60]/255); 

end 
xlabel('Number of Water Place(s)');
ylabel('Population Density at Each Water Site(s)');
%axis([1 10 0 2.5]); 
title('Population Density at Each Water Site(s) with Different Number of Water Place(s)');
hold off



%People at water site
figure(11)
hold on 
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot( t, (WaterTotalMatrix(:, i)/i ), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('People at Each Water Site(s)');
%axis(); 
title('People at Each Water Site(s) with Different Number of Water Place(s)');
legend('1 Water Place',  '2 Water Places',...
       '3 Water Places', '4 Water Places', '5 Water Places',...
       '6 Water Place',  '7 Water Places', '8 Water Places',...
       '9 Water Places', '10 Water Places'); 
hold off



%Population density at water site
figure(12)
hold on 
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(t, (WaterTotalMatrix(:, i) / (i*AW) ), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('Population density at Water Site(s)');
%axis(); 
title('Population density at Water Site(s) with Different Number of Water Place(s)');
legend('1 Water Place',  '2 Water Places',...
       '3 Water Places', '4 Water Places', '5 Water Places',...
       '6 Water Place',  '7 Water Places', '8 Water Places',...
       '9 Water Places', '10 Water Places'); 
hold off



%Overall E+I
figure(1)
hold on 
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(t, TotalEIMatrix(:, i), 'LineWidth', 1, 'color', ColorV);
%xline(tRange(OutbreakTime(i)), 'LineWidth', 1, 'color', ColorV);    

end  
xlabel('Days');
ylabel('(E+I) Population');
%axis([0 1000 0 600]); 
title('3-Setting Overall (E+I) Population with Different Number of Water Place(s)');
legend('1 Water Place',  '2 Water Places',...
       '3 Water Places', '4 Water Places', '5 Water Places',...
       '6 Water Place',  '7 Water Places', '8 Water Places',...
       '9 Water Places', '10 Water Places'); 
hold off



%Outbreak time v.s Number of Water Place(s)
figure(2)
hold on
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(nWaterVector(:, i), tRange(OutbreakTime(:, i)), 'o', 'LineWidth', 2, 'color', ColorV);
line(nWaterVector, tRange(OutbreakTime), 'color', [231 76 60]/255); 

end 
xlabel('Number of Water Place(s)');
ylabel('Outbreak Time (days)');
%axis([1 10 80 240])
title('Outbreak Time of (I+E) with Varying Number of Water Place(s)'); 
hold off



%Outbreak max scale v.s Number of Water Place(s)
TotalIvector = zeros(1, length(nWaterVector)); 
figure(3)
hold on
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(nWaterVector(:, i), TotalEIMatrix(OutbreakTime(:, i), i), 'o', 'LineWidth', 2, 'color', ColorV);
TotalIvector(:, i) = TotalEIMatrix(OutbreakTime(:, i), i);

end 
line(nWaterVector(1, :), TotalIvector(1, :), 'color', [231 76 60]/255); 
xlabel('Number of Water Place(s)');
ylabel('Peak level of outbreak');
%axis([1 10 0 600])
title('Maximum of (I+E) with Varying Number of Water Place(s)'); 
hold off
      


%Overall transmission--cumulative infection
figure(4)
hold on 
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(t, TotalCumuInfection(:, i) / N, 'LineWidth', 1, 'color', ColorV);
    
end 
xlabel('Days');
ylabel('Attack Rate%');
%axis([0 1000 0 0.9])
title('3-Setting Attack Rate with Different Number of Water Place(s)'); 
legend('1 Water Place',  '2 Water Places',...
       '3 Water Places', '4 Water Places', '5 Water Places',...
       '6 Water Place',  '7 Water Places', '8 Water Places',...
       '9 Water Places', '10 Water Places'); 
hold off



%Maximum of Overall Infection Time 
figure(5)
hold on
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(nWaterVector(:, i), tRange(MaxTotalInfectionTime(:, i)), 'o', 'LineWidth', 2, 'color', ColorV);
line(nWaterVector, tRange(MaxTotalInfectionTime), 'color', [231 76 60]/255); 

end 
xlabel('Number of Water Place(s)');
ylabel('Time for Maximum of Attack Rate (days)');
%axis([1 10 150 500])
title('Time for Maximum of Attack Rate with Different Number of Water Place(s)');
hold off



%Maximum level of Overall Infection
figure(6)
hold on
for i = 1:length(nWaterVector)
    ColorV = ColorM(i, :)/255; 
plot(nWaterVector(:, i), MaxTotalInfection(:, i)/N, 'o', 'LineWidth', 2, 'color', ColorV);
line(nWaterVector, MaxTotalInfection/N, 'color', [231 76 60]/255); 

end 
xlabel('Number of Water Place(s)');
ylabel('Maximum Attack Rate');
%axis([1 10 0.45 0.9])
title('Maximum Attack Rate by Number of Water Place(s)'); 
hold off

