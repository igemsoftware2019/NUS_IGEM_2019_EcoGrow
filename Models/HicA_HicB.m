function [time,output] = HicA_HicB(Inducer_t,Inducer_at,time_t,time_at,selection, duration)
    %to clear workspace and close all figures
%     clear all; 
%     close all;
    %---------------------edit window-----------------------
     %define the time step and time span for printing data points
    dt = 0.001; %time step
    tmin = 0;
    tmax = duration;
    
   if selection == 1
       InputParam.name = 'Inducer_t';
       InputParam.data = Inducer_t;
       InputParam.scal = 1e-6;
       InputParam.unit = 'M';
       InputParam.length = length(Inducer_t);
   elseif selection ==2 
       InputParam.name = 'Inducer_at';
       InputParam.data = Inducer_at;
       InputParam.scal = 1;
       InputParam.unit = 'M';
       InputParam.length = length(Inducer_at);
   elseif selection ==3
       InputParam.name = 'time_t';
       InputParam.data = time_t;
       InputParam.scal = 1;
       InputParam.unit = 'hr';
       InputParam.length = length(time_t);
   elseif selection ==4
       InputParam.name = 'time_at';
       InputParam.data = time_at;
       InputParam.scal = 1;
       InputParam.unit = 'hr';
       InputParam.length = length(time_at);
   end
           
Mumax = 0.24186;
ODmax = 1.0407;
Ktox   = 25.795;
ntox   = 3.21;
Kanti = 0.00026938;
syntoxin = 86.536;
Ki   = 98.919;
synantitoxin = 0.00059868;
Ka    = 0.018984; 
n_l    = 0.5076;
ODmin = 0.77606;
degt = 0.62872;
degat = 0.78863;

params = [Mumax ODmax Ktox ntox Kanti syntoxin Ki synantitoxin Ka n_l ODmin degt degat]; 

%     %What parameters are you varying? 
%     InputParam(1).name = 'Inducer_t'; %individual graph lines         
%     InputParam(1).length = 4; %number of values for the parameter
%     InputParam(1).unit = 'M';%unit of input parameter for legend
%     InputParam(1).scal = 1; %scaling factor for legend
%     InputParam(1).index = 16; %position of that parameter in the params array 
%     
%     %COMMENT OUT InputParam(2) CODE BELOW IF ONLY USING 1 INPUT PARAM
%     %second input param
%     InputParam(2).name = 'Inducer_at';%sets of graph lines, can leave it as any parameter if u r not plotting sets        
%     InputParam(2).length = 1;%number of values for the parameter
%     InputParam(2).unit = 'muM';
%     InputParam(2).scal = 5;  
%     InputParam(2).index = 17;
    
    %variables to plot 
    variables(1).name = "IPTG";
%     variables(2).init = 1e-3;    %M (1mM)
    variables(1).unit = " (\muM)";
    variables(1).scal = 1e6;  
    
    variables(2).name = "Toxin";
    variables(2).init = 0;    %M
    variables(2).unit = " (\muM)";
    variables(2).scal = 1e6;  
    
    variables(3).name = "Arabinose";
%     variables(4).init = 0.0133;    %M 2% Ara
    variables(3).unit = " (\muM)";
    variables(3).scal = 1e6;  
    
    variables(4).name = "Antitoxin";
    variables(4).init = 0;    %M
    variables(4).unit = " (\muM)";
    variables(4).scal = 1e6;  
    
    variables(5).name = "OD600";
    variables(5).init = 0.2;          %g/L, dry mass of E.coli at 0.03 OD
    variables(5).unit = " ";
    variables(5).scal = 1;%1/0.39
    
    %plot graph? (all:0, 1st variable: 1, 2nd variable: 2..., none: -1)
    graph = 0;
    %---------------------edit window-----------------------
%     if length(InputParam)==1
%         InputParam(2).length = 1; 
%     end
  
%     %Extracting all parameter values from an Excel File
%     opts = detectImportOptions('params_hicA_hicB.xlsx');
%     Allparams = readmatrix('params_hicA_hicB.xlsx',opts);
%     
%     %Extracting Parameters 1 and 2
%      opts.SelectedVariableNames = InputParam(1).name;
%      Parameter1 = readmatrix('params_hicA_hicB.xlsx',opts);
%      InputParam(1).data = Parameter1;
%      opts.SelectedVariableNames = InputParam(2).name;
%      Parameter2 = readmatrix('params_hicA_hicB.xlsx',opts);
%      InputParam(2).data = Parameter2;
    
%       %Defining the params array 
%       for i = 1:InputParam(2).length   
%        for j = 1:InputParam(1).length
%         params(j,:,i)= Allparams(j,:);
%         params(:,InputParam(2).index,i)= Parameter2(i);
%        end   
%       end

%Assigning the parameters to the initial value of the variables 
%  variables(2).init = params(1,16);%Inducer for toxin
%  variables(4).init = params(1,17); %Inducer for antitoxin
      
      
%call the numerical integrator/solver
results = solve(params,dt,tmin,tmax,graph,variables,InputParam,selection,Inducer_t,Inducer_at,time_t,time_at);

% %Reading experimental data 
% data = load('hicA_B_plotting_16.mat','ODmean');
% expdata1 = data.ODmean(:,:,1);
% expdata2 = data.ODmean(:,:,2);
% expdata3 = data.ODmean(:,:,3);
% 
% writematrix(expdata3,'Brod_prediction_experiment.xlsx');
%   
% %plot graphs
time = results(:,1,1);

growthdata = results(:,6,:);
[a b c] = size(growthdata);

for i = 1:c
    output(:,i) = results(:,6,i);
end
% hold on 
% plot(expdata3(:,1),expdata3(:,2),'Color','k');
% plot(expdata3(:,1),expdata3(:,5),'Color','r');
% plot(expdata3(:,1),expdata3(:,6),'Color','b');
% plot(expdata3(:,1),expdata3(:,7),'Color','m');%plotting experimental data 
end
function results = solve(params,dt,tmin,tmax,graph,variables,InputParam,selection,Inducer_t,Inducer_at,time_t,time_at)
    tspan = tmin:dt:tmax; %time span
    numtspan = length(tspan); %the number of time points in tspan
    numVariables = length(variables); %the number of variables
    
    %initial conditions for the variables (must follow the same sequences)
    if selection == 1
        for i = 1:InputParam.length
        Init(i,:)  = [InputParam.data(i) 0 Inducer_at(1) 0 0.1317];  
        end
    elseif selection ==2 
        for i = 1:InputParam.length
        Init(i,:)  = [Inducer_t(1) 0 Inducer_at(1) 0 0.1317];  
        end
    elseif selection ==3
        Init(1,:)  = [Inducer_t(1) 0 Inducer_at(1) 0 0.1317];  
    elseif selection ==4
        Init(1,:)  = [Inducer_t(1) 0 Inducer_at(1) 0 0.1317];    
    end
    
    
%     for i = 1:InputParam(1).length
%     Init(i,:)  = [params(i,16) 0 params(1,17) 0 0.1317];  %the initial values of IPTGout,IPTGin,hicAmRNA,hicA, OD, GFP
%     end
%     
%     for i = 1:numVariables
%         for j = 1:InputParam(1).length
%         Init(:,i) = variables(i).init; %create an initial values array to pass into ODE solver
%         Init(j,2) = params(j,16);
%         Init(j,4) = params(j,17);
%         end
%     end
    
%     for i = 1:InputParam(2).length
        for j = 1:InputParam.length
%             InputParam(1).length 
            if selection ==1
            [t,y] = euler_modified(@(t,y) SolveODE(t,y,params,Inducer_t(j),Inducer_at,time_t,time_at), tmin, Init(j,:), tmax, numtspan);
            elseif selection ==2
            [t,y] = euler_modified(@(t,y) SolveODE(t,y,params,Inducer_t,Inducer_at(j),time_t,time_at), tmin, Init(j,:), tmax, numtspan);
            elseif selection ==3
            [t,y] = euler_modified(@(t,y) SolveODE(t,y,params,Inducer_t,Inducer_at,time_t(j),time_at), tmin, Init(1,:), tmax, numtspan);
            elseif selection ==4
            [t,y] = euler_modified(@(t,y) SolveODE(t,y,params,Inducer_t,Inducer_at,time_t,time_at(j)), tmin, Init(1,:), tmax, numtspan);  
            end
            results(1:numtspan,1,j) = t;
            for k = 1:numVariables
                results(1:numtspan,k+1,j) = y(:,k);
            end
        end
%     end
end
function dydt = SolveODE(t,y,params,Inducer_t,Inducer_at,time_t,time_at)
%define the variables
IPTG = y(1);
Toxin = y(2);
Arabinose = y(3);
Antitoxin = y(4);
OD = y(5);

Mumax = params(1);
ODmax = params(2);
Ktox   = params(3);
ntox   = params(4);
Kanti = params(5);
syntoxin = params(6);
Ki   = params(7);
synantitoxin = params(8);
Ka    = params(9); 
n_l    = params(10);
ODmin = params(11);
degt = params(12);
degat = params(13);


 if t<time_t %Toxin induced at 1h
    IPTG = 0;
 end

 if t<time_at %Antitoxin induced at 3.5h
     Arabinose = 0;
 end

dIPTG = 0;  
dToxin = syntoxin* IPTG/(IPTG+Ki) - degt*Toxin;
dArabinose = 0;
dAntitoxin = synantitoxin* Arabinose/(Arabinose+Ka) - degat*Antitoxin; 

toxeff = Toxin*(Kanti)/(Kanti + Antitoxin);
mu     =  Mumax*(1-(OD/ODmax)^n_l)*Ktox^ntox /(Ktox^ntox + toxeff^ntox);
dOD    =  (OD + ODmin)*mu;


dydt = [dIPTG dToxin dArabinose dAntitoxin dOD];    
end   
function [t,y]=euler_modified(f,tinit,yinit,tfinal,n)
    % Euler approximation for ODE initial value problem
    % Euler modified method
    h=(tfinal-tinit)/n;
    % Initialization of x and y as column vectors
    t = zeros(n, 1);
    t(1,1) = tinit; 
    y = zeros(n, length(yinit));
    y(1,:) = yinit;
    % Calculation of x and y
    for i=1:n-1
        t(i+1,1)=t(i,1)+h;
        delta_y=f(t(i,1),y(i,:));
        ynew=y(i,:)+h*delta_y;
        y(i+1,:)=y(i,:)+(h/2)*(f(t(i,1),y(i,:))+f(t(i+1,1),ynew));
        %[ynew, y(i,:)+(h/2)*(f(t(i,1),y(i,:))+f(t(i+1,1),ynew))]
    end

end   

    