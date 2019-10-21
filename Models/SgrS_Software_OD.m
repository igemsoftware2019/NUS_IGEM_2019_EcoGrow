
function OD = SgrS_Software_OD(Duration, aTc_concentration, kG, ka, MewMax)
    close all;
    %---------------------edit window-----------------------
  
    %define the time step and time span for printing data points
    dt = 0.0001;                %time step 0.0001
    tmin = 0;
    tmax = Duration;

    %define parameters
    degaTc = 3.57876065e-6*60; 
%     MewMax = 0.017413640536691700*60;
%     kG = 0.31736719; 
%     ka = 700.3832460001160;
    Yield = 4.4822656250000000;
    p = 0.5154904145493300; 
    
    %variables to plot       
    variables(1).name = "aTc concentration";
    variables(1).init = 0;
    variables(1).unit = " (nM)";
    variables(1).scal = 1;
    
    variables(2).name = "OD_6_0_0";
    variables(2).init = 0.14;
    variables(2).unit = " ";
    variables(2).scal = 1;;

    variables(3).name = "Glucose concentration";
    variables(3).init = 0.2;             %g/L
    variables(3).unit = " (%)";
    variables(3).scal = 1;
    
    %---------------------edit window-----------------------
    
    
    %create parameters array
    params = [degaTc MewMax kG ka Yield p];

    
    %call the numerical integrator/solver
    results = solve(params,dt,tmin,tmax,graph,variables,aTc_concentration);
       
    %output
    for j = 1:(Duration*100)
        OD(j,1) = (j-1)/100;
    end
    for i = 1:length(aTc_concentration)
        for j = 1:(Duration*100)
            OD(j,i+1) = results(1+(j-1)*100,3,i);
        end
    end
    
%     writematrix(ODmean,'ODmean.csv');

end


function results = solve(params,dt,tmin,tmax,graph,variables,aTc_concentration)
    tspan = tmin:dt:tmax; %time span
    numtspan = length(tspan); %the number of time points in tspan
    numVariables = length(variables); %the number of variables
    numConc = length(aTc_concentration);

    
    for l = 1:numConc
        variables(1).init = aTc_concentration(l);
        for i = 1:numVariables
            Init(i) = variables(i).init; %create an initial values array to pass into ODE solver
        end
        [t,y] = euler_modified(@(t,y) SolveODE(t,y,params), tmin, Init, tmax, numtspan);
        results(1:numtspan,1,l) = t;
        for k = 1:numVariables
            results(1:numtspan,k+1,l) = y(:,k);
        end
    end
end


function dydt = SolveODE(t,y,params)
    %define the variables
    Inducer = y(1);
    OD = y(2);
    Glucose = y(3);

    degaTc = params(1);
    MewMax = params(2);
    kG = params(3); 
    ka = params(4);
    Yield = params(5);
    p = params(6);
    
    Mew = MewMax * (  Glucose / (Glucose + kG)  ) * (  ((ka).^p) / (((Inducer).^p) + ((ka).^p))  );
    %ODE expressions
    dInducer = -(degaTc + Mew) * Inducer;
    dOD =  OD * Mew;
    dGlucose   = -dOD/ Yield;
    
    %make sure the sequence is consistent with the variables sequence.
    dydt = [dInducer dOD dGlucose];     
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



