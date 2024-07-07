a = [105, 44.10, 40.60];
b = [2.45, 3.51, 3.89];
c = [0.005, 0.005, 0.005];
D = 250;

% optimization parameters
options = optimoptions(@fmincon, ...
    'MaxIterations', 100, ...
    'Display', 'iter', ...
    'TolFun', 1e-6);

P_min = [10, 20, 20];
P_max = [160, 80, 50];

% Define the objective function for optimization
objective_function = @(P) ELD_Objective(P, a, b, c);

% Run optimization
[x,fval] = fmincon(objective_function,[0,0,0],[],[],[1,1,1],D,P_min,P_max,[],options);

% Display results
fprintf('The optimal generation is: \nP1 = %.2f MW \nP2 = %.2f MW \nP3 = %.2f MW \n', x)
fprintf('Total generated power = %.2f MW\n',sum(x))
fprintf('The total cost is %s $/h\n', num2str(fval))


function total_cost = ELD_Objective(P, a, b, c)
    total_cost = 0;
    n = 3;
    for i = 1:n
        total_cost = total_cost + (c(i) * P(i)^2 + b(i) * P(i) + a(i));
    end
end