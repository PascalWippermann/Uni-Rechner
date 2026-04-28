clc
close all
clearvars

%% a)
x = logspace(0,1,100);
func = zeros(1,100);


ex = exp(x);
figure
loglog(x,func,x,ex);
legend;
%% b)
r = 100;
nx = 10;
while r >= 1
    func = zeros(1,100);
    for n=0:nx
        func = func +x.^n ./ factorial(n);
    end
    r = rmse(func,ex)
    if r < 1
        break
    else
        nx = nx+1;
    end
end

