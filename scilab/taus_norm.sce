
//clc();
clear;
stacksize('max');
path = "/media/data/evo/manipulator_dynamics/data_for_identification_Youbot/";
N = 51;
p = 5;

//path = "/media/data/evo/manipulator_dynamics/data_for_identification_2dof/";
//N = 13;
//p = 2;


Chi = zeros(N,1); 
for l = 1:1
    data_type = 'raw';
    xi = path + "bigs/'+data_type+'/big_xi" + string(l) + ".txt";
    tau = path + "bigs/'+data_type+'/big_tau" + string(l) + ".txt";
    
    
    tau = read(tau, -1, 1);
    xi = read(xi, -1, N);
    
    i = 30;
    xi = [xi(:,1:i), xi(:,(i+2):N)];
    N = N-1;
    Chi = zeros(N,1); 

    sz = size(xi);
    
    for i = 1:p
        taus(:,i) = tau(i:p:sz(1),:);
    end

    m = sz(1) / p;

    for i = 1:p
        norma(i) = norm(taus(:, i))
    end

    winH=waitbar('Нормировка данных...');
    for i = 1:sz(1)
        waitbar(i/sz(1),winH);
        tau(i) = tau(i) / norma(pmodulo(i-1, p)+1);
        xi(i, :) = xi(i, :) / norma(pmodulo(i-1, p)+1);
    end
    close(winH);
    
    Chi(:,l) = inv(xi' * xi) * xi' * tau;
    
//    scf()
//    ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'
    
    tau_calc = xi * Chi;
//    tau_calc_real = xi * ChiReal;
    for i = 1:p
//        subplot(2, p, i+p);
        subplot(2, 3, i);
        plot2d(1:m, tau(i:p:sz(1)), 3);
        plot2d(1:m, tau_calc(i:p:sz(1)), 1);
//        plot2d(1:m, tau_calc_real(i:p:sz(1)), 5);
//        legend("raw data", "ident data", 'calc data')
        a = gca()
        a.title.text = 'Calculations for normalized data. Link ' + string(i);
    end                                                           
    
end    

disp("Normalized calcs")
disp('cond(xi) = ' + string(cond(xi)))
disp('rank([xi]) = ' + string(rank([xi])))
disp('rank([xi, tau_calc]) = ' + string(rank([xi, tau_calc])))


