clc();
clear;
stacksize('max');
path = "/media/data/evo/manipulator_dynamics/data_for_identification_Youbot/";
N = 51;
p = 5;

//path = "/media/data/evo/manipulator_dynamics/data_for_identification_2dof/";
//N = 13;
//p = 2;
//

Chi = zeros(N,1); 
for l = 1:1
    data_type = 'raw';
    xi = path + "bigs/'+data_type+'/big_xi" + string(l) + ".txt";
    tau = path + "bigs/'+data_type+'/big_tau" + string(l) + ".txt";
   
    tau = read(tau, -1, 1)
    xi = read(xi, -1, N);
    
//    i = 30;
//    xi = [xi(:,1:i), xi(:,(i+2):N)];
//    N = N-1;
//    Chi = zeros(N,1); 
//    for i = 1:N-2 do
//        xi1 = [xi(:,1:i), xi(:,(i+2):N)];
//        printf("%f\t%f %d\n", cond(xi1), rank(xi1),i+1)
//    end
    
//    for i = 1:N-2 do
//        xi1 = [xi(:,1:i), xi(:,(i+2):N)];
//        for j = 1:N-3 do
//            xi2 = [xi1(:,1:j), xi1(:,(j+2):N-1)];
//            printf("%f\t%f ij = %d %d\n", cond(xi2), rank(xi2),i,j) 
//        end
//    end
//    break;

    sz = size(xi);
    
    m = sz(1) / p;

    Chi(:,l) = inv(xi' * xi) * xi' * tau;
 
    scf()

//    ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'
    
    tau_calc = xi * Chi;
    chi2 = xi \ tau_calc;
    chi3 = lsq(xi, tau_calc);
    
//    tau_calc_real = xi * ChiReal;
    for i = 1:p
        subplot(2, p, i);
        plot2d(1:m, tau(i:p:sz(1)), 3);
        plot2d(1:m, tau_calc(i:p:sz(1)), 1);
//        plot2d(1:m, tau_calc_real(i:p:sz(1)), 5);

        a = gca()
        a.title.text = 'Usual calculations. Link ' + string(i);
    end
    legend("raw data", "ident data", 'calc data')
end    

disp("Usual calcs")
disp('cond(xi) = ' + string(cond(xi)))
disp('rank([xi]) = ' + string(rank([xi])))
//disp('rank([xi, tau_calc]) = ' + string(rank([xi, tau_calc])))


