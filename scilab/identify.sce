clc();
//clear;
stacksize('max');
gstacksize('max');
data_type = 'filt';

path = get_absolute_file_path('identify.sce') + '../data_for_identification_Youbot/';
N = 51;
//N = 70;
p = 5;

//path = get_absolute_file_path('identify.sce') + "../data_for_identification_2dof/";
//N = 28; // full regressor
////N = 13;
//p = 2;
//

Chi = zeros(N,1); 
for l = 1:1
    xi = path + 'bigs/'+data_type+'/big_xi' + string(l) + ".txt";
    tau = path + 'bigs/'+data_type+'/big_tau' + string(l) + ".txt";
   
    Tau = read(tau, -1, 1)
    Xi = read(xi, -1, N);
    [r, c] = size(Xi);

    /////////////////ONE    
//    scf();
//    ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'
//    for i = 1:length(ChiReal) do
//        plot2d(1:r/p, ones(r/p,1) * ChiReal(i), i)
//    end;
    
//    [Xi, Tau] = normirovka(p, Xi, Tau)
    
    est_Chi = rls(p, Xi, Tau, %f);
    Tau_calc_rls = Xi * est_Chi;
    
    //////////////////TWO
//    [U,S,V] = svd([Xi]);
//    Splus = inv(S' * S) * S';
//    Chi_ols = V * Splus * U' * Tau;
//    Tau_calc_ols = Xi * Chi_ols;
//    

    //////////////////THREE
//    Chi(:,l) = inv(Xi' * Xi + 0.001 * eye(N,N)) * Xi' * Tau;
//    Tau_calc = Xi * Chi;
//    
    //plot taus
    scf()
    for i = 1:p
        subplot(1, p, i);

        plot2d(1:r/p, Tau(i:p:r), 3);
//        plot2d(1:r/p, Tau_calc(i:p:r), 1);
        plot2d(1:r/p, Tau_calc_rls(i:p:r), 6);
//        plot2d(1:r/p, Tau_calc_ols(i:p:r), 2);
        
        a = gca()
        a.title.text = 'Usual calculations. Link ' + string(i);
    end
    legend("raw data", 'rls')
end    

disp("Usual calcs")
disp('cond(xi) = ' + string(cond(Xi)))
disp('rank([xi]) = ' + string(rank([Xi])))
//disp('rank([xi, tau_calc]) = ' + string(rank([xi, tau_calc])))

