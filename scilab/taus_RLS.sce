clc();
clear;
stacksize('max');

p = get_absolute_file_path('taus_RLS.sce');
path = p + '../data_for_identification_Youbot/';
N = 51;
//N = 70;
p = 5;

//p = get_absolute_file_path('taus_RLS.sce')
//path = p + "../data_for_identification_2dof/";
////N = 28; // full regressor
//N = 13;
//p = 2;
//

Chi = zeros(N,1); 
for l = 1:1
    data_type = 'filt';
    xi = path + 'bigs/'+data_type+'/big_xi' + string(l) + ".txt";
    tau = path + 'bigs/'+data_type+'/big_tau' + string(l) + ".txt";
   
    Tau = read(tau, -1, 1)
    Xi = read(xi, -1, N);
    
    sz = size(Xi);
    m = sz(1) / p;

    //Chi(:,l) = inv(Xi' * Xi) * Xi' * Tau;
    Chi(:,l) = inv(Xi' * Xi + 0.001 * eye(N,N)) * Xi' * Tau;
    
    //////////////////////////////////////////////////
    ///// 1 METHOD: RLS //////////////////////////////
    //////////////////////////////////////////////////
    //init
    r = 1;
    xi_k = Xi((r*p-p+1):(r*p), :)';
    P(:,:,1) = eye(N, N);
    g(:,:,1) = P(:,:,1) * xi_k;
    // est_Chi(:,:,1) = Chi;
    est_Chi(:,:,1) = zeros(N, 1);
    eps(:,:,1) = zeros(p, 1);
   
    // compute and plot
    scf();
    // plot truth parameters. if is
//    ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'
//    for i=1:length(ChiReal) do
//        plot2d(1:m, ones(m,1) * ChiReal(i), i)
//    end;
//    for k = 2:m do
//        xi_k = Xi(r*p-p+1:r*p, :)';
//        y(:,:,k) = Tau(r*p-p+1:r*p,:);
//        r = r + 1;
//        
//        eps(:,:,k) = y(:,:,k) - xi_k' * est_Chi(:,:,k-1);
//
//        g(:,:,k) = P(:,:,k-1) * xi_k * inv(1*eye(p,p) + xi_k' * P(:,:,k-1) * xi_k);
//        P(:,:,k) = P(:,:,k-1) - g(:,:,k) * xi_k' * P(:,:,k-1);
//        est_Chi(:,:,k) = est_Chi(:,:,k-1) + g(:,:,k) * eps(:,:,k);
//
//        // plot eps
//        //disp(sum(eps(:,:,k)));
//        //plot(k-1:k, [eps(:,:,k-1), eps(:,:,k)]);
//
//        //plot parameters
////        s = est_Chi(:,:,k-1);
////        e = est_Chi(:,:,k)
////        plot(k-1:k, [s,e]);
////        disp(k);
//    end;

    //////////////////////////////////////////////////
    ///// 2 METHOD: RLS //////////////////////////////
    //////////////////////////////////////////////////
    [U,S,V] = svd([Xi]);
    Splus = inv(S' * S) * S';
    Chi_ols = V * Splus * U' * Tau;
    //Chi2 = Xi \ tau_calc;
    //Chi3 = lsq(Xi, Tau_calc);


    Tau_calc = Xi * Chi;
    Tau_calc_ols = Xi * Chi_ols;
    Tau_calc_rls = Xi * est_Chi(:,:,m);
    //Tau_calc_real = xi * ChiReal;

    //plot taus
    scf()
    for i = 1:p
        subplot(1, p, i);

        plot2d(1:m, Tau(i:p:sz(1)), 3);
        plot2d(1:m, Tau_calc(i:p:sz(1)), 1);
        //plot2d(1:m, Tau_calc_ols(i:p:sz(1)), 5);
        plot2d(1:m, Tau_calc_rls(i:p:sz(1)), 6);

        a = gca()
        a.title.text = 'Usual calculations. Link ' + string(i);
    end
    legend("raw data", "ls", 'ols_svd', 'rls')
end    

disp("Usual calcs")
disp('cond(xi) = ' + string(cond(Xi)))
disp('rank([xi]) = ' + string(rank([Xi])))
//disp('rank([xi, tau_calc]) = ' + string(rank([xi, tau_calc])))

