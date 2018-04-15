clc();
clear;
stacksize('max');
path = "/media/kirix/data/thesis/manipulator_dynamics/data_for_identification_Youbot/";
//path = "d:\thesis\manipulator_dynamics\data_for_identification_Youbot\";
N = 46;
p = 5;

R = 10; // number of data files
Chi = zeros(N, R); 
scf();
est_Chi_k_1 = zeros(N, 1);
for l = 1:R
    disp(l);
    // READING DATA
    data_type = 'filt';
    xi = path + "bigs/"+data_type+"/big_xi" + string(l) + ".txt";
    tau = path + "bigs/"+data_type+"/big_tau" + string(l) + ".txt";

    tau = read(tau, -1, 1);
    xi = read(xi, -1, N);

    // NORMALIZING DATA
    sz = size(xi);
//    for i = 1:p
//        taus(:,i) = tau(i:p:sz(1),:);
//    end
    m = sz(1) / p;
//    for i = 1:p
//        norma(i) = norm(taus(:, i))
//    end
//    winH=waitbar('Нормировка данных...');
//    for i = 1:sz(1)
//        waitbar(i/sz(1),winH);
//        tau(i) = tau(i) / norma(pmodulo(i-1, p)+1);
//        xi(i, :) = xi(i, :) / norma(pmodulo(i-1, p)+1);
//    end
//    close(winH);

    //////////////////////////////////////////////
    /////// ESTIMATION CHI PARAMETERS ////////////
    //////////////////////////////////////////////
    
    // simple methods
    // Chi(:,l) = xi \ tau;
    // Chi(:,l) = lsq(xi, tau);

    // 1.1 Ordinary LS whith 0.001 * eye(N,N)
    // Chi(:,l) = inv(xi' * xi + 0.001 * eye(N,N)) * xi' * tau;    
    // 1.2 Ordinary LS whith removed liner colomns
    // Chi(:,l) = inv(xi' * xi) * xi' * tau;    
    // 1.3 Numerical implimentation of the linear Ordinary LS
    //[U,S,V] = svd([xi]);
    //Splus = inv(S' * S) * S';
    //Chi(:,l) = V * Splus * U' * Tau;
    
    // 2. Total LS (need scilab 6 for mooore memory)
    // [U,S,V] = svd([xi, tau]);
    // Chi(:,l) = -V(1:$-1, $) / V($,$);

    // 3. Recursive Least Squares
    r = 1;
    xi_k = xi((r*p-p+1):(r*p), :)';
    P_k_1 = eye(N, N);
    g_k = P_k_1 * xi_k;
    eps_k_1 = zeros(p, 1);

 
    for k = 2:m do
        xi_k = xi(r*p-p+1:r*p, :)';
        y_k = tau(r*p-p+1:r*p,:);
        r = r + 1;
        
        eps_k = y_k - xi_k' * est_Chi_k_1;
        
        g_k = P_k_1 * xi_k * inv(1*eye(p,p) + xi_k' * P_k_1 * xi_k);
        P_k = P_k_1 - g_k * xi_k' * P_k_1;
        est_Chi_k = est_Chi_k_1 + g_k * eps_k;

        // plot eps
        // plot(k-1:k, [eps_k_1, eps_k]);

        ////plot parameters
        //s = est_Chi_k_1;
        //e = est_Chi_k;
        //plot(k-1:k, [s,e]);
        //disp(k);
        
        eps_k_1 = eps_k;
        P_k_1 = P_k;
        est_Chi_k_1 = est_Chi_k;
    end;
    Chi(:,l) = est_Chi_k;
    
    // CLEAR MEMORY
    clear tau;
    clear xi;
    clear taus;
end;

//tau_calc = xi * Chi(:,$);
//
//for i = 1:p
//    subplot(2, 3, i);
//    plot2d(1:m, tau(i:p:sz(1)), 3);
//    plot2d(1:m, tau_calc(i:p:sz(1)), 1);
//    legend("raw data", 'calc data')
//    a = gca()
//    a.title.text = 'Link ' + string(i);
//end 
//
                                        



