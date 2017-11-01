clc();
clear;
stacksize('max');
path = "/media/data/evo/robotics_report/symbolic_computation/data_for_identification/";
Chi = zeros(50,10); 
for l = 1:1
    
    xi = path + "bigs/filt/big_xi" + string(l) + ".txt";
    tau = path + "bigs/filt/big_tau" + string(l) + ".txt";
    
    //xi = '~/bigX.txt';
    //tau = '~/Y.txt';
    
    N = 50;
    
    tau = read(tau, -1, 1);
    xi = read(xi, -1, N);
    
    sz = size(xi);
    
    //xis = zeros(sz(1)/5, sz(2), 5);
//    for i = 1:5
//        xis_raw(:,:,i) = xi(i:5:sz(1),:);
//        taus(:,i) = tau(i:5:sz(1),:);
//    end
////    
//    function answer=in(num, vec)
//        for i = 1:length(vec)
//            if num == vec(i) then
//                answer = %t;
//                return;
//            end
//        end
//        answer = %f;
//    endfunction
//    
//    bad_cols = [0,0; 0, 25; 13, 25; 25,34; 37,38] ;
//    for i = 1:5
//        ind = 1;
//        for j = 1:N    
//            if ~in(j, bad_cols(i,:)) then
//                xis(:, ind, i) = xis_raw(:, j, i);
//                ind = ind + 1;
//            end
//        end
//    end
//    
    m = sz(1) / 5;
//    
//    counts = zeros(5,1)
//    for i = 1:5    
//        ind = 1;
//        for j = 1:N
//            help_matr = clean(xis(:,j,i), 10^-20);
//            if sum(help_matr == zeros(m, 1)) == m then
//                //do nothing
//            else
//                xis_lite(:,ind,i) = xis(:, j, i);
//                ind = ind + 1;
//                counts(i) = counts(i) + 1;
//            end
//        end
//    end
//    
//    Xi1 = xis_lite(:,1:counts(1), 1);
//    Xi2 = xis_lite(:,1:counts(2), 2);
//    Xi3 = xis_lite(:,1:counts(3), 3);
//    Xi4 = xis_lite(:,1:counts(4), 4);
//    Xi5 = xis_lite(:,1:counts(5), 5);
//    
//    chi1 = inv(Xi1' * Xi1) * Xi1' * taus(:, 1)
//    chi2 = inv(Xi2' * Xi2) * Xi2' * taus(:, 2)
//    chi3 = inv(Xi3' * Xi3) * Xi3' * taus(:, 3)
//    chi4 = inv(Xi4' * Xi4) * Xi4' * taus(:, 4)
//    chi5 = inv(Xi5' * Xi5) * Xi5' * taus(:, 5)
//    
//    tau_calc(:,1) = Xi1 * chi1;
//    tau_calc(:,2) = Xi2 * chi2;
//    tau_calc(:,3) = Xi3 * chi3;
//    tau_calc(:,4) = Xi4 * chi4;
//    tau_calc(:,5) = Xi5 * chi5;
//    
//    
//    for i = 1:5
//        subplot(2, 3, i);
//        plot2d(1:m, taus(:,i), 3);
//        plot2d(1:m, tau_calc(:,i), 1);
//    end
//    
//    
//    for i = 1:5
////        sigma(i) = 1/m * (taus(:,i) - tau_calc(:, i))' * (taus(:,i) - tau_calc(:, i));
//        norma(i) = norm(taus(:, i))
//    end
//    sq_sigma = sigma.^2
//    
//    
//    
//    clear taus tau_calc xis xis_lite 
//    
//    //for i = 1:sz(1)
//    //    tau(i) = tau(i) / sigma(pmodulo(i-1, 5)+1);
//    //    for j = 1:N
//    //        xi(i, j) = xi(i, j) / sigma(pmodulo(i-1, 5)+1);
//    //        printf("Now i = %d; j = %d\r", i, j)
//    //    end
//    //end
//    //
//    
//    
//    for i = 1:sz(1)
//        tau(i) = tau(i) / norma(pmodulo(i-1, 5)+1);
//        for j = 1:N
//            xi(i, j) = xi(i, j) / norma(pmodulo(i-1, 5)+1);
//            printf("Now i = %d; j = %d\n", i, j)
//        end
//    end
////    
    Chi(:,l) = inv(xi' * xi) * xi' * tau;
    
    scf()
    
    tau_calc = xi * Chi;
    for i = 1:5
        subplot(2, 3, i);
        plot2d(1:m, tau(i:5:sz(1)), 3);
        plot2d(1:m, tau_calc(i:5:sz(1)), 1);
    end
end    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
