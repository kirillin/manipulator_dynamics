/////////////////////////////////////////////////////////
////////////// RECURSIVE LEAST SQUARES ALHORITHM ////////
/////////////////////////////////////////////////////////
function [est_theta]=rls(n, Omega, Y, isplot)
    [qty_rows, qty_paramerers] = size(Omega);
    qty_mesuarments = qty_rows / n; // quantity regressors
    
    prev_P = eye(qty_paramerers, qty_paramerers);
    prev_est_theta = zeros(qty_paramerers, 1);
    prev_eps = zeros(n, 1);
    
    winRLS = waitbar('Рекуррентный МНК...');
    
    for k = 1:qty_mesuarments do
        waitbar(k / qty_mesuarments, winRLS);
        omega = Omega(k*n-n+1:k*p, :)';
        y = Y(k*n-n+1:k*p,:);
        
        eps = y - omega' * prev_est_theta;
        g = prev_P * omega * inv(eye(n,n) + omega' * prev_P * omega);
        est_theta = prev_est_theta + g * eps;
        P = prev_P - g * omega' * prev_P;        

        if isplot == %t then
            plot([k-1:k]', [prev_est_theta, est_theta]');
            printf('%d\n', k);
        end;
    
        prev_P = P;
        prev_est_theta = est_theta;
        prev_eps = eps;
    end;
    close(winRLS);
endfunction

function [Omega, Y]=normirovka(n, Omega, Y)
    [qty_rows, qty_paramerers] = size(Omega);
    qty_mesuarments = qty_rows / n; // quantity regressors

    for i = 1:n
        Ys(:,i) = Y(i:n:qty_rows,:);
        norma(i) = norm(Ys(:, i))
    end

    winH=waitbar('Нормировка данных...');
    
    for i = 1:qty_rows do
        waitbar(i / qty_rows,winH);
        Y(i) = Y(i) / norma(pmodulo(i-1, n)+1);
        Omega(i, :) = Omega(i, :) / norma(pmodulo(i-1, n)+1);
    end
    close(winH);
endfunction
