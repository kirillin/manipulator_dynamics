clc();
clear;

function ans = read_big_xi(spath, qty_files, qty_joints)
    for i = 1:qty_files
        path = sprintf(spath, string(i-1));
        Xi = read(path, -1, qty_joints);
        h_bigXi(:,:,i) = Xi;
    end
    bigXi = h_bigXi(:,:)';
    s = size(bigXi);
    printf("Succuessfully read; size of matrices: %dx%d\n", s(1), s(2))
    ans = bigXi;
endfunction

function ans=C(n, k)
    ans = factorial(n) / (factorial(k) * factorial(n - k));
endfunction

function res=combinations(n, k)
    // Returns combinations of a number sequance from 1 to n for k.
    // Example: 
    // -->combinations(3,2)
    // ans  =
    //    1.    2.  
    //    1.    3.  
    //    2.    3. 
    combs = [1:k];
    p = 1;
    for i = 1:C(n,k) do
        A = combs(p, :);
        j = k;
        while j >= 1 do
            if A(j) < n - k + j then
                A(j) = A(j) + 1;
                for l = j:k-1 do 
                    A(l + 1) = A(l) + 1;
                end
                p = p + 1;
                combs(p, :) = A;
                break;
            end
            j = j - 1;
        end
    end
    res = combs;
endfunction

function res = get_cols_combinations(cols, k)
    // Returns combinations all columns for k.
    // Example:
    // -->get_cols_combinations([3,13,42], 2)
    // ans  =
    //    3.     13.  
    //    3.     42.  
    //    13.    42.  
    qty_cols = length(cols);
    combs = combinations(qty_cols, k);
    n = length(combs(:,1))
    cols_combs = zeros(n, k);
    for i = 1:n do
        for j = 1:k do
            cols_combs(i, j) = cols(combs(i,j));
        end
    end
    res = cols_combs;
endfunction

function init = converter(cols)
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
//    no =  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17, 18, 20, 28, 31, 32, 42, 44, 45, 46, 15, 14]+1;
    no = [];// + 1;
    ind = 0;
    for i = 1:70
        in = %t;
        for j = 1:length(no)
            if no(j) == i then
                in = %f;
            end
        end
        if in then 
            ind  = ind + 1;
        end
        if ind == cols
            init = i-1
            return 
        end
    end
endfunction

path = get_absolute_file_path("fld_test.sce");
path = path + "../data_for_identification_Youbot/ee/EE%s.txt";
bigXi = read_big_xi(path, 20, 5);   // for example: 100x70

//bigXi = [1,2,3,4,5,6;7,8,9,10,11,12];
//bigXi = [   1,2,3,4;
//            7,8,9,10];
//bigXi = [1,2,3;3,5,9];
[rosw, qty_cols] = size(bigXi);

for i = 1:qty_cols
    printf("%d (%d) column depend on\n", i, converter(i))
    B = bigXi(:, i);
    select i
        case 1 then
            cols = [2:qty_cols];
        case 2 then
            cols = [1,3:qty_cols];
        case qty_cols - 1 then
            cols = [1:qty_cols-2,qty_cols];
        case qty_cols then
            cols = [1:qty_cols-1];
        else 
            cols = [1:i-1,i+1:qty_cols];
    end
  
    for j = 2 do
        cols_combs = get_cols_combinations(cols, j);
        for k = 1:length(cols_combs(:,1)) do // all combinations
            columns = cols_combs(k,:);
            Ac = bigXi(:,[columns]);
            if rank(Ac) >= rank([Ac,B]) then
                Ab=[Ac,-B];[ma,na]=size(Ab);
                [W,rk]=colcomp(Ab);
                W=W(:,1:na-rk);last=W(na,:);
                [W2,rk1]=colcomp(last);
                if rk1==0 then
                    continue;
                end
                [x0, kerA] = linsolve(Ac,-B);
                if kerA == [] & x0 <> [] then
                    disp(x0)
                    // disp('one solve')
                    res = [];
                    p = 1;
                    for l = 1:length(columns) do
                        if length(columns) > 1 then
                            if x0(l) <> 0 then
                                res(:,p) = [columns(l), x0(l)]';
                                p = p + 1;
                            end 
                        else 
                            if x0 <> 0 then
                                res = [columns, x0]';
                            end 
                        end
                    end
                    if res <> [] then
                        printf("New combination:\n");
                        disp(res);
                    end
                elseif kerA <> [] & x0 <> [] then
                    continue;
                    // disp('inf')
                end
            end 
        end                    
    end
end
