////////////////////////////////////////////////////////////////////////
//////// SCRIPT FOR FINDING THE LINEAR COLUMNS IN MATRIX ///////////////
////////////////////////////////////////////////////////////////////////
clc();
clear;
stacksize('max');

////////////////////////////////////////////////////////////////////////
////////// INITIALIZATION //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

path = get_absolute_file_path("find_lin_depend_v2.sce");
exec(path+"nchoosek.sci");
path = path + "../data_for_identification_Youbot/ee/EE%s.txt";
//path = path + "../data_for_identification_2dof/ee/EE%s.txt";
QTY_FILES = 20;
QTY_JOINTS = 5;
TOTAL_COLS = 70;

QTY_COMBINATIONS = 0;   // if 0 then all combinations
QTY_DEPENDS_TEST = 2;
LINEAR_COLS = [2,3,4,6,7,8,9,0,1,5,17,18,20,31,10,14,32,  15,16,19,21,22,23,61] + 1;
//LINEAR_COLS = [3,4,5,7,8,9,16,18,20,21,22,23,6,0,1] + 1;  // 2dof


////////////////////////////////////////////////////////////////////////
////////// FUNCTIONS FOR SOME WORKS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////
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

function res=combinations(n, k) //very slow func
    // Returns combinations of a number sequance from 1 to n for k.
    // Example: 
    // -->combinations(3,2)
    // ans  =
    //    1.    2.  
    //    1.    3.  
    //    2.    3. 
    combs = [1:k];
    p = 1;
    for i = 1:10000 do //C(n,k) do
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
//    combs = combinations(qty_cols, k);
    combs = nchoosek_re(1:qty_cols, k);
    n = length(combs(:,1))
    cols_combs = zeros(n, k);
    for i = 1:n do
        for j = 1:k do
            cols_combs(i, j) = cols(combs(i,j));
        end
    end
    res = cols_combs;
endfunction

function init = converter(cols, no)
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
    ind = 0;
    for i = 1:TOTAL_COLS do
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

function v=isZeroVec(vec)
    n = max(size(vec))
    flag = %t;
    for i = 1:n
        clear_vec = clean(vec, 10^-15);
        if sum(clear_vec == zeros(vec)) <> n then 
            flag = %f;
            break;
        end    
    end
    v = flag;
endfunction

function find_lin_depend()
    for i = 1:qty_cols
        printf("%d {py:%d} column depend on\n", converter(i,LINEAR_COLS)+1, converter(i,LINEAR_COLS))
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
//        real_cols = []
//        for i = 1:qty_cols do
//            real_cols(i) = converter(i, LINEAR_COLS)+1;
//        end
        for j = QTY_DEPENDS_TEST do
            cols_combs = get_cols_combinations(cols, j);
//            disp('cols_comb ok!!')
//            disp(cols_combs')
            if QTY_COMBINATIONS == 0 then
                qty_combs = length(cols_combs(:,1))
            else
                qty_combs = QTY_COMBINATIONS;
            end
            for k = 1:qty_combs do // all combinations
                columns = cols_combs(k,:);
                Ac = bigXi(:,[columns]);
                if rank(Ac) >= rank([Ac,B]) then
                    Ab=[Ac,-B];[ma,na]=size(Ab);
                    [W,rk]=colcomp(Ab);
                    W=W(:,1:na-rk);last=W(na,:);
                    [W2,rk1]=colcomp(last);
                    if rk1==0 then
                        continue;   // no solves
                    end
                    [x0, kerA] = linsolve(Ac,-B);
                    x0 = clean(x0);
                    if kerA == [] & x0 <> [] then
                        // disp('one solve')
                        res = [];
                        p = 1;
                        for l = 1:length(columns) do
                            if length(columns) > 1 then
                                if x0(l) <> 0 then
                                    res(:,p) = [converter(columns(l),LINEAR_COLS)+1, x0(l)]';
                                    p = p + 1;
                                end 
                            else 
                                if x0 <> 0 then
                                    res = [converter(columns,LINEAR_COLS)+1, x0]';
                                end 
                            end
                        end
                    elseif kerA <> [] & x0 <> [] then
                        //disp('inf')
                        break;
                    end
                    if res <> [] then
                        printf("  New combination:\n");
                        disp(res);
                    end
                end 
            end                    
        end
    end
endfunction

function print_zero_columns()
    printf("For !python! use these indexes:\n");
    for j = 1:qty_cols
        if isZeroVec(bigXi(:, j)) then
            printf('%d,',j-1)
        end
    end
    printf("\nFor !scilab! use these indexes:\n");
    for j = 1:qty_cols
        if isZeroVec(bigXi(:, j)) then
            printf('%d,',j)
        end
    end
endfunction

////////////////////////////////////////////////////////////////////////
////////// USING EXAMPLE ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

bigXi = read_big_xi(path, QTY_FILES, QTY_JOINTS);   // for example: 100x70
[rosw, qty_cols] = size(bigXi);
printf("cond(bigXi) = %f\n",cond(bigXi));
printf("rank(bigXi) = %f\n",rank(bigXi));
//print_zero_columns();
find_lin_depend();


