clc();
clear;
N = 20;
p = 2;
for i = 1:N
    name = "/media/data/evo/manipulator_dynamics/data_for_identification_2dof/ee/EE" + string(i-1) + ".txt";
    r = read(name, -1, p);
    EE1(:,:,i) = r'
end
s = size(EE1(:,:,N));
printf("Succuessfully read; size of matrices: %dx%d\n", s(1), s(2))

global A B pos columns flag

function init = converter(cols)
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
//    no =  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17, 18, 20, 28, 31, 32, 42, 44, 45, 46, 15, 14]+1;
    no = [3,4,5,7,8,9,16,18,20,21,22,23 ,0,1,6] + 1;
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

function make_combination(n,m, num)
    global A B pos columns flag
    pos = pos + 1
    for i = num:n
        if the_column <= i then
            columns(1,pos) = i+1;
        else
            columns(1,pos) = i;
        end
        for j = 1:n_E
            A((p*j-p+1):(p*j),pos) = Eh(:,i,j)
            B((p*j-p+1):(p*j)) = aim_vect(:,j)
        end
        if m == 1 then
            if rank(A) >= rank([A,B]) then
                [x0, kerA] = linsolve(A,-B);
                if %t then
                    k = clean(x0);
                    printf("Column No %2d (%d) is dependent on cols No\n", the_column, converter(the_column));

                    for r = 1:length(k) do
                        if k(r) <> 0.0 then
                            printf("%d (%d)\t%f\n", columns(r), converter(columns(r)), k(r));
                        end
                    end;
                    printf(' len(x0)=%d len(ker)=%d\n', length(x0), length(kerA))
//                    disp(columns);
//                    printf("with koefficients");
//                    disp(k');
//                    printf("\n")
                 end
            end 
            flag = %t
        else
            make_combination(n, m-1, i+1);
            if flag == %t then
                return;
            end;
        end
    end
    pos = pos - 1
endfunction


//values of regressor for several random values of q, dq, ddq
E = EE1;   


size_E = size(E);
n_dots = size_E(3);

nc = max(size(E(:,:,1)));

for i = 1:nc
//    printf("checking of %d column:\n", i)
    //prepare helpful variables
    for k = 1:n_dots
        aim_vect(:,k) = E(:,i,k);
        ind = 1;
        for kk = 1:nc
            if kk <> i then
                Eh(:,ind,k) = E(:,kk,k)
                ind = ind + 1;
            end 
        end
    end
    
    l = nc - 1
    //find linear dependence for i-th column
    for j = l:l
        n_E = N //ceil(j/5) 
        A = zeros(n_E * p, j);
        B = zeros(n_E * p, 1);
        the_column = i;
        columns = zeros(1,j);
        pos = 0;
        flag = %f;
        make_combination((nc-1), j, 1)
        
        printf("Checked dependency on %d/%d columns\n", j, nc);
    end
end
