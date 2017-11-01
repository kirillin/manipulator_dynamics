clc();
clear;
N = 50;
for i = 1:N
    name = "/media/data/evo/robotics_report/symbolic_computation/data_for_identification/ee/EE" + string(i-1) + ".txt";
    r = read(name, -1, 5);
    EE1(:,:,i) = r'
end
s = size(EE1(:,:,N));
printf("Succuessfully read; size of matrices: %dx%d\n", s(1), s(2))

global A B pos columns flag

function init = converter(cols)
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
//    no =  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17, 18, 20, 28, 31, 32, 42, 44, 45, 46, 15, 14]+1;
    no = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 18, 20, 28, 31, 32, 44, 45, 46, 47] + 1;
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
            A((5*j-4):(5*j),pos) = Eh(:,i,j)
            B((5*j-4):(5*j)) = aim_vect(:,j)
        end
        if m == 1 then
            if rank(A) >= rank([A,B]) then
                k = linsolve(A,-B);
                k = clean(k);
                printf("Column No %2d (%d) is dependent on cols No", the_column, converter(the_column));
                disp(columns);
                printf("with koefficients");
                disp(k');
                printf("\n")
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
    printf("checking of %d column:\n", i)
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
        A = zeros(n_E * 5, j);
        B = zeros(n_E * 5, 1);
        the_column = i;
        columns = zeros(1,j);
        pos = 0;
        flag = %f;
        make_combination((nc-1), j, 1)
        
        printf("Checked dependency on %d/%d columns\n", j, nc);
    end
end
