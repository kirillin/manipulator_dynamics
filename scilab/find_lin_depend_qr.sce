clc();
clear;
N = 21;
for i = 1:N
    name = "/media/data/evo/robotics_report/ros_packages/youbot_arm_control/sympy/EE/EE" + string(i-1) + ".txt";
    r = read(name, -1, 5);
    EE1(:,:,i) = r'
end
s = EE1(:,:,N);
printf("Succuessfully read; size of matrices: %dx%d\n", s(1), s(2))


E = EE1;
s = size(E)
nc = max(s) //number of columns
nr = s(1);             //number of rows

count = 0
total = factorial(nc) / (factorial(nc - nr) * factorial(nr));

N = 21
//must be nr nested loops below
for i = 1:nc
    for j = (i+1):nc
        printf("Checking (%2d)(%2d)(...) combination; ", i, j)
        printf("%3d%% of work done\n", floor(count/total*100))
        for k = (j+1):nc
            for x = (k+1):nc
                for y = (x+1):nc
                    count = count + 1
                    zero_count = 0;
                    for z = 1:N
                        A = [E(:,i,z), E(:,j,z), E(:,k,z), E(:,x,z), E(:,y,z)];
                        [Q,R] = qr(A);
                        dR = diag(R);
                        flag = %f;
                        for n = 1:5
                            if abs(dR(n)) < 10^(-12) then
                                flag = %t;
                                break;
                            end
                        end
                        if flag then
                            zero_count = zero_count+1;
                        end
                    end
                    if zero_count == N then
                        printf("Attention: %d, %d, %d, %d, %d!!!\n", i,j,k,x,y);
                        disp(A,Q,R)
                    end
                end
            end
        end
    end
end
