function score = scores_KCDE2(x,z,KCDE_fun, loo, fast, KCDE_fun_fast, w)

if nargin < 5
    fast = false;
end
if nargin < 7
    w = ones(1,size(z,2));
end

N = size(x,2);

if fast

    score = sum( w'.*log(KCDE_fun_fast))/N;
else
    score = 0;
    

    for i = 1:N
        if loo
            score = score + w(i)*log(KCDE_fun(x(:,i),z(:,i),i));
        else
            score = score + w(i)*log(KCDE_fun(x(:,i),z(:,i)));
        end
    end
    score =  sum( score )/N;
end

