function d = pqdistances(p, q, dim,  type, rangex, N)


if nargin < 3
    dim = 1;
end

if nargin < 4
    type = 'Hellinger';
end

if nargin < 5
    rangex = ones(dim,1)*[-20,20];
end

if nargin < 6
    if dim == 1
        N = 400;%max(4000, round(100*(max(rangex(:,2)-rangex(:,1))))) ;
    else
        N = 100;%max(100, round(10*(max(rangex(:,2)-rangex(:,1))))) ;
    end
end

if dim == 1
    x = linspace(rangex(1), rangex(2), N);
    switch type
        case 'Hellinger'
            d = 1 - trapz( x, sqrt( p(x).*q(x) ));
        case 'ISE'
            d = trapz( x, (p(x) - q(x)).^2 );
	case 'TV'
	    d = trapz(x, abs(p(x) - q(x)) );
    end
else
    x1 = linspace(rangex(1,1), rangex(1,2), N);
    x2 = linspace(rangex(2,1), rangex(2,2), N);
    [X1,X2] = meshgrid(x1,x2);
    
    switch type
        case 'Hellinger'
            d = 1- trapz(x2, trapz(x1,sqrt( reshape(p([X1(:),X2(:)]').*q([X1(:),X2(:)]'), N, N )), 2));
        case 'ISE'
            d = trapz(x2, trapz(x1,sqrt( reshape( (p([X1(:),X2(:)]') - q([X1(:),X2(:)]')).^2, N, N )), 2));
    end
end

d = max(d,0);
% d = sqrt(d);

