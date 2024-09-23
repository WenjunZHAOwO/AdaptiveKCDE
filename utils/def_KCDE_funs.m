G = @(x,h,m) 1./(sqrt(2*pi*h.^2)).*exp( -1/2*( x - m').^2./h./h);

switch size(x,1)
    case 1
        Gx = G;
    case 2
       
        if center_x
            Gx = @(x,h,m) G(x(1,:),h(1,:),m(1,:)).*G(x(2,:),h(2,:),m(2,:));
        else
            Gx = @(x,h,m) G(x(1,:),h(:,1),m(1,:)).*G(x(2,:),h(:,2),m(2,:));
        end
end


if optim_flag == 0
    switch center_x
        case true
            
            KCDE_fun = @(xstar, zstar, hx, hz) ( Gz(zstar,hz,z)./sum( Gz(zstar,hz,z)  ))'*Gx(xstar,hx(zstar),x);
            KCDE_fun_loo = @(xstar, zstar, hx, hz, idx) ( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))./sum( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))  ))'*Gx(xstar,hx(zstar),x(:,[1:idx-1,idx+1:N]));
            KCDE_fun_fast = @(hx, hz) diag( ( Gz(z,hz,z)./sum( Gz(z,hz,z)  ))'*Gx(x,hx(z),x) );
            
            KCDE_fun_loo_fast = @(hx, hz) diag( ( Gz(z,hz,z).*(ones(N)-eye(N))./sum( Gz(z,hz,z).*(ones(N)-eye(N))  ))'*(Gx(x,hx(z),x).*(ones(N)-eye(N))) );
            

       
        case false
            
            KCDE_fun = @(xstar, zstar, hx, hz) ( Gz(zstar,hz,z)./sum( Gz(zstar,hz,z)  ))'*Gx(xstar,hx(z)',x);
            KCDE_fun_loo = @(xstar, zstar, hx, hz, idx) ( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))./sum( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))  ))'*Gx(xstar,subsref(hx(z),struct('type','()','subs', {{1:D,[1:idx-1,idx+1:N]}}))',x(:,[1:idx-1,idx+1:N]));
            KCDE_fun_fast = @(hx, hz) diag( ( Gz(z,hz,z)./sum( Gz(z,hz,z)  ))'*Gx(x,hx(z)',x) );
            
            KCDE_fun_loo_fast = @(hx, hz) diag( ( Gz(z,hz,z).*(ones(N)-eye(N))./sum( Gz(z,hz,z).*(ones(N)-eye(N))  ))'*(Gx(x,hx(z)',x).*(ones(N)-eye(N))) );
            
           
    end
    
else
    switch center_x
        
        case true
%             hz = hz(z);
            KCDE_fun = @(xstar, zstar) ( Gz(zstar,hz,z)./sum( Gz(zstar,hz,z)  ))'*Gx(xstar,hx(zstar),x);
            KCDE_fun_loo = @(xstar, zstar,idx) ( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))./sum( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))  ))'*Gx(xstar,hx(zstar),x(:,[1:idx-1,idx+1:N]));
            
       
        case false
%             hz = hz(z);
            hx = hx(z);
            KCDE_fun = @(xstar, zstar) ( Gz(zstar,hz,z)./sum( Gz(zstar,hz,z)  ))'*Gx(xstar,hx',x);
            KCDE_fun_loo = @(xstar, zstar, idx) ( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))./sum( Gz(zstar,hz,z(:,[1:idx-1,idx+1:N]))  ))'*Gx(xstar,subsref(hx,struct('type','()','subs', {{1:D,[1:idx-1,idx+1:N]}}))',x(:,[1:idx-1,idx+1:N]));
            %sum( Gz(zstar,subsref(hz,struct('type','()','subs', {{[1:idx-1,idx+1:N]}})) ,z([1:idx-1,idx+1:N])).*Gx(xstar,subsref(hx(z),struct('type','()','subs', {{[1:idx-1,idx+1:N]}})) ,x([1:idx-1,idx+1:N]))./sum( Gz(zstar,subsref(hz(z),struct('type','()','subs', {{[1:idx-1,idx+1:N]}})) ,z([1:idx-1,idx+1:N]))  ));
    end

end

hz_fun = @(beta, zstar) beta + 0*zstar;
    