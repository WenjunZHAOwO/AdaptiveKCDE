function Gz = Gz_param(z, hz, mz, z_type)

L = size(z,1);



Gz = ones(size(mz,2),size(z,2));
% size(mz)
% size(z)
% size(hz)

for l = 1:L

    switch z_type(l)
        case 'd'
%             if size(hz,1) == L
                Gz = Gz.* ( (1/2+1/2*cos(2*pi*hz(l)).^2).*(z(l,:)==mz(l,:)') + (1/2-1/2*cos(2*pi*hz(l)).^2).*(z(l,:)~=mz(l,:)') );
%             else
%                 Gz = Gz.* ( cos(2*pi*hz(:,l)).^2.*(z(l,:)==mz(l,:)') + sin(2*pi*hz(:,l)).^2.*(z(l,:)~=mz(l,:)') );
%             end
        case 'c'
%             if size(hz,1) == L
                Gz = Gz.*  (1./(sqrt(2*pi*hz(l).^2)).*exp( -1/2*( z(l,:) - mz(l,:)').^2./hz(l)./hz(l)));
%             else
%                 Gz = Gz.*  (1./(sqrt(2*pi*hz(:,l).^2)).*exp( -1/2*( z(l,:) - mz(l,:)').^2./hz(:,l)./hz(:,l)));
%             end

        case 'p'
     %           Gz = Gz.*  (1./(sqrt(2*pi*hz(l).^2)).*exp( -1/2*( cos(2*pi*z(l,:)) - cos(2*pi*mz(l,:)')).^2./hz(l)./hz(l)));
     %           Gz = Gz.*  (1./(sqrt(2*pi*hz(l).^2)).*exp( -1/2*( sin(2*pi*z(l,:)) - sin(2*pi*mz(l,:)')).^2./hz(l)./hz(l)));
             Gz = Gz.*  (1./(2*pi*hz(l).^2).*exp( (cos(2*pi*(z(l,:)-mz(l,:)')) - 1  )./hz(l)./hz(l)));
                          
                
    end

end



