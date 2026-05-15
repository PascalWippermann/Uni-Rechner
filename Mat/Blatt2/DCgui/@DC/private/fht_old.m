%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fast hankel transform of type 'type'
%   
%   type:
%   'J-1/2' :   f(r)    =   int(0->infty){F(k)cos(kr)dk} 
%   'J0'    :   f(r)    =   int(0->infty){F(k)J0(kr)dk}
%   'J+1/2' :   f(r)    =   int(0->infty){F(k)sin(kr)dk} 
%   'J1'    :   f(r)    =   int(0->infty){F(k)J1(kr)dk}
%   'J1qu'  :   f(r)    =   int(0->infty){F(k)J1^2(kr)dk}
%   'J2'    :   f(r)    =   int(0->infty){F(k)J2(kr)dk}     -> not implemented!!
%
%   what    :   string 'Null' returns only coefficients
%               < usage: [bla,rfil,kfil]  =   fht(type,'Null',[],[],r) >
%               string 'fht' performs Hankel transform F
%               < usage: [f,rfil,kfil]    =   fht(type,'fht',F,kfil,r) >
%               Hankel transform is then performed on rows of input matrix F(wNxkN)
%   F       :   input data maybe a matrix
%   kfil    :   wavenumber row vector   F   =   F(kfil) [k1, ... kN]
%   r       :   distance row vector [r1, ... rN]
%               if kfil is not correct with respect to r, function will return an 
%               error message and break 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MB2003

function [f,rfil,kfil]  =   fht_old(type,what,F,k,r)

nd          =   log10(r(end))-log10(r(1));
[H, Spec]   =   Hcf_old(type,nd);

if isempty(H)
    error('Zu grosser r-bereich => Kein passender Filter !!!');
end

nrfil       =   floor(1.9+Spec.ndec*nd);
nH          =   length(H);
D           =   log(10)/Spec.ndec;
rfil        =   exp((0:nrfil-1)*D)*r(1);
kfil        =   exp((nH-Spec.mid:-1:-nrfil-Spec.mid+2)*D)/min(rfil);

if isequal(what,'NULL')
    f       =   [];
    
elseif isequal(what,'fht')

    switch type
        case 'J-1/2'
        
        if  k(1,:)   ==  kfil;
        else
            error('spectrum not identical')
            return;
        end
        F       =   sqrt(repmat(kfil,size(F,1),1)).*F;
        for ri  =   1:length(rfil)
            f(:,ri) = sqrt(pi./(2*rfil(ri)))*(F(:,ri:ri+nH-1))*H(end:-1:1)';
        end
    case 'J0'
        if  (k(1,:))   ==  kfil;
        else
            error('spectrum not identical')
            return;
        end
        for ri  =   1:length(rfil)
            f(:,ri) = 1./rfil(ri)*F(:,ri:ri+nH-1)*H(end:-1:1)';
        end
        
    case 'J+1/2'
        
        if  k(1,:)   ==  kfil;
        else
            error('spectrum not identical')
            return;
        end
        F       =   sqrt(repmat(kfil,size(F,1),1)).*F;
        for ri  =   1:length(rfil)
            f(:,ri) = sqrt(pi./(2*rfil(ri)))*(F(:,ri:ri+nH-1))*H(end:-1:1)';
        end
   
    case 'J1'
        if  k(1,:)   ==  kfil;
        else
            error('spectrum not identical')
            return;
        end
        for ri  =   1:length(rfil)
            f(:,ri) = 1./rfil(ri)*F(:,ri:ri+nH-1)*H(end:-1:1)';
        end
        
    case 'J1qu'
        if  k(1,:)   ==  kfil;
        else
            error('spectrum not identical')
            return;
        end
        for ri  =   1:length(rfil)
            f(:,ri) = 1./rfil(ri)*F(:,ri:ri+nH-1)*H(end:-1:1)';
        end
        
    case 'J2'       %   not implemented
%         if  k(1,:)   ==  kfil;
%         else
%             error('spectrum not identical')
%             return;
%         end
%         for ri  =   1:length(rfil)
%             f(:,ri) = 1./rfil(ri)*F(:,ri:ri+nH-1)*H(end:-1:1)';
%         end
        
    end
end
