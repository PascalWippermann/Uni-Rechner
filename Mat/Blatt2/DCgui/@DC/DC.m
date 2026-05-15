% Calculates DC responses on the surface (z=0) of  a layered earth
classdef DC < handle
    properties
        i = 1;
        rA
        rB
        rM
        rN
        sigma   = [];
        d       = [];
    end
    properties (Dependent = false,  SetAccess = private)
        kappa  % wavenumbers
        rh     % distances; may need interpolation to original distances r
        r      % range of distances for which a solution is sought
        rAM
        rAN
        rBM
        rBN
        k      % konfigurationsfaktor
        B
        u      % potential [uMA uMB uNA uNB; ...]
        U      % 4-point potential
        rhoa   % apparent resistivity
    end
    properties (Dependent = true,  SetAccess = private)
        
        L
        l
        M       % number of layers
    end
    methods
        % methods are functions which can operate on objects
        % first method is the constructor method, which allows to set the
        % important properties
        function obj    = DC(varargin)
            
        end
        
        % ---------------------------------
        % SET and GET methods
        % ---------------------------------
        function M = get.M(obj)
            M = numel(obj.sigma);
        end
        function rAM = get.rAM(obj)
            rAM = sqrt((obj.rA(1,:)-obj.rM(1,:)).^2+(obj.rA(2,:)-obj.rM(2,:)).^2);
        end
        function rAN = get.rAN(obj)
            rAN = sqrt((obj.rA(1,:)-obj.rN(1,:)).^2+(obj.rA(2,:)-obj.rN(2,:)).^2);
        end
        function rBM = get.rBM(obj)
            rBM = sqrt((obj.rB(1,:)-obj.rM(1,:)).^2+(obj.rB(2,:)-obj.rM(2,:)).^2);
        end
        function rBN = get.rBN(obj)
            rBN = sqrt((obj.rB(1,:)-obj.rN(1,:)).^2+(obj.rB(2,:)-obj.rN(2,:)).^2);
        end
        function k = get.k(obj)
            k = 2*pi./(1./obj.rAM-1./obj.rBM-1./obj.rAN+1./obj.rBN);
        end
        function [u] = potential(obj)
            for icfg = 1:size(obj.rA,2)
                % Potential at electrode M due to point source at A
                obj.r = obj.rAM(icfg);
                obj.getkappa('J0');
                rh      = obj.rh;
                kappa   = obj.kappa;
                obj.getB;
                A       = obj.kappa./obj.B(1,:);
                C       = A-1;
                uMA = obj.i./(obj.sigma(1)*2*pi).*(1./obj.r+fht_old('J0','fht',C,kappa,rh));
                % Potential at electrode M due to point source at A
                obj.r = obj.rBM(icfg);
                obj.getkappa('J0');
                rh      = obj.rh;
                kappa   = obj.kappa;
                obj.getB;
                A       = obj.kappa./obj.B(1,:);
                C       = A-1;
                uMB = obj.i./(obj.sigma(1)*2*pi).*(1./obj.r+fht_old('J0','fht',C,kappa,rh));
                % Potential at electrode M due to point source at A
                obj.r = obj.rAN(icfg);
                obj.getkappa('J0');
                rh      = obj.rh;
                kappa   = obj.kappa;
                obj.getB;
                A       = obj.kappa./obj.B(1,:);
                C       = A-1;
                uNA = obj.i./(obj.sigma(1)*2*pi).*(1./obj.r+fht_old('J0','fht',C,kappa,rh));
                % Potential at electrode M due to point source at A
                obj.r = obj.rBN(icfg);
                obj.getkappa('J0');
                rh      = obj.rh;
                kappa   = obj.kappa;
                obj.getB;
                A       = obj.kappa./obj.B(1,:);
                C       = A-1;
                uNB = obj.i./(obj.sigma(1)*2*pi).*(1./obj.r+fht_old('J0','fht',C,kappa,rh));
                u(icfg,:) = [uMA -uMB uNA -uNB];
                
            end
            obj.u = u;
        end
        % UTM coordinates
        function getkappa(obj,J)
            switch J
                case 'J0'
                    [~,rh,kappa]  =   fht_old('J0','NULL',[],[],sort(obj.r));
                case 'J1'
                    [~,rh,kappa]  =   fht_old('J0','NULL',[],[],sort(obj.r));
            end
            obj.kappa = kappa;
            obj.rh    = rh;
        end
        % recursion kernel
        function getB(obj)
            kappa = obj.kappa;
            B       = zeros(obj.M,numel(kappa));
            beta    = zeros(obj.M-1,numel(kappa));
            for m = 1:obj.M-1
                beta(m,:) = (obj.sigma(m))./(obj.sigma(m+1));
            end
            B(obj.M,:) = kappa;
            for m = obj.M-1:-1:1
                B(m,:) = kappa.*(B(m+1,:)+kappa.*beta(m,:).*tanh_z(kappa*obj.d(m)))./...
                    (kappa.*beta(m,:)+B(m+1,:).*tanh_z(kappa*obj.d(m)));
            end
            obj.B = B;
        end
        % four point  potential
        function getU(obj)
            u     = obj.potential;
            obj.U   = [u(:,1)+u(:,2)-(u(:,3)+u(:,4))];
        end
        % app. res.
        function appres(obj)
            obj.getU;
            R       = obj.U./obj.i;
            % apparent resistivity
            obj.rhoa     = obj.k'.*R;
        end
            
    end
    
end
