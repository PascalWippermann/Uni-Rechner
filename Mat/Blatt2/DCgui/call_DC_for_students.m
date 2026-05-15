% 1D layered earth solution for DC resistivity
obj = DC;
% model layers: d thickness; sigma conductivity (Dimensionen 1 x Nl, Nl = Number of layers)
obj.d    = [Schichtmaechtigkeiten];
obj.sigma  = [Leitfaehigkeiten];
% electrode positions; can be arbitrary, just make sure that they are of
% dim. 2 x Nd (Nd = Number of electrode configurations)
obj.rA = [x-Koordinaten; y-Koordinaten];
obj.rB = [x-Koordinaten; y-Koordinaten];
obj.rM = [x-Koordinaten; y-Koordinaten];
obj.rN = [x-Koordinaten; y-Koordinaten];

% calculate apparent resistivity
obj.appres;

% some plotting routine...