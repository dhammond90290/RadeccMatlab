function [ mb ] = LLSQ( Data,x1,x2 )
%This function will fit a linear regression y = mx + b to a data set in the
%array Data.  This array must have columns of x value, y value, weighting
%of 0 or 1. Rows to consider in Data must be x1 and x2. It returns a coefficient array of mb
%mb outputvformatting=[slope slope+/-; int int+/-]

%  Define for linear regression y = mx+b

%n = # data pts used
%Sx = sum of all x
%Sy = sum of all y
%Sxx = sum of all x^2
%Sxy = sum of all x*y
%Syy = sum of all y^2

%Eliminate points with zero wt and put others in AA
AA=[ ];
i2=0;

for i=x1:x2
    if Data(i,3)==1
    i2=i2+1;
    AA(i2,1)=Data(i,1);
    AA(i2,2)=Data(i,2);
    end 
end

n=0;
Sx = 0;
Sy = 0;
Sxx = 0;
Sxy = 0;
Syy = 0;

for i=1:length(AA(:,1))
n=n+1;
Sx = Sx+AA(i,1);
Sy = Sy +AA(i,2);
Sxx = Sxx+AA(i,1)*AA(i,1);
Sxy = Sxy+AA(i,1)*AA(i,2);
Syy = Syy+AA(i,2)*AA(i,2);
end

m = (n*Sxy - Sx*Sy)/(n*Sxx - Sx*Sx);
b = Sy/n - m*Sx/n;
A= Syy-Sy*Sy/n;
B = Sxx - Sx*Sx/n;
sigm = sqrt( ( (Syy-Sy*Sy/n) - (m^2)*(Sxx - Sx*Sx/n) ) / ( (n-2)*(Sxx - Sx*Sx/n) ) );
sigb = sqrt( ( (A - m*m*B)*(Sxx) ) / ( (n-2)*(n*Sxx - Sx*Sx) )  );

mb(1,1)=m;
mb(1,2)=sigm;
mb(2,1)=b;
mb(2,2)=sigb;




end

