
function R=cardatol(q,i,j,k)


%CARDATOL (Spacelib): Cardan angles to L matrix. 
%
% Builds L  matrix of a frame whose orientation is specified by an Euler or
% Cardanic convention.
% The parameters i, j, k  specify  the sequences of the rotation axes (their
% value must be the constant X,Y or Z). 
% j must be different from i and k, k could be equal to i.
% q:   3-element vector containing the 1st, 2nd and 3rd angle.
% Usage:
%
%                       A=cardatol(q,i,j,k)
%
% Related functions: CARDTOOM, CARDOMPT, HTOCARDA.
%
% (c) G.Legnani, C. Moiola 1998; adapted from: G.Legnani and R.Adamini 1993
%___________________________________________________________________________

% spheader   % Don't assign all these global variables
X = 1; Y = 2; Z = 3;

alfa= q(1);
beta= q(2);
gamma=q(3);

sa=sin(alfa);	sb=sin(beta);	sc=sin(gamma);
ca=cos(alfa);	cb=cos(beta);	cc=cos(gamma);

if ( i<X | i>Z | j<X | j>Z | k<X | k>Z | i==j | j==k )
	error('	**** 	Spacelib Error: Illegal rotation axis	****');
end

if ( rem(j-i+3,3)==1 )	sig=1;  % ciclic 
	else 	        sig=-1;	% anticiclic
end

R = zeros(3,3);
% Cardanic Convention 
if (i~=k)
	R(i,X)=1;
		R(i,Y)=0;
			R(i,Z)=sig*sb;

	R(j,X)= 0;
		R(j,Y)= ca;
			R(j,Z)= -sig*sa*cb;

	R(k,X)= 0;
		R(k,Y)= sig*sa;
			R(k,Z)= ca*cb;  

% Eulerian Convention
else
l=6-i-j;

R(i,X)= 1;
	R(i,Y)= 0;
		R(i,Z)= cb;

R(j,X)= 0;
	R(j,Y)= ca;
		R(j,Z)= sa*sb;

R(l,X)= 0;
	R(l,Y)= sig*sa;
		R(l,Z)= -sig*ca*sb;

end
