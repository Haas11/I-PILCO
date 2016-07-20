
function omega=cardtoom(q,qp,i,j,k)

%CARDTOOM (Spacelib): Cardan angles to angular velocity.
%
% Evaluates  the  angular  velocity  of a moving frame from three Cardan (or
% Eulerian) angles q and their time derivatives qp.
% The parameters i, j, k  specify  the sequences of the rotation axes (their
% value must be the constant X,Y or Z). 
% j must be different from i and k, k could be equal to i.
% q:     3-element row vector containing the 1st, 2nd and 3rd angle.
% qp :   3-element row vector containing the time derivative of q.
% omega: 3-element vector containing the angular velocity.
% Usage:
%			omega= cardtoom(q,qp,i,j,k);
%
% Related functions: CARDTOME, CARDOMPT 
%
% (c) G.Legnani, C. Moiola 1998; adapted from: G.Legnani and R.Adamini 1993
%___________________________________________________________________________


mat=cardatol(q,i,j,k);  

% NOTE: qp must be a row vector: qp=[ vx vy vz]

omega=mat*qp';

