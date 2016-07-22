tool_link = 'peg_link_ee_kuka';

T = 10;
dt = 0.10;

t_01 = T/5/dt;
t_12 = T/5/dt;
t_23 = T*2/5/dt;
t_33 = T/5/dt;

initRot = t2r(quat2tform([0 1 0 0]));
firstRot = t2r(quat2tform([0.1 0.9 0.1 0.1]));

T0 = rt2tr(initRot,[0.5 0 0.4]');
T1 = rt2tr(firstRot,[0.5 0 0.2]');
T2 = rt2tr(initRot,[0.6 0 0.1]');
T3 = rt2tr(initRot,[0.6 0.1 0.1]');

Tcart1 = ctraj(T0, T1, t_01);           % Cartesian trajectory generation
Tcart2 = ctraj(T1, T2, t_12);
Tcart3 = ctraj(T2, T3, t_23);
Tcart4 = ctraj(T3, T3, t_33);

Ttot = cat(3,Tcart1,Tcart2,Tcart3,Tcart4,Tcart4);

t = 0:dt:T-dt;

Href = Ttot;

Href1 = cat(3,T0,T1,T2,T3);