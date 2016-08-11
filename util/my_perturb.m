function  r2 = my_perturb(r, p)

	if nargin == 1
		p = 0.1;	% 10 percent disturb by default
	end

    r2 = SerialLink(r);

    links = r2.links;
	for i=1:r.n
% 		s = (2*rand-1)*p + 1;
		links(i).m = links(i).m * (p+1);

% 		s = (2*rand-1)*p + 1;
		links(i).I = links(i).I * (p+1);
	end

	r2.name = ['P/' r.name];