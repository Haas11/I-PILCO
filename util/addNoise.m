function [ u_out ] = addNoise(u_in, plant)
%addNoise() impose gaussian white noise upon signals passed to PILCO


if isfield(plant,'subplant')
    subi = plant.subi;
else
    plant.subplant = inline('[]',1);
    subi = []; 
end
odei = plant.odei;
idx = sort([odei subi]);

u_out = u_in;

u_out(idx) = u_in(idx) + randn(size(idx))*chol(plant.noise);

end

