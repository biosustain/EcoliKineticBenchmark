function [ T, Y, FLUX ] = runSimulation(ode_file,flux_file,span,y0,options)

[T, Y] = ode15s(ode_file,span,y0,options);
tmp    = flux_file(T(1),Y(1,:));
n_FLUX = length(tmp);
n_time = length(T);
FLUX   = zeros(n_time,n_FLUX);

for i = 1:n_time
    FLUX(i,:) = flux_file(T(i),Y(i,:));
end
