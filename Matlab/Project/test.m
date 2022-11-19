function G = test(r_max , r_min ,step , z_max, dens)
% Make a cylindrical triangular grid
    arguments
        r_max(1,1) double  = log(0.22*micro*meter);
        r_min(1,1) double  = log(0.1*nano*meter);
        step(1,1) double   = 0.25;
        z_max(1,1) double = convertFrom(ones(1,1)*(1.5), meter); %Does not work for nano*meter
        dens(1,1) double   = 20;
    end
    P = [];
    for r = exp(r_min:step:r_max)
       [x,y,~] = cylinder(r,dens); P = [P [x(1,:); y(1,:)]];          %#ok<AGROW>
    end

    P = unique(P', "rows");

    G = triangleGrid(P);
end