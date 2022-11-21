classdef ReactionDiffusionInputParamsCrossWalk < InputParams

    properties
        
        G
        
        N % parameter for diffusion component model N
        R % parameter for diffusion component model R
        C % parameter for diffusion component model C
        T % Diffusion model for component T
        Tb % Diffusion model for result Tb
        N_inac % Diffusion model for N_inac
        % Second synaptic cleft:
        N2 % parameter for diffusion component model N
        R2 % parameter for diffusion component model R
        C2 % parameter for diffusion component model C
        T2 % Diffusion model for component T
        Tb2 % Diffusion model for result Tb
        N_inac2 % Diffusion model for N_inac
        
        k1 % reaction rate of R + N -> C
        k_1 % reaction rate of C -> R + N

        k_t1 % reaction constant of T + C -> Tb
        k_t_1 % reaction constant of Tb -> T + C
        k_t2  % reaction constant of Tb -> T + N_inac
    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParamsCrossWalk(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.N = DiffusionComponentInputParams(pick('N'));
            paramobj.R = DiffusionComponentInputParams(pick('R'));
            paramobj.C = DiffusionComponentInputParams(pick('C'));
            paramobj.T = DiffusionComponentInputParams(pick('T'));
            paramobj.Tb = DiffusionComponentInputParams(pick('Tb'));
            paramobj.N_inac = DiffusionComponentInputParams(pick('N_inac'));

            paramobj.N2 = DiffusionComponentInputParams(pick('N2'));
            paramobj.R2 = DiffusionComponentInputParams(pick('R2'));
            paramobj.C2 = DiffusionComponentInputParams(pick('C2'));
            paramobj.T2 = DiffusionComponentInputParams(pick('T2'));
            paramobj.Tb2 = DiffusionComponentInputParams(pick('Tb2'));
            paramobj.N_inac2 = DiffusionComponentInputParams(pick('N_inac2'));

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            if ~isempty(paramobj.G)
                paramobj.N.G = paramobj.G;
                paramobj.R.G = paramobj.G;
                paramobj.C.G = paramobj.G;
                paramobj.T.G = paramobj.G;
                paramobj.Tb.G = paramobj.G;
                paramobj.N_inac.G = paramobj.G;

                paramobj.N2.G = paramobj.G;
                paramobj.R2.G = paramobj.G;
                paramobj.C2.G = paramobj.G;
                paramobj.T2.G = paramobj.G;
                paramobj.Tb2.G = paramobj.G;
                paramobj.N_inac2.G = paramobj.G;
            end
            
        end
        
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Cattery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
