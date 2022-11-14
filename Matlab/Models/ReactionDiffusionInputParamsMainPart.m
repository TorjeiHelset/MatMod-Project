classdef ReactionDiffusionInputParamsMainPart < InputParams

    properties
        
        G
        
        N % parameter for diffusion component model N
        R % parameter for diffusion component model R
        C % parameter for diffusion component model C
        
        k1 % reaction rate of R + N -> C
        k_1 % reaction rate of C -> R + N

    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.N = DiffusionComponentInputParams(pick('N'));
            paramobj.R = DiffusionComponentInputParams(pick('R'));
            paramobj.C = DiffusionComponentInputParams(pick('C'));

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            if ~isempty(paramobj.G)
                paramobj.N.G = paramobj.G;
                paramobj.R.G = paramobj.G;
                paramobj.C.G = paramobj.G;
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
