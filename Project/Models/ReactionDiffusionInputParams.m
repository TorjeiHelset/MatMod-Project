classdef ReactionDiffusionInputParams < InputParams

    properties
        
        G
        
        A % parameter for diffusion component model A
        B % parameter for diffusion component model B
        C % parameter for diffusion component model C
        
        k % reaction rate

    end
    
    methods
        
        function paramobj = ReactionDiffusionInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.A = DiffusionComponentInputParams(pick('A'));
            paramobj.B = DiffusionComponentInputParams(pick('B'));
            paramobj.C = DiffusionComponentInputParams(pick('C'));

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            if ~isempty(paramobj.G)
                paramobj.A.G = paramobj.G;
                paramobj.B.G = paramobj.G;
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
