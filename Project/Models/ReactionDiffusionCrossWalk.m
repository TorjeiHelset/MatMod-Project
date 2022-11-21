classdef ReactionDiffusionCrossWalk < BaseModel

    properties
        
        N % Diffision model for component N
        R % Diffision model for component R
        C % Diffision model for result
        T % Diffusion model for component T
        Tb % Diffusion model for result Tb
        N_inac % Diffusion model for N_inac
        N2 % Diffision model for component N
        R2 % Diffision model for component R
        C2 % Diffision model for result
        T2 % Diffusion model for component T
        Tb2 % Diffusion model for result Tb
        N_inac2 % Diffusion model for N_inac

        k1 % reaction constant of N + R -> C
        k_1 % reaction constant of C -> N + R

        k_t1 % reaction constant of T + C -> Tb
        k_t_1 % reaction constant of Tb -> T + C
        k_t2  % reaction constant of Tb -> T + N_inac
        
    end
    
    methods
        
        function model = ReactionDiffusionCrossWalk(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k1', 'k_1', 'k_t1', 'k_t_1', 'k_t2'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.N = DiffusionComponent(paramobj.N);
            model.R = DiffusionComponent(paramobj.R);
            model.C = DiffusionComponent(paramobj.C);
            model.T = DiffusionComponent(paramobj.T);
            model.Tb = DiffusionComponent(paramobj.Tb);
            model.N_inac = DiffusionComponent(paramobj.N_inac);
            model.N2 = DiffusionComponent(paramobj.N2);
            model.R2 = DiffusionComponent(paramobj.R2);
            model.C2 = DiffusionComponent(paramobj.C2);
            model.T2 = DiffusionComponent(paramobj.T2);
            model.Tb2 = DiffusionComponent(paramobj.Tb2);
            model.N_inac2 = DiffusionComponent(paramobj.N_inac2);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion.updateSourceTerm;
            
            inputnames = {{'N', 'c'}, ...
                          {'R', 'c'}, ...
                          {'C', 'c'}, ...
                          {'T', 'c'}, ...
                          {'Tb', 'c'}, ...
                          {'N_inac', 'c'}, ...
                          {'N2', 'c'}, ...
                          {'R2', 'c'}, ...
                          {'C2', 'c'}, ...
                          {'T2', 'c'}, ...
                          {'Tb2', 'c'}, ...
                          {'N_inac2', 'c'}};

            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'R', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'C', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'T', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'Tb', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N_inac', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N2', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'R2', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'C2', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'T2', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'Tb2', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N_inac2', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model);
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm(model, state)

            k1 = model.k1;
            k_1 = model.k_1;
            k_t1 = model.k_t1;
            k_t_1 = model.k_t_1;
            k_t2 = model.k_t2;

            vols = model.G.cells.volumes;
            nc = model.G.cells.num;
            
            
            cN = state.N.c;
            cR = state.R.c;
            cC = state.C.c;
            cT = state.T.c;
            cTb = state.Tb.c;
            cN2 = state.N2.c;
            cR2 = state.R2.c;
            cC2 = state.C2.c;
            cT2 = state.T2.c;
            cTb2 = state.Tb2.c;

            %Calculate diffusion through hole:
            hole1 = [1 + 23*50, 1 + 24*50, 1+25*50, 1+26*50];
            hole2 = [(24*50), (25*50), (26*50), (27*50)]-1;
            % Diffusion given by the Laplacian of concentration,
            % approximate this by central difference
            % What is delta x?
            dx = 1/50;
            fluxOut = cN;

            fluxIn = cN;
            fluxOut.val(1:nc) = zeros(nc,1);
            fluxOut.val(hole1) = (cN2.val(hole2) - 2.*cN.val(hole1) + cN.val(hole1+1)) / (dx^2);
            fluxIn.val(1:nc) = zeros(nc,1);
            fluxIn.val(hole2) = (cN2.val(hole2-1) - 2.*cN2.val(hole2) + cN.val(hole1)) / (dx^2);

            R1 = k1.*vols.*cN.*cR;
            R_1 = k_1.*vols.*cC;
            R_t1 = k_t1.*vols.*cT.*cN;
            R_t_1 = k_t_1.*vols.*cTb;
            R_t2 = k_t2.*vols.*cTb;

            R1_2 = k1.*vols.*cN2.*cR2;
            R_1_2 = k_1.*vols.*cC2;
            R_t1_2 = k_t1.*vols.*cT2.*cN2;
            R_t_1_2 = k_t_1.*vols.*cTb2;
            R_t2_2 = k_t2.*vols.*cTb2;


            state.N.source = -R1 + R_1 -R_t1 + R_t_1 + fluxOut;
            state.R.source = -R1 + R_1;
            state.C.source = R1 - R_1;
            state.T.source = -R_t1 + R_t_1 + R_t2;
            state.Tb.source = R_t1 - R_t_1 - R_t2;
            state.N_inac.source = R_t2;
            
            state.N2.source = -R1_2 + R_1_2 -R_t1_2 + R_t_1_2+ fluxIn;
            state.R2.source = -R1_2 + R_1_2;
            state.C2.source = R1_2 - R_1_2;
            state.T2.source = -R_t1_2 + R_t_1_2 + R_t2_2;
            state.Tb2.source = R_t1_2 - R_t_1_2 - R_t2_2;
            state.N_inac2.source = R_t2_2;
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.N = model.N.updateFlux(state.N);
            state.R = model.R.updateFlux(state.R);
            state.C = model.C.updateFlux(state.C);
            state.T = model.T.updateFlux(state.T);
            state.Tb = model.Tb.updateFlux(state.Tb);
            state.N_inac = model.N_inac.updateFlux(state.N_inac);
            state.N2 = model.N2.updateFlux(state.N2);
            state.R2 = model.R2.updateFlux(state.R2);
            state.C2 = model.C2.updateFlux(state.C2);
            state.T2 = model.T2.updateFlux(state.T2);
            state.Tb2 = model.Tb2.updateFlux(state.Tb2);
            state.N_inac2 = model.N_inac2.updateFlux(state.N_inac2);

            state = model.updateSourceTerm(state);

            state.N = model.N.updateMassAccum(state.N, state0.N, dt);
            state.R = model.R.updateMassAccum(state.R, state0.R, dt);
            state.C = model.C.updateMassAccum(state.C, state0.C, dt);
            state.T = model.T.updateMassAccum(state.T, state0.T, dt);
            state.Tb = model.Tb.updateMassAccum(state.Tb, state0.Tb, dt);
            state.N_inac = model.N_inac.updateMassAccum(state.N_inac, state0.N_inac, dt);
            state.N2 = model.N2.updateMassAccum(state.N2, state0.N2, dt);
            state.R2 = model.R2.updateMassAccum(state.R2, state0.R2, dt);
            state.C2 = model.C2.updateMassAccum(state.C2, state0.C2, dt);
            state.T2 = model.T2.updateMassAccum(state.T2, state0.T2, dt);
            state.Tb2 = model.Tb2.updateMassAccum(state.Tb2, state0.Tb2, dt);
            state.N_inac2 = model.N_inac2.updateMassAccum(state.N_inac2, state0.N_inac2, dt);

            state.N = model.N.updateMassConservation(state.N);
            state.R = model.R.updateMassConservation(state.R);
            state.C = model.C.updateMassConservation(state.C);
            state.T = model.T.updateMassConservation(state.T);
            state.Tb = model.Tb.updateMassConservation(state.Tb);
            state.N_inac = model.N_inac.updateMassConservation(state.N_inac);
            state.N2 = model.N2.updateMassConservation(state.N2);
            state.R2 = model.R2.updateMassConservation(state.R2);
            state.C2 = model.C2.updateMassConservation(state.C2);
            state.T2 = model.T2.updateMassConservation(state.T2);
            state.Tb2 = model.Tb2.updateMassConservation(state.Tb2);
            state.N_inac2 = model.N_inac2.updateMassConservation(state.N_inac2);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.N.massCons;
            names{end + 1} = 'massCons N';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.R.massCons;
            names{end + 1} = 'massCons R';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.C.massCons;
            names{end + 1} = 'massCons C';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.T.massCons;
            names{end + 1} = 'massCons T';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.Tb.massCons;
            names{end + 1} = 'massCons Tb';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.N_inac.massCons;
            names{end + 1} = 'massCons N_inac';
            types{end + 1} = 'cell';
                        
            primaryVars = model.getPrimaryVariables();

            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
        function state = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            newstate = initStateAD@BaseModel(model, state);
            newstate.time = state.time;
            state = newstate;
            
        end 

        
        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = {{'N', 'c'}, ...
                               {'R', 'c'}, ...
                               {'C', 'c'}, ...
                               {'T', 'c'}, ...
                               {'Tb', 'c'}, ...
                               {'N_inac', 'c'}, ...
                               {'N2', 'c'}, ...
                               {'R2', 'c'}, ...
                               {'C2', 'c'}, ...
                               {'T2', 'c'}, ...
                               {'Tb2', 'c'}, ...
                               {'N_inac2', 'c'}};
            
        end
        

        function model = validateModel(model, varargin)
        % nothing special to do
        end


    end

end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

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
