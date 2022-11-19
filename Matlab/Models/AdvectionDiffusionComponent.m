classdef AdvectionDiffusionComponent < BaseModel

    properties
        
        D % Diffusion coefficient
        %v % Speed of flow

    end

    methods
        
        function model = AdvectionDiffusionComponent(paramobj)

            model = model@BaseModel();
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'D'};
                       %'v'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.operators = localSetupOperators(model.G);
            

        end

        function model = registerVarAndPropfuncNames(model)
        %% Declaration of the Dynamical Variables and Function of the model
        %  Note : this is more like a documentation and is not used in assembly

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {'c'       , ...
                        'accum'   , ...
                        'source'  , ...
                        'flux'    , ...
                        'massCons'};
            
            model = model.registerVarNames(varnames);

            fn = @AdvectionDiffusionComponent.updateFlux;
            inputnames = {'c'};
            model = model.registerPropFunction({'flux', fn, inputnames});

            fn = @AdvectionDiffusionComponent.updateMassConservation;
            inputnames = {'accum', 'flux', 'source'};
            model = model.registerPropFunction({'massCons', fn, inputnames});
            
            fn = @AdvectionDiffusionComponent.updateMassAccum;
            inputnames = {'c'};
            model = model.registerPropFunction({'accum', fn, inputnames});

        end

        
        function state = updateFlux(model, state)
        % Assemble electrical current which is stored in :code:`state.j`
            % This should not be recomputed at every step
            dim = model.G.cartDims;
            N = dim(1); % x-direction
            M = dim(2); % y-direction
            
            a = 1:N-1;
            c = transpose(ones(M,1));
            A = c'*a;
            d = 0:M-1;
            e = N.*d;
            A = A + e'*transpose(ones(N-1,1));
            A_ = transpose(A);
            index1 = reshape(A_,1,[]);
            %index2 = reshape(A_+1,1,[]);
            

            D  = model.D;
            v1 = -3e-7;
            v2 = 0;
            op = model.operators;
            T  = op.T;
            
            c   = state.c;

            %replace 4900 with numbers depending on actual grid size
            % Add advection term, find better method to do this later
            % For now assume 50 nodes in each direction fixed
            tot = (N-1)*M + (N)*(M-1);
            
            j_ad = zeros(tot,1);
            %j_ad(1:2450) = v1.*c.val(index1);
            %j_ad(1:2:4899) = 1/2.*v1.*(c.val(index1) + c.val(index1+1));
            j_ad(1:(N-1)*M) = 1/2.*v1.*(c.val(index1) + c.val(index1+1));
            % rest of j_ad = 0
            %j_ad.val(2451:4900) = v2.*c.val(index1);
            %j_ad(2:2:4900) = v2.*c.val(index1);
            %j_ad((tot)/2+1:tot) = 1/2.*v2.*(c.val(index1) + c.val(index1+1));
            scaled = 1;
            flux = - D.*T.*op.Grad(c) + scaled.*j_ad;% + [v1.*c, v2.*c];
            
            state.flux = flux;
            
        end

        function state = updateMassAccum(model, state, state0, dt)

            vols = model.G.cells.volumes;
            
            c  = state.c;
            c0 = state0.c;

            state.accum = vols.*(c - c0)/dt;
            
        end
        
        
        function state = updateMassConservation(model, state)
        % Assemble residual of the charge conservation equation which is stored in :code:`state.chargeCons`

            accum    = state.accum;
            flux     = state.flux;
            source   = state.source;
            
            massCons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.massCons = massCons;
            
        end
        
    end
    
end
