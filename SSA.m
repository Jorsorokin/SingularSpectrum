classdef SSA < handle
    % A matlab class for performing Singluar Spectrum Analysis (SSA) on time series data
    %
    % SSA decomposes a time series by first embedding the time series into a phase-space trajectory
    % via time lag (see phaseSpace.m), then performing PCA on the resulting matrix. The principal
    % components represent principal "modes" found within the time series, which one can use to
    % guide reconstructon using just a subset of those modes
    %
    % USAGE:
    %   X = rand( 1000,1 );
    %   L = 250; % # of time points to embed at at time
    %   percVar = 0.8; % keep 80 % of variance with decomposition
    %
    %   ssa = SSA( X )
    %   ssa.embed( L )
    %   ssa.decompose( percVar );
    %   R = ssa.reconstruct(); % <- can provide indices of the principal components to use
    %
    % Properties:
    %   X - the original time series vector or matrix (if matrix, should be column-major)
    %   embedding - the embedded trajectory
    %   U - left singular vectors (scores)
    %   S - singular values
    %   V - right singular vectors (principal modes)
    %   groupIDs - vector of IDs, with each unique ID defining a group of similar PCs
    %   groupFeatures - the features used for clustering the PCs
    %   varExp - % variance explained from the decomposition
    %
    % Methods:
    %   embed
    %   decompose
    %   reconstruct
    %   groupPCs
    %   plotReconstruction
    %
    % Copyright Jordan Sorokin (jorsorokin@gmail.com), 5/10/18
    % You are free to use/distribute this code, but please keep a referene to the original author
    % and code base   
    
    properties(SetAccess = protected)
        X
        embedding
        U
        S
        V
        groupIDs
        groupFeatures
        varExp
    end
    
    properties(Access = private)
        N
        nChan
        M
        L
        K
        totalVar
    end
         
    methods(Access=public)
        
        function self = SSA( X )
            % self = ssa( X )
            %
            % initiate the singular-spectrum analysis instance by providing
            % column-oriented data in "X"
            
            self.X = X;
            [self.N,self.nChan] = size( X );
        end
  
        
        function embed( self,L )
            % embed( self,L )
            %
            % embed the column-time series in self.X via L-point lagging
            
            % set the embedding params
            self.M = self.N - L + 1;
            self.L = L;
            
            % embed the data
            self.embedding = phaseSpace( self.X,self.M,1,'MOD' );
        end
        
        function decompose( self,K )
            % decompose( self,K )
            %
            % decompose the embedded data using SVD, and keep only the
            % first "K" components, or if K is < 1, the # of components
            % that explain K*100 % of the variance
            
            % check if we've embedded
            if isempty( self.embedding )
                error( 'Please embed the time series!' );
            end
            
            % reshape the embedding if more than 1 channel
            emb = self.embedding;
            if self.nChan > 1
                emb = reshape( permute( self.embedding,[1,3,2] ),self.L*self.nChan,self.M );
            end 
            
            % perform SVD
            [u,s,v] = svd( emb','econ' );
            
            % check if K < 1, and if so, find the # of components to
            % explain K-% variance
            s = diag( s );
            self.totalVar = sum( s );
            if K < 1
                K = find( cumsum( s ) / self.totalVar >= K,1 );
            end
            
            inds = 1:K;
            self.S = s(inds);
            self.U = u(:,inds);
            self.V = v(:,inds);
            self.K = K;
            self.varExp = sum( self.S ) / self.totalVar;
            
            % print to screen
            fprintf( '%02f variance explained using %i components\n',self.varExp,self.K );
        end
        
        function groupPCs( self,nGroups )
            % groupPCs( self,nGroups )
            %
            % group the eigen vectors/values into nGroups distinct groups
            % by using the eigenvalues and mean frequencies of the
            % corresponding PCs as features
            
            % check if we've decomposed
            if isempty( self.S )
                error( 'Please embed and decompose the time series first!' );
            end
            
            F = meanfreq( self.U );
            features = [self.S,F'];
            self.groupIDs = cluster( linkage( [self.S,F'],'ward','euclidean' ),'maxclust',nGroups );
            self.groupFeatures = features;
        end
        
        function R = reconstruct( self,idx )
            % R = reconstruct( self,idx )
            %
            % reconstruct the approximated time series by summing the
            % rank-1 matricies U_i*S_i*V_i' for all i specified by idx,
            % followed by diagonal averaging over the result
            
            % check for decomposition
            if isempty( self.S )
                error( 'Please embed and decompose the time series first!' );
            end
            
            % pre-allocate and sum over the indicies
            R = zeros( self.N,self.nChan );
            Y = self.U(:,idx)*diag( self.S(idx) )*self.V(:,idx)';
            startDiag = -self.L;
            
            for c = 1:self.nChan
                inds = (c-1)*self.L+1:c*self.L;
                temp = fliplr( Y(:,inds)' );
                for i = 1:self.N
                    R(i,c) = mean( diag( temp,startDiag+i ) );
                end
            end
            
            R = flipud( R );
        end
        
        function plotReconstruction( self,R )
            % plotReconstruction( self,R )
            %
            % plots the raw data and the reconstruction R on top
            
            figure;
            yval = mean( range( self.X ) );
            multisignalplot( self.X,[],'k',yval );
            multisignalplot( R,[],'r',yval );
            ylabel( 'Channels' );
        end           
    end
end
            
            
            
        
        
    
    