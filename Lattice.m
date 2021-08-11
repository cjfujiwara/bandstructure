classdef Lattice
    %LATTICE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_states
        lambda
        k_lambda
        h
        hbar
        m
        ErJ
        ErHz
        depth     
    end
    
    methods
        function obj = Lattice(depth)
            %LATTICE Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin==0
                depth=10;
            end
            
            obj.depth = depth; 
            obj.n_states = 21;                                  % Number of bloch states
            obj.lambda = 1053.6E-9;                             % Lattice light wavelength [m]
            obj.k_lambda =(2*pi)/obj.lambda;                   % Lattice light wavevector
            obj.h = 6.626E-34;
            obj.hbar = obj.h/(2*pi);                            % reduced planck's constant [J s]
            obj.m = 39.96399817 * 1.66053907E-27;               % 40K mass [kg]
            obj.ErJ = obj.hbar^2*obj.k_lambda^2/(2*obj.m);       % Recoil energy [J]    
            obj.ErHz = obj.ErJ/obj.h;            
            obj.depth = 10;                                     % Depth [Er]            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

