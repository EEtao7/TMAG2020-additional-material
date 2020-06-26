classdef MotorTestEngineer  < handle
    % MotorTestEngineer who can use femm to test some effects of motor.
    
    properties
        testproject     % define the builder which will be constructed by MotorDesignEngineer.
    end
    
    methods
        function test(obj)
            % define the sequence of construct a builder
            obj.testproject.Load_motor();
            obj.testproject.FEA_analysis();
            obj.testproject.Save_Data();
        end
    end  
    
end