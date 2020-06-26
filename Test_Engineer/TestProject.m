classdef TestProject < handle
    properties
        testmotor    
    end
    methods(Abstract)
        Load_motor(obj);
        FEA_analysis(obj);
        Save_Data(obj);
    end
end