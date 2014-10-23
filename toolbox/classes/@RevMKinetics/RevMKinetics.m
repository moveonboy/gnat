classdef RevMKinetics < EnzKinetics
    %RevMKinetics Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        eleconsts=struct('kf1',[],'kr1',[],'kf2',[],'kr2',[],'kf3',[],'kr3',[],'enzconc',[]);        
        sumconsts=struct('vf',[],'vr',[],'ks',[],'kp',[]);
        form='';
    end
    
    methods
        function obj=RevMKinetics(form,rateconsts)
            if(strcmpi(form,'elemental'))
                obj.eleconsts.kf1=rateconsts.kf1;
                obj.eleconsts.kr1=rateconsts.kr1;
                obj.eleconsts.kf2=rateconsts.kf2;
                obj.eleconsts.kr2=rateconsts.kr2;
                obj.eleconsts.kf3=rateconsts.kf3;
                obj.eleconsts.kr3=rateconsts.kr3;
                obj.eleconsts.enzconc=rateconsts.enzconc; 
            elseif(strcmpi(form,'ensemble'))
                obj.sumconsts.vf=rateconsts.vf;
                obj.sumconsts.vr=rateconsts.vr;
                obj.sumconsts.ks=rateconsts.ks;
                obj.sumconsts.kp=rateconsts.kp;
            else
                error('MATLAB:GNAT:NOTSUPPORTEDFORM','This form is not supported')
            end
        end
    end
    
    methods
        function vel=simVel(obj,subsconc,prodconc)
            vel= (obj.sumconsts.vf*subsconc/obj.sumconsts.ks-...
                obj.sumconsts.vr*prodconc/obj.sumconsts.kp)/...
                (1+subsconc/obj.sumconsts.ks+prodconc/obj.sumconsts.kp);
        end
    end  
    
    methods
        function ele2ems(obj)
            obj.sumconsts.vf=obj.eleconsts.kf2*obj.eleconsts.enzconc;
            obj.sumconsts.vr=obj.eleconsts.kr2*obj.eleconsts.enzconc;
            obj.sumconsts.ks=obj.eleconsts.kr1/obj.eleconsts.kf1;
            obj.sumconsts.kp=obj.eleconsts.kf3/obj.eleconsts.kr3;
        end
    end
end