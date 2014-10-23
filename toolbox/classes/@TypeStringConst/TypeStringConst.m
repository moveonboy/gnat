classdef TypeStringConst
    %TYPESTRINGCONST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        enzRxn='Enzymatic Reaction';
        tranpRxn = 'Transport Event';
        nonenzRxn ='Nonenzymatic Reaction';        
    end
    
    properties (Constant)
        bibiSeqOrder = 'Sequential Order';
        bibiRanOrder = 'Random Order';
        bibiPingPong = 'Ping Pong';
        bibiSeqOrderwSubstInhib = 'Sequential Order with substrate Inhibition';
        mechanismdb = {TypeStringConst.bibiSeqOrder,TypeStringConst.bibiRanOrder,...
            TypeStringConst.bibiPingPong,TypeStringConst.bibiSeqOrderwSubstInhib}
    end
    
    properties (Constant)
        SMM      ='Michaelis-Menten';
        Noncomp  ='Noncompetitive';
        Uncomp   ='Uncompetitive';
        Comp     ='Competitive';
        Subst    ='Substrate Inhibition';
        Elem     ='elemental';
        Ensm     = 'ensemble';
        commechs = {TypeStringConst.SMM;TypeStringConst.Noncomp;...
            TypeStringConst.Uncomp;TypeStringConst.Comp;TypeStringConst.Subst};
        rulesbasedkinetics = {'Structure-based rule ';'Substructure-based rule'};
    end
    
    properties(Constant)
       CSTR          ='CSTR';
       PFR           ='PFR';
       BatchR        ='BatchR'; 
       Rectortypes   = {TypeStringConst.CSTR,TypeStringConst.PFR,TypeStringConst.BatchR};
    end
    
    
    
    
    
end

