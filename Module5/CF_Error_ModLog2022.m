function [calcerror] = CF_Error_ModLog2022(x,t,p)

    [h,N,Time]  = CF_Sim_ModLog2022(x(1),x(2),x(3),x(4),t,false);
    
    [calcerror] = sum((p-h).^2);    

end