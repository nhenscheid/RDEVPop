function n = ADCMAP(gADC,mask)

% Computes an estimate of tumor cell density given an ADC map gADC
% Refs: Yankeelov 2013 "Clinically relevant modeling of tumor growth" 
%       Hormuth 2015 "Predicting in vivo glioma growth with the reaction
%       diffusion equation constrained by quantitative magnetic resonance
%       imaging data" 



kappa = 1; 

ADC_water = 2.5e-3;    % ADC of water (mm^2s^-1)


g_mask  = gADC(mask); 

ADC_min = min(g_mask);  
n = zeros(size(gADC)); 


n(mask) = kappa*(ADC_water - g_mask)/(ADC_water-ADC_min); 