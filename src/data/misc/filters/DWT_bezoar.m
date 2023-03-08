classdef DWT_bezoar < TFilter
    properties (Constant)
        name                = 'DWT'
        type                = 'single'
        col_wise            = true
        parameterNames      = ["Lower Bound", "Upper Bound"]
        parameterUnits      = ["cpm", "cpm"]
        parameterDefaults   = [0.02, 12]
    end
    properties (Access = private)
        elim_freq              double
        wavelet                 string
        levels                  double
        f                       double
    end
    
    methods
        function data = filter(obj, data)
            
            %dwt of data in rows.
            dec = mdwtdec('r',data',obj.levels,obj.wavelet);
            %zero out the coefficients according to the levels (pseudofreqs) we don't
            %want to keep
            DEC2 = chgwdeccfs(dec,'cd',0,obj.elim_freq);
            %reconstruct the waveforms based on filtered dwt coeffs
            data = mdwtrec(DEC2)';        
        end
        function p = validateProperties(obj, p)

            obj.wavelet = 'db4';
            
            obj.levels = round(log2(p.samples));
            
            %figure out which scales (pseudofreqs we want to block)
            % scales are on dyadic grid (in powers of 2)
            %large level --> large scale --> low pseudofreq
            bandpass = [obj.parameterValues(1), obj.parameterValues(2)];
            f_low = bandpass(1)./60; f_high = bandpass(2)./60;  % convert cpm to Hz
            
            pseudofreqs = scal2frq(2.^[1:obj.levels], obj.wavelet, 1/p.frequency);
            pseudofreqs = fliplr(pseudofreqs); %flip the list so it is easier to work with
            reverselvl = fliplr(1:obj.levels); %reversed list of levels
            if f_low < 0
                f_low = min(pseudofreqs);
            end
            
            if f_high < 0
                f_high = max(pseudofreqs);
            end
            
            [~, I_low] = min(abs(pseudofreqs - f_low));
            [~, I_high] = min(abs(pseudofreqs - f_high));
            f_low_val = pseudofreqs(I_low);
            f_high_val = pseudofreqs(I_high);
                       
            % Set filter properties
            obj.parameterValues(1) = f_low_val.*60; % convert Hz to cpm
            obj.parameterValues(2) = f_high_val.*60; % convert Hz to cpm
            
            obj.elim_freq = [reverselvl(1:I_low-1), reverselvl(I_high+1:end)];
            obj.f = p.frequency;
            
            obj.description = [ ...
                '<html>' ...
                '<b>Discrete Wavelet Denoising</b><br/>' ...
                '<i>MATLAB Wavelet Toolbox</i><br/>' ...
                '&emsp ' ...
                obj.parameterValues(1) ' cpm lower bound <br/>' ...
                '&emsp ' ...
                obj.parameterValues(2) ' cpm upper bound <br/>' ...
                '&emsp ' ...
                obj.wavelet 'wavelet' ...
                '</html>'];
        end
    end
end