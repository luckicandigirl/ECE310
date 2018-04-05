% Daniel Nakhimovich and Sara Huang
function [ out ] = srconvertPolyPhase( in )
    UP = 320;
    DOWN = 147;
    
    % Create staged filters and store the coefficients
    % of their polyphase decomposition in a cell array
    B = {};
    filts = unique([factor(UP)]);
    for k = filts
        filt = designfilt('lowpassfir', 'PassbandFrequency', 1/k, 'StopbandFrequency', 1.2/k, 'PassbandRipple', 0.03, 'StopbandAttenuation', 85, 'DesignMethod', 'equiripple');
        B(k) = {poly1(filt.Coefficients,k)};
    end
    
    out = in;
    
    tic
    % Loop through stages determined by factoring
    % the overall upsample rate we want to achieve
    for L = factor(UP)
        sig = [];
        % Loop through polyphased components
        for k = 1:L
            s = conv(B{L}(k,:),out);
            s = [zeros(k-1,1);upsample(s,L)]; % upsample and delay
            sig = [sig;zeros(length(s)-length(sig),1)] + s; 
        end
        out = sig;
    end
    
    out = downsample(out,DOWN);

    fprintf('srconvertPolyPhase time: ')
    toc
end
