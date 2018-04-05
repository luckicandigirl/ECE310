% Daniel Nakhimovich and Sara Huang
function [ out ] = srconvertMultiRate( in )
    UP = 320;
    DOWN = 147;
    
    % Create staged filters and store their coefficients in a cell array
    B = {};
    filts = unique([factor(UP)]);
    for k = filts
        filt = designfilt('lowpassfir', 'PassbandFrequency', 1/k, 'StopbandFrequency', 1.2/k, 'PassbandRipple', 0.03, 'StopbandAttenuation', 85, 'DesignMethod', 'equiripple');
        B(k) = {filt.Coefficients};
    end
    
    out = in;
    
    tic
    % loop through stages determined by factoring
    % the overall upsample rate we want to achieve
    for L = factor(UP)
        out = upsample(out,L);
        out = conv(B{L},out);
    end
    
    out = downsample(out,DOWN);

    fprintf('srconvertMultiRate time: ')
    toc    
end
