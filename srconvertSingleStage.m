% Daniel Nakhimovich and Sara Huang
function [ out ] = srconvertSingleStage( in )
    UP = 320;
    DOWN = 147;

    [B A] = ellip(5,0.1,70,1/UP);
    
    tic
    
    out = upsample(in,UP);
    
    out = filter(B,A,out);
    
    out = downsample(out,DOWN);
    
    fprintf('srconvertSingleStage time: ')
    toc
end
