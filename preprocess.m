function [dsm2, trees]=preprocess(dsm, wSmalls, thrSmalls)

dsm2 = imopen( dsm, strel('disk', wSmalls) );

trees = abs(dsm-dsm2) > thrSmalls;
