function [ S ] = Speedup( k, fp, fo )

S = 1./(1 - fp + (fp./k) + (k - 1).*fo);

end