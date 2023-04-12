function [i,j]=vec2squ(rows,cols,x)
numvec=1:(rows*cols);
nummat=reshape(numvec,rows,cols);
[i,j]=find(nummat==x);
end