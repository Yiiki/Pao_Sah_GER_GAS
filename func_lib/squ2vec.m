function x=squ2vec(rows,cols,i,j)
numvec=1:(rows*cols);
nummat=reshape(numvec,rows,cols);
x=nummat(i,j);
end