function idx = my_sub2ind2(D,i,j)
% D = #rows, i = row subscript, j = column subscript
i = i(:); j = j(:)';
idx =  reshape(bsxfun(@plus,D*(j-1),i),1,[]);