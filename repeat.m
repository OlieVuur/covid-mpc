% output = repeat(A,num,dim)
% Repeats each row/column of the matrix A, num times, along the dimension
% specified by dim. If dim = 1 then repeat along the columns, if dim = 2 
% then repeat along the rows
function outmat = repeat(A,num,dim)
    if (nargin < 3)
        dim = 1;
    end
    if (dim == 1)
        outmat = zeros(size(A,1),num*size(A,2));
        for i = 1:(num*size(A,2))
            outmat(:,i) = A(:,ceil(i/num));
        end
    else
        outmat = zeros(num*size(A,1),size(A,2));
        for i = 1:(num*size(A,1))
            outmat(i,:) = A(ceil(i/num),:);
        end
    end
end