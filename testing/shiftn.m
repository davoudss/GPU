function shiftn(X,Y,dim)

%
% for n=1:N
%     for m=1:N
%         fprintf('%d\t',x(n + (m-1)*N))
%     end
%     fprintf('\n')
% end
N = size(X,1);

if(dim==1)
    for n=1:N
        for m=1:N/2
            temp=X(n,m+N/2);
            Y(n,m+N/2)=X(n,m);
            Y(n,m) = temp;
        end
    end
else
    for n=1:N/2
        for m=1:N
            temp=X(n+N/2,m);
            Y(n+N/2,m)=X(n,m);
            Y(n,m) = temp;
        end
    end
end

X 
Y



end