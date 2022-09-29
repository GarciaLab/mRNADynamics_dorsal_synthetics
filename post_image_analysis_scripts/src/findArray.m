function pos = findArray(A, n)

    %a is array
    %n is number of consecutive numbers
    
    k = [true;diff(A(:))~=1 ];
    s = cumsum(k);
    x =  histc(s,1:s(end));
    idx = find(k);
    
    %return position of each string of n-consecutive numbers within A
    pos = idx(x==n);

end
