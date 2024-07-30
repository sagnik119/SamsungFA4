function pSet = powerSet(S)
    n = length(S);            
    numSubsets = 2^n;         
    pSet = cell(1, numSubsets); 

    for i = 0:numSubsets-1
        binaryIndex = bitget(i, n:-1:1);  
        subset = S(logical(binaryIndex)); 
        pSet{i+1} = subset;               
    end
end