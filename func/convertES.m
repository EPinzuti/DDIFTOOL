
% need to be corrected to interact with list 
function[kseq, subset]=  convertES(d,t)
kseq=(d-1).*t;

subset=bsxfun(@plus,(round(linspace(0,d-1,d).*t)),1);



