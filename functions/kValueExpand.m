function kval = kValueExpand(mitIndexStart, mitIndexEnd, mitList, middlePart, kvalFull)
index = 1;

%%logical fix/hack
mitIndexStart = mitIndexStart*1;

for i = mitIndexStart:mitIndexEnd-1
    kvalMid(mitList(i):mitList(i+1)-1) = middlePart(index);
    index = index + 1;
end
%kvalMid(mitList(mitIndexEnd))=kvalIn(index);

kval = [kvalFull(1:mitList(mitIndexStart)-1);...
    kvalMid(mitList(mitIndexStart):mitList(mitIndexEnd)-1)';...
    kvalFull(mitList(mitIndexEnd):end)];
end