function data = pre_zeroreset(data)

[Nch,Nt] = size(data);

for ch = 1 : Nch
data(ch,:)=data(ch,:)-mean(data(ch,:),2);
end