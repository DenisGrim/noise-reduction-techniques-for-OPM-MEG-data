function [data_seg, Nseg] = pre_segmentdata(data, trigger, start_status, from, to)
% Segment Trial Data 
%
% data[Nch,Nt] --> data_seg [Nch*Nt*Ntrial]
% start_status : correspond to zero time point of single trial 
%
% 2009/10/21 OY
% * multiple start_status is supported 

[Nch,Nt] = size(data);

l = [];
for i = 1 : length(start_status)
l= [l, event_start(trigger, start_status(i))]; 
end
l = sort(l); % ascending order (timeseries order)

% Count the number of segmentations (i.e. trials)
k=1;
for i=1:length(l);
    if l(i)+from+1>0 && l(i)+to<=size(data,2);
        k=k+1;
    end
end;
Nseg = k-1;

% Segment data
Ntseg = to - from;
data_seg = zeros(Nch,Ntseg,Nseg);
k=1;
for i=1:length(l);
    if l(i)+from+1>0 && l(i)+to<=size(data,2);
        data_seg(:,:,k)=data(:,l(i)+from+1:l(i)+to);
        k=k+1;
    end
end;


        