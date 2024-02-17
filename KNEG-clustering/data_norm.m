function [data] = data_norm(data)
%数据归一化1
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;
%数据归一化2
% N=size(data,1);
% maxval=max(data);
% minval=min(data);
% range=maxval-minval;
% for i=1:size(data,2)
%     if range(i)==0
%         range(i)=0.000001;
%     end
% end
% for i=1:N
%     data(i,:)=data(i,:)-minval;
%     data(i,:)=data(i,:)./range(1, :);
% end

end