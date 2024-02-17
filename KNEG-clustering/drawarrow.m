function h = drawarrow(x,y,color,lineType,ax)
% https://blog.csdn.net/qq_41093957/article/details/115256114
 
% 这个程序最大的问题是 改变了plot的范围 从而使线段位置不对了
switch nargin
    case 3
        lineType='arrow';
        ax=gca;
    case 4
        ax=gca;
end
% 调整坐标大小以适应箭头长度
% xlimraw=get(ax,'XLim');
% ylimraw=get(ax,'YLim');
xlim=get(ax,'XLim');
ylim=get(ax,'YLim');
xlimmin=xlim(1);xlimmax=xlim(2);
ylimmin=ylim(1);ylimmax=ylim(2);
if xlimmin>min(x(1),y(1)), xlimmin=min(x(1),y(1));end
if xlimmax<max(x(1),y(1)), xlimmax=max(x(1),y(1));end
if ylimmin>min(x(2),y(2)), ylimmin=min(x(2),y(2));end
if ylimmax<max(x(2),y(2)), ylimmax=max(x(2),y(2));end
% ax.XLim = [xlimmin,xlimmax];
% ax.YLim = [ylimmin,ylimmax];
set(ax,'XLim',[xlimmin,xlimmax]);
set(ax,'YLim',[ylimmin,ylimmax]);
xlim=get(ax,'XLim');
ylim=get(ax,'YLim');
% xlim=ax.XLim;
% ylim=ax.YLim;
pos=get(ax,'Position');
x_ratio = pos(3)/(xlim(2)-xlim(1));
y_ratio = pos(4)/(ylim(2)-ylim(1)); % 缩放比例
orig_pos=[-xlim(1)*x_ratio+pos(1),-ylim(1)*y_ratio+pos(2)]; % figure坐标系中的原点坐标
x=x.*[x_ratio,y_ratio];y=y.*[x_ratio,y_ratio];
x=x+orig_pos;y=y+orig_pos;
h = annotation(lineType,[x(1),y(1)],[x(2),y(2)],'headStyle','cback2','HeadLength',5,'HeadWidth',5);
h.Color=color;
% set(ax,'XLim',xlimraw);
% set(ax,'YLim',ylimraw);
end