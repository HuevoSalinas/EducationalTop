function [] = PlotWindowLimits(NODE, d)
narginchk(1,2)
winlim = [min(min(NODE,[], 3),[], 1) , max(max(NODE,[], 3),[], 1)];
if size(NODE,2)==2
    daux = max(abs([winlim(3)- winlim(1); ...
                 winlim(4)- winlim(2)]))/20; 
else
    daux = max(abs([winlim(4)- winlim(1); ...
                 winlim(5)- winlim(2); ...
                 winlim(6)- winlim(3)]))/20; 

end
if nargin == 1
    d = daux;
end

if size(NODE,2)==2
    winlim(1:2) = winlim(1:2) - d;
    winlim(3:4) = winlim(3:4) + d;
    bool_ind = winlim(1:2) == winlim(3:4);
    ind_1 = 1:2;
    ind_2 = 3:4;
    winlim(ind_1(bool_ind)) = winlim(ind_1(bool_ind)) - daux;
    winlim(ind_2(bool_ind)) = winlim(ind_2(bool_ind)) + daux;
    xlim([winlim(1), winlim(3)])
    ylim([winlim(2), winlim(4)])
   
else
    winlim(1:3) = winlim(1:3) - d;
    winlim(4:6) = winlim(4:6) + d;
    bool_ind = winlim(1:3) == winlim(4:6);
    ind_1 = 1:3;
    ind_2 = 4:6;
    winlim(ind_1(bool_ind)) = winlim(ind_1(bool_ind)) - daux;
    winlim(ind_2(bool_ind)) = winlim(ind_2(bool_ind)) + daux;
    xlim([winlim(1), winlim(4)])
    ylim([winlim(2), winlim(5)])
    zlim([winlim(3), winlim(6)])
end
end