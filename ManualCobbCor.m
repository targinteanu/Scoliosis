fn_in = uigetfile('*.png');
fn_out = [fn_in,' -- Cobbed ',datestr(datetime, 'yyyy-mm-dd HH.MM.SS'),'.png'];

img_in = imread(fn_in);

%%
f = figure; 
imshow(img_in);
hold on;

[xs,ys] = getpts; 
xs = xs(1:2); ys = ys(1:2);
[xi,yi] = getpts;
xi = xi(1:2); yi = yi(1:2);
% consider making more robust by reordering to avoid 'z' shape


if xs(2) > xs(1)
    al = 'right';
else
    al = 'left';
end

dri = [diff(xi), diff(yi)]; drs = [diff(xs), diff(ys)];
a = (-drs*dri')/(norm(drs)*norm(dri)); a = acos(a)*180/pi

displength = 1.5;
Rs = [xs(1),ys(1)] + displength*drs; xs(2) = Rs(1); ys(2) = Rs(2);
Ri = [xi(2),yi(2)] - displength*dri; xi(1) = Ri(1); yi(1) = Ri(2);

plot([xs; xi], [ys; yi], 'r', 'LineWidth', 1); 
xytxt = [mean([xs(2), xi(1)]), mean([ys(2), yi(1)])];
text(xytxt(1), xytxt(2), [' ',num2str(a,4),'\circ '], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', al, 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');

%%
saveas(f, fn_out, 'png');