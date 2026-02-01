figure
plot([35,30,30,35],[7,7,13,13],'b')
axis equal
axis([30 50 7 13])
hold on
z_left=linspace(35,37.82842712,50);
y1=14-sqrt(9-(z_left-37.82842712).^2);
y2=sqrt(9-(z_left-37.82842712).^2)+6;
plot(z_left,y1,'b',z_left,y2,'b')
hold on
plot([37.82842712,50],[11,11],'b')
hold on
plot([50,37.82842712],[9,9],'b')
hold on
