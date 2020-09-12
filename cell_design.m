function cell_design = cell_design(xo,yo,R,design,index,color)

if(R<=0)
    cell_design = 0;
    disp('Radius of the cell must be a positive number');
    return;
else
   N = 1e+3;
   x = xo;
   y = yo;
   angle = -pi/2-pi/6;%-120 degrees
   for n = 1:1:6%6 lines per cell
     x_min = x + R*cos(angle);
     y_min = y + R*sin(angle);
     angle = angle + pi/3;
     x_max = x + R*cos(angle);
     y_max = y + R*sin(angle);
     if(design==1)
        side_x = linspace(x_min,x_max,N);
        side_y = linspace(y_min,y_max,N);
        line(side_x,side_y,'Color',color,'LineWidth',1);
     end
   end
end

%plot sectors
if(design==1)
  x_min = xo;
  x_max = xo + R;
  y_min = yo;
  y_max = yo;
  side_x = linspace(x_min,x_max,N);
  side_y = linspace(y_min,y_max,N);
  line(side_x,side_y,'Color','r','LineWidth',1);

  x_min = xo;
  x_max = xo - R*(1/2);
  y_min = yo;
  y_max = yo + R*(sqrt(3)/2);
  side_x = linspace(x_min,x_max,N);
  side_y = linspace(y_min,y_max,N);
  line(side_x,side_y,'Color','r','LineWidth',1);

  x_min = xo;
  x_max = xo - R*(1/2);
  y_min = yo;
  y_max = yo - R*(sqrt(3)/2);
  side_x = linspace(x_min,x_max,N);
  side_y = linspace(y_min,y_max,N);
  line(side_x,side_y,'Color','r','LineWidth',1);
end

if(design==1) 
    grid on;
    text(xo,yo,num2str(index),'FontSize',13)
    title('Hexagonal topology with BS in the center of each cell');
end
cell_design = 1;