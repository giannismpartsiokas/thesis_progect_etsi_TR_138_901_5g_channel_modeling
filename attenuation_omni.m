function [A,theta_BS] = attenuation_omni(x_pos,y_pos,position,number_of_cells)
   
sectors = 3;
theta_BS = zeros(1,number_of_cells);
A = zeros(1,sectors*number_of_cells);

for number_of_cell = 1:1:number_of_cells
  
  side1 = sqrt((x_pos-position(number_of_cell,1))^2 + (y_pos-position(number_of_cell,2))^2);
  side2 = abs(x_pos-position(number_of_cell,1));
  angle = acos(side2/side1);
  if((x_pos - position(number_of_cell,1)<0)&&(y_pos - position(number_of_cell,2)>0))
     angle = pi - angle;
  elseif((x_pos - position(number_of_cell,1)<0)&&(y_pos - position(number_of_cell,2)<0))
     angle = angle + pi;
  elseif((x_pos - position(number_of_cell,1)>0)&&(y_pos - position(number_of_cell,2)<0))
     angle = 2*pi - angle;
  end
  if((angle>=0)&&(angle<=(2*pi)/3))
      theta_BS(1,sectors*(number_of_cell-1)+1) = angle;
      theta_BS(1,sectors*(number_of_cell-1)+2) = angle + 2*(pi/3);
      theta_BS(1,sectors*(number_of_cell-1)+3) = angle + 4*(pi/3);
      angle = abs(pi/3 - angle);
      angle = round(angle*(180/pi)+1);
      A(1,sectors*(number_of_cell-1)+1) =  0; 
      A(1,sectors*(number_of_cell-1)+2) = 20;
      A(1,sectors*(number_of_cell-1)+3) = 20;
  elseif((angle>(2*pi)/3)&&(angle<=(4*pi)/3))
      theta_BS(1,sectors*(number_of_cell-1)+1) = angle;
      theta_BS(1,sectors*(number_of_cell-1)+2) = angle - 2*(pi/3);
      theta_BS(1,sectors*(number_of_cell-1)+3) = angle + 2*(pi/3);
      angle = angle - (2*pi)/3;
      angle = abs(pi/3 - angle);
      angle = round(angle*(180/pi)+1);
      A(1,sectors*(number_of_cell-1)+1) = 20;
      A(1,sectors*(number_of_cell-1)+2) = 0;
      A(1,sectors*(number_of_cell-1)+3) = 20;
  else
     theta_BS(1,sectors*(number_of_cell-1)+1) = angle;
     theta_BS(1,sectors*(number_of_cell-1)+2) = angle;
     theta_BS(1,sectors*(number_of_cell-1)+3) = angle - 4*(pi/3); 
     angle = angle - (4*pi)/3;
     angle = abs(pi/3 - angle);
     angle = round(angle*(180/pi)+1);   
     A(1,sectors*(number_of_cell-1)+1) = 20;
     A(1,sectors*(number_of_cell-1)+2) = 20;   
     A(1,sectors*(number_of_cell-1)+3) = 0;
  end
end