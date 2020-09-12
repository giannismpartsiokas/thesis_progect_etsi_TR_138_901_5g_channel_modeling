function cell_array_design = cell_array_design(xo,yo,R,tiers,design)

index = 1;
color = 'b';
cell = cell_design(xo,yo,R,design,index,color);
position(index,1) = xo;
position(index,2) = yo;
index = index + 1;

for m = 1:1:tiers
    for n = 1:1:m
       if(n==1) 
         phase = 0;
         side = 2*m*R*(sqrt(3)/2);
       else
          side1 = yo + 2*m*R*(sqrt(3)/2) - (n-1)*R*(sqrt(3)/2);
          vector(n,1) = xo - 3*R*(1/2)*(n-1);
          vector(n,2) = yo + 2*m*R*(sqrt(3)/2) - (n-1)*R*(sqrt(3)/2);
          side2 = sqrt((vector(n,1) - xo)^2 + ((vector(n,2) - yo)^2));
          phase = acos(side1/side2);
          side = side2;
       end
       for k = 1:1:6
          vector(n,1) = xo + side*cos(phase + pi/2 + (k-1)*(pi/3));
          vector(n,2) = yo + side*sin(phase + pi/2 + (k-1)*(pi/3));
          cell = cell_design(vector(n,1),vector(n,2),R,design,index,color); 
          position(index,1) = vector(n,1);
          position(index,2) = vector(n,2);
          index = index + 1;
       end
    end
end

cell_array_design = position;
    