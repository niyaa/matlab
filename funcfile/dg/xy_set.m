function [ node_xy ] = xy_set ( nx, ny, node_num,node_numd, xl, xr, yb, yt  )

node_xy = zeros(2,node_num );
node_dxy=zeros(2,node_numd);
for j = 1 : 2*ny-1
 for i = 1 : 2*nx - 1

      node_xy(1,i+(j-1)*(2*nx-1)) =    ...
        ( ( 2 * nx - i - 1 ) * xl   ...
        + (          i - 1 ) * xr ) ...
        / ( 2 * nx     - 2 );

      node_xy(2,i+(j-1)*(2*nx-1)) =    ...
        ( ( 2 * ny - j - 1 ) * yb   ...
        + (          j - 1 ) * yt ) ...
        / ( 2 * ny     - 2 );
    
    k=1;

if i==1||j==1   
node_dxy(1,k)=node_xy(1,i);
node_dxy(2,k)=node_xy(2,i);

else
   for k=i:i+2
node_dxy(1,k)=node_xy(1,i);
node_dxy(2,k)=node_xy(2,i);
   end
end
 end
end


  return
end
