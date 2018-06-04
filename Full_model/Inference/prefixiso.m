function [output]=prefixiso(n_grid,y)

mean_val = zeros(n_grid+1,1);
left_val = ones(n_grid+1,1);
error_val = ones(n_grid+1,1);
mean_val(1)=-Inf;

sumy=zeros(n_grid+1,1);
sumy2=zeros(n_grid+1,1);
sumw=zeros(n_grid+1,1);

for i = 2:(n_grid+1)
   mean_val(i)= y(i-1);
   left_val(i)=i;
   sumy(i)=y(i-1);
   sumy2(i)=y(i-1)^2;
   sumw(i)=1;
   
   while mean_val(i)< mean_val(left_val(i)-1)
       sumy(i) = sumy(i)+sumy(left_val(i)-1);
       sumy2(i) = sumy2(i)+sumy2(left_val(i)-1);
       sumw(i) = sumw(i)+sumw(left_val(i)-1);
       mean_val(i) = sumy(i)/sumw(i);
       left_val(i)=left_val( left_val(i)-1);
       
   end
   levelerror = sumy2(i)- sumy(i)^2/sumw(i);
   error_val(i) = levelerror+error_val(left_val(i)-1);
end
output=struct;
output.error=error_val;
output.mean=mean_val;
output.left=left_val;
