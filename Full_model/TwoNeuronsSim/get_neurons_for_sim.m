function [neurons] = get_neurons_for_sim(setups)

connections=setups.connections;
gains=setups.gains;
test=setups.test;
neurons=struct;
neurons(1)=struct;
neurons(2)=struct;

neurons(1).gain=0.03;
neurons(1).shape_index=randsample(1:10,1); %2:11

neurons(1).PR=1;
neurons(1).location=0;
neurons(1).delay_mean=30;
neurons(1).delay_var=60;


neurons(2).gain=0.04;
neurons(2).shape_index=randsample(1:10,1); %2:11
neurons(2).PR=1;
neurons(2).location=30;
neurons(2).delay_mean=50;
neurons(2).delay_var=50;

switch connections
    case 'full'
        
    case 'weak'
       neurons(1).PR=0.8;neurons(2).PR=0.2;
    case 'dis'
         neurons(1).PR=0.8;neurons(2).PR=0;

end


switch gains
    case 'balance'
        
    case 'imba'
       neurons(1).gain=0.06;neurons(2).gain=0.015;
end
switch test
    case 'far'
        neurons(2).location=60;
    case 'correct'
        neurons(1).shape_index=11;
        neurons(2).shape_index=11;
end


