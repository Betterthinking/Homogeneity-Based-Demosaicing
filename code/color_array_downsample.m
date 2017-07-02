function out = color_array_downsample(f)
%@brief: this function simulate the channle downsample procedure in bayer filter
                Usage: out = color_array_downsample(f)
%@params: f : referrence to the input image

shape = size(f);
if(length(shape)~=3)
    disp('Error! The input figure must have 3 channles');
    return;
end

%new = uint8(zeros(shape(1:2)));

%for j = 1:shape(1)
%    for i = 1:shape(2)
%        if(mod(i,2)==0&&mod(j,2)==0)
%            new(j,i) = f(j,i,1);
%        elseif(mod(i,2)==1&&mod(j,2)==1)
%            new(j,i) = f(j,i,3);
%        else
%            new(j,i) = f(j,i,2);
%        end
%    end
%end

new = uint8(zeros(shape));
for j = 1:shape(1)
    for i = 1:shape(2)
        if(mod(i,2)==0&&mod(j,2)==0)
            new(j,i,1) = f(j,i,1);
        elseif(mod(i,2)==1&&mod(j,2)==1)
            new(j,i,3) = f(j,i,3);
        else
            new(j,i,2) = f(j,i,2);
        end
    end
end

out = new;