function rec = cfa_reconstruct(f,method)
% @brief: this function offers some simple interpolation method for demosaicing 
                usage: rec = cfa_reconstruct(img,method)
% @param:
%             img: the referrence to the input filtered image
%             method: a string to define the interpolation method : "cubic" ,"linear" or "nearest"

temp = 0;
shape = size(f);
shape = shape(1:2);
buf = [];

green_map_shape = 0;

xi = [];
yi = [];

for i = 1:shape(1)
    xi = [xi;1:shape(2)];
end
for j = 1:shape(2)
    yi = [yi;1:shape(1)];
end
yi = yi';

for c = 1:3
    temp = double(f(:,:,c));
    x0 = [];
    y0 = [];
    z0 = [];
    switch(c)
        case 1%red
            z0 = temp(2:2:shape(1),2:2:shape(2));
            for i = 2:2:shape(1)
                x0 = [x0;2:2:shape(2)];
            end
            for j = 2:2:shape(2)
                y0 = [y0;2:2:shape(1)];
            end
            y0 = y0';
            temp = interp2(x0,y0,z0,xi,yi,method);
            f(:,:,c) = uint8(temp);
        case 2%green
            for i = 1:shape(1)
                if(mod(i,2)==1)
                    if(mod(shape(2),2)==0)
                        buf = 2:2:shape(2);
                    else
                        buf = [2:2:shape(2) 1];
                    end
                    x0 = [x0;buf];
                else
                    buf = 1:2:shape(2);
                    x0 = [x0;buf];
                end
                y0 = [y0;ones(1,length(buf)).*i];
            end
            green_map_shape = size(x0);
            for i = 1:green_map_shape(1)
                for j = 1:green_map_shape(2)
                    z0(i,j) = temp(y0(i,j),x0(i,j));
                end
            end
            temp = griddata(x0,y0,z0,xi,yi,method);
            f(:,:,c) = uint8(temp);
        case 3%blue
            z0 = temp(1:2:shape(1),1:2:shape(2));
            for i = 1:2:shape(1)
                x0 = [x0;1:2:shape(2)];
            end
            for j = 1:2:shape(2)
                y0 = [y0;1:2:shape(1)];
            end
            y0 = y0';
            temp = interp2(x0,y0,z0,xi,yi,method);
            f(:,:,c) = uint8(temp);
    end
end

rec = f;