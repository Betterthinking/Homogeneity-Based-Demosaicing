function  reverse_img(img_name)

img_name = [img_name,'.jpg'];

img = imread(img_name);

sizemat = size(img);

output_image = 255-img;

frame = uint8(zeros(sizemat(1)+8,sizemat(2)+8,3));

for i = 1:sizemat(1)
    for j = 1:sizemat(2)
    frame(i+4,j+4,:)=output_image(i,j,:);
    end
end


%imshow(frame);

imwrite(frame,['reverse_',img_name]);