function  Fast(img_name)

img_name = [img_name,'.jpg'];
img = imread(img_name);

processed_img = color_array_downsample(img);
imwrite(processed_img,['processed_',img_name]);

ans_img = homogeneity(processed_img);
imwrite(ans_img,['ans_',img_name]);

ans_bi_img =  cfa_reconstruct(img,'linear');
imwrite(ans_bi_img,['ans_bi_',img_name]);

img = double(img);
ans_img = double(ans_img);
ans_bi_img = double(ans_bi_img);

diff_img = uint8(abs(img-ans_img));
diff_bi_img = uint8(abs(img-ans_bi_img));
imwrite(diff_img,['diff_',img_name]);
imwrite(diff_bi_img,['diff_bi_',img_name]);

disp('normalPSNR(R):');
disp(EZPSNR(img,ans_img,1));
disp('biPSNR(R):');
disp(EZPSNR(img,ans_bi_img,1));
disp('normalPSNR(G):');
disp(EZPSNR(img,ans_img,2));
disp('biPSNR(G):');
disp(EZPSNR(img,ans_bi_img,2));
disp('normalPSNR(B):');
disp(EZPSNR(img,ans_img,3));
disp('biPSNR(B):');
disp(EZPSNR(img,ans_bi_img,3));


end
