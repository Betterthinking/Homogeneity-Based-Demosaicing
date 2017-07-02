function psnr = EZPSNR(former,current,i)

sizemat = size(former);
m=sizemat(1);
n=sizemat(2);
former_plot =former(:,:,i);
current_plot =current(:,:,i);


psnr = 10*log10(255*255*m*n/(sum(sum((former_plot-current_plot).^2))));

end