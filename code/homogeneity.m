%% homogeneity_rec: 
%@brief: This function implement the demosaic algorithm  in paper ADAPTIVE HOMOGENEITY-DIRECTED DEMOSAICING ALGORITHM
%               on IEEE TRANSACTION ON IMAGE PROCESSING  VOL. 14, NO. 3, MARCH 2005
%	    To find more information about the referrence, please visit: http://photo-lovers.org/pdf/hirakawa05mndemosaictip.pdf
%	     Usage: outputs = homogeneity(fig)

function [outputs] = homogeneity(fig)
	f_h = homo_interpolation(fig,'horizontal');
	f_v = homo_interpolation(fig,'vertical');
	H_fh  = CalHomogeneity(f_h,'horizontal');
	H_fv = CalHomogeneity(f_v,'vertical');
	final = Postprocess(f_h,f_v,H_fh,H_fv);
	outputs = uint8(final);
end
%% row_interpolation: function description
function [outputs] = homo_interpolation(fig, arg)
	shape = size(fig);
	shape = shape(1:2);
	R = double(fig(:,:,1));
	G = double(fig(:,:,2));
	B = double(fig(:,:,3));

	R1 = 0;
	G1 = 0;
	G0 = 0;
	B0 = 0;

	h0 = [-0.25,0,0.5,0,-0.25];
	h1 = [0,0.5,0,0.5,0];
	%construct guassian LPF H(u,v)
	D = sqrt(shape(1)^2+shape(2)^2)/6;
	LP = [];
	for u = 1:2*shape(1)
		for v = 1:2*shape(2)
			LP(u,v) =4*exp(-1.*((u-shape(1))^2+(v-shape(2))^2)/(2*D^2));
		end
	end
	temp = [];
	conv_1 = [];
	conv_2 = [];
	conv_3 = [];
	conv_4 = [];

	delta_BG = [];
	delta_RG = [];
	delta_GB = [];
	delta_GR = [];

	if(strcmp(arg,'horizontal'))
		for row=1:shape(1)
			%%for each line
			if(mod(row,2)==0)
				G0 = G(row,1:shape(2));
				R1 = R(row,1:shape(2));
				conv_1 = conv(G0,h1,'full');
				conv_2 = conv(R1,h0,'full');
				conv_1 = conv_1(3:(shape(2)+2));
				conv_2 = conv_2(3:(shape(2)+2));
				temp = conv_1 + conv_2;
				G(row,1:shape(2)) = G0+temp;
				%G1 = temp;
				%temp = conv(R1-G1,LP);
				%temp = temp(2:shape(2)+1) + G(row,1:shape(2));
				%R(row,1:shape(2)) = temp;
			else
				G0 = G(row,1:shape(2));
				B1 = B(row,1:shape(2));
				conv_1 = conv(G0,h1,'full');
				conv_2 = conv(B1,h0,'full');
				conv_1 = conv_1(3:(shape(2)+2));
				conv_2 = conv_2(3:(shape(2)+2));
				temp = conv_1 + conv_2;
				G(row,1:shape(2)) = G0+temp;
				%G1 = temp;
				%temp = conv(B1-G1,LP);
				%temp = temp(2:shape(2)+1) + G(row,1:shape(2));
				%B(row,1:shape(2)) = temp;
			end
		end
	elseif(strcmp(arg,'vertical'))
		for col=1:shape(2)
			%%for each col
			if(mod(col,2)==0)
				G0 = G(1:shape(1),col)';
				R1 = R(1:shape(1),col)';
				conv_1 = conv(G0,h1,'full');
				conv_2 = conv(R1,h0,'full');
				conv_1 = conv_1(3:(shape(1)+2));
				conv_2 = conv_2(3:(shape(1)+2));
				temp = conv_1 + conv_2;
				G(1:shape(1),col) = (G0+temp)';
				%G1 = temp;
				%temp = conv(R1-G1,LP);
				%temp = temp(2:shape(2)+1) + G(1:shape(1),col)';
				%R(1:shape(1),col) = temp';
			else
				G0 = G(1:shape(1),col)';
				B1 = B(1:shape(1),col)';
				conv_1 = conv(G0,h1,'full');
				conv_2 = conv(B1,h0,'full');
				conv_1 = conv_1(3:(shape(1)+2));
				conv_2 = conv_2(3:(shape(1)+2));
				temp = conv_1 + conv_2;
				G(1:shape(1),col) = (G0+temp)';
				%G1 = temp;
				%temp = conv(B1-G1,LP);
				%temp = temp(2:shape(2)+1) + G(1:shape(1),col)';
				%B(1:shape(1),col) = temp';
			end
		end
	else
		disp('Error!! The homo_interpolation could only select horizontal or vertical orientations!');
		return;
    end
    
    G_sample1 = G;
    G_sample2 = G;
    G_sample1([2:2:shape(1)],[2:2:shape(2)]) = 0;
    G_sample2([1:2:shape(1)],[1:2:shape(2)]) = 0;
    G_sample1 = G - G_sample1;
    G_sample2 = G - G_sample2;
    
    delta_RG = R - G_sample1;
    delta_BG = B - G_sample2;
    fRG = [];
    fBG = [];
	%for i = 1:shape(1)
	%	for j=1:shape(2)
	%		conv_1(i,j) = delta_RG([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]).*LP;
	%		conv_2(i,j) = delta_BG([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]).*LP;
	%	end
	%end
	%Low pass filter the differential images
    for y = 1:shape(1)
 	for x = 1:shape(2)
		fRG(y,x) = delta_RG(y,x).*(-1)^(x+y);
		fBG(y,x) = delta_BG(y,x).*(-1)^(x+y);
	end
    end
    % Here enlarge the size of picture to twice to avoid the aliasing in
    % circular convolution
    conv_1 = ifft2(fft2(fRG,2*shape(1),2*shape(2)).*LP);
    conv_2 = ifft2(fft2(fBG,2*shape(1),2*shape(2)).*LP);
    conv_1 = conv_1(1:shape(1),1:shape(2));
    conv_2 = conv_2(1:shape(1),1:shape(2));
    for y = 1:shape(1)
	for x = 1:shape(2)
		conv_1(y,x) = real(conv_1(y,x).*(-1)^(x+y));
		conv_2(y,x) = real(conv_2(y,x).*(-1)^(x+y));
	end
    end
    R = G+conv_1;
    B = G+conv_2;

    buf = zeros(size(fig));
    buf(:,:,1) = R;
    buf(:,:,2) = G;
    buf(:,:,3) = B;

    outputs = buf;
end
%% CalHomogeneity: function description
function [outputs] = CalHomogeneity(fig, arg)
	shape = size(fig);
	shape2D = shape(1:2);
	lab = rgb2lab(fig);
	L = lab(:,:,1);
	a = lab(:,:,2);
	b = lab(:,:,3);
	L_critiria = 0;
	ab_critiria = 0;
	delta = 2;
	el = 0;
	ec = 0;
	H = [];
	mask = [];
	mask_a = [];
	mask_b = [];
	total_neighbor = [];
	for i=1:shape2D(1)
		for j=1:shape2D(2)
			if(strcmp(arg,'vertical'))
			                if(j==1)
			                    el = abs(L(i,j+1)-L(i,j));
			                    ec = sqrt((b(i,j+1)-b(i,j))^2+(a(i,j+1)-a(i,j))^2);
			                elseif(j==shape2D(2))
			                    el = abs(L(i,j-1)-L(i,j));
			                    ec = sqrt((b(i,j-1)-b(i,j))^2+(a(i,j-1)-a(i,j))^2);
			                else
			                    el = max(abs(L(i,j-1)-L(i,j)),abs(L(i,j+1)-L(i,j)));
			                    ec = max(sqrt((a(i,j-1)-a(i,j))^2+(b(i,j-1)-b(i,j))^2),sqrt((a(i,j+1)-a(i,j))^2+(b(i,j+1)-b(i,j))^2));
			                end
			elseif(strcmp(arg,'horizontal'))
			                 if(i==1)
			                    el = abs(L(i+1,j)-L(i,j));
			                    ec = sqrt((b(i+1,j)-b(i,j))^2+(a(i+1,j)-a(i,j))^2);
			                elseif(i==shape2D(1))
			                    el = abs(L(i-1,j)-L(i,j));
			                    ec = sqrt((b(i-1,j)-b(i,j))^2+(a(i-1,j)-a(i,j))^2);
			                else
			                    el = max(abs(L(i-1,j)-L(i,j)),abs(L(i+1,j)-L(i,j)));
			                    ec = max(sqrt((a(i-1,j)-a(i,j))^2+(b(i-1,j)-b(i,j))^2),sqrt((a(i+1,j)-a(i,j))^2+(b(i+1,j)-b(i,j))^2));
			                 end
			else
				disp('Error!! The homo_interpolation could only select horizontal or vertical orientations!');
				return;
			end

			mask = L([max(i-delta,1):min(i+delta,shape(1))],[max(j-delta,1):min(j+delta,shape(2))]);
			total_neighbor = size(mask);
			L_critiria = (abs(mask - L(i,j))<el);
			mask_a =  a([max(i-delta,1):min(i+delta,shape(1))],[max(j-delta,1):min(j+delta,shape(2))]);
			mask_b =  b([max(i-delta,1):min(i+delta,shape(1))],[max(j-delta,1):min(j+delta,shape(2))]);
			ab_critiria = sqrt((mask_a - a(i,j)).^2+(mask_b-b(i,j)).^2)<ec;
			H(i,j) = sum(sum(L_critiria&ab_critiria))/(total_neighbor(1)*total_neighbor(2));
		end
    end
    outputs = H;
end
%% Postprecess: function description
function[outputs] =  Postprocess(fig_h, fig_v, H_fh, H_fv)
	shape = size(fig_h);
	shape2D = shape(1:2);

	final = zeros(shape);
	h_homogeneity = 0;
	v_homogeneity = 0;
	for i=1:shape2D(1)
		for j=1:shape2D(2)
			h_homogeneity = sum(sum(H_fh([max(i-1,1):min(i+1,shape(1))],[max(j-1,1):min(j+1,shape(2))])));
			v_homogeneity = sum(sum(H_fv([max(i-1,1):min(i+1,shape(1))],[max(j-1,1):min(j+1,shape(2))])));
			if(h_homogeneity>=v_homogeneity)
				final(i,j,:) = fig_h(i,j,:);
			else
				final(i,j,:) = fig_v(i,j,:);
			end
		end
	end

	R = final(:,:,1);
	G = final(:,:,2);
	B = final(:,:,3);
	%reduce discontinuity effect
	for m = 1:3
		delta_RG = R-G;
		delta_BG = B-G;
		delta_GR = -1.*delta_RG;
		delta_GB = -1.*delta_BG;
		conv_1 = [];
		conv_2 = [];
		conv_3 = [];
		conv_4 = [];
		for i = 1:shape2D(1)
			for j=1:shape2D(2)
				conv_1(i,j) = ArrayMedian(delta_RG([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]));
				conv_2(i,j) = ArrayMedian(delta_BG([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]));
				conv_3(i,j) = ArrayMedian(delta_GR([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]));
				conv_4(i,j) = ArrayMedian(delta_GB([max(i-2,1):min(i+2,shape(1))],[max(j-2,1):min(j+2,shape(2))]));
			end
		end
		R = conv_1 + G;
		B = conv_2 + G;
		G = 0.5.*(conv_3+conv_4+R+B);
	end
	final(:,:,1) = R;
	final(:,:,2) = G;
	final(:,:,3) = B;
	outputs = final;
end

function mid = ArrayMedian(X)
    shape = size(X);
    temp = [];
    for i = 1:shape(1)
            temp = [temp X(i,:)];
    end
    mid = median(temp);
end
