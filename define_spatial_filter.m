function H = define_spatial_filter(filter_width)
%written October 18, 2016 SDK since everywhere
%want to make sure it's consistent across codes

%filter_width = 6; %std of 2D guassian filter ~ 3 dva, could use 2 dva (4) as well
filter_size = filter_width*6+1;
H = fspecial('gaussian',filter_size,filter_width);

end