function yp = iodine_c(t,y);
yp = [-2.52 0 .08;.84 -.08 0;0 .08 -.1]*y + [150 0 0]';