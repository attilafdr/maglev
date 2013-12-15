% clear all; close all; clc;

low = mkfilter(10,4,'cheby',0.1);
tflow = tf(low);
bode(low)

hi = mkfilter(160/2/pi,8,'cheby',0.5);
[n d] = tfdata(hi);
[numlow denlow] = tfdata(low);
high = tf([1 0 0 0 0 0 0 0 0],d);

[numhigh denhigh] = tfdata(high)

d_high = c2d(high, Ts, 'zoh');
%d_high = c2d(high, 0.01, 'zoh');
[nd, dd, T] = tfdata(d_high)

d_low = c2d(low, 0.0001, 'zoh');
[ndl, ddl, T] = tfdata(d_low)

hi2 = mkfilter(1,8,'cheby',0.1);
tflow2 = tf(hi2);

[n2 d2] = tfdata(hi2);
high2 = tf([1 0 0 0 0 0 0 0 0],d2);
[n2 d2] = tfdata(high2);
%bode(high2)
bode(d_high)