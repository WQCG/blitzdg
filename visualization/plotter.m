# Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
# See COPYING and LICENSE files at project root for more details.

x = load('../x.dat');

files = glob('../u*.dat');

for ii=1:length(files)
    u = load(files{ii});
    plot(x(:),u(:));
    drawnow;
    pause(0.1);
endfor
