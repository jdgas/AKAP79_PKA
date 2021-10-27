function ranges=get_ranges(file_name, scale)
%  Functions for making simulated trajectories from posterior distribution
%  sample and reproducing figures in eLife 2021;10:e68164 
 
%  Copyright (C) 2021 Olivia Eriksson (olivia@kth.se) and Parul Tewatia (parul.tewatia@scilifelab.se)

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.

T=readtable(file_name);
values=T.Value;
ranges=NaN(length(values),2);
ranges(:,1)=values./scale;
ranges(:,2)=values.*scale;
%set specific ranges
ranges(7,1)=values(7)./scale;
ranges(7,2)=values(7).*scale;
ranges(10,1)=values(10)./scale;
ranges(10,2)=values(10).*scale;
ranges(21,1)=values(21)./1.9;
ranges(21,2)=values(21).*1.9;
ranges(25:27,1)=values(25:27)./1.25;
ranges(25:27,2)=values(25:27).*1.25;
ranges(28:29,1)=values(28:29)./1.5;
ranges(28:29,2)=values(28:29).*1.5;
ranges(30,1)=values(30)./2;
ranges(30,2)=values(30).*2;

ranges=log10(ranges);

end