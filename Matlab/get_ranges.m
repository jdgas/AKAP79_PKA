function ranges=get_ranges(file_name, scale)

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