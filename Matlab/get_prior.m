function draws=get_prior(ranges,nSamples, const_flag)
    draws=NaN(nSamples,length(ranges));
    for i=1:length(ranges)
        a=ranges(i,1);
        b=ranges(i,2);
        r = a + (b-a).*rand(nSamples,1);
        
        draws(:,i)=r;
       
    end
    
end