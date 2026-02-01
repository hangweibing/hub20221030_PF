
fin=fopen('Asteixspectra.txt','r');
indexDfeFind=0;
k=0;
while ~feof(fin) 
    tline = fgets(fin);
    if(~isempty(find(isspace(tline)==0, 1))&&indexDfeFind==0)    %文件中会有空行，跳过
        [PartionList,~] = CalStrPartionList(tline);% 将字符串划分成几个
        if( length(PartionList{1})==2)
            if( PartionList{1}=='88')
                indexDfeFind=1;
            end
        end
    end
    if(indexDfeFind==1)
        if(~isempty(find(isspace(tline)==0, 1)))
            [tmpNumList,tmpNum] = CalStrPartion2NumList(tline);
            for i=1:tmpNum
                k=k+1;
                normalizeSpectra(k)=tmpNumList(i);
            end
        end
    end
end
fclose(fin);
spectra=130/100*normalizeSpectra;
plot(spectra)
% nCycle=length(normalizeSpectra)/2;
nCycle=12000;
cycle=0.5:0.5:nCycle;
plot(cycle,normalizeSpectra(1:2*nCycle)/100)
xlabel('循环数')
ylabel('归一化的载荷幅值')
for i=1:12000
    ration(i)=spectra(2*i)/spectra(2*i-1);
end
maxStress=max(spectra)
minStress=min(normalizeSpectra)
plot(ration)
index=find(ration<0.7);
length(index)




