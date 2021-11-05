function S=ehcofusion(s1,s2,f1,f2,Nf)
%参数说明
[m1,n1]=size(s1);
[m2,n2]=size(s1);
    if m1==m2
        for i=1:m1
            [x1,y1]=size(f1);
            [x2,y2]=size(f2);

            f=[f1 f2];% 索引号大于y1为脉冲2的排序
            [f_sort,index]=sort(f);%对频率序列进行排序
            S1=[s1 s2];
            [x,y]=size(S1);
            S=zeros(1,y);
            for i=1:y
                S2(i)=S1(index(i));

            end

            ff=linspace(f_sort(1),f_sort(end),Nf);

            S=interp1(f_sort,S2,ff);%

        end
    end



end