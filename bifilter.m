function out=bifilter(ksize,sigmac,sigmas,in)
% 三个参数分别指定kernel size、空间域的标准差、值域的标准差和原始图像

[h,w,d]=size(in);
% 获取输入图像三个维度的大小，用于之后的处理及判断

out=zeros(h,w,d);

[X,Y]=meshgrid(1:ksize,1:ksize);
dist=(X-ceil(ksize/2)).^2+(Y-ceil(ksize/2)).^2;
Gc=exp(-dist/2/sigmac^2);
% 根据kernel size直接算出空间域的距离矩阵，即空间权值，方便之后调用，省去重复计算

if d==1
% 判断图像是rgb图像还是灰度图像
% 灰度图像直接处理即可，rgb图像需要分三次处理

    im=in;
    padim=padarray(im,[floor(ksize/2),floor(ksize/2)],0,'both');
    % 根据kernel size对图像边界做扩展，补零

    for i=1:h
        for j=1:w
            temp=double(padim(i:i+4,j:j+4));
            % 取出要处理的块
            Gs=exp(-(temp-temp(ceil(ksize/2),ceil(ksize/2))).^2/2/sigmas^2);
            % 计算相似权值
            Wp=sum(sum(Gc.*Gs));
            % 计算总权值，用于归一化
            out(i,j)=sum(sum(Gc.*Gs.*temp))/Wp;
            % 计算输出像素值
        end
    end

else
% 以下与之前内容相似，分三次分别处理三个通道的图像即可
    for k=1:d
        im=in(:,:,k);
        padim=padarray(im,[floor(ksize/2),floor(ksize/2)],0,'both');

        for i=1:h
            for j=1:w
                temp=double(padim(i:i+ksize-1,j:j+ksize-1));
                Gs=exp(-(temp-temp(ceil(ksize/2),ceil(ksize/2))).^2/2/sigmas^2);
                Wp=sum(sum(Gc.*Gs));
                out(i,j,k)=sum(sum(Gc.*Gs.*temp))/Wp;
            end
        end
    end
end

out=uint8(out);
% 类型转换
end



