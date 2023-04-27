function out=gaufilter(ksize,sigmac,in)
% 三个参数分别指定kernel size、标准差和原始图像
[h,w,d]=size(in);
[X,Y]=meshgrid(1:ksize,1:ksize);
dist=(X-ceil(ksize/2)).^2+(Y-ceil(ksize/2)).^2;
Gc=exp(-dist/2/sigmac^2);

if d==1
    im=in;
    padim=padarray(im,[floor(ksize/2),floor(ksize/2)],0,'both');

    for i=1:h
        for j=1:w
            temp=double(padim(i:i+4,j:j+4));
            Wp=sum(sum(Gc));
            out(i,j)=sum(sum(Gc.*temp))/Wp;
        end
    end

else
    for k=1:d
        im=in(:,:,k);
        padim=padarray(im,[floor(ksize/2),floor(ksize/2)],0,'both');

        for i=1:h
            for j=1:w
                temp=double(padim(i:i+ksize-1,j:j+ksize-1));
                Wp=sum(sum(Gc));
                out(i,j,k)=sum(sum(Gc.*temp))/Wp;
            end
        end
    end
end

out=uint8(out);
end



