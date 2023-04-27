# Bilateral_Filtering

matlab实现双边滤波器

### 基本思路

高斯滤波是一种常用而且简单的降噪算法。但其缺点也很明显，就是会不加区分地将噪声与图像边缘等细节一起平滑处理。

而双边滤波是基于高斯滤波进行改进的方法。在高斯滤波的基础上，双边滤波引入了像素值相似性的权值，使得滤波的结果不仅能够平滑图像，还能保留边缘。

对一幅图像，在平坦的地方像素变化程度较小，用高斯滤波降噪可以取得不错的效果。但是在边缘处，往往会出现剧烈的像素值变化，这时候考虑到边缘两侧的像素值应该是有很大不同的，所以引入像素值相似性权重。即，在一个区域S内，一条边缘将S划分为两部分A和B，A和B的像素值整体差距较大。在做降噪处理时，则把更多的权重倾斜到与中心点相似的区域。如果中心点在AB交界处，但是在A内，则A内的像素点会更多贡献给该点像素值的计算。而B内的像素点则贡献较小。

综上可知，与高斯滤波一样，其基本思想就是通过权值来控制某一个点周边像素值对该点像素值的贡献值。在双边滤波中，考虑了两项权值：距离和相似性。距离即是像素点与中心点的欧式距离，而相似性则是像素点值与中心点值的差值。基于以上思想，得到如下公式：
$$

\overline I(p)=\frac{1}{W_p}\sum_{q\in S}G_c(||p-q||)G_s(|I(p)-I(q)|)I(q)

$$
其中：
$$

W_p=\sum_{q\in S}G_c(||p-q||)G_s(|I(p)-I(q)|)

$$
p为需要计算像素值的中心点，S为卷积核覆盖的区域，$\overline I(p)$为计算出的像素值，$I(p)$为中心点原始像素值，$I(q)$为区域内某个点的像素值，$G_c$和$G_s$分别代表距离权值和相似权值。如果我们用(x, y)表示p的坐标，(i, j)表示q的坐标，基于二维高斯函数，则有：
$$

G_c=e^{(-\frac{(x-i)^2+(y-j)^2}{2\sigma_c^2})}\\
G_s=e^{(-\frac{(I(p)-I(q))^2}{2\sigma_s^2})}\\

$$
$\sigma_c$和$\sigma_s$则分别代表两个标准差。

值得注意的是，在高斯滤波中，相比于标准差小的时候，标准差增大会使得距离中心远的点的权重增加，这一性质由高斯函数决定：当标准差减小时，高斯函数的图像会变得更加平坦。而在双边滤波中，不论是空间标准差还是相似标准差，其增大时，同样在空间域或者值域上体现出相似效果，即标准差越大，与中心点距离远（不相似）的点的权重越大。

同时，方差越大，说明权重差别越小，因此表示不强调这一因素的影响，反之，则表示更强调这一因素导致的权重的不均衡。因此两个方面的某个的方差相对变小 表示这一方面相对较重要，得到强调。如$\sigma_c$变小，表示更多采用近邻的值作平滑，说明图像的空间信息更重要，即相近相似。如$\sigma_r$变小，表示和自己同一类的条件变得苛刻，从而强调值域的相似性。

### 代码实现

基于以上思想，我尝试使用matlab实现双边滤波

```matlab
function out=bifilter(ksize,sigmac,sigmas,in)
% 三个参数分别指定kernel size、空间域的标准差、值域的标准差和原始图像

[h,w,d]=size(in);
% 获取输入图像三个维度的大小，用于之后的处理及判断


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
```

### 效果演示

#### 空间域效果

原图整体偏暗，所以我先对原图用gamma函数进行处理
![image](https://user-images.githubusercontent.com/76733734/234800545-ff84cf39-398c-47cd-a8b2-508a8ee94aa3.png)

在对比度处理的基础上，我先使用matlab自带的高斯滤波和双边滤波进行了处理。

可以看到高斯滤波后的图像较为模糊，而双边滤波后的图像整体更加清晰，同时有一定的降噪效果。

![image](https://user-images.githubusercontent.com/76733734/234800593-104cf533-4156-4a62-9047-32c7033863f8.png)

对比我自己实现的双边滤波与matlab自带的双边滤波，都实现了减少噪点且保留细节的目标。因为我自己的函数，输入的标准差比matlab默认的更大，所以相对而言平滑的效果更加明显。

![image](https://user-images.githubusercontent.com/76733734/234800667-631ca6c9-637e-49a4-bb01-a4108f0190c9.png)

同时我也实现了高斯滤波，与自己实现的双边滤波对比，整体较为模糊，还是双边滤波效果比较理想。

![image](https://user-images.githubusercontent.com/76733734/234800708-c7f7e977-4c3f-423d-a9a2-f2a8a42e9195.png)

#### 频域效果

原始图像频谱取对数后如下（以下所有频谱都是对数处理后的频谱）

![image](https://user-images.githubusercontent.com/76733734/234800765-72f38b37-3f13-44a6-a1e2-e082bef2a3c9.png)

可以看到高斯滤波类似于一个低通滤波器，将高频成分大部分滤除，主要保留低频成分

![image](https://user-images.githubusercontent.com/76733734/234800806-425cd696-5754-4769-870d-fdf254129a50.png)

双边滤波则并不是简单把高频滤除保留低频，滤波后的频谱图与保留着与原始频谱相似的特征，但也滤去了一些成分，其中既有高频也有低频。

![image](https://user-images.githubusercontent.com/76733734/234800822-ec6b1116-8e09-431a-99da-bb2ec7e6701a.png)

