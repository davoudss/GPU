function test2d
tic
% 
imglist = {'flujet', ... Fluid Jet
           'spine', ... Bone
           'gatlin', ... Gatlinburg
           'durer', ... Durer
           'detail', ... Durer Detail
           'cape', ... Cape Cod
           'clown', ... Clown
           'earth', ... Earth
           'mandrill', ... Mandrill
           'spiral'};

colorlabels = {'default', 'hsv','hot','pink',...
               'cool','bone','prism','flag',...
               'gray','rand'};

load(imglist{4},'X','map');


X=X(20:40,20:40);

x = 1:size(X,1);
y = 1:size(X,2);
subplot(311)
colormap(map)
imagesc(X)

[x y] = meshgrid(x,y);
g = exp(-(x.^2 + y.^2)/(2*.1));
h  =conv2(X,g);

subplot(312)
imagesc(h)
axis([1 size(X,2) 1 size(X,1)])

h = zeros(size(X,1),size(X,2));

for ii=1:size(X,1)
    for jj=1:size(X,2)
        for c1=1:size(X,1)
            for c2=1:size(X,2)
                h(ii,jj) = h(ii,jj) + X(ii,jj)*gauss(c1-ii,c2-jj);
            end
        end
    end
end

subplot(313)
imagesc(h)
axis([1 size(X,2) 1 size(X,1)])

toc

function val = gauss(x, y)
exponent = (x^2 + y^2)./(2*.10);
val       = (exp(-exponent)); 