function [] = draw_demo(id)

%%
if id == 0
load cuprite/cupriteS1_R188.mat
load cuprite/groundTruth_Cuprite_nEnd12.mat
load cuprite.mat

x = Y;
x(find(x > 60000 )) = 0;
x = hyperNormalize(x);

w = 250; 
h = 190;
numbands = 188;

img   = reshape(x', 250, 190, numbands);
img2d = hyperConvert2d(img)';  

Mp = M(slctBnds,:);
H=12;

%%
elseif id == 1

load pavia/PaviaU
load pavia/PaviaU_gt
load pavia.mat

[w, h, c] = size(paviaU);
x    = hyperNormalize(paviaU);
x    = reshape(x, w*h, c);
paviaU_gt = reshape(paviaU_gt, w*h, 1);

img2d  = double(x);
numbands = 103;
H=7;

Mp = zeros(H, numbands);
Mp(1,:) = mean(x((paviaU_gt == 1) | (paviaU_gt == 7),:));
Mp(2,:) = mean(x((paviaU_gt == 2),:));
Mp(3,:) = mean(x((paviaU_gt == 4),:));
Mp(4,:) = mean(x((paviaU_gt == 5),:));
Mp(5,:) = mean(x((paviaU_gt == 6),:));
Mp(6,:) = mean(x((paviaU_gt == 3) | (paviaU_gt == 8),:));
Mp(7,:) = mean(x((paviaU_gt == 9),:));

Mp = Mp';

%%
elseif id==2
    
load gulfport/muufl_gulfport_campus_1_hsi_220_label
load gulfport.mat

[w, h, c] = size(hsi.Data);
x    = reshape(hsi.Data, w*h, c);
paviaU_gt = reshape(hsi.sceneLabels.labels, w*h, 1);

img2d  = double(x);
               
H = 7;
numbands = 64;

Mp = zeros(H, numbands);

Mp(1,:) = mean(x((paviaU_gt == 1),:));
Mp(2,:) = mean(x((paviaU_gt == 2),:));
Mp(3,:) = mean(x((paviaU_gt == 3),:));
Mp(4,:) = mean(x((paviaU_gt == 4),:));
Mp(5,:) = mean(x((paviaU_gt == 5),:));
Mp(6,:) = mean(x((paviaU_gt == 7),:));
Mp(7,:) = mean(x((paviaU_gt == 8),:));
 
Mp = Mp';

%%
elseif id==3
    
load jasper/jasperRidge2_R198.mat
load jasper/end4.mat
load jasper.mat

x =Y;
x(find(x > 60000 )) = 0;

x = hyperNormalize(x);

w = 100;
h = 100;
numbands = 198;

img = reshape(x', w, h, numbands);
img2d=hyperConvert2d(img)';  

Mp = M;
H = 4;

%%
elseif id==4
load urban/Urban_R162.mat
load urban/end4_groundTruth.mat
load urban.mat

x =Y;
x(find(x > 60000 )) = 0;
x = hyperNormalize(x);

w = 307;
h = 307;
numbands = 162;

img = reshape(x', w, h, numbands);
img2d=hyperConvert2d(img)';  

Mp = M;
H = 4;

%%
elseif id==5
load samson/samson.mat
load samson/end3.mat
load samson.mat

x =V;

w=95;
h=95;
numbands = 156;

img = reshape(x', w, h, numbands);
img2d=hyperConvert2d(img); 
img2d = img2d';

Mp = M;
H = 3;

end

[L, D] = correspond_end(Mp, W1');

Wact = W1';
[N,~] = size(img2d);
distf=@(x,y) hyperSam(x,y);
d = zeros(H,N);

d1 = zeros(H,N);
for e=1:H
    for s=1:N
        d1(e,s) = distf(img2d(s,:)', Wact(:,e));   
    end    
end

d2 = zeros(H,H);
for e1=1:H
    for e2=1:H
        d2(e1,e2) = distf(Wact(:,e1), Wact(:,e2));
    end
end

d1 = abs(d1);
d2 = abs(d2);
P = abs(simplex_project_dist(d1, d2));


figure;
for i=1:H
    Pact = zeros(w*h,1);
    Pact(:,1) = P(L(i),:)';
    Pact = reshape(Pact, w, h, 1);
    
    subplot(3,H,i)
    plot(Mp(:,i));
    title(['Endmember #' num2str(i)])
    
    subplot(3,H,i+H)
    plot(W1(L(i),:))
    title(['Endmember #' num2str(i)])
    
    subplot(3,H,i+2*H)
    imshow(Pact);
    title(['Endmember #' num2str(i)])
end

end
