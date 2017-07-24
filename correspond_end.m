function [L, D] = correspond_end(Wgt, W)
    [~,c] = size(Wgt);
    dist = hyperSam(Wgt,W);
    L = zeros(c, 1);
    D = zeros(c, 1);
    for i=1:c
        [v1, k1] = min(dist);
        [v2, k2] = min(v1);
        dist(k1(k2),:) = 5.0;
        dist(:,k2) = 5.0;

        L(k2) = k1(k2);
        D(k2) = v2;
    end
end