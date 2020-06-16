% Created Date: Saturday January 4th 2020                                      %
% Author: Li Hongmin (li.hongmin.xa@alumni.tsukuba.ac.jp)                      %
%##############################################################################%

% this function is used to generate N base clusterings for ensemble clustering
% where k is range of [k,sqrt(n)]
function baseCls = RPS(fea, k, N)

    n = size(fea, 1);
    %     upperK = round(sqrt(n));
    upperK = max(k, floor(sqrt(n)));
    clsNums = randi([k, upperK], N, 1);
    clsNums = randi([50, 100], N, 1);
    baseCls(n, N) = 0;

    for i = 1:N% This is a paralleled version
        baseCls(:, i) = kmeans(fea, clsNums(i));
    end

end