function [score,res] = acc(gnd, res)
    res = bestMap(gnd, res);
    score = sum(gnd == res) / length(gnd);
end
