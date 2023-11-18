

function z = constructFilter(spikeHistory,baseline)


mu    = mean(baseline);

sigma = std(baseline);

z = (spikeHistory - mu) / sigma;

return
















