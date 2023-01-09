function Algs = Estimate(Algs,Problem,Setting,indInstance,Surrogate)
% Estimate the designed algorithm's performance by surrogate.

EmbedAlgs = Surrogate.UseEmbed(Algs,Setting); 

for i = 1:length(Algs)
    for j = 1:length(Problem(indInstance))
        Algs(i).performanceApprox(indInstance(j),:) = predict(Surrogate.model,EmbedAlgs(i,:));
    end
end
end