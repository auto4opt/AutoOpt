function NewAlgs = Estimate(NewAlgs,Problem,Setting,indInstance,Surrogate)
% Estimate the designed algorithm's performance by surrogate.

EmbedAlgs = Surrogate.UseEmbed(NewAlgs,Setting); 

for i = 1:length(NewAlgs)
    for j = 1:length(Problem(indInstance))
        NewAlgs(i).performanceApprox(indInstance(j),:) = predict(Surrogate.model,EmbedAlgs(i,:));
    end
end
end