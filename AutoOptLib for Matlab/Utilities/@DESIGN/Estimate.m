function NewAlgs = Estimate(NewAlgs,Problem,Setting,indInstance,Surrogate)
% Estimate the designed algorithm's performance by a surrogate model.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

EmbedAlgs = Surrogate.UseEmbed(NewAlgs,Setting); 

for i = 1:length(NewAlgs)
    for j = 1:length(Problem(indInstance))
        NewAlgs(i).performanceApprox(indInstance(j),:) = predict(Surrogate.model,EmbedAlgs(i,:));
    end
end
end