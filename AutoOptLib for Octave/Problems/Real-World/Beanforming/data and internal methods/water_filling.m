function [vk_fill] = water_filling(PT, vk)
% filling line: miu
K=numel(vk);
miu=max(vk);
if sum(miu-vk)>PT
    [vkd,~]=sort(vk,'descend');
    t=1;
    vka_sum=sum(miu-vk);
    while vka_sum>PT && t<numel(vk)
        miu=vkd(t+1);
        vka_sum=sum(max(miu-vk,0));
        if vka_sum<=PT
            miu=miu+(PT-vka_sum)/(K-t);
            break;
        end
        t=t+1;
    end
    vk_fill=max(miu-vk,0);
    
else
    miu=miu+(PT-sum(miu-vk))/K;
    vk_fill=miu-vk;
end
end

