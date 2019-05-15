function [which_img,novel_vs_repeat,img_cnd] = get_image_numbers(cfg,itmlist,sequence_items,img_on_code)
% define which imagse are novel and which images are repeat

% Code rechecked for bugs October 17, 2016 SDK

which_img = NaN(1,96*2); %item # for image
img_cnd = NaN(1,96*2); %condition # for image
img_count = 1;
for t = 1:length(cfg.trl);
    if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)
        if any(cfg.trl(t).allval == img_on_code); %should also have codes 200 and 24
            if length (cfg.trl(t).cnd) > 1
                which_img(img_count) = itmlist(cfg.trl(t).cnd(1)-1000);
                img_cnd(img_count) = cfg.trl(t).cnd(1);
                warning('More than 1 condition for this trial')
            else
                which_img(img_count) = itmlist(cfg.trl(t).cnd-1000);
                img_cnd(img_count) = cfg.trl(t).cnd; 
            end
            img_count = img_count+1;
        end
    end
end
which_img = which_img-which_img(1)+1; %convert item # into image# by subracting first item #

novel_vs_repeat = NaN(1,96*2);
for img = 1:max(which_img)
    imgind = find(which_img == img);
    if ~isempty(imgind); 
        novel_vs_repeat(imgind(1)) = 1; %first presentation is novel
    end
    if length(imgind) == 2
        novel_vs_repeat(imgind(2)) = 2; %second presentation is repeated
    elseif length(imgind) == 3
        %emailme(['Temporal Analysis Importing Data found 3 image presentations. Img #' num2str(img) ' ' cfg.dataset(end-20:end-11)])
        %remove that image from analysis. For TO set 13 seems to be a random
        %cortex accidental bug/error not in item/cnd file. Error is found
        %in both .plx and cortex save file. Should be image 65
        which_img(imgind) = NaN; %will get removed when laundried
%     else
%         emailme(['Temporal Analysis Importing Data found 1 image presentations. Img #' num2str(img) ' ' cfg.dataset(end-20:end-11)])
%     	% For TO set 13 seems to be a random
%         %cortex accidental bug/error not in item/cnd file. Error is found
%         %in both .plx and cortex save file. Should be image 39, Will have
    end
end
end