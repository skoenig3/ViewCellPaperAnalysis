% File generates item and conditions files for the listSQ task
% Written by Seth Konig July 2014

% you can find the section you are looking for by searching [#]
% 1. Generate blank item file
% 2. Generate blank condition file
% 3.0 Generate item and condition files with 2 predictable sequences for sets 1-6
% 3.1 Generate item and condition files with 2 predictable sequences for sets 7-26
% 3.2 Generate item and condition files with 2 predictable sequences forset 27-55
% 4. Sort images into set folders
%%
%%---[1] Generate blank item file---%%
% Create Blank item file so can create a save file but it won't run unless
% you load the correct item file for the day

total_items = 115; %number of items in task

line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
line2 =' -4    1      1    0.00    0.00      0   1.00  1.00  0.00                0   0   0 l';
line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50        200 200 200 x';
line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';

fid = fopen(['blank' num2str(total_items) '.itm'],'w+');
for line = 1:6
    fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
end
for l = 1:total_items;
    if l < 10
        itmspace = '  '; %2 spaces
    else
        itmspace = ' ';%1 space
    end
    str = [itmspace num2str(l)  '\r\n'];
    fprintf(fid,str);
end
fclose(fid);
%%
%%---[2] Generate blank condition file---%%
% Create Blank condition file so can create a save file but it won't run unless
% you load the correct item file for the day

total_cond = 597; %total number of conditions in task

line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
fid = fopen(['blank' num2str(total_cond) '.cnd'],'w+');
for line = 1
    fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
end
for l = 1:total_cond;
    if l < 10
        cndspace = '  ';%2 spaces
    elseif l < 100
        cndpsace = ' ';% 1 space
    else
        cndspace = '';% 0 space
    end
    str = [cndspace num2str(l)  '\r\n'];
    fprintf(fid,str);
end
fclose(fid);
%%
%%---[3.0] Generate item and condition file with 2 predictable sequences for sets 1-6---%%
number_of_sequences = 2;
number_of_items = 4; %# of items per sequence
number_sequence_trials_btwn_images = 3;
set_nums = [1:6];%always start from 1 otherwise won't regenerate perfectly
number_of_images = 96;
image_spacing = 16; %number of images between novel and repeat presentations

%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
familiarization_block = true; %if want familrization trials before image sets else set to false
familiar_trials = 10; %# of familirization trails if want them
if familiarization_block
    max_sequence_conditions = 288+number_of_sequences*familiar_trials;
else
    max_sequence_conditions = 288;
end

item_locations = {};
overlap_index = [];
buffer = 5; %minimum distance between 2 items
maxdist = 15; %maximum distance between 2 items
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',140722); %seed randomly so we can regenerate files if lost
for i = 1:500;
    valid_points = ones(sizey,sizex);
    item_locations{i} = NaN(number_of_items*2,number_of_sequences);
    seq = 1;
    for cross = 1:number_of_items
        if cross == 1
            available_points = find(valid_points);
            choosen_ind = available_points(randi(length(available_points)));
            [y,x] = ind2sub([sizey,sizex],choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2)<= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        elseif cross == 2
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            choosen_ind = Cind(randi(length(Cind)));
            [y,x] = ind2sub(size(valid_points),choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
            %find angle of quickest scan path from previous 2 items
            dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
            dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
            angle12 = atan2d(dy12,dx12);
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            [Cy,Cx] = ind2sub(size(valid_points),Cind);
            potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
            for pc = 1:length(Cy)
                dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                angle23 = atan2d(dy23,dx23);
                potential_angles(pc) = angle23;
            end
            potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
            angle_difference = potential_angles-angle12;
            angle_difference = abs(angle_difference);
            good_angles = find(angle_difference <=  90);
            if isempty(good_angles);
                break;
            else
                choosen_ind = good_angles(randi(length(good_angles)));
                item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                x =  Cx(choosen_ind);
                y =  Cy(choosen_ind);
                Cind = find(C);
                valid_points(Cind) = 0;
            end
        end
    end
    
    if all(~isnan(item_locations{i}(:,1)))
        overlapind = rem(i,4)+1;
        overlap_index(i) = overlapind;
        seq = 2;
        switch overlapind
            case 1
                cross = 1;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_items
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        item_locations{i}(2*cross-1,seq) = x;
                        item_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 items
                        dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                        dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <=  90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
            case 2
                cross = 2;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                cross = 1;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 items
                    dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                    dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
            case 3
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 2nd is now the 3rd cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 3;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                cross = 2;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                cross  = 1;%rewrtie to 2 so that when we flip it will work out
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 items
                    dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                    dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
                xs = item_locations{i}(1:2:end,seq);
                ys = item_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                item_locations{i}(1:2:end,seq) = xs;
                item_locations{i}(2:2:end,seq) = ys;
            case 4
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 1st is now the 4th cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 4;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                cross = 1;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_items
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        item_locations{i}(2*cross-1,seq) = x;
                        item_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 items
                        dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                        dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <= 90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
                xs = item_locations{i}(1:2:end,seq);
                ys = item_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                item_locations{i}(1:2:end,seq) = xs;
                item_locations{i}(2:2:end,seq) = ys;
        end
    end
    
    %code to double check that the sequences we generated meat the desired
    %criterion
    if all(~isnan(item_locations{i}))
        item_locations{i}(1:2:end,:) = item_locations{i}(1:2:end,:)-13;
        item_locations{i}(2:2:end,:) = item_locations{i}(2:2:end,:)-10;
        for seq = 1:size(item_locations{i},2);
            xs =item_locations{i}(1:2:end,seq);
            ys =item_locations{i}(2:2:end,seq);
            
            %double check to make sure distances are good
            d = pdist([xs,ys]);
            if any(d < buffer)
                disp('error item locations too close')
                [min(d),seq]
            end
            dx = diff(xs);
            dy = diff(ys);
            if any(sqrt(dx.^2+dy.^2) > maxdist)
                disp('error locations too far appart')
                crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1;
                if length(crossnum) > 1
                    [crossnum;seq]
                else
                    [crossnum seq]
                end
            end
            if any(abs(xs) > 12) || any(abs(ys) > 9)
                disp('error locations out of bounds')
            end
            
            dy12 = ys(2)-ys(1);
            dx12 = xs(2)-xs(1);
            angle12 = atan2d(dy12,dx12);
            
            %angle from item 2 to 3
            dy23 = ys(3)-ys(2);
            dx23 = xs(3)-xs(2);
            angle23 = atan2d(dy23,dx23);
            
            %angle from item 3 to 4
            dy34 = ys(4)-ys(3);
            dx34 = xs(4)-xs(3);
            angle34 = atan2d(dy34,dx34);
            
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            angle23(angle23 < 0) = 360+angle23(angle23 < 0);
            angle34(angle34 < 0) = 360+angle34(angle34 < 0);
            
            dang123 = angle23-angle12;
            dang234 = angle34-angle23;
            dang123(dang123 < 0) =  dang123(dang123 < 0) + 360;
            dang234(dang234 < 0) =  dang123(dang234 < 0) + 360;
            dang123(dang123 > 180) = 360-dang123(dang123 > 180);
            dang234(dang234 > 180) = 360-dang234(dang234 > 180);
            if dang123 > 90 || dang234 > 90
                disp('locations not in a forward direction')
                dang123
                dang234
                seq
            end
        end
    end
end

all_nans = [];
for s = 1:length(item_locations);
    if any(any(isnan(item_locations{s})))
        all_nans = [all_nans s];
    end
end
item_locations(all_nans) = [];
overlap_index(all_nans) = [];
overlap1 = find(overlap_index == 1);
overlap2 = find(overlap_index == 2);
overlap3 = find(overlap_index == 3);
overlap4 = find(overlap_index == 4);
new_item_locations = [];
new_overlap_index = [];
%ensure that the items that overlaps in a sequences rotates uniformly but
% randomly through 1 to 4
for o = 1:min([length(overlap1),length(overlap2),length(overlap3),length(overlap4)]);
    rr = randperm(4);
    for rp = 1:length(rr);
        temp =  eval(['overlap' num2str(rr(rp))]);
        new_item_locations = [new_item_locations item_locations(temp(o))];
        new_overlap_index = [new_overlap_index overlap_index(temp(o))];
    end
end
item_locations = new_item_locations; %reassign cuz this is the variable the subsequent code uses

clr=['255 255   0 x';
    '  0 255 255 x';
    '255   0 255 x';
    '  0 255   0 x';
    '  0   0   0 x'];


itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

for itm = set_nums
    take = randperm(length(itemtype));
    itemtypeorder = itemtype(take);
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.itm'],'w+');
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ int1';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = item_locations{itm}; %first 2 are controlled sequences, 2nd 2 random
    
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for seq = 1:size(temp_matrix,2)
        rot = rotation{take(seq)}(randi(length(rotation{take(seq)})));
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,seq);
            y_pos(point) = temp_matrix(2*point,seq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            if itemtypeorder(seq) >= 10
                type_space = '   ';%3 spaces
            else
                type_space = '    ';%3 spaces
            end
            filled_space = '      ';%5 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rot < 100
                rotationspace = '  ';
                if rot < 10
                    rotation_add = '.00';
                else
                    rotation_add = '.0';
                end
            else
                rotationspace = '   ';
                rotation_add = '';
            end
            
            if any(itemtypeorder(seq) == [1,9]);
                height_width_space = '0.75  0.38';
            elseif any(itemtypeorder(seq) == [2,12]);
                height_width_space = '0.38  0.38';
            else
                height_width_space = '0.75  0.75';
            end
            
            integer_space = '                        ';
            if any(itemtypeorder(seq) == [1,9,14]);
                inneroutter_space = '              ';
            elseif any(itemtypeorder(seq) == [2,12]);
                inneroutter_space = '  0.38        ';
            else
                inneroutter_space = '  0.75        ';
            end
            if isnan(int1(take(seq)))
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) '\r\n'];
            else
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) integer_space num2str(int1(take(seq))) '\r\n'];
            end
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
            filled_space = '      ';%6 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   5.00  5.00  0.00'...
                '              ' color_order(end,:) '\r\n'];
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    
    base_image_itmnum = itmnum;  % get the starting item number for image 1
    % to use later in condition file
    for img = 1:number_of_images;
        if itmnum < 100
            itmspace = ' '; %1 spaces
        else
            itmspace = '';%0 space
        end
        if itm < 10
            setzero = '0';
        else
            setzero = '';
        end
        if img < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        str = [itmspace num2str(itmnum)...
            '    8           0.00    0.00      0                                  '...
            '75  75  75 x   C:\\LSQ' setzero num2str(itm) ...
            '\\S' setzero num2str(itm) 'I' imgzero num2str(img) '.bmp' '\r\n'];
        fprintf(fid,str,'%s');
        itmnum = itmnum+1;
    end
    fclose(fid);
    
    %----------------------------------------------------------------%
    %write unique random conditions files for each day to psueorandomize
    %sequence trials
    line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
    line2 = '  1     -3    1       1                          2     3'; %color change line
    
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.cnd'],'w+');
    
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_conditions = cell(1,number_of_sequences);
    trials_per_condition = max_sequence_conditions/number_of_sequences;
    
    %for predictable sequnces
    for seq = 1:size(item_locations{itm},2);
        for t = 1:trials_per_condition;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:19];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 20:27];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 28:35];
            end
        end
    end
    
    %organize conditions as desired
    conditions = [];
    %put into familirization block if desired
    if familiarization_block
        order = randperm(numel(all_conditions));
        for i=1:numel(all_conditions);
            conditions = [conditions;all_conditions{order(i)}(1:familiar_trials,:)];
            all_conditions{order(i)}(1:familiar_trials,:) = [];
        end
    end
    for i=1:numel(all_conditions);
        conditions = [conditions;all_conditions{i}];
    end
    if familiarization_block
        order = 1:size(all_conditions,1)*familiar_trials*number_of_sequences;
        order = [order randperm(size(conditions,1)-length(order))+length(order)];
        conditions = conditions(order,:);
    else
        conditions = conditions(randperm(size(conditions,1)),:);
    end
    
    %write to file
    cndline = 2; %1st one is devoted to clrchng
    seq_cnd = 1;
    %if there is familirization block write it first to the file
    if familiarization_block %if want familrization trials before image sets else set to false
        for seq = 1:number_of_sequences
            for ft = 1:familiar_trials
                if cndline < 10
                    cndspace = '  ';
                elseif cndline < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                btfc = '     -3    2                             '; %background timing, fixid, color palate
                teststr = [];
                for t = 1:8;
                    if conditions(seq_cnd,t) < 10;
                        testspace = '     ';
                    else
                        testspace = '    ';
                    end
                    teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                end
                fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                cndline = cndline+1;
                seq_cnd = seq_cnd + 1;
            end
        end
    end
    
    %write images with number_sequence_trials_btwn_images trials in between
    for block = 1:number_of_images/image_spacing
        for nov_rep = 1:2; %write 2x for novel then for repeat presentation
            for imgpair = 1:image_spacing/2;
                imgnums = [imgpair*2-1 imgpair*2]+(block-1)*image_spacing;
                if cndline < 10
                    cndspace = '  ';
                elseif cndspace < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(1)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(2)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
            end
        end
    end
    fclose(fid);
end
%%
%%---[3.1] Generate item and condition file with 2 predictable sequences for sets 7+---%%
set_nums = [7:26];%always start from 7 otherwise won't regenerate perfectly
number_of_sequences = 2;
number_of_items = 4; %# of items per sequence
number_sequence_trials_btwn_images = 2;
number_of_images = 96;
image_spacing = 16; %number of images between novel and repeat presentations

%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
familiarization_block = true; %if want familrization trials before image sets else set to false
familiar_trials = 10; %# of familirization trails if want them
if familiarization_block
    max_sequence_conditions = 384+number_of_sequences*familiar_trials;
else
    max_sequence_conditions = 384;
end

item_locations = {};
overlap_index = [];
buffer = 5; %minimum distance between 2 items
maxdist = 15; %maximum distance between 2 items
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',140804); %seed randomly so we can regenerate files if lost
for i = 1:500;
    valid_points = ones(sizey,sizex);
    item_locations{i} = NaN(number_of_items*2,number_of_sequences);
    seq = 1;
    for cross = 1:number_of_items
        if cross == 1
            available_points = find(valid_points);
            choosen_ind = available_points(randi(length(available_points)));
            [y,x] = ind2sub([sizey,sizex],choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2)<= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        elseif cross == 2
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            choosen_ind = Cind(randi(length(Cind)));
            [y,x] = ind2sub(size(valid_points),choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
            %find angle of quickest scan path from previous 2 items
            dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
            dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
            angle12 = atan2d(dy12,dx12);
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            [Cy,Cx] = ind2sub(size(valid_points),Cind);
            potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
            for pc = 1:length(Cy)
                dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                angle23 = atan2d(dy23,dx23);
                potential_angles(pc) = angle23;
            end
            potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
            angle_difference = potential_angles-angle12;
            angle_difference = abs(angle_difference);
            good_angles = find(angle_difference <=  90);
            if isempty(good_angles);
                break;
            else
                choosen_ind = good_angles(randi(length(good_angles)));
                item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                x =  Cx(choosen_ind);
                y =  Cy(choosen_ind);
                Cind = find(C);
                valid_points(Cind) = 0;
            end
        end
    end
    
    if all(~isnan(item_locations{i}(:,1)))
        overlapind = rem(i,4)+1;
        overlap_index(i) = overlapind;
        seq = 2;
        switch overlapind
            case 1
                cross = 1;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_items
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        item_locations{i}(2*cross-1,seq) = x;
                        item_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 items
                        dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                        dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <=  90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
            case 2
                cross = 2;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                cross = 1;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 items
                    dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                    dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
            case 3
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 2nd is now the 3rd cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 3;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                cross = 2;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                cross  = 1;%rewrtie to 2 so that when we flip it will work out
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 items
                    dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                    dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
                xs = item_locations{i}(1:2:end,seq);
                ys = item_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                item_locations{i}(1:2:end,seq) = xs;
                item_locations{i}(2:2:end,seq) = ys;
            case 4
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 1st is now the 4th cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 4;
                x = item_locations{i}(2*cross-1,1);
                y = item_locations{i}(2*cross,1);
                cross = 1;
                item_locations{i}(2*cross-1,seq) = x;
                item_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_items
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        item_locations{i}(2*cross-1,seq) = x;
                        item_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 items
                        dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
                        dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <= 90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
                xs = item_locations{i}(1:2:end,seq);
                ys = item_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                item_locations{i}(1:2:end,seq) = xs;
                item_locations{i}(2:2:end,seq) = ys;
        end
    end
    
    %code to double check that the sequences we generated meat the desired
    %criterion
    if all(~isnan(item_locations{i}))
        item_locations{i}(1:2:end,:) = item_locations{i}(1:2:end,:)-13;
        item_locations{i}(2:2:end,:) = item_locations{i}(2:2:end,:)-10;
        for seq = 1:size(item_locations{i},2);
            xs =item_locations{i}(1:2:end,seq);
            ys =item_locations{i}(2:2:end,seq);
            
            %double check to make sure distances are good
            d = pdist([xs,ys]);
            if any(d < buffer)
                disp('error item locations too close')
                [min(d),seq]
            end
            dx = diff(xs);
            dy = diff(ys);
            if any(sqrt(dx.^2+dy.^2) > maxdist)
                disp('error locations too far appart')
                crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1;
                if length(crossnum) > 1
                    [crossnum;seq]
                else
                    [crossnum seq]
                end
            end
            if any(abs(xs) > 12) || any(abs(ys) > 9)
                disp('error locations out of bounds')
            end
            
            dy12 = ys(2)-ys(1);
            dx12 = xs(2)-xs(1);
            angle12 = atan2d(dy12,dx12);
            
            %angle from item 2 to 3
            dy23 = ys(3)-ys(2);
            dx23 = xs(3)-xs(2);
            angle23 = atan2d(dy23,dx23);
            
            %angle from item 3 to 4
            dy34 = ys(4)-ys(3);
            dx34 = xs(4)-xs(3);
            angle34 = atan2d(dy34,dx34);
            
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            angle23(angle23 < 0) = 360+angle23(angle23 < 0);
            angle34(angle34 < 0) = 360+angle34(angle34 < 0);
            
            dang123 = angle23-angle12;
            dang234 = angle34-angle23;
            dang123(dang123 < 0) =  dang123(dang123 < 0) + 360;
            dang234(dang234 < 0) =  dang123(dang234 < 0) + 360;
            dang123(dang123 > 180) = 360-dang123(dang123 > 180);
            dang234(dang234 > 180) = 360-dang234(dang234 > 180);
            if dang123 > 90 || dang234 > 90
                disp('locations not in a forward direction')
                dang123
                dang234
                seq
            end
        end
    end
end

all_nans = [];
for s = 1:length(item_locations);
    if any(any(isnan(item_locations{s})))
        all_nans = [all_nans s];
    end
end
item_locations(all_nans) = [];
overlap_index(all_nans) = [];
overlap1 = find(overlap_index == 1);
overlap2 = find(overlap_index == 2);
overlap3 = find(overlap_index == 3);
overlap4 = find(overlap_index == 4);
new_item_locations = [];
new_overlap_index = [];
%ensure that the items that overlaps in a sequences rotates uniformly but
% randomly through 1 to 4
for o = 1:min([length(overlap1),length(overlap2),length(overlap3),length(overlap4)]);
    rr = randperm(4);
    for rp = 1:length(rr);
        temp =  eval(['overlap' num2str(rr(rp))]);
        new_item_locations = [new_item_locations item_locations(temp(o))];
        new_overlap_index = [new_overlap_index overlap_index(temp(o))];
    end
end
item_locations = new_item_locations; %reassign cuz this is the variable the subsequent code uses

clr=['255 255   0 x';
    '  0 255 255 x';
    '255   0 255 x';
    '  0 255   0 x';
    '  0   0   0 x'];


itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

for itm = set_nums
    take = randperm(length(itemtype));
    itemtypeorder = itemtype(take);
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.itm'],'w+');
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ int1';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = item_locations{itm}; %first 2 are controlled sequences, 2nd 2 random
    
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for seq = 1:size(temp_matrix,2)
        rot = rotation{take(seq)}(randi(length(rotation{take(seq)})));
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,seq);
            y_pos(point) = temp_matrix(2*point,seq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            if itemtypeorder(seq) >= 10
                type_space = '   ';%3 spaces
            else
                type_space = '    ';%3 spaces
            end
            filled_space = '      ';%5 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rot < 100
                rotationspace = '  ';
                if rot < 10
                    rotation_add = '.00';
                else
                    rotation_add = '.0';
                end
            else
                rotationspace = '   ';
                rotation_add = '';
            end
            
            if any(itemtypeorder(seq) == [1,9]);
                height_width_space = '0.75  0.38';
            elseif any(itemtypeorder(seq) == [2,12]);
                height_width_space = '0.38  0.38';
            else
                height_width_space = '0.75  0.75';
            end
            
            integer_space = '                        ';
            if any(itemtypeorder(seq) == [1,9,14]);
                inneroutter_space = '              ';
            elseif any(itemtypeorder(seq) == [2,12]);
                inneroutter_space = '  0.38        ';
            else
                inneroutter_space = '  0.75        ';
            end
            if isnan(int1(take(seq)))
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) '\r\n'];
            else
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) integer_space num2str(int1(take(seq))) '\r\n'];
            end
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
            filled_space = '      ';%6 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   5.00  5.00  0.00'...
                '              ' color_order(end,:) '\r\n'];
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    
    base_image_itmnum = itmnum;  % get the starting item number for image 1
    % to use later in condition file
    for img = 1:number_of_images;
        if itmnum < 100
            itmspace = ' '; %1 spaces
        else
            itmspace = '';%0 space
        end
        if itm < 10
            setzero = '0';
        else
            setzero = '';
        end
        if img < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        str = [itmspace num2str(itmnum)...
            '    8           0.00    0.00      0                                  '...
            '75  75  75 x   C:\\LSQ' setzero num2str(itm) ...
            '\\S' setzero num2str(itm) 'I' imgzero num2str(img) '.bmp' '\r\n'];
        fprintf(fid,str,'%s');
        itmnum = itmnum+1;
    end
    fclose(fid);
    
    %----------------------------------------------------------------%
    %write unique random conditions files for each day to psueorandomize
    %sequence trials
    line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
    line2 = '  1     -3    1       1                          2     3'; %color change line
    
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.cnd'],'w+');
    
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_conditions = cell(1,number_of_sequences);
    trials_per_condition = max_sequence_conditions/number_of_sequences;
    
    %for predictable sequnces
    for seq = 1:size(item_locations{itm},2);
        for t = 1:trials_per_condition;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:19];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 20:27];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 28:35];
            end
        end
    end
    
    %organize conditions as desired
    conditions = [];
    %put into familirization block if desired
    if familiarization_block
        order = randperm(numel(all_conditions));
        for i=1:numel(all_conditions);
            conditions = [conditions;all_conditions{order(i)}(1:familiar_trials,:)];
            all_conditions{order(i)}(1:familiar_trials,:) = [];
        end
    end
    for i=1:numel(all_conditions);
        conditions = [conditions;all_conditions{i}];
    end
    if familiarization_block
        order = 1:size(all_conditions,1)*familiar_trials*number_of_sequences;
        order = [order randperm(size(conditions,1)-length(order))+length(order)];
        conditions = conditions(order,:);
    else
        conditions = conditions(randperm(size(conditions,1)),:);
    end
    
    %write to file
    cndline = 2; %1st one is devoted to clrchng
    seq_cnd = 1;
    %if there is familirization block write it first to the file
    if familiarization_block %if want familrization trials before image sets else set to false
        for seq = 1:number_of_sequences
            for ft = 1:familiar_trials
                if cndline < 10
                    cndspace = '  ';
                elseif cndline < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                btfc = '     -3    2                             '; %background timing, fixid, color palate
                teststr = [];
                for t = 1:8;
                    if conditions(seq_cnd,t) < 10;
                        testspace = '     ';
                    else
                        testspace = '    ';
                    end
                    teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                end
                fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                cndline = cndline+1;
                seq_cnd = seq_cnd + 1;
            end
        end
    end
    
    %write images with number_sequence_trials_btwn_images trials in between
    for block = 1:number_of_images/image_spacing
        for nov_rep = 1:2; %write 2x for novel then for repeat presentation
            for imgpair = 1:image_spacing/2;
                imgnums = [imgpair*2-1 imgpair*2]+(block-1)*image_spacing;
                if cndline < 10
                    cndspace = '  ';
                elseif cndspace < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(1)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(2)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
            end
        end
    end
    fclose(fid);
end
%%
%%---[3.2] Generate item and condition file with 2 predictable sequences for sets 27+---%%
% modified by Seth Konig October 19-20, 2014
set_nums = [27:72];%always start from 27 otherwise won't regenerate perfectly
number_of_sequences = 2;
number_of_items = 4; %# of items per sequence
number_sequence_trials_btwn_images = 2;
number_of_images = 96;
image_spacing = 16; %number of images between novel and repeat presentations

%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
familiarization_block = true; %if want familrization trials before image sets else set to false
familiar_trials = 10; %# of familirization trails if want them
if familiarization_block
    max_sequence_conditions = 384+number_of_sequences*familiar_trials;
else
    max_sequence_conditions = 384;
end

item_locations = {};
overlap_index = [];
buffer = 5; %minimum distance between 2 items
maxdist = 15; %maximum distance between 2 items
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',141019); %seed randomly so we can regenerate files if lost
for i = 1:1000;
    valid_points = ones(sizey,sizex);
    item_locations{i} = NaN(number_of_items*2,number_of_sequences);
    seq = 1;
    for cross = 1:number_of_items
        if cross == 1
            available_points = find(valid_points);
            choosen_ind = available_points(randi(length(available_points)));
            [y,x] = ind2sub([sizey,sizex],choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2)<= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        elseif cross == 2
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            choosen_ind = Cind(randi(length(Cind)));
            [y,x] = ind2sub(size(valid_points),choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            item_locations{i}(2*cross-1,seq) = x;
            item_locations{i}(2*cross,seq) = y;
        else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
            %find angle of quickest scan path from previous 2 items
            dx12 = item_locations{i}(2*cross-3,seq)-item_locations{i}(2*cross-5,seq);
            dy12 = item_locations{i}(2*cross-2,seq)-item_locations{i}(2*cross-4,seq);
            angle12 = atan2d(dy12,dx12);
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            [Cy,Cx] = ind2sub(size(valid_points),Cind);
            potential_angles = NaN(1,length(Cy)); %angles between previous item and potential locations
            for pc = 1:length(Cy)
                dx23 = Cx(pc)-item_locations{i}(2*cross-3,seq);
                dy23 = Cy(pc)-item_locations{i}(2*cross-2,seq);
                angle23 = atan2d(dy23,dx23);
                potential_angles(pc) = angle23;
            end
            potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
            angle_difference = potential_angles-angle12;
            angle_difference = abs(angle_difference);
            good_angles = find(angle_difference <=  90);
            if isempty(good_angles);
                break;
            else
                choosen_ind = good_angles(randi(length(good_angles)));
                item_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                item_locations{i}(2*cross,seq) = Cy(choosen_ind);
                C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                x =  Cx(choosen_ind);
                y =  Cy(choosen_ind);
                Cind = find(C);
                valid_points(Cind) = 0;
            end
        end
    end
    
    if all(~isnan(item_locations{i}(:,1)))
        seq = 2;
        cross = 1; %same location as previous sequence
        item_locations{i}(2*cross-1,seq) = item_locations{i}(2*cross-1,1);
        item_locations{i}(2*cross,seq) = item_locations{i}(2*cross,1);
        cross = 2; %same location as previous sequence
        item_locations{i}(2*cross-1,seq) = item_locations{i}(2*cross-1,1);
        item_locations{i}(2*cross,seq) = item_locations{i}(2*cross,1);
        cross = 3; %same location as previous sequence
        item_locations{i}(2*cross-1,seq) = item_locations{i}(2*cross-1,1);
        item_locations{i}(2*cross,seq) = item_locations{i}(2*cross,1);
        
        %%%%%%%%%%%%%
        cross = 4;
        dx = item_locations{i}(5,1)-item_locations{i}(3,1);
        dy = item_locations{i}(6,1)-item_locations{i}(4,1);
        angle23 = atan2d(dy,dx); %angle from item 2 to 3
        
        dx = item_locations{i}(7,1)-item_locations{i}(5,1);
        dy = item_locations{i}(8,1)-item_locations{i}(6,1);
        angle34 = atan2d(dy,dx);   %angle from item 3 to 4
        dist34= sqrt(dx.^2+dy.^2); %distance from item 3 to 4 in seq 1
        angle234 = angle34-angle23;
        new_angle = angle23-angle234;
        new_angle(new_angle > 360) = new_angle-360;
        new_angle(new_angle < -360) = new_angle + 360;
        new_change_x = dist34*cosd(new_angle);
        new_change_y = dist34*sind(new_angle);
        item_locations{i}(2*cross-1,seq) = round(10*(item_locations{i}(2*(cross-1)-1,2)+new_change_x))/10; %round to nearst 10th
        item_locations{i}(2*cross,seq) = round(10*(item_locations{i}(2*(cross-1),2)+new_change_y))/10; %round to nearst 10th
        
        %if they are not at least 5 dva apart then skip to next sequence
        if sqrt((item_locations{i}(7,1)-item_locations{i}(7,2))^2+...
                (item_locations{i}(8,1)-item_locations{i}(8,2))^2) < 5
            item_locations{i} = NaN(8,2);
        end
        
        %if get a bad location then skip to next sequence
        if any(item_locations{i}(1:2:end,2) < 1 | item_locations{i}(1:2:end,2) > sizex ...
                | item_locations{i}(2:2:end,2) < 1 | item_locations{i}(2:2:end,2) > sizey)
            item_locations{i} = NaN(8,2);
        end
    end
    
    %code to double check that the sequences we generated meat the desired
    %criterion
    if all(~isnan(item_locations{i}))
        % 1st check that items 1-3 are the same for sequence 1 and sequence 2
        if ~all(item_locations{i}(1:2:end-2,1) == item_locations{i}(1:2:end-2,2))
            disp('Item locations for 1st 3 items are different between sequence 1 and sequence 2')
        end
        % 2nd check that the distance and angles between items 3 and 4 are the same
        % for the 4th item in both sequences
        if sqrt((item_locations{i}(7,1)-item_locations{i}(5,1))^2+...
                (item_locations{i}(8,1)-item_locations{i}(6,1))^2) - ...
                sqrt((item_locations{i}(7,2)-item_locations{i}(5,2))^2+...
                (item_locations{i}(8,2)-item_locations{i}(6,2))^2) > .25 %within 1 dva accuracy
            disp('Distance between 3rd and 4th items in sequences are not the same')
        end
        seq =1;
        xs =item_locations{i}(1:2:end,seq);
        ys =item_locations{i}(2:2:end,seq);
        
        %angle from item 2 to 3
        dy23 = ys(3)-ys(2);
        dx23 = xs(3)-xs(2);
        angle23 = atan2d(dy23,dx23);
        
        %angle from item 3 to 4
        dy34 = ys(4)-ys(3);
        dx34 = xs(4)-xs(3);
        angle34 = atan2d(dy34,dx34);
        
        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
        angle23(angle23 < 0) = 360+angle23(angle23 < 0);
        angle34(angle34 < 0) = 360+angle34(angle34 < 0);
        
        dang234_seq1 = angle34-angle23;
        dang234_seq1(dang234_seq1 <-180) = dang234_seq1+360;
        dang234_seq1(dang234_seq1 > 180) = dang234_seq1-360;
        
        seq =2;
        xs =item_locations{i}(1:2:end,seq);
        ys =item_locations{i}(2:2:end,seq);
        
        %angle from item 2 to 3
        dy23 = ys(3)-ys(2);
        dx23 = xs(3)-xs(2);
        angle23 = atan2d(dy23,dx23);
        
        %angle from item 3 to 4
        dy34 = ys(4)-ys(3);
        dx34 = xs(4)-xs(3);
        angle34 = atan2d(dy34,dx34);
        
        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
        angle23(angle23 < 0) = 360+angle23(angle23 < 0);
        angle34(angle34 < 0) = 360+angle34(angle34 < 0);
        
        dang234_seq2 = angle34-angle23;
        dang234_seq2(dang234_seq2 <-180) = dang234_seq2+360;
        dang234_seq2(dang234_seq2 > 180) = dang234_seq2-360;
        if abs(abs(dang234_seq1) - abs(dang234_seq2)) > 10 %so more or less similar angles
            disp('Angles are not opposite of each other')
        end
        % 3nd check that the sequences do have the right distances and angles
        % between items within a sequence and that they contain a forward bias
        item_locations{i}(1:2:end,:) = item_locations{i}(1:2:end,:)-13;
        item_locations{i}(2:2:end,:) = item_locations{i}(2:2:end,:)-10;
        for seq = 1:size(item_locations{i},2);
            xs =item_locations{i}(1:2:end,seq);
            ys =item_locations{i}(2:2:end,seq);
            
            %double check to make sure distances are good
            d = pdist([xs,ys]);
            if any(d < buffer)
                disp('error item locations too close')
                [min(d),seq]
            end
            dx = diff(xs);
            dy = diff(ys);
            if any(sqrt(dx.^2+dy.^2) > maxdist)
                disp('error locations too far appart')
                crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1;
                if length(crossnum) > 1
                    [crossnum;seq]
                else
                    [crossnum seq]
                end
            end
            if any(abs(xs) > 12) || any(abs(ys) > 9)
                disp('error locations out of bounds')
            end
            
            dy12 = ys(2)-ys(1);
            dx12 = xs(2)-xs(1);
            angle12 = atan2d(dy12,dx12);
            
            %angle from item 2 to 3
            dy23 = ys(3)-ys(2);
            dx23 = xs(3)-xs(2);
            angle23 = atan2d(dy23,dx23);
            
            %angle from item 3 to 4
            dy34 = ys(4)-ys(3);
            dx34 = xs(4)-xs(3);
            angle34 = atan2d(dy34,dx34);
            
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            angle23(angle23 < 0) = 360+angle23(angle23 < 0);
            angle34(angle34 < 0) = 360+angle34(angle34 < 0);
            
            dang123 = angle23-angle12;
            dang234 = angle34-angle23;
            dang123(dang123 < 0) =  dang123(dang123 < 0) + 360;
            dang234(dang234 < 0) =  dang123(dang234 < 0) + 360;
            dang123(dang123 > 180) = 360-dang123(dang123 > 180);
            dang234(dang234 > 180) = 360-dang234(dang234 > 180);
            if dang123 > 90 || dang234 > 90
                disp('locations not in a forward direction')
                dang123
                dang234
                seq
            end
        end
    end
end

all_nans = [];
for s = 1:length(item_locations);
    if any(any(isnan(item_locations{s})))
        all_nans = [all_nans s];
    end
end
item_locations(all_nans) = [];

clr=['255 255   0 x';
    '  0 255 255 x';
    '255   0 255 x';
    '  0 255   0 x';
    '  0   0   0 x'];

itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

for itm = set_nums
    take = randperm(length(itemtype));
    itemtypeorder = itemtype(take);
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.itm'],'w+');
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ int1';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = item_locations{itm-26}; %first 2 are controlled sequences, 2nd 2 random
    
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for seq = 1:size(temp_matrix,2)
        rot = rotation{take(seq)}(randi(length(rotation{take(seq)})));
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,seq);
            y_pos(point) = temp_matrix(2*point,seq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            if itemtypeorder(seq) >= 10
                type_space = '   ';%3 spaces
            else
                type_space = '    ';%3 spaces
            end
            filled_space = '      ';%5 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rot < 100
                rotationspace = '  ';
                if rot < 10
                    rotation_add = '.00';
                else
                    rotation_add = '.0';
                end
            else
                rotationspace = '   ';
                rotation_add = '';
            end
            
            if any(itemtypeorder(seq) == [1,9]);
                height_width_space = '0.75  0.38';
            elseif any(itemtypeorder(seq) == [2,12]);
                height_width_space = '0.38  0.38';
            else
                height_width_space = '0.75  0.75';
            end
            
            if rem(abs(x_pos(point)),ceil(abs(x_pos(point)))) > 0.01 %is a decimal number, can't use zero its a float
                decimalx = '0';
            else
                decimalx = '.00';
            end
            
            if  rem(abs(y_pos(point)),ceil(abs(y_pos(point)))) > 0.01 %is a decimal number, can't  zero its a float
                decimaly = '0';
            else
                decimaly = '.00';
            end
            
            integer_space = '                        ';
            if any(itemtypeorder(seq) == [1,9,14]);
                inneroutter_space = '              ';
            elseif any(itemtypeorder(seq) == [2,12]);
                inneroutter_space = '  0.38        ';
            else
                inneroutter_space = '  0.75        ';
            end
            
            if isnan(int1(take(seq)))
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) decimalx  y_space num2str(y_pos(point)) decimaly '      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) '\r\n'];
            else
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) decimalx y_space num2str(y_pos(point)) decimaly '      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) integer_space num2str(int1(take(seq))) '\r\n'];
            end
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
            filled_space = '      ';%6 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rem(abs(x_pos(point)),ceil(abs(x_pos(point)))) > 0.01 %is a decimal number, can't use zero its a float
                decimalx = '0';
            else
                decimalx = '.00';
            end
            
            if  rem(abs(y_pos(point)),ceil(abs(y_pos(point)))) > 0.01 %is a decimal number, can't  zero its a float
                decimaly = '0';
            else
                decimaly = '.00';
            end
            
            str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
                num2str(x_pos(point)) decimalx y_space num2str(y_pos(point)) decimaly '      0   5.00  5.00  0.00'...
                '              ' color_order(end,:) '\r\n'];
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    
    base_image_itmnum = itmnum;  % get the starting item number for image 1
    % to use later in condition file
    for img = 1:number_of_images;
        if itmnum < 100
            itmspace = ' '; %1 spaces
        else
            itmspace = '';%0 space
        end
        if itm < 10
            setzero = '0';
        else
            setzero = '';
        end
        if img < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        str = [itmspace num2str(itmnum)...
            '    8           0.00    0.00      0                                  '...
            '75  75  75 x   C:\\LSQ' setzero num2str(itm) ...
            '\\S' setzero num2str(itm) 'I' imgzero num2str(img) '.bmp' '\r\n'];
        fprintf(fid,str,'%s');
        itmnum = itmnum+1;
    end
    fclose(fid);
    
    %----------------------------------------------------------------%
    %write unique random conditions files for each day to psueorandomize
    %sequence trials
    line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
    line2 = '  1     -3    1       1                          2     3'; %color change line
    
    if itm < 10
        set = ['ListSQ0' num2str(itm)];
    else
        set = ['ListSQ' num2str(itm)];
    end
    fid = fopen([set '.cnd'],'w+');
    
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_conditions = cell(1,number_of_sequences);
    trials_per_condition = max_sequence_conditions/number_of_sequences;
    
    %for predictable sequnces
    for seq = 1:size(item_locations{itm-26},2);
        for t = 1:trials_per_condition;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:19];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 20:27];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 28:35];
            end
        end
    end
    
    %organize conditions as desired
    conditions = [];
    %put into familirization block if desired
    if familiarization_block
        order = randperm(numel(all_conditions));
        for i=1:numel(all_conditions);
            conditions = [conditions;all_conditions{order(i)}(1:familiar_trials,:)];
            all_conditions{order(i)}(1:familiar_trials,:) = [];
        end
    end
    for i=1:numel(all_conditions);
        conditions = [conditions;all_conditions{i}];
    end
    if familiarization_block
        order = 1:size(all_conditions,1)*familiar_trials*number_of_sequences;
        order = [order randperm(size(conditions,1)-length(order))+length(order)];
        conditions = conditions(order,:);
    else
        conditions = conditions(randperm(size(conditions,1)),:);
    end
    
    %write to file
    cndline = 2; %1st one is devoted to clrchng
    seq_cnd = 1;
    %if there is familirization block write it first to the file
    if familiarization_block %if want familrization trials before image sets else set to false
        for seq = 1:number_of_sequences
            for ft = 1:familiar_trials
                if cndline < 10
                    cndspace = '  ';
                elseif cndline < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                btfc = '     -3    2                             '; %background timing, fixid, color palate
                teststr = [];
                for t = 1:8;
                    if conditions(seq_cnd,t) < 10;
                        testspace = '     ';
                    else
                        testspace = '    ';
                    end
                    teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                end
                fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                cndline = cndline+1;
                seq_cnd = seq_cnd + 1;
            end
        end
    end
    
    %write images with number_sequence_trials_btwn_images trials in between
    for block = 1:number_of_images/image_spacing
        for nov_rep = 1:2; %write 2x for novel then for repeat presentation
            for imgpair = 1:image_spacing/2;
                imgnums = [imgpair*2-1 imgpair*2]+(block-1)*image_spacing;
                if cndline < 10
                    cndspace = '  ';
                elseif cndspace < 100
                    cndspace = ' ';
                else
                    cndspace = '';
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(1)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
                str = [cndspace num2str(cndline) '     -2    3      -4                         ' ...
                    num2str(imgnums(2)+base_image_itmnum-1) '\r\n'];
                fprintf(fid,str);
                cndline = cndline+1;
                for seq = 1:number_sequence_trials_btwn_images
                    if cndline < 10
                        cndspace = '  ';
                    elseif cndline < 100
                        cndspace = ' ';
                    else
                        cndspace = '';
                    end
                    btfc = '     -3    2                             '; %background timing, fixid, color palate
                    teststr = [];
                    for t = 1:8;
                        if conditions(seq_cnd,t) < 10;
                            testspace = '     ';
                        else
                            testspace = '    ';
                        end
                        teststr = [teststr testspace num2str(conditions(seq_cnd,t))];
                    end
                    fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
                    cndline = cndline+1;
                    seq_cnd = seq_cnd+1;
                end
            end
        end
    end
    fclose(fid);
end
%%
%%---[4] Sort images into Image Sets---%%
% secition sorts images from unused library into image sets and keeps track
% of the original name of files. Original files are put into used libary. 
root_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
unused_dir = 'Unused\';
used_dir ='Used\';

cd([root_dir unused_dir]);
setnums = [1:72]; %do not write over original sets!!! Start at 1 will skip ones already created if in root_dir
num_images_per_set = 96;
for set = setnums;
    d=dir([root_dir unused_dir '*.bmp']);
    if set < 10
         set_dir = [root_dir 'LSQ0' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
    else
        set_dir = [root_dir 'LSQ' num2str(set) '\']; %shortened from ListSQ or directory name too long for cortex
    end
    if exist(set_dir,'dir')
        disp('Image Set already exists') %do not let it rewrite original folders
        continue;
    else
        mkdir(set_dir)
    end
    rr = randperm(length(d));
    for img = 1:num_images_per_set;
        original = d(rr(img)).name;
        if set < 10;
            if img < 10
                new = ['S0' num2str(set) 'I0' num2str(img) '.bmp'];
            else
                new = ['S0' num2str(set) 'I' num2str(img) '.bmp'];
            end
        else
            if img < 10
                new = ['S' num2str(set) 'I0' num2str(img) '.bmp'];
            else
                new = ['S' num2str(set) 'I' num2str(img) '.bmp'];
            end
        end
        copyfile([root_dir unused_dir original],[set_dir new])
        movefile([root_dir unused_dir original],[root_dir used_dir original],'f')
    end
end
%% if name files accidentally wrong and want to correct assuming files are the correct files
% 
% root_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
% 
% setnums = [1:8]; %do not write over original sets
% num_images_per_set = 96;
% for set = setnums;
%     if set < 10
%          set_dir = [root_dir 'ListSQ0' num2str(set) '\'];
%     else
%         set_dir = [root_dir 'ListSQ' num2str(set) '\'];
%     end
%     cd(set_dir)
%     for img = 1:num_images_per_set
%         if img < 10
%             old = ['S0' num2str(set) 'I0' num2str(img)];
%         else
%             old = ['S0' num2str(set) 'I' num2str(img)];
%         end
%         if img < 10
%             new = ['S0' num2str(set) 'I0' num2str(img) '.bmp'];
%         else
%             new = ['S0' num2str(set) 'I' num2str(img) '.bmp'];
%         end
%         movefile([set_dir old],[set_dir new],'f')
%     end
%         
% end