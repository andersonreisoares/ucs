function [final] = fix_seg(temp)
classes = unique(temp);
[l,c]=size(temp);
final = zeros(l,c,length(classes));
for k=1:length(classes)
    CC = bwconncomp(logical(temp==classes(k)), 8);
    
    %Compute the area of each component:
    S = regionprops(CC, 'Perimeter');

    %Remove small objects:
    L = labelmatrix(CC);
    BW2 = ismember(L, find([S.Perimeter] >= 3));
%     BW2 = imfill(BW2,'holes');
    final(:,:,k) = double(BW2)*classes(k);
end

final = sum(final,3);