function [label_up,centroid,spectrall,spectral_v]=merge_seeds2(input,labels,centroid,spectral,similaridade)

    label_up = labels;
    [~,~,d]=size(input);

    for k=1:max(labels(:))
        for t=1:max(labels(:))
            ic=centroid(k).Centroid(2);
            jc=centroid(k).Centroid(1);

            i=centroid(t).Centroid(2);
            j=centroid(t).Centroid(1);
            dist_seed(k,t) = sqrt((i-ic)^2+(j-jc)^2);
        end
    end
    
    i=1;
    while i ~= max(labels(:))
        lista = dist_seed(i,:);
        w = 1;
        for j = lista
            if j> 0
               spectral_dif = pdist([spectral(i,:);spectral(w,:)],'euclidean');
                if spectral_dif <= similaridade
                    label_up(label_up==w) = i;
                    q(label_up==w) = i;
                end 
            end
            w=w+1;
        end
        i = i+1;
    end
    
    classes = unique(label_up);
    classes = setdiff(classes,0);

    for r=1:d
            spectrall(:,r) = struct2array(regionprops(label_up,input(:,:,r),'MeanIntensity'))';
            x = struct2cell(regionprops(label_up,input(:,:,r),'PixelValues'))';
            for i=1:length(unique(classes))
                spectral_v(i,1)=var(x{i});
            end
    end

end

