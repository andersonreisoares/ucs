function [e3,dist]=context_2(idx,labels,centroid,contexto_slic,svm_im,number,probs_c)
 
    classes = unique(number);
    %classes = setdiff(classes,999);
    dist=zeros(1,length(unique(classes)));
    e3=ones(1,length(unique(classes)));
    n=ones(1,length(unique(classes)));
    [i,j]=ind2sub(size(labels),idx);

    region = contexto_slic(i,j);
    lista = regionprops(contexto_slic,svm_im,'PixelValues');  
    %e3 = regionprops(contexto_slic,im_probs,'PixelValues'); 
    tamanho = length(lista(region).PixelValues);
    
    parfor r=1:length(unique(classes))
%         asdsa2 = regionprops(contexto_slic,im_probs(:,:,r),'MeanIntensity');
%         probs_c(:,r) = cell2mat(struct2cell(asdsa2));
        
        t = classes(r);
        
        ic=centroid(t).Centroid(2);
        jc=centroid(t).Centroid(1);
        dist(r) = sqrt((i-ic)^2+(j-jc)^2);
        
        n(r) = sum(lista(region).PixelValues==t);

    end

    e3 = probs_c(region,:)./(tamanho-n);
    e3(e3 == Inf) = 1;
    
%     parfor k =1:length(classes)
%         soma=1;    
%         t = classes(k);
%         
%         ic=centroid(t).Centroid(2);
%         jc=centroid(t).Centroid(1);
%         dist(k) = sqrt((i-ic)^2+(j-jc)^2);
%        
% 
%         n = sum(lista(region).PixelValues==t);
%         x = ((n)/tamanho);
%         soma = 0.5+(1-x);
% 
%         e3(k) = soma;
% 
%      end
end