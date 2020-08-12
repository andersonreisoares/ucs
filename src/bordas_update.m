%Entropia - #ªEnergia
function [e2]=bordas_update(idx,labels,svm_im,input,similaridade,spectral,classes)

    [l_max,c_max,~] = size(labels);
    
    [i,j]=ind2sub(size(labels),idx);
    centro = svm_im(i,j);
    svm_im(i,j) = 999;
    if (i==1 && j==1)
        w = svm_im(i:i+1,j:j+1);
        e = input(i:i+1,j:j+1,:);
    elseif (i==1 && j<c_max)
        w = svm_im(i:i+1,j-1:j+1);
        e = input(i:i+1,j-1:j+1,:);
    elseif (i==1 && j==c_max)
        w = svm_im(i:i+1,j-1:j);
        e = input(i:i+1,j-1:j,:);
    elseif (i==l_max && j==c_max)
        w = svm_im(i-1:i,j-1:j);
        e = input(i-1:i,j-1:j,:);
    elseif (i==l_max && j>1)
        w = svm_im(i-1:i,j-1:j+1);
        e = input(i-1:i,j-1:j+1,:);
    elseif (i==l_max && j>1)
        w = svm_im(i-1:i,j-1:j+1);
        e = input(i-1:i,j-1:j+1,:);
    elseif (i>1 && j==c_max)
        w = svm_im(i-1:i+1,j-1:j);
        e = input(i-1:i+1,j-1:j,:);
    elseif (i==l_max && j==1)
        w = svm_im(i-1:i,j:j+1);
        e = input(i-1:i,j:j+1,:);
    elseif (i>1 && j==1)
        w = svm_im(i-1:i+1,j:j+1);
        e = input(i-1:i+1,j:j+1,:);
    else
        w = svm_im(i-1:i+1,j-1:j+1);
        e = input(i-1:i+1,j-1:j+1,:);
    end

    [li,co]=size(w);
    classes = unique(classes);
    e2=ones(1,length(classes));
    
    for q=1:li
        for s=1:co   
            dif = norm(double(squeeze(input(i,j,:))' - squeeze(e(q,s,:))'));
            dif_m(q,s) = dif;
        end
    end
    
    dif_m = dif_m;
    pair = exp(-((dif_m.^2)/2*similaridade^2));
   
    for k = 1:length(classes)
        if ismember(classes(k),w)
            for q=1:li
                for s=1:co   
                    if(w(q,s) ~= 999)
                        if (w(q,s) == classes(k))% && classes(k) == centro)
                            e2(classes(k)) = e2(classes(k))+0;%(1+pair(q,s));
                        else
                            e2(classes(k)) = e2(classes(k))+(1-pair(q,s));
                        end
                    else
                        continue
                    end
                end            
            end
            e2(k) = e2(k)/nnz(w==k);
        else
            e2(classes(k)) = 99;
        end
    end
end