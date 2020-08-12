% files_path = 'D:\Tese\WV\dataset\';
 files_path = 'C:\Users\Anderson\Documents\Anderson Reis\Algoritmos\Matlab\Code_tese\new\Conjunto_teste\';
% files_path = 'C:\Users\Anderson\Documents\Anderson Reis\Experimentos Tese\WV_fusionada\Buffer\';
% files_path = 'C:\Users\Anderson\Desktop\Proposto\fudidos\';
files = dir(strcat(files_path,'*.tif'));
out_path = 'C:\Users\Anderson\Desktop\Cmeans\WV\Escala\';
fileID = fopen('tempos.txt','wt'); %prepara o arquivo para escrita dos tempos
fprintf(fileID,'Input Seeds Slic Watershed SVM Classificação CRF Total \n' );

clc;

for file = files'
    
    %Clean!
    clearvars -except file files fileID files_path out_path;
    
    disp(['Processing file: ',file.name]);
    
    %Check if is a valid geotiff
    info = imfinfo(strcat(files_path,file.name));
    tag = isfield(info,'GeoKeyDirectoryTag');
    geoinfo = info.GeoKeyDirectoryTag;
    epsg = geoinfo(28);

    if tag == 1
        [input, R] = geotiffread(strcat(files_path,file.name));
        depth = info.BitsPerSample(1);
    else
        input = imread(strcat(files_path,file.name));
    end
    input = double(input);
    C = strsplit(file.name,'.');
    name = C{1};
    
    redChannel = imadjust(mat2gray(input(:,:,4)));
    greenChannel = imadjust(mat2gray(input(:,:,3)));
    blueChannel = imadjust(mat2gray(input(:,:,2)));
    rgbImage = cat(3, redChannel, greenChannel, blueChannel);
    
    tic % Ler input e prepara variáveis
    [l,c,d]=size(input);
    %landsat
    escala = 1000;
    similaridade = 220;
    que = 1;

    rotulada_watershed = zeros(l*c,1);
    h = waitbar(0,'Initializing waitbar...');
    time_input=toc;
    
    tic %obtém seed centroids e imagem de bordas
    waitbar(0,h,'Computing seeds...')
    [Im_bordas,labels,centroid,spectral,edges]=seed(input);
    numberss = length(unique(labels));
    [labels,centroid,spectral]=merge_seeds2(input,labels,centroid,spectral,similaridade,escala);
    number = unique(labels);
    number = setdiff(number,0);
    
    for i=1:length(number)
       labels(labels==number(i))=i;
    end
    number = unique(labels);
    number = setdiff(number,0);
    time_seed = toc;
    
    rotulos = [1:length(number)];

    tic %Segmentação por watershed
    imgDist = imimposemin(Im_bordas,labels);
    bacias = watershed(imgDist);
    bacias = mcleanupregions(double(bacias+1),1);
    time_watershed=toc;

    input_n = reshape(input,l*c,1,d);
    labels_n = reshape(labels,l*c,1);
    
    %Treinamento com SVM
    tic
    waitbar(0,h,'Training SVM...')
    time_svm=toc;
       
    %Classificação
    tic
    waitbar(0,h,'Classifying SVM...')
    [~,probabilidades_svm]=fcm(squeeze(input_n),length(number)+1,[1.5 10000 1e-15 true]);
    time_classify=toc;
    
   
    probabilidades_svm = probabilidades_svm';
    [~,dim]= size(probabilidades_svm);
    im_probs = reshape(probabilidades_svm,l,c,dim);
    
    parfor idx=1:l*c
        region = find((probabilidades_svm(idx,:))==max(probabilidades_svm(idx,:)))
        svm_i(idx) = region;
    end
    
    svm_im = reshape(svm_i,l,c);
    number = unique(svm_im);

    escala = 5;
    while escala<=60
        tic %Gera o contexto com o SLICO
        waitbar(0,h,'Analysing context...')
        [contexto, numlabels] = slicomex(uint16(input(:,:,2:4)), (l*c)/escala);
        contexto_slic=contexto+1;
        time_slic = toc;
        
        escala
        
        %party is ON!
        tic
        temp = svm_im;
        q = 1;
        stab = 0;
        loop = 1;
        fig2 = imshow(mat2gray(temp));        % line plot
        titulo = sprintf('%s loop %d',file.name,loop);
        title(titulo)  
        %Inicio do cálculo das energias100
        while q<25
            waitbar(0,h,sprintf('Loop %d - Performing...',loop));
            idx=1;

            classes = unique(svm_im);
            parfor r=1:length(unique(classes))
                asdsa2 = regionprops(contexto_slic,im_probs(:,:,r),'MaxIntensity');
                probs_c(:,r) = cell2mat(struct2cell(asdsa2));
            end

            parfor idx=1:l*c

                e2 = bordas_update(idx,labels,temp,input,similaridade,spectral,classes) %2ªEnergia
                [e3,dist] = context_2(idx,labels,centroid,contexto_slic,temp,classes,probs_c);

                e1 = -log(probabilidades_svm(idx,:))
                
                e1 = floor(e1);
                e2 = floor(e2);
                e3 = floor(e3);
                
                if min(e2(:)) < 0
                    e2 = abs(min(e2(:)))+e2;
                end
                
                e1 = e1/sum(e1(:));
                e2 = e2/sum(e2(:));
                
                if mean(e3) ~= 0
                    e3 = e3/sum(e3(:));
                end
                
                lmap = exp(-(e1+e2+e3))
                z = sum(lmap(:));
                probs = lmap/z


                %Em caso de empate o rótulo é definido pela distância
                if size(find((probs)==max(probs(:))),2)>1
                    if size(find((dist)==min(dist(:))),2)>1
                        rotulada_watershed(idx) = find(min(dist(:)));
                    else
                        rotulada_watershed(idx)=find((dist)==min(dist(:)));%%%% distância entre pixels

                    end
                else
                    rotulada_watershed(idx) = find((probs)==max(probs(:)));
                end

                probabilidades_post(idx,:) =  probs;
            end

            if loop == 1
                q=q+1;
                loop = loop + 1;
                temp = reshape(rotulada_watershed,[l,c]);

            else
                test = reshape(rotulada_watershed,[l,c]);

                di = isequal(test,temp);
                if di == 0
                    temp = test;
                    q = q+1;
                    loop = loop + 1;
                else 
                    if stab < 2
                        temp = test;
                        q = q+1;
                        loop = loop + 1;
                        stab = stab+1;
                    else
                        q = 101;
                    end
                end
            end
            im_probs = reshape(probabilidades_post,l,c,dim);
            fig2 = imshow(mat2gray(temp));        % line plot
            titulo = sprintf('%s loop %d',file.name,loop);
            title(titulo)     

        end
        time_crf=toc;

        %Filtragem para eleminação de ruídos  
        [rotulada_watershed, ~] = makeregionsdistinct(test);
        [rotulada_watershed] = cleanupregions(rotulada_watershed, 4, 8);

        %Gera Figura
        fig=figure('position', [0, 0, 2100, 900]);
        subplot(1,7,1)       
        imshow(rgbImage)           % line plot
        title('Input')

        subplot(1,7,2)      
        Labels = label2rgb(labels, 'jet', 'w', 'shuffle');
        imshow(Labels), title('Seeds');

        subplot(1,7,3)       
        imshow(mat2gray(svm_im))
        title('SVM')

        subplot(1,7,4)       
        imshow(mat2gray(contexto_slic))
        title('Watershed - Contexto')

        subplot(1,7,5)
        imshow(mat2gray(rotulada_watershed))
        title('CRF with Watershed')

        ei = imerode(rotulada_watershed,ones(3));
        di = imdilate(rotulada_watershed,ones(3));
        boundaries = ei-di;

        r = imoverlay(rgbImage,boundaries,[1 0 0]);

        subplot(1,7,6)
        imshow(r)
        title('CRF + Watershed - Overlap')

    %     subplot(1,7,7)       
    %     imshow(mat2gray(seg))
    %     title('Watershed - Contexto')

        saveas(fig,strcat(file.name,'.png'))

        close(fig)

        %Write image
        if tag==1
            geotiffwrite(strcat(out_path,name,'_watershed_',num2str(escala),'.tif'),uint8(rotulada_watershed),R,'CoordRefSysCode',epsg);
            geotiffwrite(strcat(out_path,name,'_watershed_sementes_',num2str(escala),'.tif'),uint8(labels),R,'CoordRefSysCode',epsg);
            geotiffwrite(strcat(out_path,name,'_watershed_svm_',num2str(escala),'.tif'),uint8(svm_im),R,'CoordRefSysCode',epsg);
        else
            imwrite(cut1, strcat(out_path,name,'_watershed.tif'),'WriteMode', 'append');
    %         imwrite(mat2gray(rotulada_slic),strcat(out_path,name,'_slic.tif'))
            imwrite(mat2gray(rotulada_watershed),strcat(out_path,name,'_watershed.tif'))
            imwrite(mat2gray(labels),strcat(out_path,name,'_sementes.tif'))
            imwrite(mat2gray(svm_im),strcat(out_path,name,'_svm.tif'))
        end
        escala = escala+5;
        clearvars contexto_slic asdsa2 probs_c
    end
    %Escreve no TXT
    time_total = time_input+time_seed+time_slic+time_watershed+time_svm+time_classify+time_crf;
    fprintf(fileID,'%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %d\n' ,time_input,time_seed,time_slic, time_watershed, time_svm, time_classify, time_crf,time_total,loop-1);
    
    delete(h); %deletar barra
end
fclose(fileID);