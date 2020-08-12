files_path = 'C:\Users\Anderson\Documents\Anderson Reis\Algoritmos\Matlab\Code_tese\new\conjunto 1\';
files = dir(strcat(files_path,'*.tif'));

fileID = fopen('tempos.txt','wt'); %prepara o arquivo para escrita dos tempos
fprintf(fileID,'Input Seeds Slic Watershed SVM Classificação CRF Total \n' );


for file = files'
    
    esc = 10;
    
    while esc <= 150
        
        simi = 50;
        while simi <= 400
    
            %Clean!
            clearvars -except file files fileID files_path esc simi;
            
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

            redChannel = imadjust(mat2gray(input(:,:,4)));
            greenChannel = imadjust(mat2gray(input(:,:,3)));
            blueChannel = imadjust(mat2gray(input(:,:,2)));
            rgbImage = cat(3, redChannel, greenChannel, blueChannel);

            tic % Ler input e prepara variáveis
            [l,c,d]=size(input);
            escala = esc;
            similaridade = simi;
%             rotulada_slic = zeros(l,c);
            rotulada_watershed = zeros(l,c);
            svm_i = zeros(l,c);
            h = waitbar(0,'Initializing waitbar...');
            time_input=toc;

            tic %obtém seed centroids e imagem de bordas
            waitbar(0,h,'Computing seeds...')
            [Im_bordas,labels,centroid,spectral]=seed2(input);
            [labels]=merge_seeds2(input,labels,centroid,spectral,escala,similaridade);
            time_seed = toc;

%             tic %Gera o contexto com o SLICO
%             waitbar(0,h,'Analysing context...')
%             [contexto, numlabels] = slicomex(uint16(input(:,:,1:3)), similaridade);
%             contexto_slic=contexto+1;
%             time_slic = toc;

            tic %Segmentação por watershed
            imgDist = imimposemin(Im_bordas,labels);
            bacias = watershed(imgDist);
            bacias = mcleanupregions(double(bacias+1),1);
            seeds = unique(labels);
            %bacias = merge_bacias(labels,bacias,input);
            time_watershed=toc;
            
            if (length(seeds)-1 < 2)
                continue
            end
            

            %Extraí as informações espectrais das bandas da imagem de input para
            %regiões de contexto
            for r=1:d
                %context_spectral_slic(:,r) = regionprops(contexto_slic,input(:,:,r),'MeanIntensity');
                context_spectral_bacias(:,r) = regionprops(bacias,input(:,:,r),'MeanIntensity');
            end

            %Treinamento com SVM
            tic
            waitbar(0,h,'Training SVM...')
            [Md1]=svm(input,labels);
            
            probabilidades_svm=zeros(l,c,length(seeds)-1);
            time_svm=toc;

            pixels = l*c;
            t=0;   

            tic
            %Classificação com SVM e extração das probabilidades para as classes
            %waitbar(0,h,'Computing the posterior probability');
            for i=1:size(input,1)
                for j=1:c
                    x = input(i,j,:);
                    data = squeeze(x)';
                    [label_svm,~,~,e1] = predict(Md1,double(data));
                    probabilidades_svm(i,j,:)=e1;
                    svm_i(i,j)=label_svm;
                    waitbar((t/pixels),h,sprintf('%2.2f%% Classifying with SVM...',((t/pixels)*100)));
                    t=t+1;
                end
            end
            t = isnan(e1);
            if (max(t)==1)
                continue
                continue
                continue
                continue
            end
            time_classify=toc;
            est=1;
            tic
            temp = svm_i;
            q=1;
            %Inicio do cálculo das energias
            while est~=0
                waitbar(0,h,'Initializing CRF...');
                for type=2:2
                    t=0;
                    for i=1:l
                        for j=1:c
                            e2=bordas_update(i,j,input,labels,svm_i,centroid,escala,similaridade); %2ªEnergia

                            if type==1 %Dual para processar com os dois contextos
                                [e3,dist]=context_2(contexto_slic,i,j,input,labels,centroid,context_spectral_slic);
                            else 
                                %e3=context_2(bacias,i,j,centroid,labels,input,escala,spectral,svm_i,context_spectral_bacias,similaridade);
                                [e3,dist]=context_2(bacias,i,j,input,labels,centroid,context_spectral_bacias);
                            end


                            %Calculo das energias
                            e1=squeeze(probabilidades_svm(i,j,:))';

                            e1 = -log(e1);

                            e1 = double(e1)/sum(double(e1(:)));

                            lmap=e1*0.15+e2*0.15+e3*0.7;
                            z=sum(lmap(:));
                            probs=lmap/z;

                            %Em caso de empate o rótulo é definido pela distância
                            if size(find((probs)==min(probs(:))),2)>1
                                rotulo=find((dist)==min(dist(:)));%%%% distância entre pixels
                            else
                                rotulo = find((probs)==min(probs(:)));
                            end

                            if type==1
                                rotulada_slic(i,j) = rotulo;
                                waitbar((t/pixels),h,sprintf('Loop %d - %2.2f%% Performing CRF+SLIC...',q,((t/pixels)*100)));
                                t=t+1;
                            else 
                                rotulada_watershed(i,j) = rotulo;
                                waitbar((t/pixels),h,sprintf('Loop %d - %2.2f%% Performing CRF+Watershed...',q,((t/pixels)*100)));
                                t=t+1;
                            end


                        end
                    end
                end

                di = rotulada_watershed-temp;
                if mean(di(:)) ~= 0 
                    temp = rotulada_watershed;
                    est=1;
                else
                    est = 0;
                end

                q=1+q;

            end
            time_crf=toc;

            C = strsplit(file.name,'.');

            name = C{1};

            %Filtragem para eleminação de ruídos
            %rotulada_slic = mcleanupregions(double(rotulada_slic),1);
            rotulada_watershed = mcleanupregions(double(rotulada_watershed),1);

            %Write image
            if tag==1
                geotiffwrite(strcat(name,'_',int2str(esc),'_',int2str(simi),'_watershed.tif'),mat2gray(rotulada_watershed),R,'CoordRefSysCode',epsg);
                geotiffwrite(strcat(name,'_',int2str(esc),'_',int2str(simi),'_watershed_sementes.tif'),mat2gray(labels),R,'CoordRefSysCode',epsg);
                geotiffwrite(strcat(name,'_',int2str(esc),'_',int2str(simi),'_watershed_svm.tif'),mat2gray(svm_i),R,'CoordRefSysCode',epsg);
            else
                imwrite(cut1, strcat(name,'_watershed.tif'),'WriteMode', 'append');
                imwrite(mat2gray(rotulada_slic),strcat(name,'_slic.tif'))
                imwrite(mat2gray(rotulada_watershed),strcat(name,'_watershed.tif'))
                imwrite(mat2gray(labels),strcat(name,'_sementes.tif'))
                imwrite(mat2gray(svm_i),strcat(name,'_svm.tif'))
            end

            %Escrevendo os arquivos

            delete(h); %deletar barra

            %Gera Figura
            fig=figure('position', [0, 0, 2100, 900]);
            subplot(1,9,1)       
            imshow(rgbImage)           % line plot
            title('Input')

            subplot(1,9,2)      
            Labels = label2rgb(labels, 'jet', 'w', 'shuffle');
            imshow(Labels), title('Seeds');

            subplot(1,9,3)       
            imshow(mat2gray(svm_i))
            title('SVM')

            subplot(1,9,4)       
            imshow(mat2gray(bacias))
            title('Watershed - Contexto')

            subplot(1,9,6)       
            imshow(mat2gray(rotulada_watershed))
            title('CRF with Watershed')

            ei = imerode(rotulada_watershed,ones(3));
            di = imdilate(rotulada_watershed,ones(3));
            boundaries = ei-di;

            r = imoverlay(rgbImage,boundaries,[1 0 0]);

            subplot(1,9,8)
            imshow(r)
            title('CRF + Watershed - Overlap')

%             subplot(1,9,7)
%             imshow(mat2gray(contexto_slic))
%             title('Bordas')
% 
%             subplot(1,9,8)       
%             imshow(mat2gray(rotulada_slic))
%             title('CRF + SLIC')
% 
%             ei = imerode(rotulada_slic,ones(3));
%             di = imdilate(rotulada_slic,ones(3));
%             boundaries = ei-di;
% 
%             r = imoverlay(rgbImage,boundaries,[1 0 0]);
%             subplot(1,9,9)       
%             imshow(r)
%             title('CRF + SLIC')

            saveas(fig,strcat(name,'_',int2str(esc),'_',int2str(simi),'.png'))

            close(fig)

            simi = simi+50;
        end
        esc = esc+10; 
    end
   
end