function  [filtered] = filter_regions(label, input, areaThreshold, distThreshold) 

    label = makeregionsdistinct(label);
    
    [n e] = imRAG(label);
    
    %//Read in text and extract properties
    s  = regionprops(label, input(:,:,1), 'Centroid', 'PixelList', 'Area', 'MeanIntensity');
    
    %//Create an array that tells us whether or not we have visited this
    %//centroid
    centroidVisited = false(1,length(s));

    %//Create an array that tells us which object belongs to what cluster
    membershipList = zeros(length(s),1);

    %//List of all centroids for each object
    centroidList = [s.MeanIntensity];
    
    %//Get area segments
    areaList = [s.Area];
    
    %//Initialize cluster count
    clusterNumber = 1;

    %//Map that gives us which pixel belongs to which cluster
    map = zeros(size(label));
    
    distThreshold = ones(1,length(s))*distThreshold;
    
    %//If there are any objects we haven't visited...
    while (any(centroidVisited == false))
        %disp(centroidVisited)
        %//Find one object
        ind = find(centroidVisited == false, 1);
        %disp(ind)
        if s(ind).Area < areaThreshold
            %//Extract its centroid
            cent = s(ind).MeanIntensity;
            
            %//Grab pixels where this object is valid
            pixelLocs = s(ind).PixelList;

            %//Find Euclidean distance squared between this centroid to all the
            %//other centroids
            distCentroids = abs((cent - centroidList)/1000);

            %//Find those locations that are lower than the centroid
            %//Also ensure that we filter out those locations that we have already visited
            
            n = [e(e(:,1)==ind,2)];
            neighboors = zeros(1,length(s));
            
            for el=1:n
                neighboors(n)=1;
            end
            %disp(size(neighboors))
            belowThresh = find(min(distCentroids < distThreshold) & neighboors==1 & centroidVisited==false);
          
            %//Mark them as visited
            centroidVisited(belowThresh) = true;

            %//Assign their membership number
            membershipList(belowThresh) = clusterNumber;

            %//For each object that belongs to this cluster, mark them with this
            %//membership number
            for k = 1 : length(belowThresh)
                placesToMark = s(belowThresh(k)).PixelList;
                map(sub2ind(size(label), placesToMark(:,2), placesToMark(:,1))) = ...
                   clusterNumber;
            end

            %//For the next cluster 
           clusterNumber = clusterNumber + 1; 
        end
            %//Mark them as visited
            centroidVisited(ind) = true;
           
            %//For the next cluster 
           clusterNumber = clusterNumber + 1; 
       
    end
    %//Create a colour map that is the same size as the number of clusters
    colourMap = jet(clusterNumber);
    
    %//This colour map will contain what letters belong to what cluster (colour
    %//coded)
    colourMapRed = colourMap(:,1);

    mapColumn = map(:) + 1;
    redPlane = colourMapRed(mapColumn);

    filtered = reshape(redPlane, size(label,1), size(label,2));
end