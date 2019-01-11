%This code contains all Matlab functions needed to calculate persistence barcodes, 
%bottleneck distances for persistent homology traits, simulation for berry potential, 
%and other geometric features used in the study described in:
%Mao Li, Laura L. Klein, Keith Duncan, Ni Jiang, Daniel H. Chitwood, Jason Londo, Allison J. Miller, Christopher N. Topp (2019)
%Characterizing grapevine (Vitis spp.) inflorescence architecture using X-ray imaging: implications for understanding bunch density

%Code is implemented by Mao Li (maoli0923@gmail.com).

%A few remarks:

%(a) The study was based on mesh .ply files from X-ray imaging.

%(b) Computing persistence barcode needs to install Javaplex http://appliedtopology.github.io/javaplex/

%download all the code and add the directory folder where the code is saved
%using the commend 'addpath'. Then working in the directory where the .ply
%files are.

%Step 1. Load the vertices and faces from .ply files, compute geodesic
%distance, and save as .mat files
file=dir('*.ply');
for k=1:length(file)
    [F V]=ply_read(file(k).name,'tri');
    V=V';
    F=F';
    filename=sprintf(['GR' file(k).name(1:end-4) '.mat']);
    T=genvarname(filename(1:end-4));
    eval([T '.V= V;']);
    eval([T '.F= F;']);
    E=[[F(:,1) F(:,2)];[F(:,1) F(:,3)];[F(:,2) F(:,3)]];
    E=[E;[E(:,2) E(:,1)]];
    E=unique(E,'rows');
    a=find(E(:,1)>E(:,2));
    E(a,:)=[];
    eval([T '.E= E;']);
    a=find(V(:,2)==min(V(:,2)));
    W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2+(V(E(:,1),3)-V(E(:,2),3)).^2);
    G=graph(E(:,1),E(:,2),W);
    dist=distances(G,a);
    clear fun
    for ii=1:length(V)
        fun(ii)=min(dist(:,ii));
    end
    eval([T '.fun= fun;']);
    save(filename,filename(1:end-4))
    clear(filename(1:end-4))
end

%Step 2. Compute persistent homology and save the barcode
%download javaplex and go to their matlab_examples folder
load_javaplex
import edu.stanford.math.plex4.*;
file=dir('*.mat');
 for k=1:length(file)
     load(file(k).name(1:end-4))
     T=genvarname(file(k).name(1:end-4));
     eval(['V=' T '.V;']);
     eval(['E=' T '.E;']);
     eval(['fun=' T '.fun;']);
     stream = api.Plex4.createExplicitSimplexStream();
     for i=1:length(V)
         stream.addVertex(i,-fun(i));
     end
     for i=1:length(E)
         stream.addElement([E(i,1) E(i,2)],max(-fun(E(i,1)),-fun(E(i,2))));
     end
     stream.finalizeStream();
     persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
     intervals = persistence.computeAnnotatedIntervals(stream);
     tmp=intervals.getIntervalsAtDimension(0);
     filename=file(k).name;
     eval([T '.gd= tmp;']);
     save(filename,filename(1:end-4))
 end
 
 %Step 3. Compute pairwise bottleneck distance BD
 %You can use Javaplex function to compute the bottleneck distance by
 %>> BD(i,j) = edu.stanford.math.plex4.bottleneck.BottleneckDistance.computeBottleneckDistance(gdi,gdj);
 %where gdi is gd for i sample, gdj is gd for j sample
 %But for more efficiently compute it, we computed on our server. 
 %The instruction can be found https://github.com/danforthcenter/persistent_homology.git
 
 %Step 4. Compute branch traits and other geometric traits
 file=dir('*.mat');
for k=1:length(file)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['gd=' T '.gd;']);
    output_intervals=[file(k).name(1:end-4),'_',num2str(1),'.txt'];
    diary(output_intervals)
    gd
    diary off
    filename_old=[file(k).name(1:end-4),'_',num2str(1),'.txt'];
    fileID=fopen(filename_old,'r');
    filename_formatted=[file(k).name(1:end-4),'_',num2str(1),'_','right_format','.txt'];
    fileID_formatted=fopen(filename_formatted,'w');
    line=fgetl(fileID);
    while ischar(line)
        if length(line)> 1
            line=regexprep(line,'gd','');
            line=regexprep(line,'=', '');
            line=regexprep(line,'\[','');
            line=regexprep(line,'\), ','\n');
            line=regexprep(line,',','');
            line=regexprep(line,'\)\]','');
            line=regexprep(line,'\]','');
            line=regexprep(line,'infinity','inf');
            line=regexprep(line,'diary off','');
            if length(line)>0
                fprintf(fileID_formatted,line);
            end
        end
        line=fgetl(fileID);
    end
    fclose(fileID_formatted);
    fclose(fileID);
    
    diagram1=load(filename_formatted);
    N_1=size(diagram1,1);
    %Write result to text file
    if k<10
    fileID1=fopen(['diagram000',num2str(k),'.txt'],'w');
    elseif k>9&k<100
       fileID1=fopen(['diagram00',num2str(k),'.txt'],'w'); 
    elseif k>99&k<1000
        fileID1=fopen(['diagram0',num2str(k),'.txt'],'w');
    else
        fileID1=fopen(['diagram',num2str(k),'.txt'],'w');
    end
    for i=1:N_1
        fprintf(fileID1,'%4.4f',diagram1(i,1));
        fprintf(fileID1,' ');
        fprintf(fileID1,'%4.4f\n',diagram1(i,2));
    end
    fclose(fileID1);
end
file=dir('diagram*.txt');
cut=zeros(1,length(file));
 for k=1:length(file)
    diagram1=load(file(k).name);
    N_1=size(diagram1,1);
    [m1 m2]=sort(diagram1(:,2),'descend');
    L=diagram1(m2,2)-diagram1(m2,1);
    b=find(L(2:end)>1);
    cut(k)=m1(b(1)+1);
 end %need double check
%total length, number of tips, rachis length
file=dir('diagram*.txt');
for k=1:length(file)
    PD=dlmread(file(k).name);
    L=PD(:,2)-PD(:,1);
    a=find(L>1);
    TL(k)=sum(L(a));
    NRT(k)=length(a);
    RL(k)=max(L(a));
end
TL=TL-cut;
RL=RL-cut;
%Average length of pedicel
file=dir('diagram*.txt');
PedL=zeros(length(file),1);
for k=1:length(file)
    PD=dlmread(file(k).name);
    L=PD(:,2)-PD(:,1);
    a=find(L>1);
    L=L(a);
    [f,xi]=ksdensity(L);
    [a1 a2]=max(f);
    PedL(s)=xi(a2);
end
%convex hull
file=dir('*.mat');
load cut
for k=1:length(file)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['V=' T '.V;']);
    eval(['fun=' T '.fun;']);
    a=find(fun>=-cut(k));
    [K CHV(k)]=convhulln(V(a,:));%the volume of convex hull
    clear K
    clear(file(k).name(1:end-4))
    clear T
end
%volume and surface area
file=dir('*.mat');
load cut
for k=1:length(file)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['V=' T '.V;']);
    eval(['F=' T '.F;']);
    eval(['fun=' T '.fun;']);
    a=find(fun>=-cut(k));
    V=V(a,:);
    [Lia,Locb]=ismember(F,a);
    b=find(sum(Lia,2)==3);
    F=Locb(b,:);  
    %bounary detection which is the edges that are only used by 1 face
    E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
    [u,m,n] = unique(E,'rows');
    counts = accumarray(n(:), 1);
    O = u(counts==1,:);
    w=unique([O(:,1);O(:,2)]);   
    %count the connected components for the boundary
    [m1 sub]=ismember(O,w);
    A=edgeL2adj([sub ones(size(sub,1),1)]);
    [nComponents,sizes,members] = networkComponents(A);  
    %put a point in the center for each hole and add mesh like pizza
    for jj=1:nComponents
        if length(members{jj})>2
        V=[V;mean(V(w(members{jj}),:))];
        n=length(V);
        [m1 m2]=ismember(O,w(members{jj}));
        a=find(sum(m1,2)==2);
        F=[F;O(a,:) n*ones(length(a),1)];
        end
    end    
    [Vol(k),Area(k)] = mao_mesh2volume(V,F);
    clear(file(k).name(1:end-4))
end

%Step 4 %berry simulation
     %step 4.1 find the pedicel

file2=dir('*.txt');%original_diagrams
file=dir('*.mat');
for k=1:length(file2)
    PD=dlmread(file2(k).name);
    bind=find(PD(:,2)-PD(:,1)>=3);
    L=unique(PD(bind,1));
    
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['V=' T '.V;']);
    eval(['F=' T '.F;']);
    eval(['fun=' T '.fun;']);
    eval(['E=' T '.E;']);
    W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2+(V(E(:,1),3)-V(E(:,2),3)).^2);
    G=graph(E(:,1),E(:,2),W);
    tind=[];
    for j=1:length(L)
        a=find(fun+L(j)<1&fun+L(j)>0);
        Gsub=subgraph(G,a);
        bin=conncomp(Gsub);
        for jj=1:length(unique(bin))
            b=find(bin==jj);
            w=a(b);
            [m1 m2]=max(fun==max(fun(w)));
            c=find(E(:,1)==m2|E(:,2)==m2);
            if max(fun(unique(E(c,:))))<=max(fun(w))
                tind=[tind m2];
            end
        end
    end
    dist=distances(G,tind,tind);
    tind2=[];
    for j=1:length(tind)
        if isempty(find(dist(j,j+1:end)<4))
            tind2=[tind2 tind(j)];
        end
    end
    dist=distances(G,tind2);
    Direction=[];
    for j=1:length(tind2)
        b=find(dist(j,:)<=4.5);
        b0=find(dist(j,:)==0);
        subplot(ceil(length(tind2)/10),10,j)
        [COEFF SCORE]=pca([V(b,:);V(b0,:)]);
        if SCORE(end,1)>0
            Direction=[Direction;COEFF(:,1)'];
        else
            Direction=[Direction;-COEFF(:,1)'];
        end
        scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),5,'filled');hold on;axis equal;view([0 90]);
        scatter3(SCORE(end,1),SCORE(end,2),SCORE(end,3),50,'filled');hold on;axis equal;view([0 90]);
    end
    clear(file(k).name(1:end-4))
    clear T
end

%Step 4.2 plot, double check, and correct the tind2 and direction (especially the sign)
file=dir('*.mat');
for k=1:length(file)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['V=' T '.V;']);
    eval(['F=' T '.F;']);
    eval(['tind2=' T '.tind2;']);
    eval(['Direction=' T '.Direction;']);
    meshplot(V',F');hold on
    scatter3(V(tind2,1),V(tind2,2),V(tind2,3),30,'filled');
    quiver3(V(tind2,1),V(tind2,2),V(tind2,3),Direction(:,1),Direction(:,2),Direction(:,3))
    for j=1:length(tind2)
        text(V(tind2(j),1), V(tind2(j),2),V(tind2(j),3),num2str(j),'color','g'); 
    end   
    clear(file(k).name(1:end-4))
    clear T
end

 %Step 4.3 berry potential simulation based on the real berry diameter constrain
   file=dir('*mat');
   Range=[8 12;8 20;8 15;4 8;-1 8;12 23;8 10;8 12;8 12;8 12];
   for sp=1:10
       ind=find(Species==sp);
       for ssss=1:length(ind)
           k=ind(ssss);
           load(file(k).name(1:end-4))
           T=genvarname(file(k).name(1:end-4));
           eval(['tind2=' T '.tind2;']);
           eval(['Direction=' T '.Direction;']);
           eval(['V=' T '.V;']);
           eval(['E=' T '.E;']);
           W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2+(V(E(:,1),3)-V(E(:,2),3)).^2);
           G=graph(E(:,1),E(:,2),W);
           dist=distances(G,tind2);
           
           [x0 y0 z0]=sphere(50);
           Stop=[];
           M=0.5*ones(length(tind2),3);
           while (length(Stop)<length(tind2))
               if max(M(:,1))>=Range(sp,2)/(2*Scale(k))*100 %where Scale is the scan resolution ~100 um
                   break
               end
               M(Stop,:)=M(Stop,:)-0.1;
               Center=V(tind2,:)+(M-1.5).*Direction;
               DD=pdist(Center);
               DD=squareform(DD);
               Stop=[];
               for j=1:length(tind2)
                   if length(find(DD(j,:)-(M(j,1)+M(:,1))'<0))>1
                       Stop=[Stop j];
                   end
                   w=find(dist(j,:)<=(M(j,1)+2*pi-1));
                   DD2=pdist2(Center(j,:),V(setdiff([1:length(V)],w),:));
                   if min(DD2)<M(j,1)
                       Stop=[Stop j];
                   end
               end
               Stop=unique(Stop);
               M=M+0.1;
           end
           filename=file(k).name;
           eval([T '.berrydiameter= M(:,1)*2;']);
           save(filename,filename(1:end-4))
           clear(file(k).name(1:end-4))
       end
   end
   
%NB--number of pedicel
file=dir('*.mat');
for k=1:length(file)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['tind2=' T '.tind2;']);
    NB(k)=length(tind2);
    clear(file(k).name(1:end-4));
end

%count the touching numbers
file=dir('*.mat');
for k=1:length(length)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['tind2=' T '.tind2;']);
    eval(['Direction=' T '.Direction;']);
    eval(['V=' T '.V;']);
    eval(['F=' T '.F;']);
    eval(['E=' T '.E;']);
    eval(['M=' T '.berrydiameter;']);
    M=M/2;
    W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2+(V(E(:,1),3)-V(E(:,2),3)).^2);
    G=graph(E(:,1),E(:,2),W);
    dist=distances(G,tind2);
    Center=V(tind2,:)+(M-1).*Direction;
    DD=pdist(Center);
    DD=squareform(DD);
    DDD=zeros(size(DD));
    for ii=1:size(DD,1)-1
        for jj=ii+1:size(DD,1)
        DDD(ii,jj)=M(ii)+M(jj);
        end
    end
    N(ssss,1)=length(find(DDD>=DD-0.2))-size(DD,1);
    N(ssss,2)=0;
    for j=1:length(tind2)
    w=find(dist(j,:)<=(M(j)+2*pi-1));
            DD2=pdist2(Center(j,:),V(setdiff([1:length(V)],w),:));
            if min(DD2)<M(j)+0.2
                N(ssss,2)=N(ssss,2)+1;
            end
    end
    clear(file(k).name(1:end-4))
end

%berry potential total volume and diameter
for k=1:length(length)
    load(file(k).name(1:end-4))
    T=genvarname(file(k).name(1:end-4));
    eval(['M=' T '.berrydiameter;']);
    M=M*Scale(ssss)/200;
    berry(ssss,1)=sum(4*pi*M.^3/3);
    berry(ssss,2)=mean(M);
    clear(file(k).name(1:end-4));
end